/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "list.h"
#include "hashset.h"
#include "bitvec.h"
#include "heap.h"
#include "timer.h"
#include "file_reader.h"

#define CHIMERA_FREE_TAIL_N_BIN	0

static char *spec_read_debug = NULL;

static int clp_debug = 0;

typedef struct {
	uint64_t x:20, y:20, dir:1, legal:1, aux:22;
} pb_frg_t;

typedef struct {
	uint32_t sids[2];
	pb_frg_t pair[2];
} pb_aln_t;
define_list(pbalnv, pb_aln_t);

typedef struct {
	char     *tag;
	uint64_t len:20, clp_x:20, clp_y:20, fix:1, chg:1, closed:2;
	uint32_t obts[2];
	uint64_t alnoff:40, alncnt:24;
	//uint64_t bkpoff:40, bkpcnt:24;
	//uint64_t wgtoff:40, wgtcnt:24;
} pb_seq_t;
define_list(pbseqv, pb_seq_t);

typedef struct {
	uint32_t sid:31, dir:1;
	uint16_t cx, cy;
	uint32_t inc;
} pb_evt_t;
define_list(pbevtv, pb_evt_t);

typedef struct {
	int pos, dir, spur, dep;
} spur_t;
define_list(spurv, spur_t);

typedef struct {
	pbseqv   *seqs;
	pbalnv   *hits;
	u8v      *ptrs;
	cuhash   *tag2idx;
	//u4v      *bkps;
	//u1v      *wgts; // weights * 255 for bins along sequence

	uint32_t bin_size;
	//uint32_t norm_dep;
	int      fix_contained;
	int      min_aln_len;
	//uint64_t ol_len;
	//pbevtv   *evts;
} WTCLP;

WTCLP* init_wtclp(uint32_t bin_size, int fix_contained){
	WTCLP *wt;
	wt = calloc(1, sizeof(WTCLP));
	wt->seqs = init_pbseqv(1024);
	wt->hits = init_pbalnv(1024);
	wt->ptrs = init_u8v(1024);
	wt->tag2idx = init_cuhash(1023);
	//wt->bkps = init_u4v(1024);
	//wt->wgts = init_u1v(1024);
	wt->bin_size = bin_size;
	//wt->norm_dep  = norm_dep;
	wt->fix_contained = fix_contained;
	wt->min_aln_len = 1000;
	//wt->ol_len = 0;
	//wt->evts = init_pbevtv(1024);
	return wt;
}

void free_wtclp(WTCLP *wt){
	ref_apply_array(wt->seqs->buffer, wt->seqs->size, pb_seq_t, free(a->tag));
	free_pbseqv(wt->seqs);
	free_pbalnv(wt->hits);
	free_u8v(wt->ptrs);
	free_cuhash(wt->tag2idx);
	//free_u4v(wt->bkps);
	//free_u1v(wt->wgts);
	//free_pbevtv(wt->evts);
	free(wt);
}

void load_alignments_wtclp(WTCLP *wt, float min_sm, FileReader *fr){
	cuhash_t *h, H;
	pb_aln_t *hit;
	pb_seq_t *pb;
	uint64_t i, beg;
	uint32_t cols[2], sids[2], tmp;
	int n, exists;
	float sm;
	def_counter(num);
	memset(&H, 0, sizeof(cuhash_t));
	cols[0] = 0; cols[1] = 5;
	beg_counter(num);
	while((n = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(n < 12) continue;
		sm = atof(get_col_str(fr, 11));
		if(sm < min_sm) continue;
		hit = next_ref_pbalnv(wt->hits);
		for(i=0;i<2;i++){
			H.key = get_col_str(fr, cols[i] + 0);
			h = prepare_cuhash(wt->tag2idx, H, &exists);
			if(exists){
				hit->sids[i] = h->val;
				pb = ref_pbseqv(wt->seqs, h->val);
			} else {
				h->key = strdup(H.key);
				hit->sids[i] = h->val = wt->seqs->size;
				pb = next_ref_pbseqv(wt->seqs);
				memset(pb, 0, sizeof(pb_seq_t));
				pb->tag = h->key;
				pb->len = atoi(get_col_str(fr, cols[i] + 2));
				pb->clp_x = 0;
				pb->clp_y = pb->len;
				pb->obts[0] = 0;
				pb->obts[1] = pb->len;
				pb->chg = 1;
				pb->fix = 0;
				pb->closed = 0;
			}
			hit->pair[i].dir = (get_col_str(fr, cols[i] + 1)[0] == '-');
			hit->pair[i].x   = atoi(get_col_str(fr, cols[i] + 3));
			hit->pair[i].y   = atoi(get_col_str(fr, cols[i] + 4));
			if(hit->pair[i].dir){
				swap_tmp(hit->pair[i].x, hit->pair[i].y, tmp);
				hit->pair[i].x = pb->len - hit->pair[i].x;
				hit->pair[i].y = pb->len - hit->pair[i].y;
			}
			if(hit->pair[i].x +wt->min_aln_len > hit->pair[i].y){
				wt->hits->size --; break;
			}
		}
		if(i == 2) run_counter(num, 100000, stderr, 0);
	}
	end_counter(num, stderr);
	if(wt->hits->size == 0) return;
	clear_and_encap_u8v(wt->ptrs, 2 * wt->hits->size + 1);
	wt->ptrs->size = 2 * wt->hits->size;
	for(i=0;i<wt->ptrs->size;i++) wt->ptrs->buffer[i] = i;
	sort_array(wt->ptrs->buffer, wt->ptrs->size, uint64_t, wt->hits->buffer[a>>1].sids[a&0x1] > wt->hits->buffer[b>>1].sids[b&0x1]);
	beg = 0; sids[0] = wt->hits->buffer[wt->ptrs->buffer[beg] >> 1].sids[wt->ptrs->buffer[beg] & 0x1];
	for(i=beg=0;i<=wt->ptrs->size;i++){
		if(i == wt->ptrs->size || (sids[1] = wt->hits->buffer[wt->ptrs->buffer[i] >> 1].sids[wt->ptrs->buffer[i] & 0x1]) != sids[0]){
			wt->seqs->buffer[sids[0]].alnoff = beg;
			wt->seqs->buffer[sids[0]].alncnt = i - beg;
			if(i - beg > 1){
				sort_array(wt->ptrs->buffer + beg, i - beg, uint64_t, wt->hits->buffer[a>>1].pair[a&0x1].x > wt->hits->buffer[b>>1].pair[b&0x1].x);
			}
			beg = i;
			sids[0] = sids[1];
		}
	}
}

void set_read_clip_wtclp(WTCLP *wt, char *tag, uint32_t coff, uint32_t clen, uint32_t seqlen){
	pb_seq_t *pb;
	uint32_t sid;
	if((sid = kv_get_cuhash(wt->tag2idx, tag)) == 0xFFFFFFFFU) return;
	pb = ref_pbseqv(wt->seqs, sid);
	if(pb->len != clen){
		fprintf(stderr, " -- Error, %s in %s -- %s:%d --\n", tag, __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
		exit(1);
	}
	pb->obts[0] = coff;
	pb->obts[1] = seqlen;
}

uint64_t call_legal_overlaps_wtclp(WTCLP *wt){
	pb_aln_t *hit;
	pb_seq_t *pb1, *pb2;
	uint64_t i, ret;
	int d[4], s[2];
	ret = 0;
	for(i=0;i<wt->hits->size;i++){
		hit = ref_pbalnv(wt->hits, i);
		pb1 = ref_pbseqv(wt->seqs, hit->sids[0]);
		pb2 = ref_pbseqv(wt->seqs, hit->sids[1]);
		if(pb1->closed && pb2->closed){ hit->pair[0].legal = 0; continue; }
		d[0] = (int)pb1->clp_x - (int)hit->pair[0].x;
		d[1] = (int)hit->pair[0].y - (int)pb1->clp_y;
		if(d[0] + (int)wt->bin_size > 0 && d[1] + (int)wt->bin_size > 0) pb1->fix = 1;
		d[2] = (int)pb2->clp_x - (int)hit->pair[1].x;
		d[3] = (int)hit->pair[1].y - (int)pb2->clp_y;
		if(d[2] + (int)wt->bin_size > 0 && d[3] + (int)wt->bin_size > 0) pb2->fix = 1;
		if(hit->pair[0].dir ^ hit->pair[1].dir){
			s[0] = num_max(d[0], d[3]);
			s[1] = num_max(d[1], d[2]);
		} else {
			s[0] = num_max(d[0], d[2]);
			s[1] = num_max(d[1], d[3]);
		}
		if(s[0] + (int)wt->bin_size < 0 || s[1] + (int)wt->bin_size < 0){
			hit->pair[0].legal = 0; continue;
		}
		s[0] = s[0] > 0? s[0] : 0;
		s[1] = s[1] > 0? s[1] : 0;
		if(s[0] + s[1] + wt->min_aln_len > (int)wt->bin_size + (int)hit->pair[0].y - (int)hit->pair[0].x){
			hit->pair[0].legal = 0; continue;
		}
		hit->pair[0].legal = 1;
		ret ++;
	}
	return ret;
}

void clp_high_err_region_wtclp(WTCLP *wt, uint32_t min_dep, int whole){
	pb_seq_t *pb;
	pb_aln_t *hit;
	u4v *brks;
	uint64_t i, idx;
	uint32_t j, k, fix, dep, max, xx, yy, mx, my, pos, dir;
	brks = init_u4v(1024);
	for(i=0;i<wt->seqs->size;i++){
		pb = ref_pbseqv(wt->seqs, i);
		if(pb->closed) continue;
		clear_u4v(brks);
		fix = 0;
		for(j=0;j<pb->alncnt;j++){
			idx = wt->ptrs->buffer[pb->alnoff + j];
			k = idx & 0x01;
			idx >>= 1;
			hit = ref_pbalnv(wt->hits, idx);
			if(hit->pair[0].legal == 0) continue;
			//if(wt->seqs->buffer[hit->sids[!k]].closed) continue;
			if(wt->fix_contained && hit->pair[k].x < wt->bin_size && hit->pair[k].y + wt->bin_size > pb->len) fix = 1;
			push_u4v(brks, ((hit->pair[k].x << 1)) | 0);
			push_u4v(brks, ((hit->pair[k].y << 1)) | 1);
		}
		if(brks->size == 0){
			pb->clp_x = pb->clp_y = 0; pb->closed = 3;
		} else {
			sort_array(brks->buffer, brks->size, uint32_t, (a >> 1) > (b >> 1));
			if(fix){
				pb->fix = 1;
				if(whole){
				} else {
					pb->clp_x = brks->buffer[0] >> 1;
					pb->clp_y = brks->buffer[brks->size - 1] >> 1;
				}
			} else {
				dep = max = xx = yy = mx = my = 0;
				for(j=0;j<brks->size;j++){
					pos = brks->buffer[j] >> 1;
					dir = brks->buffer[j] & 0x01;
					if(dep >= min_dep){
						yy = pos;
						if(yy - xx > max){
							max = yy - xx; mx = xx; my = yy;
						}
					}
					if(dir){
						dep --;
					} else {
						dep ++;
						if(dep == min_dep) xx = pos;
					}
				}
				if(whole){
					if(mx > wt->bin_size || pb->len - my > wt->bin_size){
						pb->clp_x = pb->clp_y = 0; pb->closed = 3;
					}
				} else {
					pb->clp_x = mx;
					pb->clp_y = my;
				}
			}
		}
	}
	free_u4v(brks);
}

int detect_chimera_wtclp(WTCLP *wt, uint32_t sid, spurv *crss, spurv *chis, int win_size, int min_crs_dep){
	pb_seq_t *pb1, *pb2;
	pb_aln_t *hit;
	spur_t *s;
	uint64_t idx;
	uint32_t i, k, mx;
	int tot_dep, avg_dep, x, y, dep, max, d[4];
	if(min_crs_dep == 0) return 0;
	pb1 = ref_pbseqv(wt->seqs, sid);
	if(pb1->closed) return 0;
	if(pb1->fix) return 0;
	if(pb1->clp_x >= pb1->clp_y) return 0;
	clear_spurv(crss);
	tot_dep = 0;
	for(i=0;i<pb1->alncnt;i++){
		idx = wt->ptrs->buffer[pb1->alnoff + i];
		k = idx & 0x01;
		idx >>= 1;
		hit = ref_pbalnv(wt->hits, idx);
		pb2 = ref_pbseqv(wt->seqs, hit->sids[!k]);
		//if(pb2->closed) continue;
		d[0] = (int)hit->pair[k].x - (int)pb1->clp_x;
		d[1] = (int)pb1->clp_y - (int)hit->pair[k].y;
		d[2] = (int)hit->pair[!k].x - (int)pb2->clp_x;
		d[3] = (int)pb2->clp_y - (int)hit->pair[!k].y;
		if(hit->pair[0].dir ^ hit->pair[1].dir){
			d[2] = d[2] ^ d[3]; d[3] = d[2] ^ d[3]; d[2] = d[2] ^ d[3]; // swap
		}
		if(d[0] > (int)wt->bin_size){
			if(d[2] > (int)wt->bin_size){ // x_spur
				push_spurv(crss, (spur_t){hit->pair[k].x, 0, 1, 0});
				push_spurv(crss, (spur_t){hit->pair[k].x, 1, 0, 0});
				x = 2;
			} else {
				x = 1;
			}
		} else x = 0;
		if(d[1] > (int)wt->bin_size){
			if(d[3] > (int)wt->bin_size){
				push_spurv(crss, (spur_t){hit->pair[k].y, 0, 0, 0});
				push_spurv(crss, (spur_t){hit->pair[k].y, 1, 1, 0});
				y = 2;
			} else {
				y = 1;
			}
		} else y = 0;
		if(x == 2 || y == 2) continue;
		tot_dep += hit->pair[k].y - hit->pair[k].x;
		x = hit->pair[k].x + (x * win_size);
		y = hit->pair[k].y - (y * win_size);
		if(x > y) continue;
		push_spurv(crss, (spur_t){x, 0, 0, 0});
		push_spurv(crss, (spur_t){y, 1, 0, 0});
	}
	clear_spurv(chis);
	sort_array(crss->buffer, crss->size, spur_t, a.pos > b.pos);
	for(i=dep=0;i<crss->size;i++){
		s = ref_spurv(crss, i);
		if(s->dir){
			s->dep = dep;
			dep --;
		} else {
			dep ++;
			s->dep = dep;
		}
		if(s->spur){
			push_spurv(chis, (spur_t){s->pos - win_size, 0, 0, s->dep});
			push_spurv(chis, (spur_t){s->pos - 1,        1, 1, s->dep});
			push_spurv(chis, (spur_t){s->pos,            0, 1, s->dep});
			push_spurv(chis, (spur_t){s->pos + win_size, 1, 0, s->dep});
		}
	}
	avg_dep = (tot_dep + pb1->clp_y - pb1->clp_x) / (pb1->clp_y - pb1->clp_x + 1);
	if((int)chis->size < avg_dep) return 0;
	sort_array(chis->buffer, chis->size, spur_t, a.pos > b.pos);
	dep = 0; max = 0; mx = 0;
	for(i=0;i<chis->size;i++){
		s = ref_spurv(chis, i);
		if(s->dir){
			if(s->spur && (int)dep >= max && s->dep < min_crs_dep){ max = dep; mx = i; }
			dep --;
		} else {
			dep ++;
			if(s->spur && (int)dep >= max && s->dep < min_crs_dep){ max = dep; mx = i; }
		}
	}
	if(max * 2 < avg_dep) return 0;
	s = ref_spurv(chis, mx);
	if(s->dep >= avg_dep) return 0;
	//if(s->dep >= min_crs_dep) return 0;
	if(s->pos <= pb1->clp_x || s->pos >= pb1->clp_y) return 0;
	if(s->pos - pb1->clp_x > pb1->clp_y - s->pos){
		pb1->clp_y = s->pos;
	} else {
		pb1->clp_x = s->pos;
	}
	return 1;
}

uint32_t filter_chimeric_seqs_wtclp(WTCLP *wt, int win_size, int min_crs_dep, int whole){
	spurv *crss, *chis;
	uint32_t i, ret;
	ret = 0;
	crss = init_spurv(1024);
	chis = init_spurv(1024);
	for(i=0;i<wt->seqs->size;i++){
		if(wt->seqs->buffer[i].closed == 0 && detect_chimera_wtclp(wt, i, crss, chis, win_size, min_crs_dep)){
			if(whole) wt->seqs->buffer[i].closed = 1;
			ret ++;
		}
	}
	free_spurv(crss);
	free_spurv(chis);
	return ret;
}

int distinguish_chimera_wtclp(WTCLP *wt, uint32_t sid, spurv *fine,  spurv *crss, spurv *chis, int min_crs_dep){
	pb_seq_t *pb1, *pb2;
	pb_aln_t *hit;
	spur_t *s, *t, S;
	uint64_t idx;
	uint32_t i, j, m, k, mx;
	int x, y, max, d[4];
	if(min_crs_dep == 0) return 0;
	pb1 = ref_pbseqv(wt->seqs, sid);
	if(pb1->closed) return 0;
	if(pb1->fix) return 0;
	if(pb1->clp_x >= pb1->clp_y) return 0;
	memset(&S, 0, sizeof(spur_t));
	S.pos = 0x7FFFFFFF;
	clear_spurv(fine);
	clear_spurv(crss);
	for(i=0;i<pb1->alncnt;i++){
		idx = wt->ptrs->buffer[pb1->alnoff + i];
		k = idx & 0x01;
		idx >>= 1;
		hit = ref_pbalnv(wt->hits, idx);
		pb2 = ref_pbseqv(wt->seqs, hit->sids[!k]);
		//if(pb2->closed) continue;
		d[0] = (int)hit->pair[k].x - (int)pb1->clp_x;
		d[1] = (int)pb1->clp_y - (int)hit->pair[k].y;
		d[2] = (int)hit->pair[!k].x - (int)pb2->clp_x;
		d[3] = (int)pb2->clp_y - (int)hit->pair[!k].y;
		if(hit->pair[0].dir ^ hit->pair[1].dir){
			d[2] = d[2] ^ d[3]; d[3] = d[2] ^ d[3]; d[2] = d[2] ^ d[3]; // swap
		}
		x = 0;
		if(d[0] > (int)wt->bin_size){
			if(d[2] > (int)wt->bin_size){ // x spur
				push_spurv(crss, (spur_t){hit->pair[k].x / wt->bin_size, 0, num_min(hit->pair[k].y, pb1->clp_y) / wt->bin_size, 0});
				x = 1;
			}
		}
		 y = 0;
		if(d[1] > (int)wt->bin_size){
			if(d[3] > (int)wt->bin_size){
				push_spurv(crss, (spur_t){hit->pair[k].y / wt->bin_size, 1, num_max(hit->pair[k].x, pb1->clp_x) / wt->bin_size, 0});
				y = 1;
			}
		}
		if(x == 0 && y == 0) push_spurv(fine, (spur_t){hit->pair[k].x / wt->bin_size, 3, hit->pair[k].y / wt->bin_size, 1});
	}
	if((int)crss->size < min_crs_dep) return 0;
	sort_array(crss->buffer, crss->size, spur_t, num_cmpgtx(a.pos, b.pos, a.dir, b.dir));
	clear_spurv(chis);
	for(j=0,i=1;i<=crss->size;i++){
		s = (i == crss->size)? &S : ref_spurv(crss, i);
		if(s->pos == crss->buffer[j].pos && s->dir == crss->buffer[j].dir) continue;
		//if(i - j >= (uint32_t)min_crs_dep && crss->buffer[j].pos > CHIMERA_FREE_TAIL_N_BIN && crss->buffer[j].pos + CHIMERA_FREE_TAIL_N_BIN < (int)(pb1->clp_y / wt->bin_size)){
		if(i - j >= (uint32_t)min_crs_dep){
			t = next_ref_spurv(chis);
			sort_array(crss->buffer + j, i - j, spur_t, a.spur > b.spur);
			s = ref_spurv(crss, j);
			t->pos = s->pos;
			t->dir = s->dir;
			t->dep = i - j;
			t->spur = s->pos;
			if(clp_debug){
				fprintf(stderr, "DETECTED [%s] pos={%d,%d} dir=%c n_spur=%d\n", pb1->tag, t->pos * wt->bin_size, t->spur * wt->bin_size, "+-"[t->dir], t->dep);
			}
			if(clp_debug){
				for(m=j;m<i;m++){
					s = ref_spurv(crss, m);
					fprintf(stderr, "Quote pos={%d,%d}\n", s->pos * wt->bin_size, s->spur * wt->bin_size);
				}
			}
			if(t->dir){
				t->spur = crss->buffer[j + min_crs_dep].spur;
			} else {
				t->spur = crss->buffer[i - min_crs_dep].spur;
			}
		}
		j = i;
	}
	if(chis->size == 0) return 0;
	// check whether segments are syntenic, also coverage
	clear_spurv(crss);
	for(i=0;i<chis->size;i++){
		s = ref_spurv(chis, i);
		if(clp_debug){
			fprintf(stderr, "SOLVING [%s] pos={%d,%d} dir=%c n_spur=%d\n", pb1->tag, s->pos * wt->bin_size, s->spur * wt->bin_size, "+-"[s->dir], s->dep);
		}
		s->dep = 0;
		s->dep ++; // pb1 support once
		for(j=0;j<fine->size;j++){
			t = ref_spurv(fine, j);
			if(s->dir){
				if(t->pos <= s->spur && t->spur > s->pos){
					s->dep ++;
					if(clp_debug){
						fprintf(stderr, "Plea pos={%d,%d}\n", t->pos * wt->bin_size, t->spur * wt->bin_size);
					}
				}
			} else {
				if(t->pos < s->pos && t->spur >= s->spur){
					s->dep ++;
					if(clp_debug){
						fprintf(stderr, "Plea pos={%d,%d}\n", t->pos * wt->bin_size, t->spur * wt->bin_size);
					}
				}
			}
		}
		if(s->dep < min_crs_dep) push_spurv(crss, *s);
	}
	if(crss->size == 0) return 0;
	sort_array(crss->buffer, crss->size, spur_t, num_cmpgt(a.pos, b.pos));
	x = 0; y = 0; mx = 0; max = 0;
	for(i=0;i<crss->size;i++){
		s = ref_spurv(crss, i);
		if(s->pos - x > max){
			mx = i; max = s->pos - x; y = x;
		}
		x = s->pos;
	}
	if(pb1->clp_y - x * wt->bin_size > (max? max : 1) * wt->bin_size){
		y = pb1->clp_y / wt->bin_size;
	} else {
		x = y;
		y = ref_spurv(crss, mx)->pos;
	}
	pb1->clp_x = x * wt->bin_size;
	pb1->clp_y = y * wt->bin_size;
	return 1;
}

uint32_t distinguish_chimeric_seqs_wtclp(WTCLP *wt, int min_crs_dep, int whole){
	spurv *fine, *crss, *chis;
	uint32_t i, ret;
	ret = 0;
	fine = init_spurv(1024);
	crss = init_spurv(1024);
	chis = init_spurv(1024);
	for(i=0;i<wt->seqs->size;i++){
		if(wt->seqs->buffer[i].closed == 0 && distinguish_chimera_wtclp(wt, i, fine, crss, chis, min_crs_dep)){
			if(whole) wt->seqs->buffer[i].closed = 1;
			ret ++;
		}
	}
	free_spurv(fine);
	free_spurv(crss);
	free_spurv(chis);
	return ret;
}

int test_chimera_wtclp(WTCLP *wt, uint32_t sid, spurv *fine,  spurv *crss, spurv *chis, int min_crs_dep){
	pb_seq_t *pb1, *pb2;
	pb_aln_t *hit;
	spur_t *s, *t, S;
	uint64_t idx;
	uint32_t i, j, m, k, mx, my;
	int x, y, len, max, d[4], ret;
	if(min_crs_dep == 0) return 0;
	pb1 = ref_pbseqv(wt->seqs, sid);
	if(pb1->closed) return 0;
	//if(pb1->fix) return 0;
	if(pb1->clp_x >= pb1->clp_y) return 0;
	memset(&S, 0, sizeof(spur_t));
	S.pos = 0x7FFFFFFF;
	clear_spurv(fine);
	clear_spurv(crss);
	if(spec_read_debug) clp_debug = !strcmp(pb1->tag, spec_read_debug);
	for(i=0;i<pb1->alncnt;i++){
		idx = wt->ptrs->buffer[pb1->alnoff + i];
		k = idx & 0x01;
		idx >>= 1;
		hit = ref_pbalnv(wt->hits, idx);
		pb2 = ref_pbseqv(wt->seqs, hit->sids[!k]);
		//if(pb2->closed) continue;
		//d[0] = (int)hit->pair[k].x - (int)pb1->clp_x;
		//d[1] = (int)pb1->clp_y - (int)hit->pair[k].y;
		//d[2] = (int)hit->pair[!k].x - (int)pb2->clp_x;
		//d[3] = (int)pb2->clp_y - (int)hit->pair[!k].y;
		d[0] = (int)hit->pair[k].x - (int)pb1->clp_x;
		d[1] = (int)pb1->clp_y - (int)hit->pair[k].y;
		d[2] = (int)hit->pair[!k].x - (int)0;
		d[3] = (int)pb2->len - (int)hit->pair[!k].y;
		if(hit->pair[0].dir ^ hit->pair[1].dir){
			d[2] = d[2] ^ d[3]; d[3] = d[2] ^ d[3]; d[2] = d[2] ^ d[3]; // swap
		}
		x = 0;
		if(d[0] > (int)wt->bin_size){
			if(d[2] > (int)wt->bin_size){ // x spur
				push_spurv(crss, (spur_t){hit->pair[k].x / wt->bin_size, 0, num_min(hit->pair[k].y, pb1->clp_y) / wt->bin_size, 0});
				x = 1;
			}
		}
		 y = 0;
		if(d[1] > (int)wt->bin_size){
			if(d[3] > (int)wt->bin_size){
				push_spurv(crss, (spur_t){hit->pair[k].y / wt->bin_size, 1, num_max(hit->pair[k].x, pb1->clp_x) / wt->bin_size, 0});
				y = 1;
			}
		}
		if(x == 0 && y == 0) push_spurv(fine, (spur_t){hit->pair[k].x / wt->bin_size, 3, hit->pair[k].y / wt->bin_size, 1});
	}
	if((int)crss->size < min_crs_dep) return 0;
	if(clp_debug){
		fprintf(stderr, "[%s] Total overlaps:%d\n", pb1->tag, (int)pb1->alncnt);
		fprintf(stderr, "[%s] Fine  overlaps:%d\n", pb1->tag, (int)fine->size);
		fprintf(stderr, "[%s] potential spur:%d\n", pb1->tag, (int)crss->size);
	}
	sort_array(crss->buffer, crss->size, spur_t, num_cmpgt(a.pos, b.pos));
	clear_spurv(chis);
	for(j=0,i=1;i<=crss->size;i++){
		s = (i == crss->size)? &S : ref_spurv(crss, i);
		if(s->pos == crss->buffer[j].pos) continue;
		if(i - j >= (uint32_t)min_crs_dep && crss->buffer[j].pos > CHIMERA_FREE_TAIL_N_BIN && crss->buffer[j].pos + CHIMERA_FREE_TAIL_N_BIN < (int)(pb1->clp_y / wt->bin_size)){
			t = next_ref_spurv(chis);
			s = ref_spurv(crss, j);
			t->pos = s->pos;
			t->dir = 0;
			t->dep = i - j;
			t->spur = s->pos;
			if(clp_debug){
				fprintf(stderr, "DETECTED [%s] pos={%d,%d} dir=%c n_spur=%d\n", pb1->tag, t->pos * wt->bin_size, t->spur * wt->bin_size, "+-"[t->dir], t->dep);
			}
			if(0 && clp_debug){
				for(m=j;m<i;m++){
					s = ref_spurv(crss, m);
					fprintf(stderr, "Quote pos={%d,%d}\n", s->pos * wt->bin_size, s->spur * wt->bin_size);
				}
			}
		}
		j = i;
	}
	if(chis->size == 0) return 0;
	// check whether segments are syntenic, also coverage
	clear_spurv(crss);
	for(i=0;i<fine->size;i++){
		s = ref_spurv(fine, i);
		x = -1; y = 0;
		for(j=0;j<chis->size;j++){
			t = ref_spurv(chis, j);
			if(!(s->pos < t->pos && s->spur > t->pos)) continue;
			if(x == -1) x = j;
			y = j;
		}
		if(x >= 0){
			push_spurv(crss, (spur_t){x, 0, y, i});
		}
	}
	ret = 1;
	mx = 0; my = 0; max = -1;
	if(crss->size){
		sort_array(crss->buffer, crss->size, spur_t, num_cmpgtx(a.pos, b.pos, a.spur, b.spur));
		for(j=i=0;i<=crss->size;i++){
			s = (i == crss->size)? &S : ref_spurv(crss, i);
			if(s->pos == crss->buffer[j].pos && s->spur == crss->buffer[j].spur) continue;
			if(j + min_crs_dep <= i){
				s = ref_spurv(crss, j);
				x = s->pos? chis->buffer[s->pos-1].pos * wt->bin_size : pb1->clp_x;
				y = (s->spur + 1 >= (int)chis->size)? pb1->clp_y : chis->buffer[s->spur + 1].pos * wt->bin_size;
				len = y - x;
				if(clp_debug){
					fprintf(stderr, "Plea %d\t%d,%d\t%d,%d\tdep=%d\n", len, x, y, s->pos, s->spur, i - j);
				}
				if(len > max){
					if(s->pos == 0 && s->spur + 1 == (int)chis->size) ret = 0;
					max = len; mx = x; my = y;
				}
			}
			j = i;
		}
	}
	if(max == -1){
		x = num_max(chis->buffer[0].pos * wt->bin_size, pb1->clp_x);
		y = num_min(chis->buffer[chis->size-1].pos * wt->bin_size, pb1->clp_y);
		if(x >= pb1->clp_y - y){
			pb1->clp_y = x;
		} else {
			pb1->clp_x = y;
		}
		max = y - x;
	} else {
		pb1->clp_x = mx;
		pb1->clp_y = my;
	}
	if(clp_debug){
		fprintf(stderr, "clp: %d\t%d,%d\n", max, pb1->clp_x, pb1->clp_y);
	}
	return ret;
}

uint32_t test_chimeric_seqs_wtclp(WTCLP *wt, int min_crs_dep, int whole){
	spurv *fine, *crss, *chis;
	uint32_t i, ret;
	ret = 0;
	fine = init_spurv(1024);
	crss = init_spurv(1024);
	chis = init_spurv(1024);
	for(i=0;i<wt->seqs->size;i++){
		if(wt->seqs->buffer[i].closed == 0 && test_chimera_wtclp(wt, i, fine, crss, chis, min_crs_dep)){
			if(whole) wt->seqs->buffer[i].closed = 1;
			ret ++;
		}
	}
	free_spurv(fine);
	free_spurv(crss);
	free_spurv(chis);
	return ret;
}

int has_alt_path_wtclp(WTCLP *wt, uint32_t sid, UUhash *hash, u8v *stack){
	pb_seq_t *pb;
	pb_aln_t *hit, *h1, *h2;
	UUhash_t *u, U;
	uint64_t idx1, idx2, ptr1, ptr2;
	uint32_t i, k, d1, d2, pid;
	int x, y;
	{
		pb = ref_pbseqv(wt->seqs, sid);
		if(pb->closed) return 0;
		if(pb->fix) return 1;
		//pb->fix = 0;
		clear_UUhash(hash);
		clear_u8v(stack);
		idx2 = 0;
		for(i=0;i<pb->alncnt;i++){
			ptr1 = wt->ptrs->buffer[pb->alnoff + i];
			k = ptr1 & 0x01;
			hit = ref_pbalnv(wt->hits, ptr1 >> 1);
			if(!hit->pair[0].legal) continue;
			if(wt->fix_contained && hit->pair[k].x < wt->bin_size && hit->pair[k].y + wt->bin_size > pb->len){  pb->fix = 1; return 1; }
			idx1 = (ptr1 << 2);
			if(hit->pair[k].x < pb->clp_x + (int)wt->bin_size){
				idx1 |= 0;
				push_u8v(stack, idx1);
				continue;
			} else if(hit->pair[k].y + (int)wt->bin_size > pb->clp_y){
				idx1 |= 1;
				idx2 ++;
			} else {
				idx1 |= 2;
			}
			kv_put_UUhash(hash, hit->sids[!k], idx1);
		}
		if(idx2 == 0) return 0;
	}
	U.val = 0;
	while(pop_u8v(stack, &ptr1)){
		h1 = ref_pbalnv(wt->hits, ptr1 >> 3);
		d1 = (ptr1 >> 2) & 0x01;
		pid = h1->sids[d1];
		pb = ref_pbseqv(wt->seqs, pid);
		for(i=0;i<pb->alncnt;i++){
			ptr2 = wt->ptrs->buffer[pb->alnoff + i];
			k = ptr2 & 0x01;
			hit = ref_pbalnv(wt->hits, ptr2 >> 1);
			if(!hit->pair[0].legal) continue;
			U.key = hit->sids[!k];
			if((u = get_UUhash(hash, U)) == NULL) continue;
			ptr2 = u->val;
			if(ptr2 & 0x01) return 1;
			delete_UUhash(hash, u);
			push_u8v(stack, ptr2);
		}
		// check pseudo-conntection from contained reads
		reset_iter_UUhash(hash);
		while((u = ref_iter_UUhash(hash))){
			ptr2 = u->val;
			if(((ptr1 ^ ptr2) & 0x02) == 0) continue;
			h2 = ref_pbalnv(wt->hits, ptr2 >> 3);
			d2 = (ptr2 >> 2) & 0x01;
			// whether overlap
			if(h1->pair[!d1].dir ^ h2->pair[!d2].dir){
				x = num_max(pb->len - h1->pair[!d1].y, h2->pair[!d2].x);
				y = num_min(pb->len - h1->pair[!d1].x, h2->pair[!d2].y);
			} else {
				x = num_max(h1->pair[!d1].x, h2->pair[!d2].x);
				y = num_min(h1->pair[!d1].y, h2->pair[!d2].y);
			}
			if(x + wt->min_aln_len > y) continue;
			if(ptr2 & 0x01) return 1;
			delete_UUhash(hash, u);
			push_u8v(stack, ptr2);
		}
	}
	return 0;
}

uint32_t filter_lonely_seqs_wtclp(WTCLP *wt){
	UUhash *hash;
	u8v *stack;
	uint32_t i, ret;
	hash = init_UUhash(13);
	stack = init_u8v(64);
	ret = 0;
	for(i=0;i<wt->seqs->size;i++){
		if(wt->seqs->buffer[i].closed == 0 && has_alt_path_wtclp(wt, i, hash, stack) == 0){
			wt->seqs->buffer[i].closed = 2;
			ret ++;
		}
	}
	free_UUhash(hash);
	free_u8v(stack);
	return ret;
}

void genome_estimation_wtclp(WTCLP *wt){
	pb_seq_t *pb1, *pb2;
	pb_aln_t *hit;
	u64list  *lens;
	u32list  *regs;
	uint64_t tot, max;
	uint32_t node_id, i, k, idx;
	uint32_t dep, avg, x, y, x1, y1, x2, y2, max_dep;
	int d[4], s[2];
	max_dep = 100;
	lens = init_u64list(max_dep);
	zeros_u64list(lens);
	regs = init_u32list(1024);
	tot = 0;
	for(node_id=0;node_id<wt->seqs->size;node_id++){
		pb1 = ref_pbseqv(wt->seqs, node_id);
		if(pb1->closed) continue;
		if(pb1->clp_x >= pb1->clp_y) continue;
		x1 = pb1->clp_x;
		y1 = pb1->clp_y;
		tot += y1 - x1;
		if(pb1->fix) continue; // contained reads
		clear_u32list(regs);
		for(i=0;i<pb1->alncnt;i++){
			idx = wt->ptrs->buffer[pb1->alnoff + i];
			k = idx & 0x01;
			idx >>= 1;
			hit = ref_pbalnv(wt->hits, idx);
			pb2 = ref_pbseqv(wt->seqs, hit->sids[!k]);
			if(pb2->closed) continue;
			x2 = pb2->clp_x;
			y2 = pb2->clp_y;
			d[0] = (int)x1 - (int)hit->pair[k].x;
			d[1] = (int)hit->pair[k].y - (int)y1;
			d[2] = (int)x2 - (int)hit->pair[!k].x;
			d[3] = (int)hit->pair[!k].y - (int)y2;
			if(hit->pair[0].dir ^ hit->pair[1].dir){
				s[0] = num_max(d[0], d[3]);
				s[1] = num_max(d[1], d[2]);
			} else {
				s[0] = num_max(d[0], d[2]);
				s[1] = num_max(d[1], d[3]);
			}
			if(s[0] + (int)wt->bin_size < 0 || s[1] + (int)wt->bin_size < 0) continue;
			s[0] = s[0] > 0? s[0] : 0;
			s[1] = s[1] > 0? s[1] : 0;
			if(s[0] + s[1] + (int)wt->bin_size > hit->pair[k].y - hit->pair[k].x) continue;
			push_u32list(regs, ((hit->pair[k].x + s[0]) << 1) | 0);
			push_u32list(regs, ((hit->pair[k].y - s[1]) << 1) | 1);
		}
		sort_array(regs->buffer, regs->size, uint32_t, (a >> 1) > (b >> 1));
		dep = 0;
		x = 0;
		for(i=0;i<regs->size;i++){
			y = regs->buffer[i] >> 1;
			k = regs->buffer[i] & 0x01;
			if(y > x && dep < max_dep) lens->buffer[dep] += y - x;
			if(k) dep --; else dep ++;
			x = y;
		}
	}
	free_u32list(regs);
	fprintf(stderr, "== Message for debug ==\n");
	fprintf(stderr, "Sequence coverage statistic:\n");
	max = 0; avg = 1;
	for(x=0;x<max_dep;x++){
		if(x && max < lens->buffer[x]){ max = lens->buffer[x]; avg = x; }
		fprintf(stderr, "%12llu", (long long unsigned int)lens->buffer[x]);
		if((x % 10) == 9) fprintf(stderr, "\n");
	}
	if((x % 10) != 9) fprintf(stderr, "\n");
	free_u64list(lens);
	avg ++; // roundup
	fprintf(stderr, "Total aviable sequences: %llu bp\n", (unsigned long long)tot);
	fprintf(stderr, "Average Coverage(?):     %d\n", avg);
	fprintf(stderr, "Genome Size(?):          %llu bp\n", (unsigned long long)(tot / avg));
}

void output_wtclp(WTCLP *wt, FILE *clp){
	pb_seq_t *pb;
	uint32_t i, x, y;
	for(i=0;i<wt->seqs->size;i++){
		pb = ref_pbseqv(wt->seqs, i);
		if(pb->closed){
			x = y = 0;
		} else {
			x = pb->clp_x;
			y = pb->clp_y;
		}
		fprintf(clp, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", pb->tag, x + pb->obts[0], y - x, pb->obts[1], x, y, pb->closed);
	}
}

int usage(){
	printf(
	"WTCLP: Maximizing legal overlap by clipping long reads\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtclp [options]\n"
	" -i <string> Overlap file from wtzmo, +, *\n"
	"             Format: reads1\\t+/-\\tlen1\\tbeg1\\tend1\\treads2\\t+/-\\tlen2\\tbeg2\\tend2\\tscore\\tidentity<float>\\tmat\\tmis\\tins\\tdel\\tcigar\n"
	" -b <string> Long reads retained region, often from wtobt/wtcyc, +\n"
	"             Format: read_name\\toffset\\tlength\\toriginal_len\n"
	" -o          Ouput of reads' regions after clipping, -:stdout, *\n"
	"             Format: read_name\\toffset\\tlength\n"
	" -f          Force overwrite output file\n"
	" -F          Keep full length or clip all\n"
	" -s <int>    Minimum length of alignment, [1000]\n"
	" -m <float>  Minimum identity of alignment, [0.6]\n"
	" -C          Trun off specical reservation for reads contained by others\n"
	"             Default: one read (A) will not be trimmed when it is contained by another read (B).\n"
	//" -c <int>    Minimum depth of overlap along read, [2]\n"
	//" -r <int>    Max margin size of legal overlap, [200]\n"
	" -k <int>    Bin size, [50]\n"
	" -w <int>    Window size used in chimera detection, [1000]\n"
	" -d <int>    Min number of solid overlaps in a suspecting region to reject chimeric, [3]\n"
	" -n <int>    Max turns of iterations, [5]\n"
	" -T          Treat read as a path of many blocks broken by possible chimeric sites, and test whether the path is valid\n"
	"             will disable iteration, connection checking\n"
	" -x <int>    For debug. 1: chimera checking; 2: conntection checking; 4: clip high error ending [7]\n"
	" -v          Verbose\n"
	" -8 <string> Print message for special read\n"
	);
	return 1;
}

int main(int argc, char **argv){
	WTCLP *wt;
	FileReader *fr;
	cplist *ovls, *obts;
	char *outf;
	FILE *out;
	int c, iter, max_iter, force, min_alen, min_dep, bin_size, fix_contained, whole;
	int win_size, min_crs_dep, block_test, debug_x;
	float min_sm;
	uint64_t tol;
	uint32_t nflt, nclp;
	ovls = init_cplist(4);
	obts = init_cplist(4);
	outf = NULL;
	force = 0;
	min_alen = 1000;
	min_sm = 0.6;
	win_size = 1000;
	min_crs_dep = 3;
	bin_size = 50;
	max_iter = 5;
	fix_contained = 1;
	whole = 0;
	block_test = 0;
	debug_x = 7;
	while((c = getopt(argc, argv, "hi:b:o:fFs:m:Cc:r:w:k:d:n:Tx:v8:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'f': force = 1; break;
			case 'F': whole = 1; break;
			case 'i': push_cplist(ovls, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			case 'o': outf = optarg; break;
			case 's': min_alen = atoi(optarg); break;
			case 'm': min_sm = atof(optarg); break;
			case 'C': fix_contained = 0; break;
			case 'c': min_dep = atoi(optarg); break;
			//case 'r': max_margin = atoi(optarg); break;
			case 'w': win_size = atoi(optarg); break;
			case 'd': min_crs_dep = atoi(optarg); break;
			case 'k': bin_size = atoi(optarg); break;
			case 'n': max_iter = atoi(optarg); break;
			case 'T': max_iter = 1; block_test = 1; break;
			case 'x': debug_x = atoi(optarg); break;
			case 'v': clp_debug = 1; break;
			case '8': spec_read_debug = optarg; break;
			default: return usage();
		}
	}
	if(outf == NULL) return usage();
	if(ovls->size == 0) return usage();
	if(!force && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	wt = init_wtclp(bin_size, fix_contained);
	wt->min_aln_len = min_alen;
	fr = fopen_m_filereader(ovls->size, ovls->buffer);
	fprintf(stderr, "[%s] loading alignments\n", date());
	load_alignments_wtclp(wt, min_sm, fr);
	fprintf(stderr, "[%s] Done, %u reads, %llu overlaps\n", date(), (unsigned int)wt->seqs->size, (unsigned long long)wt->hits->size);
	fclose_filereader(fr);
	if(obts->size){
		fprintf(stderr, "[%s] loading reads clips information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(fr->line->string[0] == '#') continue;
			if(c < 4) continue;
			set_read_clip_wtclp(wt, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)), atoi(get_col_str(fr, 3)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date());
	}
	fprintf(stderr, "[%s] clipping based on overlap depth\n", date());
	tol = call_legal_overlaps_wtclp(wt);
	fprintf(stderr, "Before: legal overlaps = %llu\n", (unsigned long long)tol);
	if(debug_x & 0x04) clp_high_err_region_wtclp(wt, min_crs_dep, whole);
	tol = call_legal_overlaps_wtclp(wt);
	fprintf(stderr, "After:  legal overlaps = %llu\n", (unsigned long long)tol);
	fprintf(stderr, "[%s] Done\n", date());
	for(iter=0;iter<max_iter;iter++){
		if(block_test){
			nflt = (debug_x & 0x02)? filter_lonely_seqs_wtclp(wt) : 0;
			fprintf(stderr, "%u reads were filtered by connection-checking\n", nflt);
			nclp = (debug_x & 0x01)? test_chimeric_seqs_wtclp(wt, min_crs_dep, whole) : 0;
		} else {
			fprintf(stderr, "[%s] iteration %d\n", date(), iter + 1);
			nflt = (debug_x & 0x02)? filter_lonely_seqs_wtclp(wt) : 0;
			fprintf(stderr, "%u reads were filtered by connection-checking\n", nflt);
			nclp = (debug_x & 0x01)? filter_chimeric_seqs_wtclp(wt, win_size, min_crs_dep, whole) : 0;
			//nclp = (debug_x & 0x01)? distinguish_chimeric_seqs_wtclp(wt, min_crs_dep, whole) : 0;
		}
		fprintf(stderr, "%u reads were truncated by chimera-checking\n", nclp);
		tol = call_legal_overlaps_wtclp(wt);
		fprintf(stderr, "legal overlaps = %llu\n", (unsigned long long)tol);
		if(nflt + nclp == 0) break;
	}
	fprintf(stderr, "[%s] Done\n", date());
	fprintf(stderr, "--------------------------------\n");
	genome_estimation_wtclp(wt);
	fprintf(stderr, "--------------------------------\n");
	fprintf(stderr, "[%s] output\n", date());
	out = strcmp(outf, "-")? fopen(outf, "w") : stdout;
	output_wtclp(wt, out);
	if(out != stdout) fclose(out);
	fprintf(stderr, "[%s] Done\n", date());
	free_wtclp(wt);
	free_cplist(ovls);
	free_cplist(obts);
	return 0;
}


