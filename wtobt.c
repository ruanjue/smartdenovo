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

/**
 * Update 2015-04-21: trimming high error ends by alignments
 */

#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "kswx.h"
#include "file_reader.h"
#include "bitvec.h"
#include "heap.h"
#include "thread.h"

#define WT_EDGE_BITS	32
#define WT_MAX_EDGE	((1LLU << (WT_EDGE_BITS)) - 1U)
#define WTOBT_CHIMERA_WIN	500

//#define TRIM_ERR

#ifdef TRIM_ERR
#define WTOBT_ERR_WIN	200
#endif

static int WTOBT_MAX_MARGIN = 200;

typedef struct {
	uint64_t edge_offs[2];
	uint32_t edge_cnts[2];
	int32_t  clips[2];
#ifdef TRIM_ERR
	uint64_t win_off:40, win_cnt:24;
#endif
} wt_frg_t;
define_list(wtfrgv, wt_frg_t);

typedef struct {
	uint32_t node_id[2];
	int8_t dir[2];
	int beg[2], end[2];
} wt_hit_t;
define_list(wthitv, wt_hit_t);

typedef struct {
	uint64_t hit_id:63, idx:1;
} wt_ovl_t;
define_list(wtovlv, wt_ovl_t);

typedef struct {
	uint32_t n_rd;
	u32list  *rdlens;
	u32list  *clips[3];
	u32list  *atts; // indicate one read is contained in which overlap
	cplist   *rdnames;
	cuhash   *rdname2id;
	wtfrgv   *frgs;
	wtovlv   *ovls;
	wthitv   *hits;
	u8list   *wins;
	uint32_t min_cov;
	int      min_score;
	int      spur_len, spur_cut;
	float    min_sm;
} WTOBT;

WTOBT* init_wtobt(uint32_t min_cov, int min_score, float min_sm){
	WTOBT *g;
	g = malloc(sizeof(WTOBT));
	g->n_rd = 0;
	g->rdlens = init_u32list(1024);
	g->clips[0] = init_u32list(1024);
	g->clips[1] = init_u32list(1024);
	g->clips[2] = init_u32list(1024);
	g->atts = init_u32list(1024);
	g->rdnames = init_cplist(1024);
	g->rdname2id = init_cuhash(1023);
	g->frgs = init_wtfrgv(1024);
	g->ovls = init_wtovlv(1024);
	g->hits = init_wthitv(1024);
	g->wins = init_u8list(1024);
	g->min_cov = min_cov;
	g->min_score = min_score;
	g->spur_len = 200;
	g->spur_cut = 500;
	g->min_sm = min_sm;
	return g;
}

void free_wtobt(WTOBT *g){
	uint32_t i;
	free_u32list(g->rdlens);
	free_u32list(g->clips[0]);
	free_u32list(g->clips[1]);
	free_u32list(g->clips[2]);
	for(i=0;i<g->rdnames->size;i++) free(get_cplist(g->rdnames, i));
	free_cplist(g->rdnames);
	free_cuhash(g->rdname2id);
	free_wtfrgv(g->frgs);
	free_wtovlv(g->ovls);
	free_wthitv(g->hits);
	free_u8list(g->wins);
	free(g);
}

void push_read_wtobt(WTOBT *g, char *name, int name_len, int seq_len){
	char *ptr;
	push_u32list(g->rdlens, seq_len);
	push_u32list(g->clips[0], 0);
	push_u32list(g->clips[1], seq_len);
	push_u32list(g->clips[2], seq_len);
	push_u32list(g->atts, 0xFFFFFFFFU);
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_cplist(g->rdnames, ptr);
	kv_put_cuhash(g->rdname2id, ptr, g->n_rd);
	g->n_rd ++;
}

void set_read_clip_wtobt(WTOBT *wt, char *name, int coff, int clen){
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->rdname2id, name)) == 0xFFFFFFFFU) return;
	if(coff < 0 || coff + clen > (int)wt->rdlens->buffer[pbid]) return;
	wt->rdlens->buffer[pbid]  = clen;
	wt->clips[0]->buffer[pbid] = coff;
	wt->clips[1]->buffer[pbid] = coff + clen;
}

void generate_frgs_wtobt(WTOBT *g){
#ifdef TRIM_ERR
	wt_frg_t *n;
	uint32_t i, j, len;
#endif
	clear_wtfrgv(g->frgs);
	encap_wtfrgv(g->frgs, g->n_rd);
	g->frgs->size = g->n_rd;
	zeros_wtfrgv(g->frgs);
#ifdef TRIM_ERR
	clear_u8list(g->wins);
	for(i=0;i<g->n_rd;i++){
		n = ref_wtfrgv(g->frgs, i);
		n->win_off = g->wins->size;
		len = g->rdlens->buffer[i];
		len = (len + WTOBT_ERR_WIN - 1) / WTOBT_ERR_WIN;
		n->win_cnt = len;
		for(j=0;j<len;j++) push_u8list(g->wins, 0);
	}
#endif
}

/**
 * Format of overlap file
 * tab-delimited
 * 1, rd1
 * 3, strand1: +/-
 * 2, len1
 * 4, beg1: 0-base position
 * 5, end1: exclude end1 site
 * 6, rd2
 * 8, strand2
 * 7,len2
 * 9, beg2
 * 10, end2
 * 11, score
 * 12, identity
 * 13, match
 * 14, mismatch
 * 15, insertion
 * 16, deletion
 * 17, cigar
 */

int parse_overlap_item_wtobt(WTOBT *g, FileReader *fr, wt_hit_t *dat, u32list *cigars){
#ifdef TRIM_ERR
	wt_frg_t *r;
	int i, j, dir, trm, off, idx, brk, f, op, ln, alns[3];
#endif
	int n, score, lens[2];
	float identity;
	while((n = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(n < 17) continue;
		score = atoi(get_col_str(fr, 10));
		identity = atof(get_col_str(fr, 11));
		if(score < g->min_score) continue;
		if(identity < g->min_sm) continue;
		if((dat->node_id[0] = kv_get_cuhash(g->rdname2id, get_col_str(fr, 0))) == 0xFFFFFFFFU) continue;
		dat->dir[0] = (get_col_str(fr, 1)[0] == '-');
		lens[0] = atoi(get_col_str(fr, 2));
		if(lens[0] != (int)g->rdlens->buffer[dat->node_id[0]]){
			fprintf(stderr, " -- Error: %s length %d != %d in %s -- %s:%d --\n", g->rdnames->buffer[dat->node_id[0]], lens[0], g->rdlens->buffer[dat->node_id[0]],  __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		dat->beg[0] = atoi(get_col_str(fr, 3));
		dat->end[0] = atoi(get_col_str(fr, 4));
		if((dat->node_id[1] = kv_get_cuhash(g->rdname2id, get_col_str(fr, 5))) == 0xFFFFFFFFU) continue;
		if(dat->node_id[0] == dat->node_id[1]) continue;
		dat->dir[1] = (get_col_str(fr, 6)[0] == '-');
		lens[1] = atoi(get_col_str(fr, 7));
		if(lens[1] != (int)g->rdlens->buffer[dat->node_id[1]]){
			fprintf(stderr, " -- Error: %s length %d != %d in %s -- %s:%d --\n", g->rdnames->buffer[dat->node_id[1]], lens[1], g->rdlens->buffer[dat->node_id[1]],  __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		dat->beg[1] = atoi(get_col_str(fr, 8));
		dat->end[1] = atoi(get_col_str(fr, 9));
		clear_u32list(cigars);
#ifdef TRIM_ERR
		kswx_string2cigar(cigars, get_col_str(fr, 16));
		for(i=0;i<2;i++){
			f = i? 2 : 1;
			r = ref_wtfrgv(g->frgs, dat->node_id[i]);
			dir = dat->dir[i];
			off = dir? lens[i] - dat->end[i] : dat->beg[i];
			idx = off / WTOBT_ERR_WIN;
			brk = (idx + 1) * WTOBT_ERR_WIN - off;
			alns[0] = alns[1] = alns[2] = 0;
			for(j=0;j<(int)cigars->size;j++){
				op = cigars->buffer[dir? (int)cigars->size - 1 - j : j] & 0x0f;
				ln = cigars->buffer[dir? (int)cigars->size - 1 - j : j] >> 4;
				while(ln > 0){
					if((op == 0 || op == f) && alns[0] + alns[f] + ln >= brk){
						ln -= brk - (alns[0] + alns[f]);
						alns[op] += brk - (alns[0] + alns[f]);
						if(brk >= 0.8 * WTOBT_ERR_WIN){
							if(alns[0] == 0) alns[0] = 1;
							trm = (alns[0] * 1.0 / (alns[0] + alns[1] + alns[2])) * 255;
							if(g->wins->buffer[r->win_off + idx] < trm) g->wins->buffer[r->win_off + idx] = trm;
						}
						alns[0] = alns[1] = alns[2] = 0;
						brk = WTOBT_ERR_WIN;
						idx ++;
					} else {
						alns[op] += ln;
						ln = 0;
					}
				}
			}
		}
#endif
		return 1;
	}
	return 0;
}

void load_overlaps_wtobt(WTOBT *g, FileReader *fr){
	u32list *cigars;
	wt_hit_t *hit, *b, HIT;
	wt_ovl_t *e1, *e2;
	wt_frg_t *n;
	uint64_t off, ret;
	uint32_t i, k, idx[2];
	idx[0] = 0;
	idx[1] = 0;
	uint32_t wushigang0 = idx[0];
	wushigang0 = idx[1];
	uint32_t tmp0 = wushigang0;
	wushigang0 = tmp0;
	int len1, len2;
	len1 = 0;
	int wushigang = len1;
	len2 = 0;
	wushigang = len2;
	int tmp = wushigang;
	wushigang = tmp;
	hit = &HIT;
	ret = 0;
	cigars = init_u32list(1024);
	while(parse_overlap_item_wtobt(g, fr, hit, cigars)){
		ret ++;
		if((ret % 10000) == 0){ fprintf(stderr, "\r%llu overlaps", (unsigned long long)ret); fflush(stderr); }
		if(ref_wtfrgv(g->frgs, hit->node_id[0])->edge_cnts[hit->dir[0]] >= WT_MAX_EDGE) continue;
		if(ref_wtfrgv(g->frgs, hit->node_id[1])->edge_cnts[!hit->dir[1]] >= WT_MAX_EDGE) continue;
		ref_wtfrgv(g->frgs, hit->node_id[0])->edge_cnts[hit->dir[0]] ++;
		ref_wtfrgv(g->frgs, hit->node_id[1])->edge_cnts[!hit->dir[1]] ++;
		push_wthitv(g->hits, HIT);
	}
	fprintf(stderr, "\r%llu overlaps\n", (unsigned long long)ret);
	free_u32list(cigars);
	off = 0;
	for(k=0;k<2;k++){
		for(i=0;i<g->frgs->size;i++){
			n = ref_wtfrgv(g->frgs, i);
			n->edge_offs[k] = off;
			off += n->edge_cnts[k];
			n->edge_cnts[k] = 0;
		}
	}
	clear_wtovlv(g->ovls);
	encap_wtovlv(g->ovls, off);
	g->ovls->size = off;
	for(ret=0;ret<g->hits->size;ret++){
		if((ret % 10000) == 0){ fprintf(stderr, "\r%llu overlaps", (unsigned long long)ret); fflush(stderr); }
		b = ref_wthitv(g->hits, ret);
		len1 = g->rdlens->buffer[b->node_id[0]];
		len2 = g->rdlens->buffer[b->node_id[1]];
		n = ref_wtfrgv(g->frgs, b->node_id[0]);
		k = b->dir[0];
		idx[0] = n->edge_cnts[k];
		e1 = ref_wtovlv(g->ovls, n->edge_offs[k] + n->edge_cnts[k]);
		e1->hit_id = ret;
		e1->idx    = 0;
		n->edge_cnts[k] ++;

		n = ref_wtfrgv(g->frgs, b->node_id[1]);
		k = !b->dir[1];
		e2 = ref_wtovlv(g->ovls, n->edge_offs[k] + n->edge_cnts[k]);
		e2->hit_id = ret;
		e2->idx    = 1;
		n->edge_cnts[k] ++;
	}
	fprintf(stderr, "\r%llu overlaps\n", (unsigned long long)ret);
}

wt_ovl_t* edge_wtobt(WTOBT *g, uint32_t node_id, int dir, uint32_t eidx){
	wt_frg_t *n;
	wt_ovl_t *e;
	n = ref_wtfrgv(g->frgs, node_id);
	e = ref_wtovlv(g->ovls, n->edge_offs[dir] + eidx);
	return e;
}

thread_beg_def(mobt);
WTOBT *wt;
thread_end_def(mobt);

typedef struct {
	uint32_t x:31, x_spur:1, y:31, y_spur:1;
} wt_reg_t;
define_list(wtregv, wt_reg_t);

typedef struct {
	int pos;
	uint32_t dir:1, spur:1, dep:30;
} wt_brk_t;
define_list(wtbrkv, wt_brk_t);

thread_beg_func(mobt);
WTOBT *wt;
uint64_t tot_dep;
uint32_t tidx, ncpu, node_id, k, i, j, dep, perc, contained;
int x_spur, y_spur, alen, blen, dir, max, mx, my, xx, yy, ol, avg_dep;
wtregv *regs;
wtbrkv *brks, *chis;
wt_frg_t *n;
wt_ovl_t *e;
wt_hit_t *hit;
wt_reg_t *r;
wt_brk_t *b;
regs = init_wtregv(1024);
brks = init_wtbrkv(1024);
chis = init_wtbrkv(1024);
wt = mobt->wt;
thread_beg_loop(mobt);
tidx = thread_index(mobt);
ncpu = thread_n_cpus(mobt);
perc = 0;
for(node_id=tidx;node_id<wt->frgs->size;node_id+=ncpu){
	if(tidx == 0 && ((10000 * ((uint64_t)node_id)) / wt->frgs->size) != perc){
		perc = ((10000 * ((uint64_t)node_id)) / wt->frgs->size);
		fprintf(stderr, "\rprogress: %0.2f%%", perc / 100.f); fflush(stderr);
	}
	n = ref_wtfrgv(wt->frgs, node_id);
	alen = wt->rdlens->buffer[node_id];
	n->clips[0] = 0; n->clips[1] = alen;
	clear_wtregv(regs);
	clear_wtbrkv(brks);
	contained = 0xFFFFFFFFU;
	tot_dep = 0;
	for(k=0;k<2;k++){
		for(j=0;j<n->edge_cnts[k];j++){
			e = ref_wtovlv(wt->ovls, n->edge_offs[k] + j);
			hit = ref_wthitv(wt->hits, e->hit_id);
			blen = wt->rdlens->buffer[hit->node_id[!e->idx]];
			r = next_ref_wtregv(regs);
			dir = hit->dir[e->idx];
			r->x = dir? alen - hit->end[e->idx] : hit->beg[e->idx];
			r->y = dir? alen - hit->beg[e->idx] : hit->end[e->idx];
			ol = r->y - r->x;
			if(ol + 100 >= alen){ contained = (j << 1) | k; break; }
			x_spur = (hit->beg[0] > WTOBT_MAX_MARGIN && hit->beg[1] > WTOBT_MAX_MARGIN);
			y_spur = (hit->end[e->idx] + WTOBT_MAX_MARGIN < alen && hit->end[!e->idx] + WTOBT_MAX_MARGIN < blen);
			if((x_spur || y_spur) && ol < 1000) continue;
			if(x_spur && y_spur) continue;
			r->x_spur = dir? y_spur : x_spur;
			r->y_spur = dir? x_spur : y_spur;
			if(r->x_spur){
				push_wtbrkv(brks, (wt_brk_t){r->x, 0, 1, 0});
				push_wtbrkv(brks, (wt_brk_t){r->x, 1, 0, 0});
			} else if(r->y_spur){
				push_wtbrkv(brks, (wt_brk_t){r->y, 0, 0, 0});
				push_wtbrkv(brks, (wt_brk_t){r->y, 1, 1, 0});
			} else {
				tot_dep += ol;
				push_wtbrkv(brks, (wt_brk_t){r->x, 0, 0, 0});
				push_wtbrkv(brks, (wt_brk_t){r->y, 1, 0, 0});
			}
		}
	}
	if(contained != 0xFFFFFFFFU){
		set_u32list(wt->atts, node_id, contained);
		continue;
	}
	clear_wtbrkv(chis);
	avg_dep = (tot_dep + alen) / (alen + 1);
	sort_array(brks->buffer, brks->size, wt_brk_t, ((a.pos << 1) | a.dir) > ((b.pos << 1) | b.dir));
	xx = yy = 0; dep = 0; max = 0; mx = my = 0;
	for(i=0;i<brks->size;i++){
		b = ref_wtbrkv(brks, i);
		if(dep >= wt->min_cov){
			yy = b->pos;
			if(yy - xx > max){
				mx = xx; my = yy; max = yy - xx;
			}
		}
		if(b->dir){ // is end of a overlap
			b->dep = dep;
			dep --;
		} else { // is beg of a overlap
			dep ++;
			b->dep = dep;
			if(dep == wt->min_cov) xx = b->pos;
		}
		if(b->spur){
			push_wtbrkv(chis, (wt_brk_t){b->pos - WTOBT_CHIMERA_WIN, 0, 0, b->dep});
			push_wtbrkv(chis, (wt_brk_t){b->pos - 1, 1, 0, b->dep});
			push_wtbrkv(chis, (wt_brk_t){b->pos, 0, 1, b->dep});
			push_wtbrkv(chis, (wt_brk_t){b->pos + WTOBT_CHIMERA_WIN, 1, 0, b->dep});
		}
	}
	n->clips[0] = mx; n->clips[1] = my;
	if((int)chis->size < avg_dep) continue;
	sort_array(chis->buffer, chis->size, wt_brk_t, a.pos > b.pos);
	dep = 0; max = 0; mx = 0;
	for(i=0;i<chis->size;i++){
		b = ref_wtbrkv(chis, i);
		if(b->dir){ // is end of a overlap
			if(b->spur && (int)dep >= max){ max = dep; mx = i; }
			dep --;
		} else { // is beg of a overlap
			dep ++;
			if(b->spur && (int)dep >= max){ max = dep; mx = i; }
		}
	}
	if(max * 2 < avg_dep) continue;
	b = ref_wtbrkv(chis, mx);
	if(b->dep >= avg_dep) continue;
	if(2 * b->dep > max + 1) continue;
	if(b->pos <= n->clips[0] || b->pos >= n->clips[1]) continue;
	if(strcmp("pb000000014672", wt->rdnames->buffer[node_id]) == 0){
		fprintf(stderr, " -- spur=%d in %s -- %s:%d --\n", b->pos, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	if(b->pos - n->clips[0] > n->clips[1] - b->pos){
		n->clips[1] = b->pos;
	} else {
		n->clips[0] = b->pos;
	}
}
if(tidx == 0){ fprintf(stderr, "\rprogress: 100.00%%\n");} fflush(stderr);
thread_end_loop(mobt);
free_wtregv(regs);
free_wtbrkv(brks);
free_wtbrkv(chis);
thread_end_func(mobt);

void genome_estimation_wtobt(WTOBT *g){
	wt_frg_t *n;
	wt_ovl_t *e;
	wt_hit_t *hit;
	u64list  *lens;
	u32list  *regs;
	uint64_t tot, max;
	uint32_t node_id, i, k, j;
	int len, dep, avg, x, y, max_dep;
	max_dep = 200;
	lens = init_u64list(max_dep);
	zeros_u64list(lens);
	regs = init_u32list(1024);
	tot = 0;
	for(node_id=0;node_id<g->n_rd;node_id++){
		n = ref_wtfrgv(g->frgs, node_id);
		if(n->clips[0] >= n->clips[1]) continue;
		tot += n->clips[1] - n->clips[0];
		if(g->atts->buffer[node_id] != 0xFFFFFFFFU) continue;
		len = g->rdlens->buffer[node_id];
		clear_u32list(regs);
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_wtovlv(g->ovls, n->edge_offs[k] + j);
				hit = ref_wthitv(g->hits, e->hit_id);
				x = hit->dir[e->idx]? len - hit->end[e->idx] : hit->beg[e->idx];
				y = hit->dir[e->idx]? len - hit->beg[e->idx] : hit->end[e->idx];
				if(x < n->clips[0]) x = n->clips[0];
				if(y > n->clips[1]) y = n->clips[1];
				if(y <= x) continue;
				push_u32list(regs, (x << 1) | 0);
				push_u32list(regs, (y << 1) | 1);
			}
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
		y = len;
		if(y > x) lens->buffer[0] += y - x;
	}
	free_u32list(regs);
	fprintf(stderr, "== Message for debug ==\n");
	fprintf(stderr, "Sequence coverage statistic:\n");
	max = 0; avg = 1;
	for(x=0;x<max_dep;x++){
		if(x && max < lens->buffer[x]){ max = lens->buffer[x]; avg = x; }
		fprintf(stderr, "%20llu", (long long unsigned int)lens->buffer[x]);
		if((x % 5) == 4) fprintf(stderr, "\n");
	}
	if((x % 5) != 4) fprintf(stderr, "\n");
	free_u64list(lens);
	avg ++;
	fprintf(stderr, "Average Coverage(?): %d\n", avg);
	fprintf(stderr, "Genome Size(?):      %llu\n", (unsigned long long)(tot / avg));
}

void process_contained_wtobt(WTOBT *g){
	wt_frg_t *n1, *n2;
	wt_ovl_t *e;
	wt_hit_t *hit;
	uint32_t node_id, j, k;
	int x, y, dx, dy, len, dir;
	for(node_id=0;node_id<g->n_rd;node_id++){
		if(g->atts->buffer[node_id] == 0xFFFFFFFFU) continue;
		n1 = ref_wtfrgv(g->frgs, node_id);
		j = g->atts->buffer[node_id] >> 1;
		k = g->atts->buffer[node_id] & 0x1;
		e = ref_wtovlv(g->ovls, n1->edge_offs[k] + j);
		hit = ref_wthitv(g->hits, e->hit_id);
		len = g->rdlens->buffer[hit->node_id[!e->idx]];
		x = hit->dir[!e->idx]? len - hit->end[!e->idx] : hit->beg[!e->idx];
		y = hit->dir[!e->idx]? len - hit->beg[!e->idx] : hit->end[!e->idx];
		n2 = ref_wtfrgv(g->frgs, hit->node_id[!e->idx]);
		dx = x < n2->clips[0]? n2->clips[0] - x : 0;
		dy = y > n2->clips[1]? y - n2->clips[1] : 0;
		dir = hit->dir[0] ^ hit->dir[1];
		if(dir){
			n1->clips[0] += dy;
			n1->clips[1] -= dx;
		} else {
			n1->clips[0] += dx;
			n1->clips[1] -= dy;
		}
		if(n1->clips[0] >= n1->clips[1]) n1->clips[0] = n1->clips[1] = 0;
	}
}

void print_obt_wtobt(WTOBT *g, FILE *out){
	uint32_t i;
	wt_frg_t *n;
#ifdef TRIM_ERR
	uint32_t j;
#endif
	for(i=0;i<g->frgs->size;i++){
		n = ref_wtfrgv(g->frgs, i);
		fprintf(out, "%s\t%d\t%d\t%d\t%d\t%d\n", get_cplist(g->rdnames, i), n->clips[0] + g->clips[0]->buffer[i], n->clips[1] - n->clips[0], g->clips[2]->buffer[i], n->clips[0], n->clips[1]);
#ifdef TRIM_ERR
		fprintf(out, "#");
		for(j=0;j<n->win_cnt;j++){
			fprintf(out, "\t%0.1f", 100.0 * g->wins->buffer[n->win_off + j] / 255);
		}
		fprintf(out, "\n");
#endif
	}
}

int usage(){
	printf(
	"WTOBT: Overlap based trimming\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtobt [options]\n"
	" -i <string> Long reads sequences file, + *\n"
	" -b <string> Long reads retained region, often from wtobt/wtcyc, +\n"
	"             Format: read_name\\toffset\\tlength\\toriginal_len\n"
	" -j <string> Overlap file(s), + *\n"
	"             Format: reads1\\t+/-\\tlen1\\tbeg1\\tend1\\treads2\\t+/-\\tlen2\\tbeg2\\tend2\\tscore\\tidentity<float>\\tmat\\tmis\\tins\\tdel\\tcigar\n"
	" -o          Ouput of reads' regions after trimming, -:stdout, *\n"
	"             Format: read_name\\toffset\\tlength\n"
	" -f          Force overwrite output file\n"
	" -C          Trun off specical trim for reads contained by others\n"
	"             One read (A) will not be trimmed when it is contained by another read (B).\n"
	"             When trun on special trim (by default), if the B read is trimmed, program will accordingly trim A read\n"
	" -s <int>    Minimum score of alignment, [200]\n"
	" -m <float>  Minimum identity of alignment , [0.5]\n"
	" -w <int>    Maximum margin of alignment, [200]\n"
	" -c <int>    Minimum depth of overlap between anchored reads along reference read, to detect chimeric reads, [1]\n"
	"\n"
	"Example:\n"
	"$> wtobt -i wt.fa -j wt.zmo.ovl -o wt.zmo.obt -m 0.6 -c 2\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	obj_desc_t wsg = wtfrgv_obj_desc;
	wsg = wthitv_obj_desc;
	wsg = wtovlv_obj_desc;
	wsg = wtregv_obj_desc;
	wsg = wtbrkv_obj_desc;
	obj_desc_t ttt = wsg;
	wsg = ttt;
	WTOBT *g;
	FileReader *fr;
	Sequence *seq;
	cplist *pbs, *ovls, *obts;
	char *outf;
	FILE *out;
	uint64_t clips[2];
	uint32_t node_id;
	int c, ncpu, force_overwrite, sw_extend, min_cov, min_score, contained_trim;
	sw_extend = 0;
	int wushigang = sw_extend;
	float min_sm;
	thread_preprocess(mobt);
	wushigang = mobt_j;
	int tmp = wushigang;
	wushigang = tmp;
	pbs = init_cplist(4);
	ovls = init_cplist(4);
	obts = init_cplist(4);
	sw_extend = 0;
	min_cov = 1;
	min_score = 200;
	min_sm = 0.50;
	ncpu = 1;
	force_overwrite = 0;
	contained_trim = 1;
	outf = NULL;
	while((c = getopt(argc, argv, "hft:i:b:j:o:c:s:m:w:C")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'f': force_overwrite = 1; break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			case 'j': push_cplist(ovls, optarg); break;
			case 'o': outf = optarg; break;
			case 's': min_score = atoi(optarg); break;
			case 'm': min_sm = atof(optarg); break;
			case 'c': min_cov = atoi(optarg); break;
			case 'w': WTOBT_MAX_MARGIN = atoi(optarg); break;
			case 'C': contained_trim = 0; break;
			default: return usage();
		}
	}
	if(outf == NULL) return usage();
	if(pbs->size == 0) return usage();
	if(ovls->size == 0) return usage();
	if(!force_overwrite && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	g = init_wtobt(min_cov, min_score, min_sm);
	if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fprintf(stderr, "[%s] loading reads\n", date());
	seq = NULL;
	while(fread_seq(&seq, fr)){
		push_read_wtobt(g, seq->name.string, seq->name.size, seq->seq.size);
		if((g->n_rd % 1000) == 0){
			fprintf(stderr, "\r%u", g->n_rd); fflush(stderr);
		}
	}
	fclose_filereader(fr);
	fprintf(stderr, "\r[%s] Done, %u reads\n", date(), (unsigned)g->n_rd);
	if(obts->size){
		fprintf(stderr, "[%s] loading reads clips information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(fr->line->string[0] == '#') continue;
			if(c < 3) continue;
			set_read_clip_wtobt(g, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date());
	}
	generate_frgs_wtobt(g);
	if((fr = fopen_m_filereader(ovls->size, ovls->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", ovls->buffer[0], __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fprintf(stderr, "[%s] loading alignments\n", date());
	load_overlaps_wtobt(g, fr);
	fclose_filereader(fr);
	fprintf(stderr, "[%s] Done\n", date());
	fprintf(stderr, "[%s] Overlap based trimming\n", date());
	thread_beg_init(mobt, ncpu);
	mobt->wt = g;
	thread_end_init(mobt);
	thread_wake_all(mobt);
	thread_waitfor_all_idle(mobt);
	thread_beg_close(mobt);
	thread_end_close(mobt);
	if(contained_trim){
		fprintf(stderr, "[%s] Specical trimming for contained reads\n", date());
		process_contained_wtobt(g);
	}
	fprintf(stderr, "[%s] Done\n", date());
	genome_estimation_wtobt(g);
	clips[0] = clips[1] = 0;
	for(node_id=0;node_id<g->n_rd;node_id++){
		clips[0] += g->clips[2]->buffer[node_id];
		clips[1] += g->frgs->buffer[node_id].clips[1] - g->frgs->buffer[node_id].clips[0];
	}
	fprintf(stderr, "%llu\t%llu\n", (unsigned long long)clips[0], (unsigned long long)clips[1]);
	out = strcmp(outf, "-")? fopen(outf, "w") : stdout;
	print_obt_wtobt(g, out);
	if(out != stdout) fclose(out);
	free_wtobt(g);
	free_cplist(pbs);
	free_cplist(ovls);
	free_cplist(obts);
	fprintf(stderr, "[%s] Done\n", date());
	return 0;
}

