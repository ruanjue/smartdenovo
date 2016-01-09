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

#ifndef __HZM_ALN_RJ_H
#define __HZM_ALN_RJ_H

#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "bitvec.h"
#include "kswx.h"

#define HZM_MAX_SEED_ZMER	16
#define HZM_MAX_SEED_LEN	0xFFFF
#define HZMP_MAX_SEED_LEN	0xFFFF
#define HZMP_MAX_SEED_OFF	0x7FFFFFFFU
#define hzm_debug_out	stderr

static int hzm_debug = 0;
static int HZM_FAST_WINDOW_KMER_CHAINING = 1;

typedef struct {
	uint32_t mer;
	uint32_t dir:1, off:31;
	uint32_t len;
} hzm_t;
define_list(hzmv, hzm_t);

typedef struct {
	uint64_t mer;
	uint64_t off:47, cnt:16, flt:1;
} hzmh_t;
define_list(hzmhv, hzmh_t);
#define hzmh_hashcode(E) u32hashcode((E).mer)
#define hzmh_hashequals(E1, E2) ((E1).mer == (E2).mer)
define_hashset(hzmhash, hzmh_t, hzmh_hashcode, hzmh_hashequals);

typedef struct {
	uint32_t dir1:1, off1:31;
	uint32_t dir2:1, off2:31;
	uint32_t len1:16, len2:16;
	uint32_t gid; // group id
} hzmp_t;
define_list(hzmpv, hzmp_t);

typedef struct {
	uint32_t pb2;
	uint32_t ovl:29, dir:1, closed:2;
	int beg[2], end[2];
	uint32_t anchors[2];
} wt_seed_t;
define_list(wtseedv, wt_seed_t);

static inline void index_single_read_seeds(uint8_t *pbseq, uint32_t pblen, uint32_t zsize, int hz, uint32_t max_kcnt, hzmhv *hash, BitVec *bits, hzmv *seeds, u32list *hzoff){
	hzm_t   *m;
	uint64_t kmer, krev, kmask, off;
	uint32_t i, j, kcnt, idx;
	uint8_t b, c, dir;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - zsize) << 1);
	clear_hzmv(seeds);
	clear_hzmhv(hash);
	recap_bitvec(bits, kmask);
	zeros_bitvec(bits);
	b = 4;
	kmer = 0;
	clear_u32list(hzoff);
	for(i=j=0;j<pblen;j++){
		c = pbseq[j];
		if(hz && c == b) continue;
		b = c;
		i ++;
		push_u32list(hzoff, j);
		kmer = ((kmer << 2) | b) & kmask;
		if(i < zsize) continue;
		krev = dna_rev_seq(kmer, zsize);
		if(krev == kmer) continue;
		dir  = krev > kmer? 0 : 1;
		krev = krev > kmer? kmer : krev;
		m = next_ref_hzmv(seeds);
		m->mer   = krev;
		m->dir   = dir;
		m->off   = hzoff->buffer[i - zsize];
		m->len   = (j + 1 - m->off > HZM_MAX_SEED_LEN)? HZM_MAX_SEED_LEN : j + 1 - m->off;
	}
	sort_array(seeds->buffer, seeds->size, hzm_t, (a.mer > b.mer)? 1 : ((a.mer == b.mer)? (a.off > b.off) : 0));
	push_hzmv(seeds, (hzm_t){0xFFFFFFFFLLU, 0, 0, 0});
	kmer = 0;
	kcnt = 0;
	for(idx=off=0;idx<seeds->size;idx++){
		if(seeds->buffer[idx].mer != kmer){
			if(kcnt && kcnt < max_kcnt){ push_hzmhv(hash, (hzmh_t){kmer, off, kcnt, 0}); one_bitvec(bits, kmer); }
			kmer = seeds->buffer[idx].mer;
			kcnt = 1;
			off = idx;
		} else kcnt ++;
	}
	push_hzmhv(hash, (hzmh_t){0xFFFFFFFFLLU, off, 0, 0});
	index_bitvec(bits);
}

static inline void query_single_read_seeds_by_region(uint8_t *pbseq, uint32_t pblen, uint32_t zsize, int hz, uint32_t max_kmer_var, hzmhv *hash, BitVec *bits, hzmv *seeds, u32list *hzoff, hzmpv *rs, int qb, int qe, int tb, int te){
	hzm_t   *p1, *p2, M;
	hzmh_t  *h;
	uint64_t kmer, krev, kmask;
	uint32_t i, j, k;
	uint8_t b, c, dir;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - zsize) << 1);
	b = 4;
	kmer = 0;
	clear_u32list(hzoff);
	p2 = &M;
	clear_hzmpv(rs);
	for(i=j=0;j<pblen;j++){
		c = pbseq[j];
		if(hz && c == b) continue;
		b = c;
		i ++;
		push_u32list(hzoff, j);
		kmer = ((kmer << 2) | b) & kmask;
		if(i < zsize) continue;
		if((int)j < qb) continue;
		if((int)j >= qe) break;
		krev = dna_rev_seq(kmer, zsize);
		if(krev == kmer) continue;
		dir  = krev > kmer? 0 : 1;
		krev = krev > kmer? kmer : krev;
		if(get_bitvec(bits, krev) == 0) continue;
		M.mer   = krev;
		M.dir   = dir;
		M.off   = hzoff->buffer[i - zsize];
		M.len   = (j + 1 - M.off > HZM_MAX_SEED_LEN)? HZM_MAX_SEED_LEN : j + 1 - M.off;
		h = ref_hzmhv(hash, rank_bitvec(bits, M.mer + 1) - 1);
		for(k=0;k<h->cnt;k++){
			p1 = ref_hzmv(seeds, h->off + k);
			if(p1->off < tb) continue;
			if(p1->off + (int)p1->len > te) break;
			if((p1->len > p2->len? p1->len - p2->len : p2->len - p1->len) > max_kmer_var) continue;
			if(p1->dir ^ p2->dir){
				push_hzmpv(rs, (hzmp_t){p1->dir, p1->off, p2->dir, pblen - (p2->off + p2->len), p1->len, p2->len, 0});
			} else {
				push_hzmpv(rs, (hzmp_t){p1->dir, p1->off, p2->dir, p2->off, p1->len, p2->len, 0});
			}
		}
	}
	sort_array(rs->buffer, rs->size, hzmp_t, ((((int64_t)a.off1) << 32) | a.off2) > ((((int64_t)b.off1) << 32) | b.off2));
}

static inline void query_single_read_seeds(uint8_t *pbseq, uint32_t pblen, uint32_t zsize, int hz, uint32_t max_kmer_var, hzmhv *hash, BitVec *bits, hzmv *seeds, u32list *hzoff, hzmpv *rs){
	hzm_t   *p1, *p2, M;
	hzmh_t  *h;
	uint64_t kmer, krev, kmask;
	uint32_t i, j, k;
	uint8_t b, c, dir;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - zsize) << 1);
	b = 4;
	kmer = 0;
	clear_u32list(hzoff);
	p2 = &M;
	clear_hzmpv(rs);
	for(i=j=0;j<pblen;j++){
		c = pbseq[j];
		if(hz && c == b) continue;
		b = c;
		i ++;
		push_u32list(hzoff, j);
		kmer = ((kmer << 2) | b) & kmask;
		if(i < zsize) continue;
		krev = dna_rev_seq(kmer, zsize);
		if(krev == kmer) continue;
		dir  = krev > kmer? 0 : 1;
		krev = krev > kmer? kmer : krev;
		if(get_bitvec(bits, krev) == 0) continue;
		M.mer   = krev;
		M.dir   = dir;
		M.off   = hzoff->buffer[i - zsize];
		M.len   = (j + 1 - M.off > HZM_MAX_SEED_LEN)? HZM_MAX_SEED_LEN : j + 1 - M.off;
		h = ref_hzmhv(hash, rank_bitvec(bits, M.mer + 1) - 1);
		for(k=0;k<h->cnt;k++){
			p1 = ref_hzmv(seeds, h->off + k);
			if((p1->len > p2->len? p1->len - p2->len : p2->len - p1->len) > max_kmer_var) continue;
			//push_hzmpv(rs, (hzmp_t){p1->dir, p1->off, p2->dir, p2->off, p1->len, p2->len});
			if(p1->dir ^ p2->dir){
				push_hzmpv(rs, (hzmp_t){p1->dir, p1->off, p2->dir, pblen - (p2->off + p2->len), p1->len, p2->len, 0});
			} else {
				push_hzmpv(rs, (hzmp_t){p1->dir, p1->off, p2->dir, p2->off, p1->len, p2->len, 0});
			}
		}
	}
	sort_array(rs->buffer, rs->size, hzmp_t, ((((int64_t)a.off1) << 32) | a.off2) > ((((int64_t)b.off1) << 32) | b.off2));
}

static inline void check_kswx_cigar(kswx_t x, u32list *cigar, uint8_t *pb1, uint8_t *pb2){
	uint32_t i, j, op, len;
	int x1, x2, mat, mis, ins, del, aln;
	mat = mis = ins = del = aln = 0;
	x1 = x.tb;
	x2 = x.qb;
	for(i=0;i<cigar->size;i++){
		len = cigar->buffer[i] >> 4;
		op  = cigar->buffer[i] & 0xF;
		aln += len;
		switch(op){
			case 1: ins += len; x2 += len; break;
			case 2: del += len; x1 += len; break;
			default:
				for(j=0;j<len;j++){
					if(pb1[x1 + j] == pb2[x2 + j]) mat ++;
					else mis ++;
				}
				x1 += len;
				x2 += len;
		}
	}
	if(x1 != x.te){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(x2 != x.qe){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(aln != x.aln){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(mat != x.mat){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(mis != x.mis){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(ins != x.ins){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(del != x.del){
		fprintf(hzm_debug_out, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
}

static inline kswx_t hz_align_hzmo(uint8_t *pb1, uint32_t len1, uint8_t *pb2, uint32_t len2, int M, int I, int D, int E, u32list *cigars){
	kswx_t x;
	uint32_t s[2], e[2], l[2];
	s[0] = s[1] = 0;
	x = KSWX_NULL;
	while(s[0] < len1 || s[1] < len2){
		if(pb1[s[0]] != pb2[s[1]]) return KSWX_NULL;
		e[0] = s[0] + 1; while(e[0] < len1 && pb1[e[0]] == pb1[s[0]]) e[0] ++;
		e[1] = s[1] + 1; while(e[1] < len2 && pb2[e[1]] == pb2[s[1]]) e[1] ++;
		l[0] = e[0] - s[0];
		l[1] = e[1] - s[1];
		if(l[0] < l[1]){
			x.aln += l[1];
			x.mat += l[0];
			x.ins += l[1] - l[0];
			x.score += l[0] * M + I + (l[1] - l[0]) * E;
			kswx_push_cigar(cigars, 0, l[0]);
			kswx_push_cigar(cigars, 1, l[1] - l[0]);
		} else if(l[0] == l[1]){
			x.aln += l[0];
			x.mat += l[0];
			x.score += l[0] * M;
			kswx_push_cigar(cigars, 0, l[0]);
		} else {
			x.aln += l[0];
			x.mat += l[1];
			x.del += l[0] - l[1];
			x.score += l[1] * M + D + (l[0] - l[1]) * E;
			kswx_push_cigar(cigars, 0, l[1]);
			kswx_push_cigar(cigars, 2, l[0] - l[1]);
		}
		s[0] = e[0]; s[1] = e[1];
	}
	x.te = x.mat + x.del;
	x.qe = x.mat + x.ins;
	return x;
}

static inline int32_t calculate_median_value(int32_t *rs, int32_t size){
	int32_t i, j, key, mid, beg, end, tmp;
	if(size == 0) return 0;
	beg = 0;
	end = size - 1;
	while(beg < end){
		mid = beg + (end - beg) / 2;
		if(rs[beg] > rs[mid]){ tmp = rs[beg]; rs[beg] = rs[mid]; rs[mid] = tmp; }
		if(rs[mid] > rs[end]){
			tmp = rs[end]; rs[end] = rs[mid]; rs[mid] = tmp;
			if(rs[beg] > rs[mid]){ tmp = rs[beg]; rs[beg] = rs[mid]; rs[mid] = tmp; }
		}
		key = rs[mid];
		i = beg + 1; j = end - 1;
		while(1){
			while(key > rs[i]) i ++;
			while(rs[j] > key) j --;
			if(i < j){
				tmp = rs[i]; rs[i] = rs[j]; rs[j] = tmp;
				i ++; j --;
			} else break;
		}
		if(i == j){ i ++; j --; }
		if(i <= size / 2) beg = i;
		else end = j;
	}
	return rs[size/2];
}

#define hzmo_roundup8x(n) (((n) + 0x7LLU) & 0xFFFFFFFFFFFFFFF8LLU)

#define KWIN_MAX_HOMOLOG	10

#define KWIN_MAX_OFFSET_DEV	50

int chaining_hzmps(hzmpv *regs, uint32_t *trans, uint32_t beg, uint32_t end, hzmpv *rs, int W, u8list *mem){
	typedef struct { int weight, bt; } node_t;
	node_t *nodes;
	hzmp_t *r1, *r2;
	uint32_t i, j, size;
	int mw, bt, band, weight, ol, lst;
	float band_penalty = 0.5;
	clear_and_encap_u8list(mem, kswx_roundup8x(sizeof(node_t) * (end - beg)));
	nodes = (node_t*)mem->buffer;
	for(i=beg;i<end;i++){
		nodes[i-beg].weight = regs->buffer[trans[i]].len2;
		nodes[i-beg].bt = -1;
	}
	mw = -1000000;
	bt = -1;
	for(i=beg;i<end;i++){
		r1 = ref_hzmpv(regs, trans[i]);
		//nodes[i-beg].weight += num_min(r1->len1, r1->len2);
		if(nodes[i-beg].weight > mw){ mw = nodes[i-beg].weight; bt = i; }
		for(j=i+1;j<end;j++){
			r2 = ref_hzmpv(regs, trans[j]);
			if(r2->off2 > r1->off2 + W) break;
			//if(r2->off2 < r1->off2 + r1->len2) continue;
			//if(r2->off1 < r1->off1 + r1->len1) continue;
			//band = num_diff(r2->off1 - (r1->off1 + r1->len1), r2->off2 - (r1->off2 + r1->len2));
			if(r2->off2 <= r1->off2) continue;
			if(r2->off1 <= r1->off1) continue;
			band = num_diff(r2->off1 - r1->off1, r2->off2 - r1->off2);
			//if(band > w) continue;
			weight = (r2->off2 + r2->len2) - (r1->off2 + r1->len2);
			if(weight > r2->len2) weight = r2->len2;
			weight = weight + nodes[i-beg].weight - band * band_penalty;
			if(nodes[j-beg].weight < weight){ nodes[j-beg].weight = weight; nodes[j-beg].bt = i; }
		}
	}
	size = rs->size;
	while(bt >= 0){
		r1 = ref_hzmpv(regs, trans[bt]);
		push_hzmpv(rs, *r1);
		bt = nodes[bt-beg].bt;
	}
	reverse_array(rs->buffer + size, rs->size - size, hzmp_t);
	ol = lst = 0;
	for(i=size;i<rs->size;i++){
		r1 = ref_hzmpv(rs, i);
		if(r1->off2 >= lst) ol += r1->len2;
		else ol += r1->off2 + r1->len2 - lst;
		lst = r1->off2 + r1->len2;
	}
	if(hzm_debug > 3 && rs->size - size){
		for(i=size;i<rs->size;i++){
			r1 = ref_hzmpv(rs, i);
			fprintf(hzm_debug_out, "CHAIN[%d]\t%d\t%d\t%d\t%d", i - size, r1->off1, r1->off1 + r1->len1, r1->off2, r1->off2 + r1->len2);
			fprintf(hzm_debug_out, "\t%d\n", ((int)r1->off1) - ((int)r1->off2));
		}
	}
	return ol;
}

static inline uint32_t potential_paired_kmers_windows(hzmpv *rs, int dir, uint32_t beg, uint32_t end, int bound, wtseedv *seeds, hzmpv *anchors, u8list *mem[2], uint32_t zsize, uint32_t kwin, uint32_t zovl){
	typedef struct { uint32_t b, e, ovl; } wreg_t;
	wreg_t *ws;
	hzmp_t *p, *p0, *p1;
	wt_seed_t *seed;
	uint32_t *ts;
	int32_t *as;
	uint32_t i, j, n, n2, ol, ol2, s, t, lst, size, ret;
	int offset, offn, off, offdev;
	offdev = KWIN_MAX_OFFSET_DEV;
	if(hzm_debug > 2){
		fprintf(hzm_debug_out, "SCAN\t%c\t%d\t%d\t%d\n", "+-"[dir], beg, end, bound);
	}
	n = 0;
	while(beg < end){
		p = ref_hzmpv(rs, beg);
		//if(p->closed || p->dir1 ^ p->dir2 ^ dir || p->off1 < bound){ beg ++; }
		if(p->dir1 ^ p->dir2 ^ dir || p->off1 < bound){ beg ++; }
		else break;
	}
	for(i=beg;i<end;i++){
		p = ref_hzmpv(rs, i);
		//if(p->closed) continue;
		if(p->dir1 ^ p->dir2 ^ dir){ continue; }
		n ++;
	}
	if(n * zsize < zovl) return 0;
	clear_u8list(mem[0]);
	encap_u8list(mem[0], hzmo_roundup8x(n * sizeof(uint32_t)) * 2 + hzmo_roundup8x(n * sizeof(wreg_t)));
	ts   = (uint32_t*)mem[0]->buffer;
	as   =  (int32_t*)(((void*)ts) + hzmo_roundup8x(n * sizeof(uint32_t)));
	ws   =   (wreg_t*)(((void*)as) + hzmo_roundup8x(n * sizeof(int32_t)));
	n = 0;
	for(i=beg;i<end;i++){
		p = ref_hzmpv(rs, i);
		//if(p->closed) continue;
		if(p->dir1 ^ p->dir2 ^ dir){ continue; }
		ts[n++] = i;
	}
	sort_array(ts, n, uint32_t, (rs->buffer[a].off2) > (rs->buffer[b].off2));
	ol = 0; lst = 0; n2 = 0;
	for(i=j=0;i<n;i++){
		p  = ref_hzmpv(rs, ts[i]);
		while((uint32_t)p->off2 + p->len2 > rs->buffer[ts[j]].off2 + kwin){
			p0 = ref_hzmpv(rs, ts[j++]);
			p1 = ref_hzmpv(rs, ts[j]);
			s = p1->off2;
			t = p0->off2 + p0->len2;
			ol2 = s < t? t - s : 0;
			ol = ol + ol2 - p0->len2;
		}
		ol += (p->off2 > lst)? p->len2 : p->off2 + p->len2 - lst;
		lst = p->off2 + p->len2;
		if(ol >= zovl){
			if(n2 &&
				( rs->buffer[ts[i]].off2 <= (uint32_t)(rs->buffer[ts[ws[n2-1].e]].off2) + kwin / 3 ||
				  rs->buffer[ts[j]].off2 <= (uint32_t)(rs->buffer[ts[ws[n2-1].b]].off2) + kwin / 3 ) ){
				if(ol > ws[n2-1].ovl){
					ws[n2-1].b = j; ws[n2-1].e = i; ws[n2-1].ovl = ol;
					if(hzm_debug > 2){
						fprintf(hzm_debug_out, "WSU\t%d\t%d\t%d\t%d\t%d\n", j, i, rs->buffer[ts[j]].off2, rs->buffer[ts[i]].off2 + rs->buffer[ts[i]].len2, ol);
					}
				}
			} else {
				if(hzm_debug > 2){
					fprintf(hzm_debug_out, "WSA\t%d\t%d\t%d\t%d\t%d\n", j, i, rs->buffer[ts[j]].off2, rs->buffer[ts[i]].off2 + rs->buffer[ts[i]].len2, ol);
				}
				ws[n2].b = j; ws[n2].e = i; ws[n2].ovl = ol;
				n2 ++;
			}
		}
	}
	ret = 0;
	for(i=0;i<n2;i++){
		size = anchors->size;
		if(hzm_debug > 2){
			fprintf(hzm_debug_out, "pair\t%d\t%d\t%d\n", rs->buffer[ts[ws[i].b]].off2, rs->buffer[ts[ws[i].e]].off2 + rs->buffer[ts[ws[i].e]].len2, ws[i].ovl);
		}
		if(HZM_FAST_WINDOW_KMER_CHAINING){
			offset = 0;
			offn = 0;
			for(j=ws[i].b;j<=ws[i].e;j++){
				p = ref_hzmpv(rs, ts[j]);
				off = ((int)p->off1) - ((int)p->off2);
				offset += off;
				as[offn] = off;
				offn ++;
				if(hzm_debug > 3){
					fprintf(hzm_debug_out, "LINEAR\t%d\t%d\t%d\t%d\t%c\t%d\n", p->off1, p->len1, p->off2, p->len2, "+-"[p->dir1^p->dir2], (int)p->off1 - (int)p->off2);
				}
			}
			offset = calculate_median_value(as, offn);
			//offset = offset / offn;
			if(hzm_debug > 2){
				fprintf(hzm_debug_out, "offset = %d\n", offset);
			}
			ol = lst = 0;
			for(j=ws[i].b;j<=ws[i].e;j++){
				p = ref_hzmpv(rs, ts[j]);
				off = ((int)p->off1) - ((int)p->off2);
				if(off < offset - offdev || off > offset + offdev) continue;
				push_hzmpv(anchors, *p);
				ol += (p->off2 > lst)? p->len2 : p->off2 + p->len2 - lst;
				lst = p->off2 + p->len2;
			}
			if(hzm_debug > 2){
				fprintf(hzm_debug_out, "ovl2 = %d\n", ol);
			}
			if(anchors->size == size){ continue; }
			sort_array(anchors->buffer + size, anchors->size - size, hzmp_t, a.off1 > b.off1);
			{
				seed = next_ref_wtseedv(seeds);
				seed->closed = 0;
				seed->dir = dir;
				seed->anchors[0] = size;
				seed->beg[0] = 0x7FFFFFFF;
				seed->beg[1] = 0x7FFFFFFF;
				seed->end[0] = 0;
				seed->end[1] = 0;
				seed->ovl = ol;
			}
			ol = lst = 0;
			for(j=size;j<anchors->size;j++){
				p = ref_hzmpv(anchors, j);
				ol += (p->off1 > lst)? p->len1 : p->off1 + p->len1 - lst;
				lst = p->off1 + p->len1;
				if(p->off1 < seed->beg[0]) seed->beg[0] = p->off1;
				if(p->off1 + p->len1 > seed->end[0]) seed->end[0] = p->off1 + p->len1;
				if(p->off2 < seed->beg[1]) seed->beg[1] = p->off2;
				if(p->off2 + p->len2 > seed->end[1]) seed->end[1] = p->off2 + p->len2;
			}
			if(hzm_debug > 2){
				fprintf(hzm_debug_out, "ovl1 = %d\n", ol);
			}
		} else {
			ol = chaining_hzmps(rs, ts, ws[i].b, ws[i].e, anchors, num_max(kwin / 8, 64), mem[1]);
			if(anchors->size == size){ continue; }
			{
				seed = next_ref_wtseedv(seeds);
				seed->closed = 0;
				seed->dir = dir;
				seed->anchors[0] = size;
				p = ref_hzmpv(anchors, size);
				seed->beg[0] = p->off1;
				seed->beg[1] = p->off2;
				offset = ((int)p->off1) - ((int)p->off2);
				p = ref_hzmpv(anchors, anchors->size - 1);
				seed->end[0] = p->off1 + p->len1;
				seed->end[1] = p->off2 + p->len2;
				seed->ovl = ol;
			}
		}
		if(ol * 2 < zovl){
			anchors->size = size;
			seeds->size --;
		} else if(ret && (seed->end[1] <= (int)(seeds->buffer[seeds->size-2].end[1] + kwin / 3) && ol <= seeds->buffer[seeds->size-2].ovl)){
			anchors->size = size;
			seeds->size --;
		} else {
			ret ++;
			seed->ovl = ol;
			seed->anchors[1] = anchors->size;
			if(hzm_debug){
				fprintf(hzm_debug_out, "WINDOW\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", "+-"[dir], seed->beg[0], seed->end[0], seed->beg[1], seed->end[1], ol, offset);
			}
		}
	}
	return ret;
}

static inline uint32_t merge_paired_kmers_window(hzmpv *rs, int dir, wtseedv *seeds, hzmpv *anchors, u8list *mem[2], uint32_t zsize, uint32_t kwin, uint32_t kstep, uint32_t zovl){
	hzmp_t *p, *p0, *p1, *p2, P;
	wt_seed_t SEED;
	uint32_t i, j, n, a, ol, ol2, lst, wlst, s, t, ret;
	int nxt;
	P.off1 = 0x1FFFFFU; P.len1 = 0x3FFU;
	memset(&SEED, 0, sizeof(wt_seed_t));
	p0 = NULL;
	ol = 0;
	lst = 0;
	wlst = 0;
	ret = 0;
	for(j=0;j<rs->size;j++){
		p0 = ref_hzmpv(rs, j);
		//if(p0->closed) continue;
		if(p0->dir1 ^ p0->dir2 ^ dir){ continue; }
		break;
	}
	if(j == rs->size) return 0;
	p = p0;
	for(i=j;i<=rs->size;i++){
		p2 = p;
		if(i < rs->size){
			p = ref_hzmpv(rs, i);
			//if(p->closed) continue;
			if(p->dir1 ^ p->dir2 ^ dir){ continue; }
		} else p = &P;
		if(p->off1 > p0->off1 + kwin){
			if(ol >= zovl){
				// check whether enough hits in cur window
				if((n = potential_paired_kmers_windows(rs, dir, j, i, wlst, seeds, anchors, mem, zsize, kwin, zovl))){
					for(a=0;a<n;a++){
						if((int)wlst < seeds->buffer[seeds->size + a - n].end[0] + 20) wlst = seeds->buffer[seeds->size + a - n].end[0] + 20;
					}
					ret += n;
					p0 = p;
					ol = p->len1;
					lst = p->off1 + p->len1;
					j = i;
				} else {
					nxt = p0->off1 + kstep;
					// move p0 by kstep
					while(p0->off1 < nxt && j < i){
						p1 = ref_hzmpv(rs, ++j);
						// revise ol, I am sure of that one kmer cannot be contained by other
						s = num_max(p0->off1, p1->off1);
						t = num_min(p0->off1 + p0->len1, p1->off1 + p1->len1);
						ol2 = s < t? t - s : 0;
						ol = ol + ol2 - p0->len1;
						p0 = p1;
					}
				}
			}
			if(p->off1 == P.off1) break;
			// move the start of window
			while(p->off1 > p0->off1 + kwin){
				p1 = ref_hzmpv(rs, ++j);
				// revise ol, I am sure of that one kmer cannot be contained by other
				s = num_max(p0->off1, p1->off1);
				t = num_min(p0->off1 + p0->len1, p1->off1 + p1->len1);
				ol2 = s < t? t - s : 0;
				ol = ol + ol2 - p0->len1;
				p0 = p1;
			}
		} else {
			if(p->off1 >= lst) ol += p->len1;
			else if(p->off1 + p->len1 > (int)lst) ol += p->off1 + p->len1 - lst;
			else continue; // should never happen, we assume that end of p always <= end of (p+1)
			lst = p->off1 + p->len1;
		}
	}
	return ret;
}

static inline int chaining_wtseedv(uint32_t pb1, uint32_t pb2, int dir, wtseedv *regs, uint32_t beg, uint32_t end, u8list *mem, int W){
	typedef struct { int weight, bt; } node_t;
	node_t *nodes;
	wt_seed_t *r1, *r2;
	uint32_t i, j;
	int mw, bt, band;
	int max_overhang = 0;
	float band_penalty = 0.05;
	if(hzm_debug > 2){
		fprintf(hzm_debug_out, "CHAINING\t%u\n", pb2);
	}
	clear_and_encap_u8list(mem, kswx_roundup8x(sizeof(node_t) * (end - beg)));
	nodes = (node_t*)mem->buffer;
	nodes = nodes - beg;
	for(i=beg;i<end;i++){
		nodes[i].weight =  0; //- band_penalty * 1000;;
		nodes[i].bt = - 1;
	}
	mw = -1000000;
	bt = -1;
	for(i=beg;i<end;i++){
		r1 = ref_wtseedv(regs, i);
		r1->closed = 1;
		nodes[i].weight += r1->ovl;
		if(nodes[i].weight > mw){ mw = nodes[i].weight; bt = i; }
		for(j=i+1;j<end;j++){
			r2 = ref_wtseedv(regs, j);
			if(r2->beg[1] + max_overhang < r1->end[1]) continue;
			if(r2->beg[0] + max_overhang < r1->end[0]) continue;
			if(r2->beg[0] - r1->end[0] > W && r2->beg[1] - r1->end[1] > W) break;
			band = num_diff(r2->beg[0] - r1->end[0], r2->beg[1] - r1->end[1]);
			if(band > W) continue;
			band = band * band_penalty;
			if(nodes[j].weight < nodes[i].weight - band){ nodes[j].weight = nodes[i].weight - band; nodes[j].bt = i; }
		}
	}
	mw = 0;
	while(bt >= 0){
		r1 = ref_wtseedv(regs, bt);
		r1->closed = 0;
		mw += r1->end[0] - r1->beg[0];
		bt = nodes[bt].bt;
	}
	if(hzm_debug){
		fprintf(hzm_debug_out, "CHAIN-HZMP\t%u\t%u\t%c\tweight=%d\n", pb1, pb2, "+-"[dir], mw);
		j = 0;
		for(i=beg;i<end;i++){
			r1 = ref_wtseedv(regs, i);
			if(r1->closed) continue;
			j ++;
			fprintf(hzm_debug_out, "CHAIN[%d]\t%d\t%d\t%d\t%d", j, r1->beg[0], r1->end[0], r1->beg[1], r1->end[1]);
			fprintf(hzm_debug_out, "\t%d\n", r1->ovl);
		}
	}
	return mw;
}

typedef struct {
	int offset;
	uint32_t off, cnt;
} diag_t;
define_list(diagv, diag_t);

static inline void denoising_hzmps(hzmpv *rs, uint32_t pb2, hzmpv *dst[2], wtseedv *regs[2], int xvar, int yvar, int min_linear_len, diagv *diags, u4v *block, u4v *grps){
	diag_t *d;
	wt_seed_t *seed;
	hzmp_t *p, *p0, P;
	uint32_t i, j, k, doff, dcnt, gid;
	int dir, len, mat, lst;
	int lst_offset, end_offset;
	sort_array(rs->buffer, rs->size, hzmp_t, ((((int64_t)a.off1 - (int64_t)a.off2) << 32) | a.off1) > ((((int64_t)b.off1 - (int64_t)b.off2) << 32) | b.off1));
	memset(&P, 0, sizeof(hzmp_t));
	P.off1 = HZMP_MAX_SEED_OFF;
	//denoising
	for(dir=0;dir<2;dir++){
		clear_diagv(diags);
		clear_hzmpv(dst[dir]);
		clear_wtseedv(regs[dir]);
		d = NULL;
		for(i=0;i<rs->size;i++){
			p = ref_hzmpv(rs, i);
			if(p->dir1 ^ p->dir2 ^ dir) continue;
			if(d && d->offset == (int)p->off1 - (int)p->off2){
				d->cnt ++;
			} else {
				d = next_ref_diagv(diags);
				d->offset = (int)p->off1 - (int)p->off2;
				d->off = i;
				d->cnt = 1;
			}
		}
		doff = 0;
		end_offset = - 0x7FFFFFFF;
		clear_u4v(grps);
		push_u4v(grps, 0);
		while(doff < rs->size){
			// find yvar offsets region
			lst_offset = diags->buffer[doff].offset;
			dcnt = 0;
			while(1){
				if(diags->buffer[dcnt+doff].offset > lst_offset + yvar) break;
				if(dcnt + doff + 1 >= diags->size) break;
				dcnt ++;
			}
			if(dcnt == 0) break;
			if(diags->buffer[doff + dcnt].offset == end_offset){
				doff += dcnt;
				continue;
			}
			end_offset = diags->buffer[doff + dcnt].offset;
			// sort local yvar offsets region by p->off1
			clear_u4v(block);
			for(i=0;i<dcnt;i++){
				d = ref_diagv(diags, i + doff);
				for(j=0;j<d->cnt;j++){
					p = ref_hzmpv(rs, d->off + j);
					if(p->dir1 ^ p->dir2 ^ dir) continue;
					push_u4v(block, d->off + j);
				}
			}
			sort_array(block->buffer, block->size, uint32_t, rs->buffer[a].off1 > rs->buffer[b].off1);
			// scan linear hzmps
			if(block->size){
				p0 = block->size? ref_hzmpv(rs, block->buffer[0]) : NULL;
				len = mat = p0->len1;
			} else {
				len = mat = 0;
				p0 = NULL;
			}
			j = 0;
			for(i=1;i<=block->size;i++){
				p = (i == block->size)? &P : ref_hzmpv(rs, block->buffer[i]);
				if(p->off1 <= p0->off1 + p0->len1){
					mat += p->off1 + p->len1 - (p0->off1 + p0->len1);
					len += p->off1 + p->len1 - (p0->off1 + p0->len1);
				} else if(p->off1 <= p0->off1 + p0->len1 + xvar){
					mat += p->len1;
					len += p->off1 + p->len1 - (p0->off1 + p0->len1);
				} else {
					if(len >= min_linear_len){
						gid = 0;
						for(k=j;k<i;k++){
							if(rs->buffer[block->buffer[k]].gid){
								if(gid == 0){
									gid = grps->buffer[rs->buffer[block->buffer[k]].gid];
								} else {
									grps->buffer[rs->buffer[block->buffer[k]].gid] = gid;
								}
							}
						}
						if(gid == 0){
							gid = grps->size;
							push_u4v(grps, gid);
						}
						for(;j<i;j++){ ref_hzmpv(rs, block->buffer[j])->gid = gid; }
					}
					j = i;
					p0 = p;
					len = mat = p0->len1;
				}
			}
			// move offset to next
			for(i=doff;i<doff+dcnt;i++){
				if(diags->buffer[i].offset > lst_offset + yvar / 2) break;
			}
			doff = i;
		}
		// tidy gid map
		for(i=1;i<grps->size;i++){
			if(grps->buffer[i] < i) continue;
			for(j=i+1;j<grps->size;j++){
				if(grps->buffer[j] != i) continue;
				for(k=j+1;k<grps->size;k++){
					if(grps->buffer[k] == j){
						grps->buffer[k] = i;
					}
				}
			}
		}
		for(i=0;i<rs->size;i++){
			p = ref_hzmpv(rs, i);
			if(p->dir1 ^ p->dir2 ^ dir) continue;
			if(p->gid == 0) continue;
			// re-map gid
			p->gid = grps->buffer[p->gid];
			push_hzmpv(dst[dir], *p);
		}
		sort_array(dst[dir]->buffer, dst[dir]->size, hzmp_t, num_cmpgtx(a.gid, b.gid, a.off1, b.off1));
		// generate hzmp block: wt_seed_t
		j = 0;
		for(i=1;i<=dst[dir]->size;i++){
			if(i < dst[dir]->size && dst[dir]->buffer[i].gid == dst[dir]->buffer[j].gid) continue;
			seed = next_ref_wtseedv(regs[dir]);
			seed->pb2 = pb2;
			seed->closed = 0;
			seed->dir = dir;
			seed->anchors[0] = j;
			seed->anchors[1] = i;
			seed->beg[0] = 0x7FFFFFFF;
			seed->beg[1] = 0x7FFFFFFF;
			seed->end[0] = 0;
			seed->end[1] = 0;
			seed->ovl = 0;
			lst = 0;
			for(k=j;k<i;k++){
				p = ref_hzmpv(dst[dir], k);
				if(p->off1 < seed->beg[0]) seed->beg[0] = p->off1;
				if(p->off1 + p->len1 > seed->end[0]) seed->end[0] = p->off1 + p->len1;
				if(p->off2 < seed->beg[1]) seed->beg[1] = p->off2;
				if(p->off2 + p->len2 > seed->end[1]) seed->end[1] = p->off2 + p->len2;
				seed->ovl += (p->off1 > lst)? p->len1 : p->off1 + p->len1 - lst;
				lst = p->off1 + p->len1;
			}
			j = i;
		}
		sort_array(regs[dir]->buffer, regs[dir]->size, wt_seed_t, a.beg[0] > b.beg[0]);
	}
}

static inline void dot_plot_hzmps(hzmpv *rs, int dir, int pass, FILE *out){
	hzmp_t *p;
	uint32_t i;
	for(i=0;i<rs->size;i++){
		p = ref_hzmpv(rs, i);
		if(p->dir1 ^ p->dir2 ^ dir) continue;
		if((p->gid > 0) < pass) continue;
		fprintf(out, "%d\t%d\t%d\n", p->off1, p->off2, p->gid);
	}
}

static inline void debug_dot_plot_hzmps(hzmpv *rs){
	{
		FILE *out;
		out = fopen("dot_plot.fwd.src.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 0, 0, out);
		fclose(out);
	}
	{
		FILE *out;
		out = fopen("dot_plot.fwd.dst.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 0, 1, out);
		fclose(out);
	}
	{
		FILE *out;
		out = fopen("dot_plot.rev.src.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 1, 0, out);
		fclose(out);
	}
	{
		FILE *out;
		out = fopen("dot_plot.rev.dst.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 1, 1, out);
		fclose(out);
	}
}

static inline int chaining_overhang_wtseedv(uint32_t pb1, uint32_t pb2, int dir, wtseedv *regs, uint32_t beg, uint32_t end, u8list *mem, int max_overhang, float band_penalty, float gap_penalty){
	typedef struct { int weight, bt; } node_t;
	node_t *nodes;
	wt_seed_t *r1, *r2;
	uint32_t i, j;
	int mw, bt, band, gap, W, score;
	if(hzm_debug > 2){
		fprintf(hzm_debug_out, "CHAINING\t%u\n", pb2);
	}
	clear_and_encap_u8list(mem, kswx_roundup8x(sizeof(node_t) * (end - beg)));
	nodes = (node_t*)mem->buffer;
	nodes = nodes - beg;
	for(i=beg;i<end;i++){
		nodes[i].weight =  0; //- band_penalty * 1000;;
		nodes[i].bt = - 1;
	}
	mw = -1000000;
	bt = -1;
	for(i=beg;i<end;i++){
		r1 = ref_wtseedv(regs, i);
		r1->closed = 1;
		nodes[i].weight += r1->ovl;
		if(nodes[i].weight > mw){ mw = nodes[i].weight; bt = i; }
		W = nodes[i].weight / gap_penalty;
		//fprintf(hzm_debug_out, "beg[0]=%d\ti=%d\tW=%d\n", r1->beg[0], i, W);
		for(j=i+1;j<end;j++){
			r2 = ref_wtseedv(regs, j);
			if(r2->beg[1] + max_overhang < r1->end[1]) continue;
			if(r2->beg[0] + max_overhang < r1->end[0]) continue;
			if(r2->beg[0] - r1->end[0] > W) break;
			band = num_diff(r2->beg[0] - r1->end[0], r2->beg[1] - r1->end[1]);
			gap  = num_max(r2->beg[0] - r1->end[0], r2->beg[1] - r1->end[1]);
			if(gap < 0) gap = 0;
			score = band * band_penalty + gap * gap_penalty;
			//fprintf(hzm_debug_out, "beg[0]=%d\ti=%d\tj=%d\t%d\t%d\t%d\t%d\t%f\t%f\n", r1->beg[0], i, j, band, gap, score, nodes[i].weight - score, band_penalty, gap_penalty);
			if(nodes[j].weight < nodes[i].weight - score){
				nodes[j].weight = nodes[i].weight - score;
				nodes[j].bt = i;
			}
			if(score < nodes[i].weight && 2 * band * band_penalty <= regs->buffer[j].ovl) break; // j can replace i in further chaining, heuristically
		}
	}
	mw = 0;
	while(bt >= 0){
		r1 = ref_wtseedv(regs, bt);
		r1->closed = 0;
		//mw += r1->end[0] - r1->beg[0];
		mw += r1->ovl;
		bt = nodes[bt].bt;
	}
	if(0 && hzm_debug){
		fprintf(hzm_debug_out, "CHAIN-HZMP\t%u\t%u\t%c\tweight=%d\n", pb1, pb2, "+-"[dir], mw);
		j = 0;
		for(i=beg;i<end;i++){
			r1 = ref_wtseedv(regs, i);
			if(r1->closed) continue;
			j ++;
			fprintf(hzm_debug_out, "CHAIN[%d]\t%d\t%d\t%d\t%d", j, r1->beg[0], r1->end[0], r1->beg[1], r1->end[1]);
			fprintf(hzm_debug_out, "\t%d\n", r1->ovl);
		}
	}
	return mw;
}

static inline kswr_t dot_matrix_align_hzmps(hzmpv *rs, hzmpv *dst[2], wtseedv *regs[2], diagv *diags, u4v *block, u4v *grps, u8list *mem, int xvar, int yvar, int min_block_len, int max_overhang, float deviation_penalty, float gap_penalty){
	wt_seed_t *seed;
	uint32_t i, d;
	int weight[2];
	kswr_t ret;
	denoising_hzmps(rs, 0, dst, regs, xvar, yvar, min_block_len, diags, block, grps);
	weight[0] = chaining_overhang_wtseedv(0, 0, 0, regs[0], 0, regs[0]->size, mem, max_overhang, deviation_penalty, gap_penalty);
	weight[1] = chaining_overhang_wtseedv(0, 0, 1, regs[1], 0, regs[1]->size, mem, max_overhang, deviation_penalty, gap_penalty);
	d = (weight[0] < weight[1]);
	ret.score = weight[d];
	ret.qb = ret.tb = 0x7FFFFFFF;
	ret.qe = ret.te = 0;
	for(i=0;i<regs[d]->size;i++){
		seed = ref_wtseedv(regs[d], i);
		if(seed->closed == 0){
			if(ret.qb > seed->beg[1]) ret.qb = seed->beg[1];
			if(ret.tb > seed->beg[0]) ret.tb = seed->beg[0];
			if(ret.qe < seed->end[1]) ret.qe = seed->end[1];
			if(ret.te < seed->end[0]) ret.te = seed->end[0];
		}
	}
	ret.score2 = d;
	ret.te2 = -1;
	if(hzm_debug){
		int qb, qe, tb, te, aln;
		for(d=0;d<2;d++){
			fprintf(hzm_debug_out, "\nBLOCKS\t%c\t%d\t%d\n", "+-"[d], (int)regs[d]->size, weight[d]);
			qb = tb = 0x7FFFFFFF;
			qe = te = 0;
			for(i=0;i<regs[d]->size;i++){
				seed = ref_wtseedv(regs[d], i);
				if(seed->closed == 0){
					if(qb > seed->beg[0]) qb = seed->beg[0];
					if(tb > seed->beg[1]) tb = seed->beg[1];
					if(qe < seed->end[0]) qe = seed->end[0];
					if(te < seed->end[1]) te = seed->end[1];
				}
				fprintf(hzm_debug_out, "WINDOW\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\n", "+-"[d], seed->beg[0], seed->end[0], seed->beg[1], seed->end[1], seed->ovl, seed->end[0] - seed->beg[0], seed->beg[0] - seed->beg[1], "*x--"[seed->closed]);
			}
			if(weight[d] == 0) continue;
			aln = num_max(qe - qb, te - tb);
			fprintf(hzm_debug_out, "OVL\t%c\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%0.3f\n", '+', qb, qe, "+-"[d], tb, te, weight[d], aln, 1.0 * weight[d] / aln);
		}
	}
	return ret;
}

// empty function
static inline void process_hzmps(hzmpv *rs){
	if(rs->size) return;
}

static inline void filter_by_region_hzmps(hzmpv *dst, hzmpv *src, int dir, uint32_t beg1, uint32_t end1, uint32_t beg2, uint32_t end2){
	hzmp_t *p;
	uint32_t i;
	for(i=0;i<src->size;i++){
		p = ref_hzmpv(src, i);
		if(p->dir1 ^ p->dir2 ^ dir) continue;
		if((p->off1 < beg1 || p->off1 + p->len1 > (int)end1) || (p->off2 < beg2 || p->off2 + p->len2 > (int)end2)) continue;
		push_hzmpv(dst, *p);
	}
}

static inline uint32_t rough_scoring_hzmps(hzmpv *rs, uint32_t beg, uint32_t end){
	hzmp_t *p;
	uint32_t i, ol, l;
	ol = 0;
	l = 0;
	for(i=beg;i<end;i++){
		p = ref_hzmpv(rs, i);
		//if(p->closed) continue;
		if(p->off1 >= l) ol += p->len1;
		else if(p->off1 + p->len1 > (int)l) ol += p->off1 + p->len1 - l;
		else continue;
		l = p->off1 + p->len1;
	}
	return ol;
}

static inline uint32_t estimate_overlap_hzmps(hzmpv *rs, uint32_t beg, uint32_t end, u8list *mem){
	uint64_t *ts;
	hzmp_t *p;
	uint32_t i, n, ol, b, e, l;
	clear_u8list(mem);
	encap_u8list(mem, hzmo_roundup8x(rs->size * sizeof(uint64_t)));
	ts = (uint64_t*)mem->buffer;
	n = 0;
	for(i=beg;i<end;i++){
		p = ref_hzmpv(rs, i);
		//if(p->closed) continue;
		ts[n++] = (((uint64_t)p->off1) << 32) | (p->off1 + p->len1);
	}
	if(n == 0) return 0;
	ol = 0; l = 0;
	for(i=0;i<n;i++){
		b = ts[i] >> 32;
		e = ts[i] & 0xFFFFFFFFU;
		if(b >= l) ol += e - b;
		else if(e > l) ol += e - l;
		else continue;
		l = e;
	}
	return ol;
}

typedef struct {
	kswx_t x;
	uint32_t cigar_off, cigar_len;
} aln_reg_t;
define_list(alnregv, aln_reg_t);

static inline kswx_t fast_seeds_align_hzmo(uint8_t *pb1, uint8_t *pb2, wt_seed_t *seed, hzmpv *anchors, u32list *cigar, u8list *mem_cache, u32list *mem_cigar, int w, int M, int X, int I, int D, int E, int T){
	kswx_t x, y;
	hzmp_t *p;
	uint32_t i;
	x = KSWX_NULL;
	if(hzm_debug > 1){
		fprintf(hzm_debug_out, "KWIN\t%d\t%d\t%d\t%d\t%d\n", seed->beg[0], seed->end[0], seed->beg[1], seed->end[1], seed->ovl);
	}
	for(i=seed->anchors[0];i<seed->anchors[1];i++){
		p = ref_hzmpv(anchors, i);
		//if(p->closed) continue;
		if(x.aln == 0){
			x.tb = x.te = p->off1;
			x.qb = x.qe = p->off2;
		}
		if(p->off1 < x.te) continue;
		if(p->off2 < x.qe) continue;
		if(hzm_debug > 1){
			fprintf(hzm_debug_out, "%s:%d\tkswx_t[%d, %d->%d, %d->%d] => [%d+%d, %d+%d]\n", __FUNCTION__, __LINE__, x.score, x.tb, x.te, x.qb, x.qe, p->off1, p->len1, p->off2, p->len2);
		}
		clear_u32list(mem_cigar);
		y = kswx_extend_align_core(p->off2 - x.qe, pb2 + x.qe, p->off1 - x.te, pb1 + x.te, 1, x.score, w, M, X, I, D, E, T, mem_cache, mem_cigar);
		x.score = y.score;
		x.aln += y.aln; x.mat += y.mat; x.mis += y.mis; x.ins += y.ins; x.del += y.del;
		x.te  += y.te;
		x.qe  += y.qe;
		if(x.te < p->off1){
			x.del += p->off1 - x.te;
			x.aln += p->off1 - x.te;
			kswx_push_cigar(mem_cigar, 2, p->off1 - x.te);
			x.te = p->off1;
		}
		if(x.qe < p->off2){
			x.ins += p->off2 - x.qe;
			x.aln += p->off2 - x.qe;
			kswx_push_cigar(mem_cigar, 1, p->off2 - x.qe);
			x.qe = p->off2;
		}
		kswx_push_cigars(cigar, mem_cigar->buffer, mem_cigar->size);
		clear_u32list(mem_cigar);
		y = hz_align_hzmo(pb1 + p->off1, p->len1, pb2 + p->off2, p->len2, M, I, D, E, mem_cigar);
		if(y.aln == 0){
			if(hzm_debug > 1) fprintf(hzm_debug_out, "%s:%d\tRET\n", __FUNCTION__, __LINE__);
			return x; // should never happen
		}
		x.score += y.score;
		x.aln += y.aln; x.mat += y.mat; x.mis += y.mis; x.ins += y.ins; x.del += y.del;
		x.te  += y.te;
		x.qe  += y.qe;
		kswx_push_cigars(cigar, mem_cigar->buffer, mem_cigar->size);
	}
	if(hzm_debug > 1){
		fprintf(hzm_debug_out, "LOCAL\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", x.tb, x.te, x.qb, x.qe, x.score, x.mat, x.mis, x.ins, x.del); fflush(hzm_debug_out);
	}
	return x;
}

static inline void chaining_regs_hzmo(alnregv *regs, alnregv *rs, u8list *mem){
	typedef struct { int weight, bt; } node_t;
	node_t *nodes;
	aln_reg_t *r1, *r2;
	uint32_t i, j, size;
	int mw, bt;
	clear_and_encap_u8list(mem, kswx_roundup8x(sizeof(node_t) * regs->size));
	nodes = (node_t*)mem->buffer;
	for(i=0;i<regs->size;i++){
		nodes[i].weight = 0;
		nodes[i].bt = -1;
	}
	mw = -1;
	bt = -1;
	for(i=0;i<regs->size;i++){
		r1 = ref_alnregv(regs, i);
		nodes[i].weight += r1->x.score;
		if(nodes[i].weight > mw){ mw = nodes[i].weight; bt = i; }
		for(j=i+1;j<regs->size;j++){
			r2 = ref_alnregv(regs, j);
			if(r2->x.qb < r1->x.qe) continue;
			if(nodes[j].weight < nodes[i].weight){ nodes[j].weight = nodes[i].weight; nodes[j].bt = i; }
		}
	}
	size = rs->size;
	while(bt >= 0){
		r1 = ref_alnregv(regs, bt);
		push_alnregv(rs, *r1);
		bt = nodes[bt].bt;
	}
	reverse_array(rs->buffer + size, rs->size - size, aln_reg_t);
	if(hzm_debug){
		fprintf(hzm_debug_out, "CHAIN\tweight=%d\tn_hsp=%d\n", mw, (int)rs->size);
		for(i=0;i<rs->size;i++){
			r1 = ref_alnregv(rs, i);
			fprintf(hzm_debug_out, "CHAIN[%d]\t%d\t%d\t%d\t%d", i, r1->x.tb, r1->x.te, r1->x.qb, r1->x.qe);
			fprintf(hzm_debug_out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", r1->x.score, 1.0 * r1->x.mat / (r1->x.aln + 1), r1->x.mat, r1->x.mis, r1->x.ins, r1->x.del);
		}
	}
}

static inline kswx_t global_align_regs_hzmo(int len1, int len2, alnregv *regs, int esti_regs[2], uint8_t *pbseqs[2], uint32_t *cigar_cache, u32list *cigar, u8list *mem_cache, u32list *mem_cigar, int W, int ew, int _w, int M, int X, int I, int D, int E, int T){
	kswx_t x, y;
	aln_reg_t *reg1, *reg2;
	uint32_t i, *_cigar;
	uint8_t *q, *t;
	int8_t matrix[16];
	int n_cigar, idx, j, x1, x2, op, len, w, max_gap, init_score, score;
	x = KSWX_NULL;
	clear_u32list(cigar);
	if(regs->size == 0) return x;
	for(i=0;i<16;i++) matrix[i] = ((i % 4) == (i / 4))? M : X;
	init_score = 100 * M;
	reg2 = ref_alnregv(regs, 0);
	x = reg2->x;
	if(x.qb && x.tb){
		// left extending
		w = ew;
		max_gap = ((num_min(x.qb, x.tb) * M + x.score + init_score + (- T)) + (I < D? D : I)) / (- E) + 1;
		if(max_gap < w) max_gap = w;
		while(1){
			clear_u32list(mem_cigar);
			y = kswx_extend_align_shift_core(x.qb, pbseqs[1] + x.qb - 1, x.tb, pbseqs[0] + x.tb - 1, -1, x.score + init_score, - w, M, X, I, D, E, T, mem_cache, mem_cigar);
			if(hzm_debug){
				fprintf(hzm_debug_out, "score=%d w=%d\n", y.score - init_score, w);
			}
			if(y.qe == x.qb || y.te == x.tb) break;
			if(x.tb - y.te <= esti_regs[0]) break;
			if(w >= ew || w >= max_gap) break;
			w <<= 1;
		}
		x.score =  y.score - init_score;
		x.aln   += y.aln;
		x.mat   += y.mat;
		x.mis   += y.mis;
		x.ins   += y.ins;
		x.del   += y.del;
		x.qb    -= y.qe;
		x.tb    -= y.te;
		reverse_u32list(mem_cigar);
		kswx_push_cigars(cigar, mem_cigar->buffer, mem_cigar->size);
		if(hzm_debug){
			fprintf(hzm_debug_out, "GAP0\t%d\t%d\t%d\t%d\t%d\t%d", len1, x.tb, x.te, len2, x.qb, x.qe);
			fprintf(hzm_debug_out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
		}
	}
	kswx_push_cigars(cigar, cigar_cache + reg2->cigar_off, reg2->cigar_len);
	reg1 = reg2;
	for(i=1;i<regs->size;i++){
		reg2 = ref_alnregv(regs, i);
		// fill gap
		q = pbseqs[1] + reg1->x.qe;
		t = pbseqs[0] + reg1->x.te;
		if(hzm_debug){
			fprintf(hzm_debug_out, "TRY\t%d\t%d\t%d\t%d\n", reg1->x.te, reg2->x.tb, reg1->x.qe, reg2->x.qb);
		}
		w = _w;
		{
			while(1){
				if(w < num_diff(reg2->x.qb - reg1->x.qe, reg2->x.tb - reg1->x.te)){
					w <<= 1; continue;
				}
				n_cigar = 0; _cigar = NULL;
				score = ksw_global2(reg2->x.qb - reg1->x.qe, q, reg2->x.tb - reg1->x.te, t, 4, matrix, - I, - E, - D, - E, w, &n_cigar, &_cigar);
				if(hzm_debug){
					fprintf(hzm_debug_out, "score=%d w=%d\n", score, w);
				}
				if(score < 0 && w < W && w < num_max(reg2->x.qb - reg1->x.qe, reg2->x.tb - reg1->x.te)){
					w <<= 1;
					free(_cigar);
				} else {
					break;
				}
			}
		}
		//if(score < 0) break;
		x.score += score;
		x.qe    =  reg2->x.qb;
		x.te    =  reg2->x.tb;
		for(idx=x1=x2=0;idx<n_cigar;idx++){
			op  = _cigar[idx] & 0xF;
			len = _cigar[idx] >> 4;
			x.aln += len;
			switch(op){
				case 0: for(j=0;j<len;j++){ if(q[x1 + j] == t[x2 + j]) x.mat ++; else x.mis ++; } x1 += len; x2 += len; break;
				case 1: x1 += len; x.ins += len; break;
				case 2: x2 += len; x.del += len; break;
			}
		}
		kswx_push_cigars(cigar, _cigar, n_cigar);
		free(_cigar);
		// append reg2
		x.score += reg2->x.score;
		x.aln   += reg2->x.aln;
		x.mat   += reg2->x.mat;
		x.mis   += reg2->x.mis;
		x.ins   += reg2->x.ins;
		x.del   += reg2->x.del;
		x.qe    =  reg2->x.qe;
		x.te    =  reg2->x.te;
		kswx_push_cigars(cigar, cigar_cache + reg2->cigar_off, reg2->cigar_len);
		if(hzm_debug){
			fprintf(hzm_debug_out, "GAP1\t%d\t%d\t%d\t%d\t%d\t%d", len1, x.tb, x.te, len2, x.qb, x.qe);
			fprintf(hzm_debug_out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
		}
		reg1 = reg2;
	}
	if(x.te < len1 && x.qe < len2){
		// right extending
		w = ew;
		max_gap = ((num_min(len2 - x.qe, len1 - x.te) * M + x.score + (- T)) + (I < D? D : I)) / (- E) + 1;
		if(max_gap < w) max_gap = w;
		while(1){
			clear_u32list(mem_cigar);
			y = kswx_extend_align_shift_core(len2 - x.qe, pbseqs[1] + x.qe, len1 - x.te, pbseqs[0] + x.te, 1, x.score, - w, M, X, I, D, E, T, mem_cache, mem_cigar);
			if(hzm_debug){
				fprintf(hzm_debug_out, "score=%d w=%d\n", y.score - x.score, w);
			}
			if(y.qe == len2 - x.qe || y.te == len1 - x.te) break;
			if(x.te + y.te >= esti_regs[1]) break;
			if(w >= ew || w >= max_gap) break;
			w <<= 1;
		}
		x.score =  y.score;
		x.aln   += y.aln;
		x.mat   += y.mat;
		x.mis   += y.mis;
		x.ins   += y.ins;
		x.del   += y.del;
		x.qe    += y.qe;
		x.te    += y.te;
		kswx_push_cigars(cigar, mem_cigar->buffer, mem_cigar->size);
		if(hzm_debug){
			fprintf(hzm_debug_out, "GAP2\t%d\t%d\t%d\t%d\t%d\t%d", len1, x.tb, x.te, len2, x.qb, x.qe);
			fprintf(hzm_debug_out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
		}
	}
	if(hzm_debug){
		fprintf(hzm_debug_out, "ALIGN\t%s\t%c\t%d\t%d\t%d\t%s\t%c\t%d\t%d\t%d", "TARGET", '+', len1, x.tb, x.te, "QUERY", '?', len2, x.qb, x.qe);
		fprintf(hzm_debug_out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
	}
	return x;
}

typedef struct {
	uint32_t zsize, hz, zwin, zstep, zovl, zmax, zvar;
	int w, W, ew, rw;
	int M, X, O, I, D, E, T;
	int QMIS, QDEL, QEXT, QCLP;

	u8list *tseq;
	hzmhv *hash;
	BitVec *bits;
	hzmv *seeds;
	int linked_index;

	kswx_t hit;
	int has_alignment;
	String *alns[3];
	u32list *cigars;

	wtseedv *windows;
	hzmpv *rs, *anchors;
	u8list *mem_cache[2];
	u32list *mem_cigar, *cigar_cache;
	alnregv *regs;
} HZMAux;

static inline HZMAux* init_hzmaux(){
	HZMAux *aux;
	aux = calloc(1, sizeof(HZMAux));
	aux->zsize = 10;
	aux->hz    = 1;
	aux->zmax  = 100;
	aux->zvar  = 2;
	aux->zwin  = 800;
	aux->zstep = 400;
	aux->zovl  = 200;
	aux->w  = 50;
	aux->W  = 3200;
	aux->ew = 800;
	aux->rw = 8;
	aux->M = 2;
	aux->X = -5;
	aux->O = -3;
	aux->I = -2;
	aux->D = -3;
	aux->E = -1;
	aux->T = -100;
	aux->QCLP = -5;
	aux->QMIS = -20;
	aux->QDEL = -15;
	aux->QEXT = -5;
	aux->has_alignment = 1;

	aux->tseq = init_u8list(1024);
	aux->hash = init_hzmhv(1024);
	aux->bits = init_bitvec(1024);
	aux->seeds = init_hzmv(1024);
	aux->linked_index = 0;

	aux->alns[0] = init_string(1024);
	aux->alns[1] = init_string(1024);
	aux->alns[2] = init_string(1024);
	aux->cigars  = init_u32list(1024);

	aux->windows = init_wtseedv(1024);
	aux->rs = init_hzmpv(1024);
	aux->anchors = init_hzmpv(1024);
	aux->mem_cache[0] = init_u8list(1024);
	aux->mem_cache[1] = init_u8list(1024);
	aux->mem_cigar = init_u32list(1024);
	aux->cigar_cache = init_u32list(1024);
	aux->regs = init_alnregv(1024);
	return aux;
}

static inline void copy_param_hzmaux(HZMAux *dst, HZMAux *src){
	dst->zsize = src->zsize;
	dst->hz = src->hz;
	dst->zmax = src->zmax;
	dst->zvar = src->zvar;
	dst->zwin = src->zwin;
	dst->zstep = src->zstep;
	dst->zovl = src->zovl;
	dst->w = src->w;
	dst->W = src->W;
	dst->ew = src->ew;
	dst->rw = src->rw;
	dst->M = src->M;
	dst->X = src->X;
	dst->I = src->I;
	dst->D = src->D;
	dst->E = src->E;
	dst->T = src->T;
	dst->QCLP = src->QCLP;
	dst->QMIS = src->QMIS;
	dst->QDEL = src->QDEL;
	dst->QEXT = src->QEXT;
	dst->has_alignment = src->has_alignment;
}

static inline void copy_index_hzmaux(HZMAux *dst, HZMAux *src){
	dst->zsize = src->zsize;
	dst->hz = src->hz;
	dst->zmax = src->zmax;
	dst->zvar = src->zvar;
	dst->zwin = src->zwin;
	dst->zstep = src->zstep;
	dst->zovl = src->zovl;
	dst->w = src->w;
	dst->W = src->W;
	dst->ew = src->ew;
	dst->rw = src->rw;
	dst->M = src->M;
	dst->X = src->X;
	dst->I = src->I;
	dst->D = src->D;
	dst->E = src->E;
	dst->T = src->T;
	dst->QCLP = src->QCLP;
	dst->QMIS = src->QMIS;
	dst->QDEL = src->QDEL;
	dst->QEXT = src->QEXT;
	dst->has_alignment = src->has_alignment;
	clear_u8list(dst->tseq); append_u8list(dst->tseq, src->tseq);
	clear_hzmhv(dst->hash); append_hzmhv(dst->hash, src->hash);
	clear_hzmv(dst->seeds); append_hzmv(dst->seeds, src->seeds);
	recap_bitvec(dst->bits, src->bits->n_cap); memcpy(dst->bits->bits, src->bits->bits, src->bits->n_cap / 8);
	if(dst->bits->sum_size != src->bits->sum_size) dst->bits->sums = realloc(dst->bits->sums, (src->bits->sum_size * 2 + 1) * sizeof(uint64_t));
	memcpy(dst->bits->sums, src->bits->sums, (src->bits->sum_size * 2 + 1) * sizeof(uint64_t)); dst->bits->sum_size = src->bits->sum_size;
	dst->bits->n_ones = src->bits->n_ones;
	if(dst->bits->hash_size != src->bits->hash_size) dst->bits->hash = realloc(dst->bits->hash, src->bits->hash_size * sizeof(uint64_t));
	memcpy(dst->bits->hash, src->bits->hash, src->bits->hash_size * sizeof(uint64_t)); dst->bits->hash_size = src->bits->hash_size;
	dst->bits->hash_mod = src->bits->hash_mod;
}

static inline void share_index_hzmaux(HZMAux *dst, HZMAux *src){
	if(dst->linked_index == 0){
		free_u8list(dst->tseq);
		free_hzmhv(dst->hash);
		free_bitvec(dst->bits);
		free_hzmv(dst->seeds);
	}
	dst->linked_index = 1;
	dst->tseq  = src->tseq;
	dst->hash  = src->hash;
	dst->bits  = src->bits;
	dst->seeds = src->seeds;
}

static inline void free_hzmaux(HZMAux *aux){
	if(aux->linked_index){
	} else {
		free_u8list(aux->tseq);
		free_hzmhv(aux->hash);
		free_bitvec(aux->bits);
		free_hzmv(aux->seeds);
	}
	free_string(aux->alns[0]);
	free_string(aux->alns[1]);
	free_string(aux->alns[2]);
	free_u32list(aux->cigars);
	free_wtseedv(aux->windows);
	free_hzmpv(aux->rs);
	free_hzmpv(aux->anchors);
	free_u8list(aux->mem_cache[0]);
	free_u8list(aux->mem_cache[1]);
	free_u32list(aux->mem_cigar);
	free_u32list(aux->cigar_cache);
	free_alnregv(aux->regs);
	free(aux);
}

static inline void reset_hzmaux(HZMAux *aux){
	if(aux->linked_index){
		aux->tseq = init_u8list(1024);
		aux->hash = init_hzmhv(1024);
		aux->bits = init_bitvec(1024);
		aux->seeds = init_hzmv(1024);
		aux->linked_index = 0;
	}
	clear_u8list(aux->tseq);
}

static inline void add_tseq_hzmaux(HZMAux *aux, uint8_t t){ push_u8list(aux->tseq, t); }

static inline void app_tseq_hzmaux(HZMAux *aux, uint8_t *ts, uint32_t tlen){ append_array_u8list(aux->tseq, ts, tlen); }

static inline void ready_hzmaux(HZMAux *aux){
	index_single_read_seeds(aux->tseq->buffer, aux->tseq->size, aux->zsize, aux->hz, aux->zmax, aux->hash, aux->bits, aux->seeds, aux->mem_cigar);
}

// beg: the min position of aux->tseq
// end: the max position of aux->tseq
// beg-end defines the region of tseq to be queried by rdseq
// rdqvs is for f5q, otherwise set it to NULL
static inline int align_hzmaux(HZMAux *aux, uint32_t rid, uint8_t *rdseq, uint8_t **rdqvs, int rdlen, int beg, int end, int refine_align, float min_sm){
	wt_seed_t *seed;
	aln_reg_t REG;
	kswx_t x, y;
	uint32_t *cigar, i;
	int n_cigar, j, k, ovl, op, len, x1, x2;
	if(beg < 0) beg = 0;
	if(end <= 0) end = aux->tseq->size;
	clear_hzmpv(aux->rs);
	query_single_read_seeds_by_region(rdseq, rdlen, aux->zsize, aux->hz, aux->zvar, aux->hash, aux->bits, aux->seeds, aux->mem_cigar, aux->anchors, 0, rdlen, beg, end);
	process_hzmps(aux->anchors);
	filter_by_region_hzmps(aux->rs, aux->anchors, 0, beg, end, 0, rdlen);
	clear_wtseedv(aux->windows);
	clear_hzmpv(aux->anchors);
	if(merge_paired_kmers_window(aux->rs, 0, aux->windows, aux->anchors, aux->mem_cache, aux->zsize, aux->zwin, aux->zstep, aux->zovl) == 0) return 0;
	if(chaining_wtseedv(0, rid, 0, aux->windows, 0, aux->windows->size, aux->mem_cache[0], aux->W) < (int)aux->zovl) return 0;
	clear_alnregv(aux->regs);
	clear_u32list(aux->cigar_cache);
	for(i=0;i<aux->windows->size;i++){
		seed = ref_wtseedv(aux->windows, i);
		if(seed->closed) continue;
		push_u32list(aux->cigar_cache, (0 << 4) | 0xFU);
		REG.cigar_off = aux->cigar_cache->size;
		REG.x = fast_seeds_align_hzmo(aux->tseq->buffer, rdseq, seed, aux->anchors, aux->cigar_cache, aux->mem_cache[0], aux->mem_cigar, aux->w, aux->M, aux->X, aux->I, aux->D, aux->E, aux->T);
		REG.cigar_len = aux->cigar_cache->size - REG.cigar_off;
		if(REG.x.aln * 2 < (int)aux->zovl || REG.x.mat < REG.x.aln * min_sm) continue;
		push_alnregv(aux->regs, REG);
	}
	if(aux->regs->size == 0) return 0;
	clear_u32list(aux->cigars);
	x = global_align_regs_hzmo(aux->tseq->size, rdlen, aux->regs, (int[2]){0, aux->tseq->size}, (uint8_t*[2]){aux->tseq->buffer, rdseq}, aux->cigar_cache->buffer, aux->cigars, aux->mem_cache[0], aux->mem_cigar, aux->W, aux->ew, aux->w, aux->M, aux->X, aux->I, aux->D, aux->E, aux->T);
	beg = x.qb - x.tb; if(beg < 0) beg = 0;
	end = x.qe + aux->tseq->size - x.te; if(end > rdlen) end = rdlen;
	ovl = end - beg;
	if(x.score < 0 || x.mat < x.aln * min_sm || x.mat < ovl * min_sm) return 0;
	n_cigar = aux->cigars->size;
	cigar   = aux->cigars->buffer;
	if(refine_align){
		clear_u32list(aux->cigar_cache); append_u32list(aux->cigar_cache, aux->cigars);
		if(rdqvs){
			y = kswx_refine_affine_alignment_5q(rdseq, x.qb, aux->tseq->buffer, x.tb, aux->rw, aux->QCLP, aux->QMIS, aux->QDEL, aux->QEXT, rdqvs, aux->cigar_cache, aux->mem_cache[0], aux->cigars);
		} else y = kswx_refine_alignment(rdseq, x.qb, aux->tseq->buffer, x.tb, aux->rw, aux->M, aux->X, aux->I, aux->D, aux->E, aux->cigar_cache, aux->mem_cache[0], aux->cigars);
		x = y;
		n_cigar = aux->cigars->size;
		cigar   = aux->cigars->buffer;
	}
	aux->hit = x;
	x1 = x.qb;
	x2 = x.tb;
	clear_string(aux->alns[0]);
	clear_string(aux->alns[1]);
	clear_string(aux->alns[2]);
	if(!aux->has_alignment) return 1;
	for(j=0;j<n_cigar;j++){
		op = cigar[j] & 0xF;
		len = cigar[j] >> 4;
		switch(op){
			case 0:
				for(k=0;k<len;k++){
					if(rdseq[x1+k] == aux->tseq->buffer[x2 + k]){
						add_char_string(aux->alns[0], bit_base_table[aux->tseq->buffer[x2 + k]]);
						add_char_string(aux->alns[1], bit_base_table[rdseq[x1 + k]]);
						add_char_string(aux->alns[2], '|');
					} else {
						add_char_string(aux->alns[0], bit_base_table[aux->tseq->buffer[x2 + k]]);
						add_char_string(aux->alns[1], bit_base_table[rdseq[x1 + k]]);
						add_char_string(aux->alns[2], 'X');
					}
				}
				x1 += len;
				x2 += len;
				break;
			case 1:
				for(k=0;k<len;k++){
					add_char_string(aux->alns[0], '-');
					add_char_string(aux->alns[1], bit_base_table[rdseq[x1 + k]]);
					add_char_string(aux->alns[2], '-');
				}
				x1 += len;
				break;
			case 2:
				for(k=0;k<len;k++){
					add_char_string(aux->alns[0], bit_base_table[aux->tseq->buffer[x2 + k]]);
					add_char_string(aux->alns[1], '-');
					add_char_string(aux->alns[2], '-');
				}
				x2 += len;
				break;
		}
	}
	return 1;
}

#endif
