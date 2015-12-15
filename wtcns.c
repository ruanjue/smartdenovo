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

#include "dagcns.h"
#include "bit2vec.h"
#include "kswx.h"
#include "hzm_aln.h"
#include "file_reader.h"
#include "timer.h"
#include "thread.h"

define_list(kswxv, kswx_t);
define_list(strpv, String*);

typedef struct {
	String *name;
	uint32_t n_rd;
	int offset_delta;
	u8list  *rdseqs;
	u8list  *rdqvs;
	u64list *seqoffs;
	u32list *seqlens;
	u64list *rdqoffs;
	u32list *rdqlens;
	u8list  *namestr;
	u32list *seqtags;
	Bit2Vec *flags;
	kswxv   *kswxs;

	int fast_align[2];
	int refine;
	uint32_t zsize, zwin, zstep, zovl, zvar;
	int hz;
	hzmhv   *hash;
	BitVec  *bits;
	hzmv    *seeds;

	int W, ew, w, rw, M, X, O, I, D, E, T;
	uint8_t QCLP, QMIS, QDEL, QEXT;
	int win_size;
	int8_t matrix[16];
	float min_id;

	DAGCNS  *dag;

	int     aln_score;
} WTCNS;

WTCNS* init_wtcns(int W, int M, int X, int O, int I, int D, int E, int T, float min_id){
	WTCNS *g;
	int i;
	g = malloc(sizeof(WTCNS));
	g->name = init_string(1024);
	g->n_rd = 0;
	g->offset_delta = 0;
	g->rdseqs  = init_u8list(1024);
	g->rdqvs   = init_u8list(1024);
	g->seqoffs = init_u64list(1024);
	g->seqlens = init_u32list(1024);
	g->rdqoffs = init_u64list(1024);
	g->rdqlens = init_u32list(1024);
	g->namestr = init_u8list(1024);
	g->seqtags = init_u32list(1024);
	g->flags   = init_bit2vec(1024);
	g->kswxs   = init_kswxv(1024);
	g->hash = init_hzmhv(1024);
	g->bits = init_bitvec(1024);
	g->seeds = init_hzmv(1024);
	g->fast_align[0] = 1;
	g->fast_align[1] = 1;
	g->refine = 1;
	g->zsize = 10;
	g->zwin  = 500;
	g->zstep = 200;
	g->zovl = 200;
	g->hz = 1;
	g->ew = W;
	g->w  = 50;
	g->rw = 8;
	g->W = W;
	g->M = M;
	g->X = X;
	g->O = O;
	g->I = I;
	g->D = D;
	g->E = E;
	g->T = T;
	g->QCLP = -5;
	g->QMIS = -20;
	g->QDEL = -15;
	g->QEXT = -5;
	g->win_size = 1000;
	g->min_id = min_id;
	for(i=0;i<4*4;i++) g->matrix[i] = ((i % 4) == (i / 4))? M : X;
	g->dag = init_dagcns();
	g->aln_score = 0;
	return g;
}

void free_wtcns(WTCNS *g){
	free_string(g->name);
	free_u8list(g->rdseqs);
	free_u8list(g->rdqvs);
	free_u64list(g->seqoffs);
	free_u32list(g->seqlens);
	free_u64list(g->rdqoffs);
	free_u32list(g->rdqlens);
	free_u8list(g->namestr);
	free_u32list(g->seqtags);
	free_bit2vec(g->flags);
	free_kswxv(g->kswxs);
	free_hzmhv(g->hash);
	free_bitvec(g->bits);
	free_hzmv(g->seeds);
	free_dagcns(g->dag);
	free(g);
}

void reset_wtcns(WTCNS *g){
	clear_string(g->name);
	g->n_rd = 0;
	clear_u8list(g->rdseqs);
	clear_u8list(g->rdqvs);
	clear_u64list(g->seqoffs);
	clear_u32list(g->seqlens);
	clear_u64list(g->rdqoffs);
	clear_u32list(g->rdqlens);
	clear_u8list(g->namestr);
	clear_u32list(g->seqtags);
	clear_bit2vec(g->flags);
	clear_kswxv(g->kswxs);
	clear_hzmhv(g->hash);
	clear_bitvec(g->bits);
	clear_hzmv(g->seeds);
	reset_dagcns(g->dag);
	g->aln_score = 0;
	g->offset_delta = 0;
}

void push_wtcns(WTCNS *g, int flag, char *tag, int tag_len, char *seq, int len, int offset){
	int i;
	push_bit2vec(g->flags, flag);
	push_u32list(g->seqlens, len);
	push_u64list(g->seqoffs, g->rdseqs->size);
	push_u32list(g->rdqlens, 0);
	push_u64list(g->rdqoffs, g->rdqvs->size);
	encap_u8list(g->rdseqs, len);
	for(i=0;i<len;i++) lazy_push_u8list(g->rdseqs, base_bit_table[(int)seq[i]]);
	push_kswxv(g->kswxs, (kswx_t){0, offset, offset + len, 0, len, 0, 0, 0, 0, 0});
	push_u32list(g->seqtags, g->namestr->size);
	append_array_u8list(g->namestr, (uint8_t*)tag, tag_len);
	push_u8list(g->namestr, 0);
	g->n_rd ++;
}

void push5q_wtcns(WTCNS *g, int flag, char *tag, int tag_len, char *seq, int len, int offset, char *qvs){
	int i;
	push_bit2vec(g->flags, flag);
	push_u32list(g->seqlens, len);
	push_u64list(g->seqoffs, g->rdseqs->size);
	push_u32list(g->rdqlens, len * 7);
	push_u64list(g->rdqoffs, g->rdqvs->size);
	encap_u8list(g->rdseqs, len);
	for(i=0;i<len;i++) lazy_push_u8list(g->rdseqs, base_bit_table[(int)seq[i]]);
	encap_u8list(g->rdqvs, len * 7);
	for(i=0;i<5*len;i++) lazy_push_u8list(g->rdqvs, qvs[i] - 33);
	for(;i<7*len;i++) lazy_push_u8list(g->rdqvs, base_bit_table[(int)qvs[i]]);
	push_kswxv(g->kswxs, (kswx_t){0, offset, offset + len, 0, len, 0, 0, 0, 0, 0});
	push_u32list(g->seqtags, g->namestr->size);
	append_array_u8list(g->namestr, (uint8_t*)tag, tag_len);
	push_u8list(g->namestr, 0);
	g->n_rd ++;
}

kswx_t standard_aln_rd_vs_cns_wtcns(WTCNS *g, int qlen, uint8_t *q, int tlen, uint8_t *t, u8list *mem_cache, u32list *cigar_cache, u32list *cigars, int iter){
	kswr_t r;
	kswx_t x, y;
	uint32_t *cigar;
	uint8_t *t2, *q2;
	int qe, te, w, nc;
	clear_u32list(cigars);
	//qe = qlen > 3 * g->win_size? 3 * g->win_size : qlen;
	//te = tlen > 9 * g->win_size? 9 * g->win_size : tlen;
	qe = qlen;
	te = tlen;
	encap_u8list(mem_cache, kswx_roundup8x(mem_cache->size) - mem_cache->size + kswx_roundup8x(te) + kswx_roundup8x(qe));
	// ksw_align2 will reverse t[0,te] and q[0,qe] to find KSW_XSTART, it is thread-unsafe
	// t2 and q2 are used to avoid that
	t2 = (mem_cache->buffer + kswx_roundup8x(mem_cache->size));
	q2 = t2 + kswx_roundup8x(te);
	memcpy(t2, t, te);
	memcpy(q2, q, qe);
	r = ksw_align2(qe, q2, te, t2, 4, g->matrix, - (iter? g->I : g->O), - g->E, - (iter? g->D : g->O), - g->E, KSW_XSTART, NULL);
	if(r.qb == -1 || r.tb == -1 || r.qe == -1 || r.te == -1){
		return KSWX_NULL;
	}
	r.qe ++; r.te ++;
	// extend left
	y = kswx_extend_align_shift_core(r.qb, q + r.qb - 1, r.tb, t + r.tb - 1, -1, r.score, g->W, g->M, g->X, iter? g->I : g->O, iter? g->D : g->O, g->E, g->T, mem_cache, cigar_cache);
	reverse_u32list(cigar_cache);
	append_u32list(cigars, cigar_cache);
	w = (r.qe - r.qb < r.te - r.tb)? ((r.te - r.tb) - (r.qe - r.qb)) : ((r.qe - r.qb) - (r.te - r.tb));
	if(w > g->W) w = g->W;
	nc = 0; cigar = NULL;
	x = kswx_gen_cigar_core2(qe, q, te, t, r, 4, g->matrix, w, iter? g->I : g->O, iter? g->D : g->O, g->E, &nc, &cigar);
	append_array_u32list(cigars, cigar, nc);
	if(cigar) free(cigar);
	x.score = y.score;
	x.aln   += y.aln;
	x.mat   += y.mat;
	x.mis   += y.mis;
	x.ins   += y.ins;
	x.del   += y.del;
	x.qb    -= y.qe;
	x.tb    -= y.te;
	// right
	y = kswx_extend_align_shift_core(qlen - x.qe, q + x.qe, tlen - x.te, t + x.te, 1, r.score, g->W, g->M, g->X, iter? g->I : g->O, iter? g->D : g->O, g->E, g->T, mem_cache, cigar_cache);
	append_u32list(cigars, cigar_cache);
	x.score = y.score;
	x.aln   += y.aln;
	x.mat   += y.mat;
	x.mis   += y.mis;
	x.ins   += y.ins;
	x.del   += y.del;
	x.qe    += y.qe;
	x.te    += y.te;
	return x;
}

int gen_backbone_wtcns(WTCNS *g){
	kswr_t r;
	kswx_t *p;
	uint32_t i;
	uint8_t *t, *q;
	int tlen, qlen, tb, te, qb, qe;
	clear_u8list(g->dag->cns);
	for(i=0;i<g->n_rd;i++){
		p = ref_kswxv(g->kswxs, i);
		if(get_bit2vec(g->flags, i) == 0){
			t = g->dag->cns->buffer;
			tlen = g->dag->cns->size;
			q = g->rdseqs->buffer + g->seqoffs->buffer[i];
			qlen = g->seqlens->buffer[i];
			if(tlen == 0){
				append_array_u8list(g->dag->cns, q, qlen);
				g->offset_delta = - p->tb;
			} else {
				tb = p->tb + g->offset_delta;
				tb = tb > g->win_size? tb - g->win_size : 0;
				te = tb + 3 * g->win_size;
				te = te > tlen? tlen : te;
				qb = 0; qe = (qb + g->win_size > qlen)? qlen : qb + g->win_size;
				r = ksw_align2(qe - qb, q + qb, te - tb, t + tb, 4, g->matrix, - g->O, - g->E, - g->O, - g->E, KSW_XSTART, NULL);
				if(r.qb == -1 || r.tb == -1 || r.qe == -1 || r.te == -1){
					r.qb = 0; r.tb = p->tb + g->offset_delta; // fail to get correct offset, use suggested offset
				}
				r.tb += tb; r.te += tb + 1;
				r.qb += qb; r.qe += qb + 1;
				g->offset_delta = (r.tb - r.qb) - (p->tb - p->qb);
				g->dag->cns->size = r.tb;
				append_array_u8list(g->dag->cns, q + r.qb, qlen - r.qb);
			}
		}
		p->tb += g->offset_delta;
		p->te += g->offset_delta;
	}
	return 1;
}

int aln_read_wtcns(WTCNS *g, uint32_t rid, String *alns[2], int iter, wtseedv *windows, hzmpv *rs, hzmpv *anchors, u8list *mem_cache[2], u32list *mem_cigar, u32list *cigar_cache, u32list *cigars, alnregv *regs, int has_mismatch, int fast_align, int refine){
	wt_seed_t *seed;
	aln_reg_t REG;
	kswx_t x, y, *s;
	uint8_t *rd, *qvs[7];
	uint32_t *cigar, i;
	int rdlen, beg, end, n_cigar, j, k, ovl, op, len, x1, x2;
	if(get_bit2vec(g->flags, rid) == 3) return 0;
	rdlen = g->seqlens->buffer[rid];
	rd = g->rdseqs->buffer + g->seqoffs->buffer[rid];
	if(hzm_debug){
		x = get_kswxv(g->kswxs, rid);
		fprintf(stdout, "LAY\t%s\t%c\t%d\t%d\t%d\t%s\t%c\t%d\t%d\t%d", "BACKBONE", '+', (int)g->dag->cns->size, x.tb, x.te, g->namestr->buffer + g->seqtags->buffer[rid], '+', rdlen, x.qb, x.qe);
		fprintf(stdout, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
	}
	s = ref_kswxv(g->kswxs, rid);
	s->tb -= s->qb; s->qb = 0;
	s->te += rdlen - s->qe; s->qe = rdlen;
	beg = s->tb;
	if(beg < 0) beg = 0;
	end = s->te;
	if(end > (int)g->dag->cns->size) end = g->dag->cns->size;
	if(fast_align){
		clear_hzmpv(rs);
		query_single_read_seeds(rd, rdlen, g->zsize, g->hz, g->zvar, g->hash, g->bits, g->seeds, mem_cigar, anchors);
		process_hzmps(anchors, rdlen);
		filter_by_region_hzmps(rs, anchors, 0, beg, end, 0, rdlen);
		clear_wtseedv(windows);
		clear_hzmpv(anchors);
		if(merge_paired_kmers_window(rs, 0, windows, anchors, mem_cache, g->zsize, g->zwin, g->zstep, g->zovl) == 0){
			set_bit2vec(g->flags, rid, 3);
			return 0;
		}
		if(chaining_wtseedv(0, rid, 0, windows, 0, windows->size, mem_cache[0], g->W) < (int)g->zovl){
			set_bit2vec(g->flags, rid, 3);
			return 0;
		}
		clear_alnregv(regs);
		clear_u32list(cigar_cache);
		for(i=0;i<windows->size;i++){
			seed = ref_wtseedv(windows, i);
			if(seed->closed) continue;
			push_u32list(cigar_cache, (0 << 4) | 0xFU);
			REG.cigar_off = cigar_cache->size;
			REG.x = fast_seeds_align_hzmo(g->dag->cns->buffer, rd, seed, anchors, cigar_cache, mem_cache[0], mem_cigar, g->w, g->M, g->X, iter? g->I : g->O, iter? g->D : g->O, g->E, g->T);
			REG.cigar_len = cigar_cache->size - REG.cigar_off;
			if(REG.x.aln * 2 < (int)g->zovl || REG.x.mat < REG.x.aln * g->min_id) continue;
			push_alnregv(regs, REG);
		}
		if(regs->size == 0){
			set_bit2vec(g->flags, rid, 3);
			return 0;
		}
		clear_u32list(cigars);
		x = global_align_regs_hzmo(g->dag->cns->size, rdlen, regs, (int[2]){0, g->dag->cns->size}, (uint8_t*[2]){g->dag->cns->buffer, rd}, cigar_cache->buffer, cigars, mem_cache[0], mem_cigar, g->W, g->ew, g->w, g->M, g->X, iter? g->I : g->O, iter? g->D : g->O, g->E, g->T);
		beg = x.qb - x.tb; if(beg < 0) beg = 0;
		end = x.qe + g->dag->cns->size - x.te; if(end > rdlen) end = rdlen;
		ovl = end - beg;
		if(x.score < 0 || x.mat < x.aln * g->min_id || x.mat < ovl * g->min_id){
			if(hzm_debug){
				fprintf(stdout, "DISCARD %s\n", g->namestr->buffer + g->seqtags->buffer[rid]);
			}
			set_bit2vec(g->flags, rid, 3);
			return 0;
		}
		set_kswxv(g->kswxs, rid, x);
		n_cigar = cigars->size;
		cigar   = cigars->buffer;
	} else {
		x = standard_aln_rd_vs_cns_wtcns(g, rdlen, rd, g->dag->cns->size, g->dag->cns->buffer, mem_cache[0], cigar_cache, cigars, iter);
		beg = x.qb - x.tb; if(beg < 0) beg = 0;
		end = x.qe + g->dag->cns->size - x.te; if(end > rdlen) end = rdlen;
		ovl = end - beg;
		if(x.score < 0 || x.mat < x.aln * g->min_id || x.mat < ovl * g->min_id){
			if(hzm_debug){
				fprintf(stdout, "DISCARD %s\n", g->namestr->buffer + g->seqtags->buffer[rid]);
			}
			set_bit2vec(g->flags, rid, 3);
			return 0;
		}
		set_kswxv(g->kswxs, rid, x);
		n_cigar = cigars->size;
		cigar   = cigars->buffer;
	}
	if(refine){
		clear_u32list(cigar_cache); append_u32list(cigar_cache, cigars);
		if((int)g->rdqlens->buffer[rid] == 7 * rdlen){
			qvs[0] = g->rdqvs->buffer + g->rdqoffs->buffer[rid];
			qvs[1] = qvs[0] + rdlen;
			qvs[2] = qvs[1] + rdlen;
			qvs[3] = qvs[2] + rdlen;
			qvs[4] = qvs[3] + rdlen;
			qvs[5] = qvs[4] + rdlen;
			qvs[6] = qvs[5] + rdlen;
			y = kswx_refine_affine_alignment_5q(rd, x.qb, g->dag->cns->buffer, x.tb, g->rw, g->QCLP, g->QMIS, g->QDEL, g->QEXT, qvs, cigar_cache, mem_cache[0], cigars);
		} else y = kswx_refine_alignment(rd, x.qb, g->dag->cns->buffer, x.tb, g->rw, g->M, g->X, iter? g->I : g->O, iter? g->D : g->O, g->E, cigar_cache, mem_cache[0], cigars);
		if(hzm_debug){
			fprintf(stdout, "refine %s %d[%d,%d,%d,%d] -> %d[%d,%d,%d,%d]\n", g->namestr->buffer + g->seqtags->buffer[rid],
				x.score, x.mat, x.mis, x.ins, x.del,
				y.score, y.mat, y.mis, y.ins, y.del);
		}
		x = y;
		set_kswxv(g->kswxs, rid, x);
		n_cigar = cigars->size;
		cigar   = cigars->buffer;
	}
	x1 = x.qb;
	x2 = x.tb;
	clear_string(alns[0]);
	clear_string(alns[1]);
	for(j=0;j<n_cigar;j++){
		op = cigar[j] & 0xF;
		len = cigar[j] >> 4;
		switch(op){
			case 0:
				for(k=0;k<len;k++){
					if(rd[x1+k] == g->dag->cns->buffer[x2 + k]){
						add_char_string(alns[0], bit_base_table[g->dag->cns->buffer[x2 + k]]);
						add_char_string(alns[1], bit_base_table[rd[x1 + k]]);
					} else {
						add_char_string(alns[0], bit_base_table[g->dag->cns->buffer[x2 + k]]);
						if(!has_mismatch){
							add_char_string(alns[1], '-');
							add_char_string(alns[0], '-');
						}
						add_char_string(alns[1], bit_base_table[rd[x1 + k]]);
					}
				}
				x1 += len;
				x2 += len;
				break;
			case 1:
				for(k=0;k<len;k++){
					add_char_string(alns[0], '-');
					add_char_string(alns[1], bit_base_table[rd[x1 + k]]);
				}
				x1 += len;
				break;
			case 2:
				for(k=0;k<len;k++){
					add_char_string(alns[0], bit_base_table[g->dag->cns->buffer[x2 + k]]);
					add_char_string(alns[1], '-');
				}
				x2 += len;
				break;
		}
	}
	return 1;
}

thread_beg_def(maln);
WTCNS *g;
uint32_t rid;
int ret, iter, has_mismatch, fast_align, refine;
String *alns[2];
thread_end_def(maln);

thread_beg_func(maln);
u8list *mem_cache[2];
u32list *mem_cigar, *cigar_cache, *cigars;
wtseedv *windows;
hzmpv *rs, *anchors;
alnregv *regs;
mem_cache[0] = init_u8list(1024);
mem_cache[1] = init_u8list(1024);
mem_cigar = init_u32list(1024);
cigar_cache = init_u32list(1024);
cigars = init_u32list(1024);
windows = init_wtseedv(64);
rs = init_hzmpv(1024);
anchors = init_hzmpv(1024);
regs = init_alnregv(64);
thread_beg_loop(maln);
maln->ret = 0;
if(maln->rid != 0xFFFFFFFFU){
	maln->ret = aln_read_wtcns(maln->g, maln->rid, (String**)maln->alns, maln->iter, windows, rs, anchors, mem_cache, mem_cigar, cigar_cache, cigars, regs, maln->has_mismatch, maln->fast_align, maln->refine);
}
thread_end_loop(maln);
free_u8list(mem_cache[0]);
free_u8list(mem_cache[1]);
free_u32list(mem_cigar);
free_u32list(cigar_cache);
free_u32list(cigars);
free_wtseedv(windows);
free_hzmpv(rs);
free_hzmpv(anchors);
free_alnregv(regs);
thread_end_func(maln);

void run_wtcns(WTCNS *g, int ncpu, int max_iter, FILE *aln, int vmsa, uint32_t min_cnt, float min_freq){
	kswx_t x, *s;
	u32list *map;
	u32list *aux;
	String *cnsid;
	strpv  *alns[2];
	BitVec *keys;
	uint16_t *bases[4];
	uint32_t I, i, j, k, n_key, a, b;
	int size, iter;
	char c;
	thread_preprocess(maln);
	gen_backbone_wtcns(g);
	fprintf(stderr, "[%s]\"%s\" generate backbone length=%d\n", date(), g->name->string, (int)g->dag->cns->size);
	//if(max_iter < 1) max_iter = 1;
	map = init_u32list(1024);
	aux = init_u32list(1024);
	alns[0] = init_strpv(g->n_rd);
	alns[1] = init_strpv(g->n_rd);
	for(i=0;i<g->n_rd;i++){
		push_strpv(alns[0], init_string(64));
		push_strpv(alns[1], init_string(64));
	}
	for(iter=0;iter<max_iter;iter++){
		if(g->fast_align[0]){
			index_single_read_seeds(g->dag->cns->buffer, g->dag->cns->size, g->zsize, g->hz, 0xFFFFFFFFU, g->hash, g->bits, g->seeds, aux);
		}
		// prepare DAG
		gen_pregraph_dagcns(g->dag);
		// Align, update DAG
		thread_beg_init(maln, ncpu);
		maln->g = g;
		maln->rid = 0xFFFFFFFFU;
		maln->ret = 0;
		maln->iter = iter;
		maln->alns[0] = NULL;
		maln->alns[1] = NULL;
		maln->has_mismatch = 0;
		maln->fast_align = g->fast_align[0];
		maln->refine = g->refine;
		thread_end_init(maln);
		g->aln_score = 0;
		if(1){
			for(i=0;i<g->n_rd;i++){
				thread_wait_one(maln);
				if(maln->ret){
					x = get_kswxv(g->kswxs, maln->rid);
					g->aln_score += x.score;
					//alignment2dagcns(g->dag, x.tb, x.te, (char*[2]){maln->alns[0]->string, maln->alns[1]->string}, maln->alns[0]->size);
					maln->ret = 0;
				}
				maln->rid = i;
				maln->iter = iter;
				maln->alns[0] = get_strpv(alns[0], maln->rid);
				maln->alns[1] = get_strpv(alns[1], maln->rid);
				thread_wake(maln);
			}
			thread_wait_all(maln);
			thread_beg_iter(maln);
			if(maln->ret){
				x = get_kswxv(g->kswxs, maln->rid);
				g->aln_score += x.score;
				//alignment2dagcns(g->dag, x.tb, x.te, (char*[2]){maln->alns[0]->string, maln->alns[1]->string}, maln->alns[0]->size);
				maln->ret = 0;
			}
			thread_end_iter(maln);
		}
		thread_beg_close(maln);
		//free_string(maln->alns[0]);
		//free_string(maln->alns[1]);
		thread_end_close(maln);
		clear_u32list(aux);
		for(i=0;i<g->n_rd;i++) push_u32list(aux, i);
		sort_array(aux->buffer, aux->size, uint32_t, g->kswxs->buffer[b].score > g->kswxs->buffer[a].score);
		for(i=0;i<g->n_rd;i++){
			j = aux->buffer[i];
			x = get_kswxv(g->kswxs, j);
			if(get_bit2vec(g->flags, j) == 3) continue;
			alignment2dagcns(g->dag, x.tb, x.te, (char*[2]){alns[0]->buffer[j]->string, alns[1]->buffer[j]->string}, alns[0]->buffer[j]->size);
		}
		fprintf(stderr, "[%s]\"%s\" align aln_score=%d\n", date(), g->name->string, g->aln_score);
		// Merge nodes
		merge_nodes_dagcns(g->dag);
		fprintf(stderr, "[%s]\"%s\" merge nodes\n", date(), g->name->string);
		// Call CNS
		size = g->dag->cns->size;
		gen_consensus_dagcns(g->dag, map);
		for(i=0;i<g->n_rd;i++){
			if(get_bit2vec(g->flags, i) == 3) continue;
			s = ref_kswxv(g->kswxs, i);
			if(s->tb < 0) s->tb = 0;
			if(s->tb > size) s->tb = size;
			s->tb = map->buffer[s->tb];
			if(s->te < 0) s->te = 0;
			if(s->te > size) s->te = size;
			s->te = map->buffer[s->te];
		}
		fprintf(stderr, "[%s]\"%s\" iter%d length=%d aln_score=%d cns_score=%f\n", date(), g->name->string, iter + 1, (int)g->dag->cns->size, g->aln_score, g->dag->cns_score);
		//fprintf(stderr, ">iter%d score=%d\n", iter + 1, tot_score);
		//print_seq_wtcns(g, stderr);
		//fflush(stderr);
	}
	{
		cnsid = clone_string(g->name);
		for(i=0;(int)i<cnsid->size;i++){
			if(cnsid->string[i] == ' ' || cnsid->string[i] == '\t'){
				cnsid->size = i;
				cnsid->string[i] = '\0';
			}
		}
	}
	if(aln){
		if(g->fast_align[1]){
			index_single_read_seeds(g->dag->cns->buffer, g->dag->cns->size, g->zsize, g->hz, 0xFFFFFFFFU, g->hash, g->bits, g->seeds, aux);
		}
		//if(1) ncpu = 1;
		thread_beg_init(maln, ncpu);
		maln->g = g;
		maln->rid = 0xFFFFFFFFU;
		maln->ret = 0;
		maln->iter = iter;
		maln->alns[0] = NULL;
		maln->alns[1] = NULL;
		maln->has_mismatch = 1;
		maln->fast_align = g->fast_align[1];
		maln->refine = g->refine;
		thread_end_init(maln);
		for(i=0;i<g->n_rd;i++){
			thread_wait_one(maln);
			maln->rid = i;
			maln->iter = iter;
			maln->alns[0] = get_strpv(alns[0], maln->rid);
			maln->alns[1] = get_strpv(alns[1], maln->rid);
			thread_wake(maln);
		}
		thread_wait_all(maln);
		thread_beg_close(maln);
		thread_end_close(maln);
		if(vmsa){
			bases[0] = calloc(g->dag->cns->size, sizeof(uint16_t));
			bases[1] = calloc(g->dag->cns->size, sizeof(uint16_t));
			bases[2] = calloc(g->dag->cns->size, sizeof(uint16_t));
			bases[3] = calloc(g->dag->cns->size, sizeof(uint16_t));
		}
		for(i=0;i<g->n_rd;i++){
			if(get_bit2vec(g->flags, i) == 3) continue;
			x = get_kswxv(g->kswxs, i);
			fprintf(aln, "%s\t+\t%d\t%d\t%d\t%s\t+\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", g->namestr->buffer + g->seqtags->buffer[i], g->seqlens->buffer[i], x.qb, x.qe, cnsid->string, (int)g->dag->cns->size, x.tb, x.te, x.score, 1.0 * x.mat / (x.aln + 1), x.mat, x.mis, x.ins, x.del);
			fprintf(aln, "Q\t%s\n", alns[1]->buffer[i]->string);
			fprintf(aln, "T\t%s\n", alns[0]->buffer[i]->string);
			fprintf(aln, "M\t");
			k = x.tb;
			c = ' ';
			int m, n, margin;
			margin = 3;
			n = 0;
			for(j=0;(int)j<alns[0]->buffer[i]->size;j++){
				if(alns[0]->buffer[i]->string[j] == '-'){
					if(alns[1]->buffer[i]->string[j] == '-') c = ' ';
					else c = '-';
					if(vmsa){
						for(m=margin+1;m+margin<=n;m++){
							bases[base_bit_table[(int)alns[1]->buffer[i]->string[j-m]]][k-m] ++;
							alns[1]->buffer[i]->string[j-m] = lc(alns[1]->buffer[i]->string[j-m]);
						}
					}
					n = 0;
				} else {
					if(alns[1]->buffer[i]->string[j] == '-'){
						c = '-';
						if(vmsa){
							for(m=margin+1;m+margin<=n;m++){
								bases[base_bit_table[(int)alns[1]->buffer[i]->string[j-m]]][k-m] ++;
								alns[1]->buffer[i]->string[j-m] = lc(alns[1]->buffer[i]->string[j-m]);
							}
						}
						n = 0;
					} else if(alns[0]->buffer[i]->string[j] == alns[1]->buffer[i]->string[j]){
						n ++;
						c = ' ';
					} else {
						n ++;
						c = '*';
					}
					//b = base_bit_table[(int)alns[1]->buffer[i]->string[j]];
					//if(b < 4) bases[b][k] ++;
					k ++;
				}
				fputc(c, aln);
			}
			if(vmsa){
				for(m=margin+1;m+margin<=n;m++){
					bases[base_bit_table[(int)alns[1]->buffer[i]->string[j-m]]][k-m] ++;
					alns[1]->buffer[i]->string[j-m] = lc(alns[1]->buffer[i]->string[j-m]);
				}
			}
			fprintf(aln, "\n\n");
		}
		if(vmsa){
			keys = init_bitvec(g->dag->cns->size);
			n_key = 0;
			for(i=0;i<g->dag->cns->size;i++){
				a = b = 0;
				for(j=1;j<4;j++){
					if(bases[j][i] >= bases[a][i]){
						b = a;
						a = j;
					} else if(bases[j][i] > bases[b][i]){
						b = j;
					}
				}
				if(a != b && bases[b][i] >= min_cnt && bases[b][i] >= min_freq * bases[a][i]){
					n_key ++;
					one_bitvec(keys, i);
				}
			}
			free(bases[0]);
			free(bases[1]);
			free(bases[2]);
			free(bases[3]);
			index_bitvec(keys);
			clear_u32list(aux);
			for(i=0;i<g->n_rd;i++) push_u32list(aux, i);
			sort_array(aux->buffer, aux->size, uint32_t, g->kswxs->buffer[a].tb > g->kswxs->buffer[b].tb);
			for(I=0;I<g->n_rd;I++){
				i = aux->buffer[I];
				if(get_bit2vec(g->flags, i) == 3) continue;
				x = get_kswxv(g->kswxs, i);
				fprintf(aln, "MATRIX\t%s\t", g->namestr->buffer + g->seqtags->buffer[i]);
				k = rank_bitvec(keys, x.tb);
				for(j=0;j<k;j++) fputc('-', aln);
				k = x.tb;
				for(j=0;(int)j<alns[0]->buffer[i]->size;j++){
					if(alns[0]->buffer[i]->string[j] == '-'){
					} else {
						if(get_bitvec(keys, k)){
							if(!(alns[1]->buffer[i]->string[j] >= 'a' && alns[1]->buffer[i]->string[j] <= 'z')) fputc('-', aln);
							else if(base_bit_table[(int)alns[0]->buffer[i]->string[j]] == base_bit_table[(int)alns[1]->buffer[i]->string[j]]) fputc('.', aln);
							else fputc(alns[1]->buffer[i]->string[j], aln);
						}
						k ++;
					}
				}
				for(j=rank_bitvec(keys, x.te);j<n_key;j++) fputc('-', aln);
				fprintf(aln, "\n");
			}
			free_bitvec(keys);
		}
		fprintf(stderr, "[%s]\"%s\" finish alignments against final consensus\n", date(), g->name->string);
	}
	for(i=0;i<g->n_rd;i++) free_string(get_strpv(alns[0], i));
	for(i=0;i<g->n_rd;i++) free_string(get_strpv(alns[1], i));
	free_strpv(alns[0]);
	free_strpv(alns[1]);
	free_string(cnsid);
	free_u32list(map);
	free_u32list(aux);
}

int usage(){
	printf(
	"WTCNS: Consensus caller\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtcns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [16]\n"
	" -P <int>    Total parallel jobs, [1]\n"
	" -p <int>    Index of current job (0-based), [0]\n"
	"             Suppose to run wtcns for the same layout file parallelly in 60 cpu. For cpu1, -P 60 -p 0; cpu2, -P 60 -p 1, ...\n"
	" -i <string> Input file, layout from wtlay, +\n"
	" -o <string> Output file, consensus sequences, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -H          Trun on homopolymer compression\n"
	" -z <int>    Zmer size, 5 <= <-z> <= 16, [10]\n"
	" -y <int>    Zmer window, [800]\n"
	" -R <int>    Minimum size of seeding region within zmer window, [200]\n"
	" -l <int>    Maximum variant of uncompressed sizes between two matched zmer, [2]\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: insertion or deletion, used in first round [-2]\n"
	" -I <int>    Alignment penalty: insertion, used in rounds after first, [-2]\n"
	" -D <int>    Alignment penalty: deletion, used in rounds after first, [-3]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping [-10]\n"
	//" -Q <X,D>    Default PhredQV, X(int): default mismatch PhredQV, D(int): default deletion PhredQV [20,15]\n"
	//"             Those parameters will be used in refine-alignment when provided the layout file with 5quality+2tag\n"
	//"             To use base quality, you need to invoke wtlay with f5q input\n"
	" -F          Disable PhreadQV in refine-alignment\n"
	" -w <int>    Minimum bandwidth, iteratively doubled to maximum [50]\n"
	" -W <int>    Maximum bandwidth, [3200]\n"
	" -e <int>    Maximum bandwidth at ending extension, [800]\n"
	" -r <int>    Basic bandwidth in refine-alignment, [8]\n"
	" -m <float>  Minimum alignment identity, [0.5]\n"
	" -Y <float>  Penalty of backbone edge in calling consensus, [0.5]\n"
	" -N <float>  Penalty of alternative edge in calling consensus, [0.2]\n"
	"             The above two options control whether the consensus look like backbone or alternative\n"
	"             Default 0.5 and 0.2, will let the consensus don't look like backbone\n"
	" -n <int>    Number of iterations for consensus calling, the larger, the accurater, the slower [6]\n"
	" -a <string> Align reads against final consensus, and output to <-a>\n"
	" -A          Disable fast zmer align in final aligning (see -a), use standard smith-waterman\n"
	"             More than once -A, will disable fast zmer align in all process\n"
	" -V <float> Ouput call variants and print to <-a>, -V 2.05 mean: min_allele_count>=2,min_allele_freq>=0.05\n"
	" -v          Verbose, BE careful, HUGEEEEEEEE output on STDOUT\n"
	"\n"
	"Example: \n"
	"$> wtcns wt.lay > wt.lay.cns.fa 2>log.cns\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	cplist *lays;
	char *outf, *alnf;
	WTCNS *g;
	FILE *aln, *out;
	long job;
	int M, X, I, D, O, E, T, W, ew, rw, w, zsize, hz, zwin, zovl, zvar, n, ncpu, c, flag, n_job, i_job, work, refine, vmsa, fast_align[2], overwrite, min_cnt;
	int use_qv;
	float min_id, min_freq, ref_penalty, alt_penalty;
	HZM_FAST_WINDOW_KMER_CHAINING = 0;
	ncpu = 16;
	zsize = 10;
	zvar = 2;
	hz = 0;
	zwin = 800;
	zovl = 200;
	refine = 1;
	use_qv = 1;
	w = 50;
	ew = 800;
	rw = 8; // basic bandwidth for refine-alignment
	W = 3200;
	M = 2;
	X = -5;
	O = -2;
	I = -2;
	D = -3;
	E = -1;
	T = -10;
	min_id = 0.5;
	n = 6;
	n_job = 1;
	i_job = 0;
	fast_align[0] = 1;
	fast_align[1] = 1;
	outf = NULL;
	alnf = NULL;
	vmsa = 0;
	min_cnt = 2;
	min_freq = 0.05;
	ref_penalty = 0.5;
	alt_penalty = 0.2;
	overwrite = 0;
	lays = init_cplist(4);
	while((c = getopt(argc, argv, "ht:P:p:i:o:fAHz:y:R:l:M:X:O:I:D:E:T:Fw:W:e:r:m:Y:N:n:a:V:v")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'P': n_job = atoi(optarg); break;
			case 'p': i_job = atoi(optarg); break;
			case 'i': push_cplist(lays, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'A': if(fast_align[1]) fast_align[1] = 0; else fast_align[0] = 0; break;
			case 'H': hz = 1; break;
			case 'z': zsize = atoi(optarg); break;
			case 'y': zwin  = atoi(optarg); break;
			case 'R': zovl  = atoi(optarg); break;
			case 'l': zvar  = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			case 'F': use_qv = 0; break;
			case 'w': w  = atoi(optarg); break;
			case 'W': W  = atoi(optarg); break;
			case 'e': ew = atoi(optarg); break;
			case 'r': rw = atoi(optarg); break;
			case 'm': min_id = atof(optarg); break;
			case 'Y': ref_penalty = atof(optarg); break;
			case 'N': alt_penalty = atof(optarg); break;
			case 'n': n = atoi(optarg); break;
			case 'a': alnf = optarg; break;
			case 'V': vmsa = 1; min_freq = atof(optarg); min_cnt = min_freq; min_freq -= min_cnt; break;
			case 'v': hzm_debug ++; break;
			default: return usage();
		}
	}
	if(outf && !overwrite && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	for(c=optind;c<argc;c++) push_cplist(lays, argv[c]);
	if(lays->size){
		fr = fopen_m_filereader(lays->size, lays->buffer);
	} else fr = stdin_filereader();
	if(outf) out = fopen(outf, "w");
	else out = stdout;
	if(alnf){
		aln = fopen(alnf, "w");
	} else aln = NULL;
	g = init_wtcns(W, M, X, O, I, D, E, T, min_id);
	g->win_size = zwin;
	g->zsize = zsize;
	g->zwin  = zwin;
	g->zstep = zwin / 2;
	g->zovl = zovl;
	g->zvar = zvar;
	g->refine = refine;
	g->ew = ew;
	g->rw = rw;
	g->w = w;
	g->fast_align[0] = fast_align[0];
	g->fast_align[1] = fast_align[1];
	g->dag->ref_penalty = ref_penalty;
	g->dag->alt_penalty = alt_penalty;
	work = 0;
	job = 0;
	while(1){
		c = fread_table(fr);
		if(c == -1 || fr->line->string[0] == '>'){
			if(g->n_rd){
				run_wtcns(g, ncpu, n, aln, vmsa, min_cnt, min_freq);
				fprintf(out, ">%s\n", g->name->string);
				print_seq_dagcns(g->dag, out);
				fflush(out);
				reset_wtcns(g);
			}
			if((job % n_job) == i_job) work = 1;
			else work = 0;
			job ++;
			if(c == -1) break;
			if(work){
				append_string(g->name, fr->line->string + 1, fr->line->size - 1);
				trim_string(g->name);
				fprintf(stderr, "[%s]\"%s\" init\n", date(), g->name->string);
			}
			continue;
		}
		if(work == 0) continue;
		if(fr->line->string[0] == '#') continue;
		if(c < 6) continue;
		switch(get_col_str(fr, 0)[0]){
			case 'Y': flag = 0; break;
			case 'N': flag = 1; break;
			case 'n': flag = 2; break;
			default: flag = 3;
		}
		if(use_qv && c > 6 && get_col_len(fr, 6) == 7 * get_col_len(fr, 5)){ // f5q
			push5q_wtcns(g, flag, get_col_str(fr, 1), get_col_len(fr, 1), get_col_str(fr, 5), get_col_len(fr, 5), atoi(get_col_str(fr, 3)), get_col_str(fr, 6));
		} else {
			push_wtcns(g, flag, get_col_str(fr, 1), get_col_len(fr, 1), get_col_str(fr, 5), get_col_len(fr, 5), atoi(get_col_str(fr, 3)));
		}
	}
	fclose_filereader(fr);
	free_cplist(lays);
	if(outf) fclose(out);
	if(aln) fclose(aln);
	free_wtcns(g);
	return 0;
}
