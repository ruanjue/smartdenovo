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
#include "dna.h"
#include "kswx.h"
#include "hzm_aln.h"
#include "file_reader.h"
#include "bitvec.h"
#include "thread.h"
#include <getopt.h>

#define HZMH_KMER_MOD	1024
#define hzmh_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7) ^ ((K) >> 12))

#define HZM_MAX_SEED_KMER	32

typedef struct {
	uint32_t rd_id:31, dir:1;
} hzmin_t;
define_list(hzminv, hzmin_t);

typedef struct {
	hzm_t p1;
	hzmin_t *p2;
	uint64_t beg, end;
} hz_ref_t;
define_list(hzrefv, hz_ref_t);
#define tmp_util1(r) (((r)->p2->rd_id << 1) | (r)->p2->dir)
#define heap_cmp_hz_ref_macro(refs, a, b)	\
({	\
	hz_ref_t *r1, *r2;	\
	r1 = ref_hzrefv(refs, a);	\
	r2 = ref_hzrefv(refs, b);	\
	(tmp_util1(r1) > tmp_util1(r2))? 1 : ((tmp_util1(r1) < tmp_util1(r2))? -1 : ((r1->p1.off > r2->p1.off)? 1 : ((r1->p1.off < r2->p1.off)? -1 : 0)));	\
})

#define WTZMO_SR_SEED_MAX_OVL	((1U<<20)-1)
#define WTZMO_SR_SEED_MAX_CNT	((1U<<12)-1)
typedef struct {
	uint32_t pb1:31, dir1:1;
	uint32_t pb2:30, dir2:1, closed:1;
	uint32_t ovl:20, cnt:12;
	uint32_t anchors[2];
	int beg[2];
	int end[2];
	int avg, std, pseudo_ovl;
} sr_seed_t;
static const sr_seed_t SR_SEED_NULL = (sr_seed_t){0, 0, 0, 0, 0, 0, 0, {0, 0}, {0, 0}, {0, 0}, 0, 0, 0};
define_list(srseedv, sr_seed_t);
#define srseed_hashcode(E) (((E).pb2 << 1) | (E).dir2)
#define srseed_hashequals(E1, E2) ((E1).pb2 == (E2).pb2 && (E1).dir2 == (E2).dir2)
define_hashset(srseedhash, sr_seed_t, srseed_hashcode, srseed_hashequals);

typedef struct {
	uint32_t pb1:31, dir1:1, pb2:31, dir2:1;
	int qb, qe, tb, te;
	int score, mat, mis, ins, del, aln;
	char *cigar;
} wt_ovl_t;
define_list(wtovlv, wt_ovl_t);

#define _ovl_uniq_long_id(id1, id2, dir) ((((uint64_t)(id1)) << 33) | (((uint64_t)(id2)) << 1) | (dir))
#define ovl_uniq_long_id(id1, id2, dir) (((id1) < (id2))? _ovl_uniq_long_id(id1, id2, dir) : _ovl_uniq_long_id(id2, id1, dir))

typedef struct {
	uint64_t rdoff:40, rdlen:24;
	char *rdname;
} pbread_t;
define_list(pbreadv, pbread_t);

typedef struct {
	uint64_t link:62, closed:1, dir:1;
} wtlink_t;

typedef struct {
	uint32_t nodes[2];
	wtlink_t links[2];
} wtedge_t;

typedef struct {
	uint64_t seqoff:40, seqlen:23, closed:1;
	struct { uint64_t bt_visit:19, bt_dir:1, bt_idx:10, bt_off:34; };
	wtlink_t edges[2];
} wtnode_t;

typedef struct {
	uint32_t n_rd;
	BaseBank *rdseqs;
	pbreadv  *reads;
	cuhash   *rdname2id;
	vplist   *rdhits;
	hzminv   *seeds;
	hzmhash  *hashs[HZMH_KMER_MOD];
	BitVec   *masked;
	BitVec   *needed;
	u64hash  *closed_alns;
	u32list  *rdcovs;
	int hk, hz;
	uint32_t ksize, zsize, n_idx, kwin, kstep, kovl, ksave, ztot, zovl, max_kmer_freq, max_zmer_freq, max_kmer_var, win_rep_cutoff;
	uint32_t ncand, nbest;
	int w, ew, W, M, X, O, E, T;
	int8_t matrix[4*4];
	int off_esti;
	int min_score;
	float min_id;
	uint32_t max_unalign_in_contained;
	uint32_t max_unalign_in_dovetail;
	float skip_contained;
	int best_overlap;
	int debug;
} WTZMO;

WTZMO* init_wtzmo(uint32_t ksize, uint32_t zsize, int w, int W, int M, int X, int O, int E, int T, int min_score, float min_id){
	WTZMO *wt;
	int i;
	wt = malloc(sizeof(WTZMO));
	wt->rdseqs = init_basebank();
	wt->reads  = init_pbreadv(1024);
	wt->rdname2id = init_cuhash(1023);
	wt->seeds = init_hzminv(1024);
	for(i=0;i<HZMH_KMER_MOD;i++){
		wt->hashs[i] = init_hzmhash(1023);
	}
	wt->masked = init_bitvec(1024);
	wt->needed = init_bitvec(1024);
	wt->closed_alns = init_u64hash(1023);
	wt->rdcovs = init_u32list(1024);
	wt->rdhits = init_vplist(1024);
	wt->n_rd = 0;
	wt->n_idx = 1;
	wt->ksize = ksize;
	wt->zsize = zsize;
	wt->kwin  = 500;
	wt->kstep = 0;
	wt->kovl  = 200;
	wt->ksave = 0;
	wt->ztot  = 200;
	wt->zovl  = 100;
	wt->win_rep_cutoff = 100;
	wt->max_kmer_freq = 10000;
	wt->max_zmer_freq = 100;
	wt->max_kmer_var = 5;
	wt->w = w;
	wt->W = W;
	wt->ew = (w + W) / 2;
	wt->M = M;
	wt->X = X;
	wt->O = O;
	wt->E = E;
	wt->T = T;
	for(i=0;i<4*4;i++) wt->matrix[i] = ((i % 4) == (i / 4))? M : X;
	wt->off_esti = 1;
	wt->skip_contained = 1;
	wt->min_score = min_score;
	wt->min_id = min_id;
	wt->ncand = 1000;
	wt->nbest = 50;
	wt->max_unalign_in_contained = 0;
	wt->max_unalign_in_dovetail  = 200;
	wt->best_overlap = 0;
	wt->debug = 0;
	return wt;
}

void free_wtzmo(WTZMO *wt){
	uint32_t i;
	free_basebank(wt->rdseqs);
	for(i=0;i<wt->n_rd;i++) free(wt->reads->buffer[i].rdname);
	free_pbreadv(wt->reads);
	free_cuhash(wt->rdname2id);
	free_u64hash(wt->closed_alns);
	free_bitvec(wt->masked);
	free_bitvec(wt->needed);
	free_u32list(wt->rdcovs);
	for(i=0;i<wt->rdhits->size;i++) free_u64list((u64list*)get_vplist(wt->rdhits, i));
	free_vplist(wt->rdhits);
	for(i=0;i<HZMH_KMER_MOD;i++){
		free_hzmhash(wt->hashs[i]);
	}
	free_hzminv(wt->seeds);
	free(wt);
}

void push_long_read_wtzmo(WTZMO *wt, char *name, int name_len, char *seq, int seq_len){
	char *ptr;
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_pbreadv(wt->reads, (pbread_t){wt->rdseqs->size, seq_len, ptr});
	seq2basebank(wt->rdseqs, seq, seq_len);
	wt->n_rd ++;
}

void set_read_clip_wtzmo(WTZMO *wt, char *name, int coff, int clen){
	pbread_t *rd;
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->rdname2id, name)) == 0xFFFFFFFFU) return;
	rd = ref_pbreadv(wt->reads, pbid);
	if(coff < 0 || coff + clen > (int)rd->rdlen) return;
	rd->rdoff += coff;
	rd->rdlen  = clen;
}

thread_beg_def(midx);
WTZMO *wt;
uint32_t beg, end;
int task;
thread_end_def(midx);

thread_beg_func(midx);
WTZMO *wt;
hzmh_t *u, U;
pbread_t *rd;
uint64_t kmer, krev, kmask, off, kidx;
uint32_t beg, end, pbid, pblen, i, j, ncpu, tidx;
uint8_t b, c, dir;
int exists;
wt = midx->wt;
beg = midx->beg;
end = midx->end;
ncpu = midx->n_cpu;
tidx = midx->t_idx;
kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - wt->ksize) << 1);
memset(&U, 0, sizeof(hzmh_t));
thread_beg_loop(midx);
if(midx->task == 1){
	for(pbid=beg;pbid<end;pbid++){
		if(tidx == 0 && ((pbid - beg) % 1000) == 0){ fprintf(stderr, "\r%u", pbid - beg); fflush(stderr); }
		rd = ref_pbreadv(wt->reads, pbid);
		pblen = rd->rdlen;
		b = 4;
		kmer = 0;
		off = rd->rdoff;
		for(i=j=0;j<pblen;j++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			if(wt->hk && c == b) continue;
			b = c;
			i ++;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->ksize) continue;
			krev = dna_rev_seq(kmer, wt->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			if(wt->ksave && (krev & 0x03) > 1) continue; // skip G and T
			kidx = hzmh_kmer_smear(krev) % HZMH_KMER_MOD;
			if((kidx % ncpu) != tidx) continue;
			U.mer = krev;
			u = prepare_hzmhash(wt->hashs[kidx], U, &exists);
			if(exists){
				if(u->cnt < 0xFFFFU) u->cnt ++;
			} else {
				u->mer = krev;
				u->off = 0;
				u->cnt = 1;
				u->flt = 0;
			}
		}
	}
	if(tidx == 0) fprintf(stderr, "\r%u reads\n", end - beg); fflush(stderr);
} else if(midx->task == 2){
	for(pbid=beg;pbid<end;pbid++){
		if(tidx == 0 && ((pbid - beg) % 1000) == 0){ fprintf(stderr, "\r%u", pbid - beg); fflush(stderr); }
		rd = ref_pbreadv(wt->reads, pbid);
		pblen = rd->rdlen;
		b = 4;
		kmer = 0;
		off = rd->rdoff;
		for(i=j=0;j<pblen;j++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			if(wt->hk && c == b) continue;
			b = c;
			i ++;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->ksize) continue;
			krev = dna_rev_seq(kmer, wt->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			if(wt->ksave && (krev & 0x03) > 1) continue; // skip G and T
			kidx = hzmh_kmer_smear(krev) % HZMH_KMER_MOD;
			if((kidx % ncpu) != tidx) continue;
			U.mer = krev;
			u = get_hzmhash(wt->hashs[kidx], U);
			if(u == NULL || u->flt || u->cnt >= wt->max_kmer_freq) continue; // too high or singleton
			wt->seeds->buffer[u->off + u->cnt] = (hzmin_t){pbid, dir};
			u->cnt ++;
		}
	}
	if(tidx == 0) fprintf(stderr, "\r%u reads\n", pbid - beg); fflush(stderr);
}
thread_end_loop(midx);
thread_end_func(midx);

void index_wtzmo(WTZMO *wt, uint32_t beg, uint32_t end, uint32_t ncpu){
	hzmh_t  *u;
	uint64_t nflt, nrem, none, ktot, off;
	uint32_t ktyp, kavg, i;
	thread_preprocess(midx);
	clear_hzminv(wt->seeds);
	for(i=0;i<HZMH_KMER_MOD;i++){
		clear_hzmhash(wt->hashs[i]);
	}
	thread_beg_init(midx, ncpu);
	midx->wt = wt;
	midx->beg = beg;
	midx->end = end;
	midx->task = 0;
	thread_end_init(midx);
	fprintf(stderr, "[%s] - scanning kmers (%d bp)\n", date(), wt->ksize);
	thread_apply_all(midx, midx->task = 1); // finish task1
	// estimate kmer_freq_cutoff
	if(wt->max_kmer_freq < 2){
		ktot = ktyp = 0;
		for(i=0;i<HZMH_KMER_MOD;i++){
			reset_iter_hzmhash(wt->hashs[i]);
			while((u = ref_iter_hzmhash(wt->hashs[i]))){
				ktot += u->cnt;
			}
			ktyp += wt->hashs[i]->count;
		}
		kavg = ktot / (ktyp + 1);
		if(kavg < 20) kavg = 20;
		wt->max_kmer_freq = kavg * 5;
		fprintf(stderr, "[%s] - high frequency kmer depth is set to %d\n", date(), wt->max_kmer_freq);
	}
	ktot = nrem = none = nflt = ktyp = 0;
	off = 0;
	for(i=0;i<HZMH_KMER_MOD;i++){
		reset_iter_hzmhash(wt->hashs[i]);
		while((u = ref_iter_hzmhash(wt->hashs[i]))){
			ktot += u->cnt;
			if(u->cnt > wt->max_kmer_freq){ u->cnt = wt->max_kmer_freq; nflt ++; }
			{
				u->off = off;
				off += u->cnt;
				u->cnt = 0;
				nrem ++;
			}
		}
		ktyp += wt->hashs[i]->count;
	}
	clear_and_encap_hzminv(wt->seeds, off);
	wt->seeds->size = off;
	kavg = ktot / (ktyp + 1);
	fprintf(stderr, "[%s] - average kmer depth = %d\n", date(), kavg);
	//fprintf(stderr, "[%s] - %llu single kmers (==1)\n", date(), (unsigned long long)none);
	fprintf(stderr, "[%s] - %llu high frequency kmers (>=%d)\n", date(), (unsigned long long)nflt, wt->max_kmer_freq);
	fprintf(stderr, "[%s] - indexing %llu kmers\n", date(), (unsigned long long)nrem);
	thread_apply_all(midx, midx->task = 2); // finish task2
	thread_beg_close(midx);
	thread_end_close(midx);
}


void query_wtzmo(WTZMO *wt, uint32_t pbid, u64list *candidates, hzrefv *refs, u32list *heap, u32list *hzoff){
	hz_ref_t *r, R;
	pbread_t *rd;
	hzmh_t *h;
	hzm_t *p1, P1;
	hzmin_t *p2, P2;
	uint64_t kmer, krev, kmask, off, x1, x2;
	uint32_t i, j, pblen, ol, lst, kidx, id1, id2, cnt, dir;
	uint8_t b, c;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - wt->ksize) << 1);
	rd = ref_pbreadv(wt->reads, pbid);
	pblen = rd->rdlen;
	R.p1 = (hzm_t){0xFFFFFFFFU, 0, 0, 0};
	P2 = (hzmin_t){0x7FFFFFFFU, 0};
	R.p2 = &P2;
	R.beg = R.end = 0;
	p1 = &P1;
	b = 4;
	kmer = 0;
	clear_hzrefv(refs);
	clear_u32list(heap);
	clear_u32list(hzoff);
	off = rd->rdoff;
	for(i=j=0;j<pblen;j++){
		c = bits2bit(wt->rdseqs->bits, off); off ++;
		if(wt->hk && c == b) continue;
		b = c;
		i ++;
		push_u32list(hzoff, j);
		kmer = ((kmer << 2) | b) & kmask;
		if(i < wt->ksize) continue;
		krev = dna_rev_seq(kmer, wt->ksize);
		if(krev == kmer) continue;
		p1->dir  = krev > kmer? 0 : 1;
		krev = krev > kmer? kmer : krev;
		if(wt->ksave && (krev & 0x03) > 1) continue; // skip G and T
		kidx = hzmh_kmer_smear(krev) % HZMH_KMER_MOD;
		h = get_hzmhash(wt->hashs[kidx], (hzmh_t){krev, 0, 0, 0});
		if(h == NULL) continue;
		if(wt->debug > 3){
			fprintf(stdout, "KMER\t%u\t%u\t%u\n", hzoff->buffer[i - wt->ksize], i - wt->ksize, h? h->cnt : 0);
		}
		if(h->flt) continue;
		p1->mer = krev;
		p1->off = hzoff->buffer[i - wt->ksize];
		p1->len = (j + 1 - p1->off > HZM_MAX_SEED_LEN)? HZM_MAX_SEED_LEN : j + 1 - p1->off;
		r = next_ref_hzrefv(refs);
		r->p1 = P1;
		r->beg = h->off;
		r->end = h->off + h->cnt;
		while(r->beg < r->end){
			r->p2 = ref_hzminv(wt->seeds, r->beg ++);
			if(r->p2->rd_id == pbid) continue;
			if(wt->reads->buffer[r->p2->rd_id].rdlen > pblen) continue;
			array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, refs->size - 1, heap_cmp_hz_ref_macro(refs, a, b));
			break;
		}
	}
	id2 = 0xFFFFFFFFU;
	ol = 0;
	lst = 0;
	cnt = 0;
	dir = 0;
	x1 = x2 = 0xFFFFFFFF00000000LLU;
	while(1){
		id1 = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, heap_cmp_hz_ref_macro(refs, a, b));
		if(id1 != 0xFFFFFFFFU){
			r = ref_hzrefv(refs, id1);
		} else r = &R;
		p1 = &r->p1;
		p2 = r->p2;
		if(id2 != p2->rd_id || dir != p2->dir){
			if(id2 != 0xFFFFFFFFU){
				if(ol >= wt->kovl){
					x2 = (((uint64_t)id2) << 32) | ol;
					if(((x1 >> 32) == (x2 >> 32))){ x1 = (x1 & 0xFFFFFFFFU) > (x2 & 0xFFFFFFFFU)? x1 : x2;}
					else if(x1 == 0xFFFFFFFF00000000LLU){ x1 = x2; }
					else {
						if(candidates->size >= wt->ncand) array_heap_remove(candidates->buffer, candidates->size, candidates->cap, uint64_t, 0,
							((a & 0xFFFFFFFFU) > (b & 0xFFFFFFFFU))? 1 : (((a & 0xFFFFFFFFU) < (b & 0xFFFFFFFFU))? -1 : 0));
						array_heap_push(candidates->buffer, candidates->size, candidates->cap, uint64_t, x1,
							((a & 0xFFFFFFFFU) > (b & 0xFFFFFFFFU))? 1 : (((a & 0xFFFFFFFFU) < (b & 0xFFFFFFFFU))? -1 : 0));
						x1 = x2;
					}
				}
				if(wt->debug > 2){
					fprintf(stdout, "FAXIAN\t%s\t%d\t%d\t%c\n", wt->reads->buffer[id2].rdname, ol, cnt, "+-"[dir]);
				}
			}
			id2 = p2->rd_id;
			dir = p2->dir;
			ol = 0;
			lst = 0;
			cnt = 0;
		}
		if(id1 == 0xFFFFFFFFU) break;
		if(p1->off >= lst) ol += p1->len;
		else ol += p1->off + p1->len - lst;
		lst = p1->off + p1->len;
		cnt ++;
		while(r->beg < r->end){
			r->p2 = ref_hzminv(wt->seeds, r->beg ++);
			if(r->p2->rd_id == pbid) continue;
			if(wt->reads->buffer[r->p2->rd_id].rdlen > pblen) continue;
			array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, id1, heap_cmp_hz_ref_macro(refs, a, b));
			break;
		}
	}
	if(x1 != 0xFFFFFFFF00000000LLU){
		if(candidates->size >= wt->ncand) array_heap_remove(candidates->buffer, candidates->size, candidates->cap, uint64_t, 0,
			((a & 0xFFFFFFFFU) > (b & 0xFFFFFFFFU))? 1 : (((a & 0xFFFFFFFFU) < (b & 0xFFFFFFFFU))? -1 : 0));
		array_heap_push(candidates->buffer, candidates->size, candidates->cap, uint64_t, x1,
			((a & 0xFFFFFFFFU) > (b & 0xFFFFFFFFU))? 1 : (((a & 0xFFFFFFFFU) < (b & 0xFFFFFFFFU))? -1 : 0));
	}
}

thread_beg_def(mzmo);
WTZMO *wt;
uint32_t rd_id, bcov;
int do_align, fast_align, just_query;
wtovlv *hits;
wtseedv *seeds, *windows;
u32list *masks;
thread_end_def(mzmo);

thread_beg_func(mzmo);
WTZMO *wt;
hzrefv *refs;
hzmhv *zhash;
hzmv *zseeds;
BitVec *zbits;
wtseedv *seeds, *windows, *windows2;
wt_seed_t *seed, SEED[2];
wt_seed_t *zp;
hzmpv  *cache, *anchors, *anchors2;
wtovlv *hits;
u32hash *masked;
String *cigar_str;
wt_ovl_t HIT, *hit;
kswx_t x;
u32list *heap;
u64list *candidates;
u16list *windeps;
u32list *hzoff;
alnregv *regs;
aln_reg_t REG;
u8list *rdseq, *hzseq, *pb1, *pb2;
u8list *mem_pbseq, *mem_cache[2];
u32list *mem_cigar, *cigar_cache, *cigars;
uint64_t val, *u;
uint32_t i, j, k, dir, id2, x1, x2, x3, x4, ol, pbid, avg, *ptr, bcov, ncand;
int alen, blen, exists;
wt = mzmo->wt;
hits = mzmo->hits;
masked = init_u32hash(1023);
seeds = mzmo->seeds;
windows = mzmo->windows;
windows2 = init_wtseedv(1024);
cache = init_hzmpv(1024);
anchors = init_hzmpv(1024);
anchors2 = init_hzmpv(1024);
hit = &HIT;
pb1 = init_u8list(1024);
pb2 = init_u8list(1024);
cigar_str = init_string(1024);
rdseq = init_u8list(1024);
hzseq = init_u8list(1024);
hzoff = init_u32list(1024);
heap = init_u32list(1024);
refs = init_hzrefv(1024);
zhash = init_hzmhv(1023);
zseeds = init_hzmv(1024);
mem_pbseq = init_u8list(1024);
mem_cache[0] = init_u8list(1024);
mem_cache[1] = init_u8list(1024);
cigars = init_u32list(1024);
mem_cigar = init_u32list(1024);
cigar_cache = init_u32list(1024);
memset(&SEED[0], 0, sizeof(wt_seed_t));
memset(&SEED[1], 0, sizeof(wt_seed_t));
regs = init_alnregv(1024);

thread_beg_loop(mzmo);
pbid = mzmo->rd_id;
alen = wt->reads->buffer[pbid].rdlen;
if(wt->rdhits->size) candidates = (u64list*)get_vplist(wt->rdhits, pbid);
else candidates = init_u64list(1024);
query_wtzmo(wt, pbid, candidates, refs, heap, hzoff);
if(mzmo->just_query) continue;
sort_array(candidates->buffer, candidates->size, uint64_t, (b & 0xFFFFFFFFU) > (a & 0xFFFFFFFFU));
if(wt->debug){
	fprintf(stdout, "=%s\n", wt->reads->buffer[pbid].rdname);
	for(i=0;i<candidates->size;i++){
		id2 = get_u64list(candidates, i) >> 32;
		fprintf(stdout, "CANDIDATE\t%s\t%d\n", wt->reads->buffer[id2].rdname, (uint32_t)get_u64list(candidates, i));
	}
	fprintf(stdout, "Total candidates: %u\n", (unsigned)candidates->size);
}
clear_wtseedv(seeds);
clear_wtseedv(windows);
clear_wtovlv(hits);
clear_hzmpv(cache);
clear_hzmpv(anchors);
windeps = init_u16list(wt->reads->buffer[pbid].rdlen);
zeros_u16list(windeps);
zbits  = init_bitvec(0xFFFFFFFFFFFFFFFFLLU >> ((32 - wt->zsize) << 1));
clear_and_encap_u8list(pb1, alen);
bitseq_basebank(wt->rdseqs, wt->reads->buffer[pbid].rdoff, alen, pb1->buffer);
index_single_read_seeds(pb1->buffer, alen, wt->zsize, wt->hz, 64, zhash, zbits, zseeds, hzoff);
for(i=0;i<candidates->size;i++){
	id2 = get_u64list(candidates, i) >> 32;
	clear_hzmpv(cache);
	blen = wt->reads->buffer[id2].rdlen;
	clear_and_encap_u8list(pb2, blen);
	bitseq_basebank(wt->rdseqs, wt->reads->buffer[id2].rdoff, blen, pb2->buffer);
	query_single_read_seeds(pb2->buffer, blen, wt->zsize, wt->hz, wt->max_kmer_var, zhash, zbits, zseeds, hzoff, cache);
	if(wt->debug > 2){
		fprintf(stdout, "Query\t%s\t%d\n", wt->reads->buffer[id2].rdname, (int)cache->size);
	}
	if(cache->size * wt->zsize < wt->ztot) continue;
	process_hzmps(cache, wt->reads->buffer[id2].rdlen);
	for(dir=0;dir<2;dir++){
		//SEED[dir].pb1 = pbid;
		SEED[dir].pb2 = id2;
		SEED[dir].dir = dir;
		SEED[dir].ovl = 0;
		SEED[dir].closed = 1;
		clear_wtseedv(windows2);
		clear_hzmpv(anchors2);
		if(merge_paired_kmers_window(cache, dir, windows2, anchors2, mem_cache, wt->zsize, wt->kwin, wt->kstep, wt->zovl) == 0) continue;
		SEED[dir].ovl = chaining_wtseedv(pbid, id2, dir, windows2, 0, windows2->size, mem_cache[0], wt->W);
		if(SEED[dir].ovl < wt->ztot) continue;
		SEED[dir].anchors[0] = windows->size;
		for(j=0;j<windows2->size;j++){
			zp = ref_wtseedv(windows2, j);
			if(zp->closed) continue;
			append_array_hzmpv(anchors, anchors2->buffer + zp->anchors[0], zp->anchors[1] - zp->anchors[0]);
			zp->anchors[1] = zp->anchors[1] - zp->anchors[0];
			zp->anchors[0] = anchors->size - zp->anchors[1];
			zp->anchors[1] = anchors->size;
			push_wtseedv(windows, *zp);
			for(k=zp->beg[0];(int)k<zp->end[0];k++) windeps->buffer[k] ++;
		}
		SEED[dir].anchors[1] = windows->size;
		SEED[dir].closed = 0;
	}
	dir = (SEED[0].ovl < SEED[1].ovl);
	if(SEED[dir].ovl >= wt->ztot) push_wtseedv(seeds, SEED[dir]);
}
free_bitvec(zbits);
if(wt->rdhits->size == 0) free_u64list(candidates);
if(wt->debug){
	fprintf(stdout, "TOTAL SEEDS: %u\n", (uint32_t)seeds->size);
}
// Filter repeat-matching only seeds
if(wt->debug > 3){
	for(i=0;(int)i<alen;i++){
		fprintf(stdout, "REP[%d]\t%d\t%c\n", i, windeps->buffer[i], windeps->buffer[i] >= wt->win_rep_cutoff? 'R' : 'U');
	}
}
if(wt->skip_contained) clear_u32hash(masked);
ncand = 0;
for(i=0;i<seeds->size;i++){
	seed = ref_wtseedv(seeds, i);
	if(seed->closed) continue;
	ol = 0;
	for(j=seed->anchors[0];j<seed->anchors[1];j++){
		zp = ref_wtseedv(windows, j);
		if(zp->closed) continue;
		avg = 0;
		for(k=zp->beg[0];(int)k<zp->end[0];k++) avg += windeps->buffer[k];
		avg = avg / (zp->end[0] - zp->beg[0] + 1);
		if(wt->debug > 1){
			fprintf(stdout, "WINREP\t%s\t%d\t%d\t%d\t%d\t%d\n", wt->reads->buffer[seed->pb2].rdname, zp->beg[0], zp->end[0], zp->beg[1], zp->end[1], avg);
		}
		if(avg < wt->win_rep_cutoff){ ol += zp->end[0] - zp->beg[0] + 1; }
	}
	if(ol < wt->ztot){
		seed->closed = 1; ncand ++;
		if(wt->debug){
			fprintf(stdout, "Filter %s, non-repetitive seed size %d\n", wt->reads->buffer[seed->pb2].rdname, ol);
		}
	} else if(wt->debug > 1){
		fprintf(stdout, "Keep %s, non-repetitive seed size %d\n", wt->reads->buffer[seed->pb2].rdname, ol);
	}
}
free_u16list(windeps);
if(wt->debug){
	fprintf(stdout, "Filtered %u repetitive-seq-caused seeds\n", ncand);
}
sort_array(seeds->buffer, seeds->size, wt_seed_t, b.ovl > a.ovl);
thread_beg_syn(mzmo);
for(i=0;i<seeds->size;i++){
	seed = ref_wtseedv(seeds, i);
	if(seed->closed) continue;
	//val = ovl_uniq_long_id(seed->pb2, seed->pb1, seed->dir);
	val = ovl_uniq_long_id(seed->pb2, pbid, seed->dir);
	u = prepare_u64hash(wt->closed_alns, val, &exists);
	if(!exists) *u = val;
	else seed->closed = 2;
}
thread_end_syn(mzmo);
//sort_array(seeds->buffer, seeds->size, wt_seed_t, b.ovl > a.ovl);
//sort_array(seeds->buffer, seeds->size, wt_seed_t, (a.closed > b.closed)? 1 : ((a.closed == b.closed)? (b.ovl > a.ovl) : 0));
//for(i=0;i<seeds->size;i++) if(seeds->buffer[i].closed == 2){ seeds->size = i; break; }
if(mzmo->do_align){
	bcov = mzmo->bcov;
	ncand = wt->ncand;
	for(i=0;i<seeds->size&&i<ncand;i++){
		seed = ref_wtseedv(seeds, i);
		if(seed->closed){ ncand ++; continue; }
		blen = wt->reads->buffer[seed->pb2].rdlen;
		clear_and_encap_u8list(pb2, blen);
		if(seed->dir) revbitseq_basebank(wt->rdseqs, wt->reads->buffer[seed->pb2].rdoff, blen, pb2->buffer);
		else             bitseq_basebank(wt->rdseqs, wt->reads->buffer[seed->pb2].rdoff, blen, pb2->buffer);
		if(wt->debug){
			fprintf(stdout, "DO_ALIGN\t%04d\t%s\t%s\t%c\t%d\n", i, wt->reads->buffer[pbid].rdname, wt->reads->buffer[seed->pb2].rdname, "+-"[seed->dir], seed->ovl);
		}
		clear_u32list(cigar_cache);
		clear_alnregv(regs);
		for(j=seed->anchors[0];j<seed->anchors[1];j++){
			zp = ref_wtseedv(windows, j);
			if(zp->closed) continue;
			REG.cigar_off = cigar_cache->size;
			REG.x = fast_seeds_align_hzmo(pb1->buffer, pb2->buffer, zp, anchors, cigar_cache, mem_cache[0], mem_cigar, wt->w, wt->M, wt->X, wt->O, wt->O, wt->E, wt->T);
			REG.cigar_len = cigar_cache->size - REG.cigar_off;
			if(REG.x.aln * 2 < (int)wt->zovl || REG.x.mat < REG.x.aln * wt->min_id) continue;
			push_alnregv(regs, REG);
		}
		if(regs->size == 0){ seed->closed = 1; ncand ++; continue; }
		x = global_align_regs_hzmo(alen, blen, regs, (int[2]){0, alen}, (uint8_t*[2]){pb1->buffer, pb2->buffer}, cigar_cache->buffer, cigars, mem_cache[0], mem_cigar, wt->W, wt->ew, wt->w, wt->M, wt->X, wt->O, wt->O, wt->E, wt->T);
		clear_alnregv(regs);
		if(wt->debug){
			fprintf(stdout, "OVL\t%s\t%c\t%d\t%d\t%d\t%s\t%c\t%d\t%d\t%d", wt->reads->buffer[pbid].rdname, '+', alen, x.tb, x.te, wt->reads->buffer[seed->pb2].rdname, "+-"[seed->dir], blen, x.qb, x.qe);
			fprintf(stdout, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
		}
		if(x.score < wt->min_score || x.mat < x.aln * wt->min_id) continue;
		clear_string(cigar_str);
		kswx_cigar2string(cigar_str, cigars->size, cigars->buffer);
		//HIT.pb1  = seed->pb1;
		HIT.pb1  = pbid;
		HIT.pb2  = seed->pb2;
		HIT.dir1 = 0;
		HIT.dir2 = seed->dir;
		HIT.score = x.score;
		HIT.tb   = x.tb;
		HIT.te   = x.te;
		HIT.qb   = x.qb;
		HIT.qe   = x.qe;
		HIT.mat  = x.mat;
		HIT.mis  = x.mis;
		HIT.ins  = x.ins;
		HIT.del  = x.del;
		HIT.aln  = x.aln;
		HIT.cigar = strdup(cigar_str->string);
		push_wtovlv(mzmo->hits, HIT);
		{
			x1 = num_min(hit->tb, hit->qb);
			x2 = num_min(((int)wt->reads->buffer[hit->pb1].rdlen) - hit->te, ((int)wt->reads->buffer[hit->pb2].rdlen) - hit->qe);
			if(x1 + x2 <= wt->max_unalign_in_dovetail){
				if(wt->skip_contained){
					hit = &HIT;
					x3 = ((hit->tb == 0 && hit->qb) || (hit->te == (int)wt->reads->buffer[hit->pb1].rdlen && hit->qe < (int)wt->reads->buffer[hit->pb2].rdlen));
					x4 = ((hit->qb == 0 && hit->tb) || (hit->qe == (int)wt->reads->buffer[hit->pb2].rdlen && hit->te < (int)wt->reads->buffer[hit->pb1].rdlen));
					x1 = wt->reads->buffer[hit->pb2].rdlen + hit->qb - hit->qe;
					x2 = wt->reads->buffer[hit->pb1].rdlen + hit->tb - hit->te;
					if(x1 <= wt->max_unalign_in_contained && x3 == 0){
						if(x2 <= wt->max_unalign_in_contained && x4 == 0){
							if(wt->reads->buffer[hit->pb1].rdlen > wt->reads->buffer[hit->pb2].rdlen){
								put_u32hash(masked, hit->pb2);
							} else if(wt->reads->buffer[hit->pb1].rdlen < wt->reads->buffer[hit->pb2].rdlen){
								put_u32hash(masked, hit->pb1); break;
							} else if(hit->pb2 > hit->pb1){
								put_u32hash(masked, hit->pb2); continue;
							} else {
								put_u32hash(masked, hit->pb1); break;
							}
						} else {
							put_u32hash(masked, hit->pb2); continue;
						}
						ncand ++;
					} else if(x2 <= wt->max_unalign_in_contained && x4 == 0){
						put_u32hash(masked, hit->pb1); break;
						ncand ++;
					}
				}
				x1 = num_max(0, HIT.tb);
				x2 = num_min((int)wt->kwin, HIT.te);
				if(x1 < x2) bcov += x2 - x1;
				x3 = num_max(alen - (int)wt->kwin, HIT.tb);
				x4 = num_min(alen, HIT.te);
				if(x3 < x4) bcov += x4 - x3;
				if(wt->debug){
					fprintf(stdout, "COV\t%d\t%d\t%d\t%d\t%d\n", bcov, x1, x2, x3, x4);
				}
				if(bcov >= 2 * wt->kwin * wt->nbest){
					if(wt->debug){
						fprintf(stdout, "COVBREAK\t%d\t%d\n", bcov, 2 * wt->kwin * wt->nbest);
					}
					break;
				}
			}
		}
	}
	//wt->rdcovs->buffer[pbid] = bcov; // thread-unsafe
	mzmo->bcov = bcov; // thread-unsafe
}
//else if(1){ seeds->size = 0; }
if(wt->skip_contained){
	reset_iter_u32hash(masked);
	while((ptr = ref_iter_u32hash(masked))){
		push_u32list(mzmo->masks, *ptr);
	}
}
thread_end_loop(mzmo);
free_u8list(rdseq);
free_u8list(hzseq);
free_u32list(hzoff);
free_string(cigar_str);
free_u8list(pb1);
free_u8list(pb2);
free_u32hash(masked);
free_wtseedv(windows2);
free_hzmpv(anchors);
free_hzmpv(anchors2);
free_hzmpv(cache);
free_u32list(heap);
free_hzrefv(refs);
free_hzmhv(zhash);
free_hzmv(zseeds);
free_u8list(mem_pbseq);
free_u8list(mem_cache[0]);
free_u8list(mem_cache[1]);
free_u32list(cigars);
free_u32list(mem_cigar);
free_u32list(cigar_cache);
free_alnregv(regs);
thread_end_func(mzmo);

//uint32_t print_hits_wtzmo(WTZMO *wt, wtseedv *seeds, wtovlv *hits, FILE *out){
uint32_t print_hits_wtzmo(WTZMO *wt, uint32_t pbid, wtseedv *seeds, wtovlv *hits, FILE *out){
	wt_ovl_t *hit;
	wt_seed_t *seed;
	uint32_t i, j, ret, x1, x2;
	ret = 0;
	if(seeds){
		for(i=j=0;i<seeds->size;i++){
			seed = ref_wtseedv(seeds, i);
			if(seed->closed) continue;
			fprintf(out, "# %s\t%c\t%d", wt->reads->buffer[pbid].rdname, "+-"[0], wt->reads->buffer[pbid].rdlen);
			fprintf(out, "\t%s\t%c\t%d", wt->reads->buffer[seed->pb2].rdname, "+-"[seed->dir], wt->reads->buffer[seed->pb2].rdlen);
			fprintf(out, "\t%d\n", seed->ovl);
			if(j == hits->size) continue;
			hit = ref_wtovlv(hits, j);
			if(hit->score == -19830203){ if(hit->cigar){ free(hit->cigar); hit->cigar = NULL; } j ++; continue; }
			if(hit->pb2 != seed->pb2) continue;
			j ++;
			ret ++;
			if(hit->aln == 0) hit->aln = 1;
			x1 = num_min(hit->tb, hit->qb);
			x2 = num_min(((int)wt->reads->buffer[hit->pb1].rdlen) - hit->te, ((int)wt->reads->buffer[hit->pb2].rdlen) - hit->qe);
			if(x1 + x2 <= wt->max_unalign_in_dovetail){
				x1 = num_max(0, hit->qb);
				x2 = num_min((int)wt->kwin, hit->qe);
				if(x1 < x2) wt->rdcovs->buffer[hit->pb2] += x2 - x1;
				x1 = num_max(((int)(wt->reads->buffer[hit->pb2].rdlen)) - ((int)wt->kwin), hit->qb);
				x2 = num_min(((int)(wt->reads->buffer[hit->pb2].rdlen)), hit->qe);
				if(x1 < x2) wt->rdcovs->buffer[hit->pb2] += x2 - x1;
			}
			fprintf(out, "%s\t%c\t%d\t%d\t%d", wt->reads->buffer[hit->pb1].rdname, "+-"[hit->dir1], wt->reads->buffer[hit->pb1].rdlen, hit->tb, hit->te);
			fprintf(out, "\t%s\t%c\t%d\t%d\t%d", wt->reads->buffer[hit->pb2].rdname, "+-"[hit->dir2], wt->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe);
			fprintf(out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d", hit->score, 1.0 * hit->mat / hit->aln, hit->mat, hit->mis, hit->ins, hit->del);
			if(hit->cigar){
				fprintf(out, "\t%s\n", hit->cigar);
				free(hit->cigar);
				hit->cigar = NULL;
			} else {
				fprintf(out, "\t0M\n");
			}
		}
		clear_wtseedv(seeds);
	} else {
		for(j=0;j<hits->size;j++){
			hit = ref_wtovlv(hits, j);
			if(hit->score == -19830203){ if(hit->cigar){ free(hit->cigar); hit->cigar = NULL; } j ++; continue; }
			ret ++;
			if(hit->aln == 0) hit->aln = 1;
			x1 = num_min(hit->tb, hit->qb);
			x2 = num_min(((int)wt->reads->buffer[hit->pb1].rdlen) - hit->te, ((int)wt->reads->buffer[hit->pb2].rdlen) - hit->qe);
			if(x1 + x2 <= wt->max_unalign_in_dovetail){
				x1 = num_max(0, hit->qb);
				x2 = num_min((int)wt->kwin, hit->qe);
				if(x1 < x2) wt->rdcovs->buffer[hit->pb2] += x2 - x1;
				x1 = num_max(((int)(wt->reads->buffer[hit->pb2].rdlen)) - ((int)wt->kwin), hit->qb);
				x2 = num_min(((int)(wt->reads->buffer[hit->pb2].rdlen)), hit->qe);
				if(x1 < x2) wt->rdcovs->buffer[hit->pb2] += x2 - x1;
			}
			fprintf(out, "%s\t%c\t%d\t%d\t%d", wt->reads->buffer[hit->pb1].rdname, "+-"[hit->dir1], wt->reads->buffer[hit->pb1].rdlen, hit->tb, hit->te);
			fprintf(out, "\t%s\t%c\t%d\t%d\t%d", wt->reads->buffer[hit->pb2].rdname, "+-"[hit->dir2], wt->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe);
			fprintf(out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d", hit->score, 1.0 * hit->mat / hit->aln, hit->mat, hit->mis, hit->ins, hit->del);
			if(hit->cigar){
				fprintf(out, "\t%s\n", hit->cigar);
				free(hit->cigar);
				hit->cigar = NULL;
			} else {
				fprintf(out, "\t0M\n");
			}
		}
	}
	clear_wtovlv(hits);
	return ret;
}

uint64_t overlap_wtzmo(WTZMO *wt, int ncpu, int do_align, int fast_align, uint32_t n_job, uint32_t i_job, FILE *out){
	uint64_t ret;
	uint32_t i, j, rd_id, pbbeg, pbend, i_idx;
	int n_cpu;
	thread_preprocess(mzmo);
	if(wt->debug) n_cpu = 1;
	else n_cpu = ncpu;
	thread_beg_init(mzmo, n_cpu);
	mzmo->wt = wt;
	mzmo->rd_id = 0xFFFFFFFFU;
	mzmo->bcov = 0;
	mzmo->do_align = do_align;
	mzmo->fast_align = fast_align;
	mzmo->just_query = 0;
	mzmo->seeds = init_wtseedv(1024);
	mzmo->windows = init_wtseedv(1024);
	mzmo->hits = init_wtovlv(1024);
	mzmo->masks = init_u32list(16);
	thread_end_init(mzmo);
	ret = 0;

	if(wt->n_idx > 1){
		for(j=0;j<wt->n_rd;j++) push_vplist(wt->rdhits, init_u64list(2));
	}

	pbbeg = pbend = 0;
	for(i_idx=0;i_idx<wt->n_idx;i_idx++){
		pbbeg = pbend;
		pbend = pbbeg + (wt->n_rd + wt->n_idx - 1) / wt->n_idx;
		fprintf(stderr, "[%s] indexing %u/%u\n", date(), i_idx + 1, wt->n_idx);
		index_wtzmo(wt, pbbeg, pbend, ncpu);
		fprintf(stderr, "[%s] Done\n", date());
		fprintf(stderr, "[%s] querying %u/%u\n", date(), i_idx + 1, wt->n_idx);
		if(i_idx + 1 >= wt->n_idx) break;
		for(j=0;j<wt->n_rd;j++){
			rd_id = j;
			if((rd_id % n_job) != i_job) continue;
			if(get_bitvec(wt->masked, rd_id)) continue;
			thread_waitfor_one_idle(mzmo);
			mzmo->just_query = 1;
			mzmo->rd_id = rd_id;
			thread_wake(mzmo);
			if((j % 100) == 0){
				fprintf(stderr, "\r%u", j); fflush(stderr);
			}
		}
		thread_waitfor_all_idle(mzmo);
		fprintf(stderr, "\rprogress: %u 100.00%%\n", wt->n_rd); fflush(stderr);
	}

	for(j=0;j<wt->n_rd;j++){
		rd_id = j;
		if((rd_id % n_job) != i_job) continue;
		if(get_bitvec(wt->masked, rd_id)) continue;
		thread_waitfor_one_idle(mzmo);
		mzmo->just_query = 0;
		if(mzmo->rd_id != 0xFFFFFFFFU) wt->rdcovs->buffer[mzmo->rd_id] = mzmo->bcov;
		ret += print_hits_wtzmo(wt, mzmo->rd_id, do_align? NULL : mzmo->seeds, mzmo->hits, out);
		if(wt->skip_contained){
			thread_beg_syn_write(mzmo);
			for(i=0;i<mzmo->masks->size;i++) one_bitvec(wt->masked, mzmo->masks->buffer[i]);
			thread_end_syn_write(mzmo);
		}
		clear_u32list(mzmo->masks);
		mzmo->rd_id = rd_id;
		mzmo->bcov = wt->rdcovs->buffer[rd_id];
		//mzmo->bcov = 0;
		thread_wake(mzmo);
		if((j % 100) == 0){
			fprintf(stderr, "\r%012u\t%llu", j, (unsigned long long)ret); fflush(stderr);
		}
	}
	thread_waitfor_all_idle(mzmo);
	thread_beg_close(mzmo);
	if(mzmo->rd_id != 0xFFFFFFFFU) wt->rdcovs->buffer[mzmo->rd_id] = mzmo->bcov;
	ret += print_hits_wtzmo(wt, mzmo->rd_id, do_align? NULL : mzmo->seeds, mzmo->hits, out);
	if(wt->skip_contained){
		thread_beg_syn_write(mzmo);
		for(i=0;i<mzmo->masks->size;i++) one_bitvec(wt->masked, mzmo->masks->buffer[i]);
		thread_end_syn_write(mzmo);
	}
	free_wtovlv(mzmo->hits);
	free_wtseedv(mzmo->seeds);
	free_wtseedv(mzmo->windows);
	free_u32list(mzmo->masks);
	thread_end_close(mzmo);
	fprintf(stderr, "\rprogress: %u %llu 100.00%%\n", wt->n_rd, (unsigned long long)ret); fflush(stderr);
	return ret;
}

uint64_t online_overlap_wtzmo(WTZMO *wt, int ncpu, int do_align, int fast_align, uint32_t n_job, uint32_t i_job, FileReader *fr, FILE *out){
	uint64_t ret;
	uint32_t i, j, rd_id;
	thread_preprocess(mzmo);
	fprintf(stderr, "[%s] indexing whole reads in online mode\n", date());
	index_wtzmo(wt, 0, wt->n_rd, ncpu);
	fprintf(stderr, "[%s] Done\n", date());
	if(wt->debug) ncpu = 1;
	thread_beg_init(mzmo, ncpu);
	mzmo->wt = wt;
	mzmo->rd_id = 0xFFFFFFFFU;
	mzmo->bcov = 0;
	mzmo->do_align = do_align;
	mzmo->fast_align = fast_align;
	mzmo->just_query = 0;
	mzmo->seeds = init_wtseedv(1024);
	mzmo->windows = init_wtseedv(1024);
	mzmo->hits = init_wtovlv(1024);
	mzmo->masks = init_u32list(16);
	thread_end_init(mzmo);
	ret = 0;
	j = 0;
	while(fread_line(fr->line, fr) != -1){
		if(fr->line->string[0] == '#') continue;
		rd_id = kv_get_cuhash(wt->rdname2id, fr->line->string);
		if(rd_id == 0xFFFFFFFFU) continue;
		j ++;
		if((rd_id % n_job) != i_job) continue;
		if(get_bitvec(wt->masked, rd_id)) continue;
		thread_waitfor_one_idle(mzmo);
		if(mzmo->rd_id != 0xFFFFFFFFU) wt->rdcovs->buffer[mzmo->rd_id] = mzmo->bcov;
		ret += print_hits_wtzmo(wt, mzmo->rd_id, mzmo->seeds, mzmo->hits, out);
		if(wt->skip_contained){
			thread_beg_syn_write(mzmo);
			for(i=0;i<mzmo->masks->size;i++) one_bitvec(wt->masked, mzmo->masks->buffer[i]);
			thread_end_syn_write(mzmo);
		}
		clear_u32list(mzmo->masks);
		mzmo->rd_id = rd_id;
		mzmo->bcov = wt->rdcovs->buffer[rd_id];
		//mzmo->bcov = 0;
		thread_wake(mzmo);
	}
	thread_waitfor_all_idle(mzmo);
	thread_beg_close(mzmo);
	if(mzmo->rd_id != 0xFFFFFFFFU) wt->rdcovs->buffer[mzmo->rd_id] = mzmo->bcov;
	ret += print_hits_wtzmo(wt, mzmo->rd_id, mzmo->seeds, mzmo->hits, out);
	if(wt->skip_contained){
		thread_beg_syn_write(mzmo);
		for(i=0;i<mzmo->masks->size;i++) one_bitvec(wt->masked, mzmo->masks->buffer[i]);
		thread_end_syn_write(mzmo);
	}
	free_wtovlv(mzmo->hits);
	free_wtseedv(mzmo->seeds);
	free_wtseedv(mzmo->windows);
	free_u32list(mzmo->masks);
	thread_end_close(mzmo);
	fprintf(stderr, "\rprogress: %u %llu 100.00%%\n", j, (unsigned long long)ret); fflush(stderr);
	return ret;
}

int usage(){
	printf(
	"WTZMO: Overlaper of long reads using homopolymer compressed k-mer seeding\n"
	"WatchTower: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtzmo [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -P <int>    Total parallel jobs, [1]\n"
	" -p <int>    Index of current job (0-based), [0]\n"
	"             Suppose to run wtzmo parallelly in 60 nodes. For node1, -P 60 -p 0; node2, -P 60 -p 1, ...\n"
	//" -N          DONOT perform pairwise smith-waterman alignment, just output potential pairs\n"
	" -i <string> Long reads sequences file, + *\n"
	" -b <string> Long reads retained region, often from wtobt/wtcyc, +\n"
	"             Format: read_name\\toffset\\tlength\\toriginal_len\n"
	" -J <int>    Jack knife of original read length, [0]\n"
	//" -x <string> Long reads index file from wtidx.\n"
	//"             If provided, will load the index and skip the reading and indexing process of long reads\n"
	//"             Index file will be shared (mmap MAP_SHARED) by all kind of those processes\n"
	" -I <string> Only find the overlaps of reads in -I file, -: stdin, +\n"
	" -o <string> Output file of alignments, *\n"
	" -f          Force overwrite\n"
	" -H <int>    Option of homopolymer compression, [3]\n"
	"             1: trun on compression on kmer\n"
	"             2: trun on compression on small-kmer(zmer)\n"
	" -k <int>    Kmer size, 5 <= <-k> <= %d, [16]\n"
	" -K <int>    Filter high frequency kmers, maybe repetitive, [0]\n"
	"             0: set K to 5 * <average_kmer_depth>, but no less than 100\n"
	" -d <int>    Minimum size of total seeding region for kmer windows, [300]\n"
	" -S          Trun off memory saving of kmer index.\n"
	"             Default, skips kmers ending with 'G' and 'T', halve the memory\n"
	" -G <int>    Build kmer index in multiple iterations to save memory, 1: once, [1]\n"
	"             Given 10M reads having 100G bases, about 100/(4)=25G used in seq storage, about 100*(6)G=600G\n"
	"             used in kmer-index. If -G = 10, kmer-index is divided into 10 pieces, thus taking 60G. But we need additional\n"
	"             10M / <tot_jobs: -P> * 8 * <num_of_cand: -A> memory to store candidates to be aligned.\n"
	" -z <int>    Smaller kmer size (z-mer), 5 <= <-z> <= %d, [10]\n"
	" -Z <int>    Filter high frequency z-mers, maybe repetitive, [100]\n"
	" -y <int>    Zmer window, [800]\n"
	" -R <int>    Minimum size of seeding region within zmer window, [200]\n"
	" -r <int>    Minimum size of total seeding region for zmer windows, [300]\n"
	" -l <int>    Maximum variant of uncompressed sizes between two matched hz-kmer, [2]\n"
	" -q <int>    Threshold of seeds-window coverage of repetitve region in long reads, [100]\n"
	"             One repetitive block in long read has too many zmer-window matched, say <win_tot_len> >= 100 * <block_len>\n"
	" -A <int>    Limit number of best candidates per read, [500]\n"
	" -B <int>    Limit coverage(total_overlap_length / read_length) of best alignment per read, [100]\n"
	"             It is tested on read ends with kmer window size\n"
	" -C          Don't skip calculation of its overlaps even when the read was contained by others\n"
	" -L <string> Already existing alignments, will avoid to calculate them again, + [NULL]\n"
	" -F <string> Reads from this file(s) are to be exclued, one line for one read name, + [NULL]\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: insertion or deletion, [-2]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping, 0: distable HSP extension, otherwise set to -50 or other [-50]\n"
	" -w <int>    Minimum bandwidth, iteratively doubled to maximum [50]\n"
	" -W <int>    Maximum bandwidth, [3200]\n"
	" -e <int>    Maximum bandwidth at ending extension, [800]\n"
	" -s <int>    Minimum alignment score, [200]\n"
	" -m <float>  Minimum alignment identity, [0.5]\n"
	" -v          Verbose, BE careful, HUGEEEEEEEE output on STDOUT\n"
	"\n"
	"Example:\n"
	"$> wtzmo -t 32 -i wt.fa -o wt.zmo.ovl -m 0.6\n"
	"\n", HZM_MAX_SEED_KMER, HZM_MAX_SEED_ZMER
	);
	return 1;
}

int main(int argc, char **argv){
	WTZMO *wt;
	cplist *pbs, *flts, *ovls, *obts, *tbas;
	FileReader *fr;
	Sequence *seq;
	char *output;
	FILE *out;
	uint64_t val, pb1, pb2, naln;
	uint32_t pbid;
	int c, ncpu, min_rdlen, w, W, ew, M, X, O, E, T, hk, hz, ksize, zsize, kwin, kstep, ksave, ztot, zovl, kovl, kcut, zcut, kvar, wrep, min_score, ncand, nbest, overwrite, debug, n_job, i_job, n_idx, best_overlap;
	float min_id, skip_contained, do_align, fast_align;
	HZM_FAST_WINDOW_KMER_CHAINING = 1;
	output = NULL;
	ncpu = 1;
	w = 50;
	ew = 800;
	W = 3200;
	M = 2;
	X = -5;
	O = -2;
	E = -1;
	T = -50;
	min_rdlen = 0;
	min_score = 200;
	min_id = 0.5;
	hk = 1;
	hz = 1;
	ksize = 16;
	zsize = 10;
	kwin = 800;
	kstep = 0;
	kovl = 300;
	ksave = 1;
	n_idx = 1;
	wrep = 100;
	ncand = 500;
	nbest = 100;
	ztot = 300;
	zovl = 200;
	kcut = 0;
	zcut = 100;
	kvar = 2;
	skip_contained = 1;
	overwrite = 0;
	do_align = 1;
	fast_align = 0;
	best_overlap = 0;
	debug = 0;
	n_job = 1;
	i_job = 0;
	pbs = init_cplist(4);
	flts = init_cplist(4);
	ovls = init_cplist(4);
	obts = init_cplist(4);
	tbas = init_cplist(4);
	while((c = getopt(argc, argv, "ht:P:p:Ni:b:J:I:o:SfCH:k:G:z:Z:y:d:r:q:l:K:A:B:r:R:L:F:W:w:e:M:X:O:E:T:s:m:v")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'P': n_job = atoi(optarg); break;
			case 'p': i_job = atoi(optarg); break;
			case 'N': do_align = 0; break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			case 'J': min_rdlen = atoi(optarg); break;
			case 'I': push_cplist(tbas, optarg); break;
			case 'o': output = optarg; break;
			case 'S': ksave = 0; break;
			case 'f': overwrite = 1; break;
			case 'C': skip_contained = 0; break;
			case 'H': hk = atoi(optarg); hz = (hk >> 1) & 0x01; hk = hk & 0x01; break;
			case 'k': ksize = atoi(optarg); break;
			case 'K': kcut = atoi(optarg); break;
			case 'z': zsize = atoi(optarg); break;
			case 'Z': zcut = atoi(optarg); break;
			case 'y': kwin = atoi(optarg); break;
			case 'l': kvar = atoi(optarg); break;
			case 'd': kovl = atof(optarg); break;
			case 'G': n_idx = atoi(optarg); break;
			case 'r': ztot = atof(optarg); break;
			case 'R': zovl = atof(optarg); break;
			case 'q': wrep = atoi(optarg); break;
			case 'A': ncand = atoi(optarg); break;
			case 'B': nbest = atoi(optarg); break;
			case 'w': w = atoi(optarg); break;
			case 'e': ew = atoi(optarg); break;
			case 'W': W = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			case 'L': push_cplist(ovls, optarg); break;
			case 'F': push_cplist(flts, optarg); break;
			case 's': min_score = atoi(optarg); break;
			case 'm': min_id    = atof(optarg); break;
			case 'v': debug ++; break;
			default: return usage();
		}
	}
	if(output == NULL) return usage();
	if(!overwrite && strcmp(output, "-") && file_exists(output)){
		fprintf(stderr, "File exists! '%s'\n\n", output);
		return usage();
	}
	if(pbs->size == 0) return usage();
	if(ksize > HZM_MAX_SEED_KMER || ksize < 5) return usage();
	if(zsize > HZM_MAX_SEED_ZMER || zsize < 5) return usage();
	wt = init_wtzmo(ksize, zsize, w, W, M, X, O, E, T, min_score, min_id);
	wt->hk = hk;
	wt->hz = hz;
	wt->kwin = kwin;
	wt->kstep = kstep = kwin / 2;
	wt->ew = ew;
	wt->kovl = kovl;
	wt->ztot = ztot;
	wt->zovl = zovl;
	wt->ksave = ksave;
	wt->n_idx = n_idx;
	wt->win_rep_cutoff = wrep;
	wt->max_kmer_freq = kcut;
	wt->max_zmer_freq = zcut;
	wt->max_kmer_var = kvar;
	wt->ncand = ncand;
	wt->nbest = nbest;
	wt->debug = debug;
	hzm_debug = debug;
	if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	fprintf(stderr, "[%s] loading long reads\n", date());
	seq = NULL;
	while(fread_seq(&seq, fr)){
		if(seq->seq.size < min_rdlen) continue;
		push_long_read_wtzmo(wt, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
		if((wt->n_rd % 1000) == 0){
			fprintf(stderr, "\r%u", wt->n_rd); fflush(stderr);
		}
	}
	fclose_filereader(fr);
	encap_bitvec(wt->masked, wt->n_rd); zeros_bitvec(wt->masked);
	encap_bitvec(wt->needed, wt->n_rd); zeros_bitvec(wt->needed);
	fprintf(stderr, "\r[%s] Done, %u reads (length >= %d)\n", date(), (unsigned)wt->n_rd, min_rdlen);
	{
		sort_array(wt->reads->buffer, wt->reads->size, pbread_t, b.rdlen > a.rdlen);
		for(pbid=0;pbid<wt->n_rd;pbid++){
			kv_put_cuhash(wt->rdname2id, wt->reads->buffer[pbid].rdname, pbid);
		}
		fprintf(stderr, "[%s] sorted sequences by length dsc\n", date());
	}
	if(obts->size){
		fprintf(stderr, "[%s] loading reads obt information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(fr->line->string[0] == '#') continue;
			if(c < 3) continue;
			set_read_clip_wtzmo(wt, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date());
	}
	if(flts->size){
		if((fr = fopen_m_filereader(flts->size, flts->buffer)) == NULL) exit(1);
		fprintf(stderr, "[%s] loading read filter list\n", date());
		c = 0;
		while(fread_line(fr->line, fr) != -1){
			if(fr->line->string[0] == '#') continue;
			pbid = kv_get_cuhash(wt->rdname2id, fr->line->string);
			if(pbid == 0xFFFFFFFFU) continue;
			if(get_bitvec(wt->masked, pbid)) continue;
			one_bitvec(wt->masked, pbid);
			c ++;
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] masked %d reads\n", date(), c);
	}
	if(ovls->size){
		if((fr = fopen_m_filereader(ovls->size, ovls->buffer)) == NULL) exit(1);
		fprintf(stderr, "[%s] loading already existing alignments\n", date());
		naln = 0;
		while(fread_table(fr) != -1){
			if(fr->line->string[0] == '#') continue;
			if((naln % 10000) == 0){ fprintf(stderr, "\r%llu", (unsigned long long)naln); fflush(stderr); }
			naln ++;
			if((pb1 = kv_get_cuhash(wt->rdname2id, get_col_str(fr, 0))) == 0xFFFFFFFFU) continue;
			if((pb2 = kv_get_cuhash(wt->rdname2id, get_col_str(fr, 5))) == 0xFFFFFFFFU) continue;
			int dir, tb, te, qb, qe, tlen, qlen;
			float iden;
			dir = (get_col_str(fr, 1)[0] != get_col_str(fr, 6)[0]);
			//val = ovl_uniq_long_id(pb1, pb2, dir);
			val = ovl_uniq_long_id(pb1, pb2, 0);
			put_u64hash(wt->closed_alns, val);
			if(skip_contained){
				iden = atof(get_col_str(fr, 11));
				if(iden < skip_contained) continue;
				tlen = atoi(get_col_str(fr, 2));
				tb = atoi(get_col_str(fr, 3));
				te = atoi(get_col_str(fr, 4));
				qlen = atoi(get_col_str(fr, 7));
				qb = atoi(get_col_str(fr, 8));
				qe = atoi(get_col_str(fr, 9));
				if(tb + tlen - te < 100){
					if(qb + qlen - qe < 100){
						one_bitvec(wt->masked, pb1 > pb2? pb1 : pb2);
					} else {
						one_bitvec(wt->masked, pb1);
					}
				} else if(qb + qlen - qe < 100){
					one_bitvec(wt->masked, pb2);
				}
			}
		}
		fclose_filereader(fr);
		fprintf(stderr, "\r[%s] there were %llu existing alignments\n", date(), (unsigned long long)wt->closed_alns->count);
	}
	clear_and_encap_u32list(wt->rdcovs, wt->n_rd);
	zeros_u32list(wt->rdcovs);
	out = strcmp(output, "-")? fopen(output, "w") : stdout;
	fprintf(stderr, "[%s] calculating overlaps, %d threads\n", date(), ncpu);
	if(tbas->size){
		fprintf(stderr, "[%s] loading reads required to align against all from specified file(s)\n", date());
		if((fr = fopen_m_filereader(tbas->size, tbas->buffer)) == NULL) exit(1);
		online_overlap_wtzmo(wt, ncpu, do_align, fast_align, n_job, i_job, fr, out);
		fclose_filereader(fr);
	} else overlap_wtzmo(wt, ncpu, do_align, fast_align, n_job, i_job, out);
	if(strcmp(output, "-")) fclose(out);
	fprintf(stderr, "[%s] Done\n", date());
	if(skip_contained && strcmp(output, "-")){
		char *maskf;
		maskf = catstr(2, output, ".contained");
		out = fopen(maskf, "w");
		uint32_t i;
		for(i=0;i<wt->n_rd;i++){
			if(get_bitvec(wt->masked, i) == 0) continue;
			fprintf(out, "%s\n", (char*)wt->reads->buffer[i].rdname);
		}
		fclose(out);
		free(maskf);
	}
	free_wtzmo(wt);
	free_cplist(pbs);
	free_cplist(flts);
	free_cplist(ovls);
	free_cplist(obts);
	free_cplist(tbas);
	return 0;
}

