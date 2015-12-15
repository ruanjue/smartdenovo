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
#include "string.h"
#include "dna.h"
#include "bitvec.h"
#include "hashset.h"
#include "counting_bloom_filter.h"
#include "heap.h"
#include "kswx.h"
#include "file_reader.h"
#include "thread.h"

#define DBG_MAX_KSIZE	31
#define DBG_MAX_BF_KSIZE	128
#define DBG_MAX_KLINK	0xFF

static int verbose = 0;
static char *debug_opt1 = NULL;
static char *debug_opt2 = NULL;

typedef struct {
	uint64_t val;
	uint8_t links[4+4];
} kmer_t;
#define dbg_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7) ^ ((K) >> 12))
#define dbg_kmer_hc(E) u64hashcode((E).val)
#define dbg_kmer_he(E1, E2) (E1).val == (E2).val
define_hashset(kmerhash, kmer_t, dbg_kmer_hc, dbg_kmer_he);

typedef struct {
	kmerhash **hashs;
	CBF      **bfs; // max count is 3
	uint32_t ksize, bf_ksize, nh; // nh = n_hash
	uint64_t kmask;
	uint8_t  min_link_cov;
} DBG;
size_t _dbg_n_hash_count(void *obj, int idx){ if(idx == 0) return ((DBG*)obj)->nh; else if(idx == 1) return ((DBG*)obj)->nh; else return 1; }
static const struct obj_desc_t DBG_obj_desc = (obj_desc_t){sizeof(DBG), 2, {3, 3}, {offsetof(DBG, hashs), offsetof(DBG, bfs)}, {&kmerhash_obj_desc, &cbf_obj_desc}, _dbg_n_hash_count};

DBG* init_dbg(uint32_t ksize, uint32_t bf_ksize, uint32_t ncpu, uint32_t min_link_cov){
	DBG *g;
	uint32_t i;
	if(ncpu < 1) ncpu = 1;
	g = calloc(1, sizeof(DBG));
	g->nh = ncpu;
	g->hashs = malloc(sizeof(kmerhash*) * ncpu);
	for(i=0;i<ncpu;i++) g->hashs[i] = init_kmerhash(1023);
	g->bfs = NULL;
	if((ksize & 0x01) == 0){
		fprintf(stderr, " -- Automaticly change kmer_size from %u to %u in %s -- %s:%d --\n", ksize, ksize + 1, __FUNCTION__, __FILE__, __LINE__);
		ksize ++;
	}
	if(ksize > DBG_MAX_KSIZE){
		fprintf(stderr, " -- Automaticly change kmer_size from %u to %u in %s -- %s:%d --\n", ksize, DBG_MAX_KSIZE, __FUNCTION__, __FILE__, __LINE__);
		ksize = DBG_MAX_KSIZE;
	}
	g->ksize = ksize;
	if(bf_ksize){
		if((bf_ksize & 0x01) == 0){
			fprintf(stderr, " -- Automaticly change bf_kmer_size from %u to %u in %s -- %s:%d --\n", bf_ksize, bf_ksize + 1, __FUNCTION__, __FILE__, __LINE__);
			bf_ksize ++;
		}
		if(bf_ksize > DBG_MAX_BF_KSIZE){
			fprintf(stderr, " -- Automaticly change bf_kmer_size from %u to %u in %s -- %s:%d --\n", bf_ksize, DBG_MAX_BF_KSIZE - 1, __FUNCTION__, __FILE__, __LINE__);
			bf_ksize = DBG_MAX_BF_KSIZE - 1;
		}
		if(bf_ksize < ksize){
			fprintf(stderr, " -- bf_kmer_size (%u) is no larger than kmer_size (%u), skip to build bloom-filter in %s -- %s:%d --\n", bf_ksize, ksize, __FUNCTION__, __FILE__, __LINE__);
			bf_ksize = 0;
		}
	}
	g->bf_ksize = bf_ksize;
	g->kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
	g->min_link_cov = min_link_cov;
	return g;
}

void free_dbg(DBG *g){
	uint32_t i;
	for(i=0;i<g->nh;i++) free_kmerhash(g->hashs[i]);
	free(g->hashs);
	if(g->bfs){
		for(i=0;i<g->nh;i++) free_cbf(g->bfs[i]);
		free(g->bfs);
	}
	free(g);
}

// Thread: reading sequence
thread_beg_def(mseq);
BaseBank *rdseqs;
u32list  *rdlens;
u8list *seqs;
uint64_t max;
FileReader *fr;
Sequence *seq;
thread_end_def(mseq);

thread_beg_func(mseq);
unsigned long long nseq, tseq, toff;
uint32_t nlen;
int i;
tseq = toff = 0;
mseq->seq = NULL;
thread_beg_loop(mseq);
clear_u8list(mseq->seqs);
nseq = 0;
if(mseq->fr){
	while(fread_seq_adv((Sequence**)&mseq->seq, mseq->fr, SEQ_FLAG_NO_NAME | SEQ_FLAG_NO_QUAL)){
		nseq ++;
		push_u8list(mseq->seqs, 4);
		for(i=0;i<mseq->seq->seq.size;i++) push_u8list(mseq->seqs, base_bit_table[(int)mseq->seq->seq.string[i]]);
		if(mseq->rdseqs){
			seq2basebank(mseq->rdseqs, mseq->seq->seq.string, mseq->seq->seq.size);
			push_u32list(mseq->rdlens, mseq->seq->seq.size);
		}
		if(mseq->seqs->size >= mseq->max) break;
	}
} else if(mseq->rdseqs){
	while(tseq + nseq < mseq->rdlens->size){
		push_u8list(mseq->seqs, 4);
		nlen = mseq->rdlens->buffer[tseq + (nseq ++)];
		encap_u8list(mseq->seqs, nlen);
		bitseq_basebank(mseq->rdseqs, toff, nlen, mseq->seqs->buffer + mseq->seqs->size);
		mseq->seqs->size += nlen;
		toff += nlen;
		if(mseq->seqs->size >= mseq->max) break;
	}
}
if(mseq->seqs->size) push_u8list(mseq->seqs, 4);
tseq += nseq;
if(verbose){ fprintf(stderr, "[%s] + %llu = %llu reads.\n", date(), nseq, tseq); fflush(stderr); }
thread_end_loop(mseq);
if(mseq->seq) free_sequence(mseq->seq);
thread_end_func(mseq);

// Thread: hashing kmers
thread_beg_def(midx);
DBG *g;
u8list *seqs;
thread_end_def(midx);

thread_beg_func(midx);
kmerhash *hash;
kmer_t *k, KMER;
uint64_t i, j, kmer;
uint32_t tidx;
int exists;
uint8_t cur, lst, nxt;
tidx = midx->t_idx;
hash = midx->g->hashs[tidx];
memset(&KMER, 0, sizeof(kmer_t));
thread_beg_loop(midx);
kmer = 0;
for(i=j=0;i<midx->seqs->size;i++){
	cur = midx->seqs->buffer[i];
	if(cur == 4){ j = 0; kmer = 0; continue; }
	kmer = ((kmer << 2) | cur) & midx->g->kmask;
	j ++;
	if(j < midx->g->ksize) continue;
	lst = j == midx->g->ksize? 8 : 3 - midx->seqs->buffer[i - midx->g->ksize];
	nxt = midx->seqs->buffer[i + 1];
	if(nxt == 4) nxt = 8;
	KMER.val = dna_rev_seq(kmer, midx->g->ksize);
	if(KMER.val < kmer){
		nxt = nxt == 8? 8 : 4 + nxt;
	} else {
		KMER.val = kmer;
		lst = lst == 8? 8 : 4 + lst;
	}
	if((dbg_kmer_smear(KMER.val) % midx->g->nh) != tidx) continue;
	k = prepare_kmerhash(hash, KMER, &exists);
	if(!exists) *k = KMER;
	if(lst != 8 && k->links[lst] < DBG_MAX_KLINK) k->links[lst] ++;
	if(nxt != 8 && k->links[nxt] < DBG_MAX_KLINK) k->links[nxt] ++;
}
thread_end_loop(midx);
thread_end_func(midx);

// Thread: bf-hashing kmers
thread_beg_def(mbf);
DBG *g;
u8list *seqs;
thread_end_def(mbf);

thread_beg_func(mbf);
CBF *bf;
uint64_t i, j, *k, kmer[DBG_MAX_BF_KSIZE / 32], krev[DBG_MAX_BF_KSIZE / 32];
uint32_t tidx, bfn;
uint8_t cur;
tidx = mbf->t_idx;
bf = mbf->g->bfs[tidx];
bfn = (mbf->g->bf_ksize + 31) / 32 * 8;
thread_beg_loop(mbf);
memset(kmer, 0, bfn);
for(i=j=0;i<mbf->seqs->size;i++){
	cur = mbf->seqs->buffer[i];
	if(cur == 4){ j = 0; memset(kmer, 0, bfn); continue; }
	dna_shl_seqs(kmer, mbf->g->bf_ksize, cur);
	j ++;
	if(j < mbf->g->ksize) continue;
	memcpy(krev, kmer, bfn);
	dna_rev_seqs(krev, mbf->g->bf_ksize);
	if(dna_cmp_seqs(krev, kmer, mbf->g->bf_ksize) == -1){
		k = krev;
	} else {
		k = kmer;
	}
	if((dbg_kmer_smear(k[0]) % mbf->g->nh) != tidx) continue;
	put_cbf(bf, k, bfn);
}
thread_end_loop(mbf);
thread_end_func(mbf);

void build_dbg(DBG *g, uint64_t seq_batch_size, FileReader *fr){
	BaseBank *rdseqs;
	u32list  *rdlens;
	u8list *seqs[2];
	kmer_t KMER;
	size_t size, totk;
	uint32_t i, seq_idx;
	thread_preprocess(mseq);
	thread_preprocess(midx);
	thread_preprocess(mbf);
	rdseqs = init_basebank();
	rdlens = init_u32list(1024);
	memset(&KMER, 0, sizeof(kmer_t));
	seqs[0] = init_u8list(1024);
	seqs[1] = init_u8list(1024);
	seq_idx = 0;
	thread_beg_init(mseq, 1);
	mseq->rdseqs = g->bf_ksize? rdseqs : NULL;
	mseq->rdlens = g->bf_ksize? rdlens : NULL;
	mseq->seqs = seqs[seq_idx];
	mseq->max  = seq_batch_size;
	mseq->fr   = fr;
	thread_end_init(mseq);
	thread_beg_operate(mseq, 0);
	thread_wake(mseq);
	thread_beg_init(midx, g->nh);
	midx->g = g;
	midx->seqs = seqs[seq_idx];
	thread_end_init(midx);
	while(1){
		thread_wait_all(midx);
		thread_wait(mseq);
		if(mseq->seqs->size == 0) break;
		thread_beg_iter(midx);
		midx->seqs = seqs[seq_idx];
		thread_wake(midx);
		thread_end_iter(midx);
		seq_idx = !seq_idx;
		mseq->seqs = seqs[seq_idx];
		thread_wake(mseq);
	}
	thread_beg_close(midx);
	thread_end_close(midx);
	thread_beg_close(mseq);
	thread_end_close(mseq);
	if(g->bf_ksize == 0) return;
	size = mem_size_obj(g->hashs, 3, &kmerhash_obj_desc, 0, g->nh);
	for(i=totk=0;i<g->nh;i++) totk += g->hashs[i]->count;
	if(verbose){
		fprintf(stdout, "[%s] Total size of kmer hash is %llu bytes\n", date(), (unsigned long long)size); fflush(stdout);
		fprintf(stdout, "[%s] Total number of kmers is %llu\n", date(), (unsigned long long)totk); fflush(stdout);
	}
	// bits = 1.44 * log2(1/FP), here we set FP = 0.001, bits ~ 14
	// The total size of CBF is 14 * kmer_count, here we assume the number of small-kmers as number of big-kmers
	size = 14 * (totk / g->nh);
	size = _rj_hashset_find_prime(size);
	if(verbose){
		fprintf(stdout, "[%s] Total memory of counting bloom filter is %llu bytes\n", date(), (unsigned long long)size * g->nh / 4); fflush(stdout);
	}
	g->bfs = malloc(sizeof(CBF*) * g->nh);
	for(i=0;i<g->nh;i++) g->bfs[i] = init_cbf(size, 2, 3);

	thread_beg_init(mseq, 1);
	mseq->rdseqs = g->bf_ksize? rdseqs : NULL;
	mseq->rdlens = g->bf_ksize? rdlens : NULL;
	mseq->seqs = seqs[seq_idx];
	mseq->max  = seq_batch_size;
	mseq->fr   = NULL;
	thread_end_init(mseq);
	thread_beg_operate(mseq, 0);
	thread_wake(mseq);
	thread_beg_init(mbf, g->nh);
	mbf->g = g;
	mbf->seqs = seqs[seq_idx];
	thread_end_init(mbf);
	while(1){
		thread_wait_all(mbf);
		thread_wait(mseq);
		if(mseq->seqs->size == 0) break;
		thread_beg_iter(mbf);
		mbf->seqs = seqs[seq_idx];
		thread_wake(mbf);
		thread_end_iter(mbf);
		seq_idx = !seq_idx;
		mseq->seqs = seqs[seq_idx];
		thread_wake(mseq);
	}
	thread_beg_close(mbf);
	thread_end_close(mbf);
	thread_beg_close(mseq);
	thread_end_close(mseq);
	free_u8list(seqs[0]);
	free_u8list(seqs[1]);
	free_basebank(rdseqs);
	free_u32list(rdlens);
}

void chk_sequences_dbg(DBG *g, char *tag, char *seq, uint32_t seqlen, uint8_t min_link_cov, uint8_t min_bf_cov, FILE *out){
	kmer_t *k, KMER;
	uint64_t kmer, krev;
	uint64_t *bk, bkmer[DBG_MAX_BF_KSIZE / 32], bkrev[DBG_MAX_BF_KSIZE / 32];
	uint32_t i, j, bfn, hidx;
	uint8_t cur, lst;
	kmer = 0;
	bfn = (g->bf_ksize + 31) / 32 * 8;
	memset(&KMER, 0, sizeof(kmer_t));
	if(bfn) memset(bkmer, 0, bfn);
	for(i=j=0;i<=seqlen;i++){
		if((cur = base_bit_table[(int)seq[i]]) == 4){
			if(j > g->ksize) fprintf(out, "%s\t%d\t%d\t%d\n", tag, i - j, i - 1, j);
			j = 0; kmer = 0; if(bfn) memset(bkmer, 0, bfn);
			continue;
		}
		kmer = ((kmer << 2) | cur) & g->kmask;
		dna_shl_seqs(bkmer, g->bf_ksize, cur);
		j ++;
		if(j < g->ksize) continue;
		lst = (j == g->ksize)? 8 : 3 - base_bit_table[(int)seq[i - g->ksize]];
		krev = dna_rev_seq(kmer, g->ksize);
		if(krev > kmer){ krev = kmer; lst = lst == 8? 8 : 4 + lst; }
		hidx = dbg_kmer_smear(krev) % g->nh;
		KMER.val = krev;
		if(((k = get_kmerhash(g->hashs[hidx], KMER)) == NULL || (lst != 8 && k->links[lst] < min_link_cov))){
			if(j > g->ksize) fprintf(out, "%s\t%d\t%d\t%d\n", tag, i - j + 1, i - 1, j - 2);
			j = 0; kmer = 0; if(bfn) memset(bkmer, 0, bfn);
			continue;
		}
		if(bfn && j >= g->bf_ksize){
			//uint32_t l;
			//for(l=0;l<g->bf_ksize;l++) fputc(bit_base_table[bits2bit(bkmer, l)], stdout);
			memcpy(bkrev, bkmer, bfn);
			dna_rev_seqs(bkrev, g->bf_ksize);
			if(dna_cmp_seqs(bkrev, bkmer, g->bf_ksize) == -1) bk = bkrev;
			else bk = bkmer;
			hidx = dbg_kmer_smear(bk[0]) % g->nh;
			//fprintf(stdout, " -- %d in %s -- %s:%d --\n", get_cbf(g->bfs[hidx], bk, bfn), __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
			if(get_cbf(g->bfs[hidx], bk, bfn) < min_bf_cov){
				if(j > g->ksize) fprintf(out, "%s\t%d\t%d\t%d\n", tag, i - j + 1, i - 1, j - 2);
				j = 0; kmer = 0; memset(bkmer, 0, bfn);
			}
		}
	}
}

typedef struct { uint32_t kidx:31, dir:1; } klnk_t;

typedef struct {
	kmer_t *k;
	klnk_t links[8]; // kidx = 0 : NULL
} kmer_cache_t;
define_list(kcachev, kmer_cache_t);
// need to call set_userdata before use kcachehash
#define E(idx) ((kcachev*)set->userdata)->buffer[idx]
#define kmer_cache_hc(idx) 64hashcode(E(idx).k->val)
#define kmer_cache_he(idx1, idx2) (E(idx1).k->val == E(idx2).k->val)
define_hashset(kcachehash, uint32_t, kmer_cache_hc, kmer_cache_he);

#define DBG_MAX_BTIDX	0x7FFFFFFFU
#define DBG_MAX_SEED	0x00000FFFU

typedef struct {
	uint32_t qidx:12;
	int score:20;
	uint32_t link:1; // qdp is stacked togethter in qdpv, link == 1 means existing next linked qdp.
	uint32_t bt_idx:31;
} q_dp_t;
define_list(qdpv, q_dp_t);

typedef struct {
	union {
		uint32_t kidx:28, path:2, qdir:1, dir:1;
		uint32_t identifier;
	};
	union {
		int score:21, inc:10;
		uint32_t cached:1;
		uint32_t scoreinfo;
	};
	uint32_t bt_idx;
} dbg_dp_t;
define_list(dbgdpv, dbg_dp_t);
int dbg_dp_heap_cmp(uint32_t idx1, uint32_t idx2, void *ref){
	dbgdpv *dps = (dbgdpv*)ref;
	dbg_dp_t *dp1 = ref_dbgdpv(dps, idx1);
	dbg_dp_t *dp2 = ref_dbgdpv(dps, idx2);
	cmp_2nums_proc(dp2->score, dp1->score);
	return 0;
}
// set->userdata is bulit-in field of hashset
#define E(idx) ((gdpv*)set->userdata)->buffer[idx]
#define dbg_dp_hc(idx) (u64hashcode(E(idx).identifier))
#define dbg_dp_he(idx1, idx2) ((E(idx1)).identifier == (E(idx2)).identifier)
define_hashset(dphash, uint32_t, dbg_dp_hc, dbg_dp_he);
#undef E

//Speedup strategies:
//TODO: 1, kmer_cache can cache traveled kmers and their linked kmers
//TODO: 2, dphash can be divided by (qpos, path, qdir, dir), to reduce search time, aln->qtop - aln->W determinates the life-time of qpos's dphash
//TODO: 3, maintianing qtop, qmaxs(length eq qlen, max score for each pos) for all qmers. First start the 1th qmer, when its qtop come to the next qmer position, start the next qmer, and so on
//TODO:    the dbg_dp_t are shared between all qmers, the max score dp will live. How to find multiple best local alignments? analysis the qmaxs! We need 
//TODO:    a complicated qmaxs, should be {int score; uint32_t dp_idx;}.

typedef struct {
	uint32_t idx;
	int score;
	int x, y;
} dbg_hsp_t;
define_list(dbghspv, dbg_hsp_t);

typedef struct {
	DBG      *g;
	uint8_t  min_link_cov, link_cov_cutoff, kmer_cov_cutoff; // link_cov is for edges of small-kmer, kmer_cov is for big-kmer's coverage
	float    min_mat_ratio; // in any move, if total mat is less than this ratio of qry run, discard it
	int      lowcov_bonus; // penalty for path with cov < min_link_cov
	float    cov_bonus; // add coverage * cov_bonus to alignment score, usually set to FLT_MIN
	int      Z, W, M, X, I, D, E, T; // M > 0, Z is somehow like score bandwidth, W is like qpos bandwidth. T * (X or QMIS) is end-clip penalty
	int      QMAT, QMIS, QDEL, QEXT; // global mat/mis/del/gap_ext PhredQV for 5q alignment
	int      aln_mode; // 0: global; 1: local
	int      all_kmer; // whether to try all qmers
	size_t   max_nodes;

	uint8_t  *qry;
	int      qlen;
	u8list   **qvs;
	int      has_5q; // for quality values guided alignment

	dbgdpv   *qmers; // solid k-mer in qry
	dphash   *qhash; // hash of qmers
	dbghspv  *hsps; // local alignment for each qmer
	kcachev    *kcs; // cache searched kmer and links
	kcachehash *kch;
	dbgdpv    *nodes;
	u64list  *bfks; // recording bf-kmers for each node
	uint32_t bfk_n; // how many uint64_t for each record

	Heap     *heap; // sort living (to be extended) nodes by score
	dphash   *hash; // hash of nodes, avoid to search traved kmer-pos
	int      qbeg;
	int      qmax;
	int      local_max; // max score, like local alignment
	uint32_t local_idx; // node idx of local_max
	//int      qtop; // The longest move on query
	//int      smin; // smin is the low bound of dp_t->score
	//int      last_cached_pos;

	dbgdpv   *paths; // best paths of extension for each solid k-mer, the first element will be used for querying dphash
	dphash   *phash; // hash of best paths
	int      max_score;
	uint32_t best_idx;
} dbg_aligner;

dbg_aligner* init_dbgaln(DBG *g, int Z, int W, int T, int M, int X, int I, int D, int E, int QMAT, int QMIS, int QDEL, int QEXT){
	//uint32_t i;
	dbg_aligner *aln;
	aln = malloc(sizeof(dbg_aligner));
	aln->g = g;
	aln->min_mat_ratio = 0.5;
	aln->min_link_cov = g->min_link_cov;
	aln->link_cov_cutoff = 0;
	aln->kmer_cov_cutoff = 0;
	aln->lowcov_bonus = -2;
	aln->cov_bonus = FLT_MIN;
	aln->Z = Z;
	aln->W = W;
	aln->T = T;
	aln->M = M;
	aln->X = X;
	aln->I = I;
	aln->D = D;
	aln->E = E;
	aln->QMAT = QMAT;
	aln->QMIS = QMIS;
	aln->QDEL = QDEL;
	aln->QEXT = QEXT;
	aln->aln_mode = 1;
	aln->all_kmer = 0;
	aln->max_nodes = 0xFFFFFFFFFFFFFFFFLLU;
	aln->qry = NULL;
	aln->qlen = 0;
	aln->qvs = NULL;
	aln->has_5q = 0;
	aln->qmers = init_dbgdpv(1024);
	aln->qhash = init_dphash(1023);
	set_userdata_dphash(aln->qhash, aln->qmers);
	aln->kcs   = init_kcachev(1024);
	aln->kch   = init_kcachehash(1023);
	aln->hsps  = init_dbghspv(1024);
	aln->nodes = init_dbgdpv(1024);
	aln->bfks  = init_u64list(1024);
	aln->bfk_n = g->bf_ksize? (g->bf_ksize + 31) / 32 : 0;
	aln->heap  = init_heap(dbg_dp_heap_cmp, aln->nodes);
	aln->hash  = init_dphash(1023);
	set_userdata_dphash(aln->hash, aln->nodes);
	aln->qmax = - 0x7FFFFFFF;
	aln->qbeg = 0;
	aln->local_max = 0;
	aln->local_idx = 1;
	aln->paths = init_dbgdpv(1024);
	aln->phash = init2_dphash(1023, 0.33);
	set_userdata_dphash(aln->phash, aln->paths);
	aln->max_score = -0x0FFFFFFF;
	aln->best_idx = 0;
	return aln;
}

void free_dbgaln(dbg_aligner *aln){
	free_kcachev(aln->kcs);
	free_kcachehash(aln->kch);
	free_dbgdpv(aln->qmers);
	free_dphash(aln->qhash);
	free_dbghspv(aln->hsps);
	free_dbgdpv(aln->nodes);
	free_u64list(aln->bfks);
	free_heap(aln->heap);
	free_dphash(aln->hash);
	free_dbgdpv(aln->paths);
	free_dphash(aln->phash);
	free(aln);
}

inline uint8_t last_base_dp_kmer(dbg_dp_t *d, uint32_t ksize){
	return d->dir? ((~(d->k->val >> ((ksize - 1) << 1))) & 0x03) : (d->k->val & 0x03);
}

void trace_aln_paths(dbg_aligner *aln, uint32_t idx, uint32_t limit, FILE *out){
	dbg_dp_t *d;
	uint32_t cnt;
	cnt = 0;
	while(idx){
		cnt ++;
		if(cnt == limit) break;
		d = ref_dbgdpv(aln->paths, idx);
		fprintf(out, "%d[%d,%d,%d,%p];\n", idx, d->bt_idx, d->qpos, d->path, d->k);
		idx = d->bt_idx;
	}
	fflush(out);
}

int count_covered_qmers(dbg_aligner *aln, uint32_t sidx){
	dbg_dp_t *d;
	int cnt;
	cnt = 0;
	while(sidx){
		d = ref_dbgdpv(aln->nodes, sidx);
		sidx = d->bt_idx;
		if(d->path != 0) continue;
		aln->qmers->buffer[aln->qmers->size] = *d;
		d = ref_dbgdpv(aln->qmers, aln->qmers->size);
		if(d->qdir){
			d->qdir = 0;
			d->dir  = d->dir ^ 1;
			d->qpos = aln->qlen - (d->qpos - aln->g->ksize);
		}
		if(exists_dphash(aln->qhash, aln->qmers->size)) cnt ++;
	}
	return cnt;
}

typedef struct {
	int score;
	int qb, qe;
	int alns[4];
} dbg_hit_t;

dbg_hit_t stat_aln_paths(dbg_aligner *aln, uint32_t idx){
	dbg_hit_t R;
	dbg_dp_t *d;
	memset(&R, 0, sizeof(dbg_hit_t));
	while(idx){
		d = ref_dbgdpv(aln->nodes, idx);
		R.alns[d->path] ++;
		R.score += d->inc;
		idx = d->bt_idx;
	}
	return R;
}

static const int dbg_path2cigar[] = {0, 2, 1, 0};

dbg_hit_t call_correct_seq(dbg_aligner *aln, uint32_t qmer_idx, u8list *seq, u32list *cigars){
	dbg_hit_t R;
	dbg_dp_t *d;
	uint64_t kmer;
	uint32_t *ptr, idx, i;
	uint8_t v;
	if(seq) clear_u8list(seq);
	if(cigars) clear_u32list(cigars);
	d = ref_dbgdpv(aln->qmers, qmer_idx);
	set_dbgdpv(aln->paths, 0, *d);
	aln->paths->buffer[0].qdir = 1; // first pass of the reverse strand
	aln->paths->buffer[0].dir  = d->dir ^ 1;
	aln->paths->buffer[0].qpos = aln->qlen - (d->qpos - aln->g->ksize);
	ptr = get_dphash(aln->phash, 0);
	if(ptr == NULL){
		fprintf(stderr, " -- Unexpected error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	idx = *ptr;
	memset(&R, 0, sizeof(dbg_hit_t));
	while(idx){
		d = ref_dbgdpv(aln->paths, idx);
		if(seq && d->path != 1){
			// v = dir? ((~(d->k->val >> ((aln->g->ksize - 1) << 1))) & 0x03) : (d->k->val & 0x03);
			// We want (~v) & 0x03, thre complementary
			v = d->dir? (d->k->val >> ((aln->g->ksize - 1) << 1)) : ((~d->k->val) & 0x03);
			if(d->path) v += 4;
			push_u8list(seq, v);
		}
		if(cigars) kswx_push_cigar(cigars, dbg_path2cigar[d->path], 1);
		idx = d->bt_idx;
		R.alns[d->path] ++;
	}
	R.qb = aln->qlen - d->qpos;
	if(seq) reverse_u8list(seq);
	if(cigars) reverse_u32list(cigars);
	d = ref_dbgdpv(aln->qmers, qmer_idx);
	set_dbgdpv(aln->paths, 0, *d);
	aln->paths->buffer[0].qdir = 0; // the forward strand
	aln->paths->buffer[0].dir  = d->dir ^ 0;
	aln->paths->buffer[0].qpos = d->qpos;
	ptr = get_dphash(aln->phash, 0);
	if(ptr == NULL){
		fprintf(stderr, " -- Unexpected error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	idx = *ptr;
	if(seq){
		d = ref_dbgdpv(aln->paths, idx);
		kmer = d->dir? dna_rev_seq(d->k->val, aln->g->ksize) : d->k->val;
		//fputc('\n', stdout);
		for(i=1;i+1<aln->g->ksize;i++){
			v = (kmer >> ((aln->g->ksize - 1 - i) << 1)) & 0x03;
			push_u8list(seq, v);
			//fputc(bit_base_table[v], stdout);
		}
		//fputc('\n', stdout);
	}
	if(cigars) kswx_push_cigar(cigars, 0, aln->g->ksize - 2);
	R.alns[0] += aln->g->ksize - 2;
	while(idx){
		d = ref_dbgdpv(aln->paths, idx);
		if(seq && d->path != 1){
			v = last_base_dp_kmer(d, aln->g->ksize);
			if(d->path) v += 4;
			push_u8list(seq, v);
		}
		if(cigars) kswx_push_cigar(cigars, dbg_path2cigar[d->path], 1);
		idx = d->bt_idx;
		R.alns[d->path] ++;
	}
	R.qe = d->qpos;
	return R;
}

inline void init_bf_kmer_dbgaln(dbg_aligner *aln, dbg_dp_t *d){
	uint64_t *k, kmer;
	uint32_t i;
	encap_u64list(aln->bfks, aln->bfk_n);
	kmer = d->dir? dna_rev_seq(d->k->val, aln->g->ksize) : d->k->val;
	k = aln->bfks->buffer + (d - aln->nodes->buffer) * aln->bfk_n;
	memset(k, 0, aln->bfk_n * sizeof(uint64_t));
	for(i=0;i<aln->g->ksize;i++) dna_shl_seqs(k, aln->g->bf_ksize, (kmer >> ((aln->g->ksize - 1 - i) << 1)) & 0x03U);
}

inline int prepare_bf_kmer_dbgaln(dbg_aligner *aln, uint32_t dpidx, uint8_t link){
	uint64_t *k1, *k2, *k3;
	uint32_t hidx;
	encap_u64list(aln->bfks, 2 * aln->bfk_n);
	k1 = aln->bfks->buffer + ((size_t)dpidx) * aln->bfk_n;
	k2 = aln->bfks->buffer + aln->nodes->size * aln->bfk_n;
	k3 = aln->bfks->buffer + (aln->nodes->size + 1) * aln->bfk_n;
	memcpy(k2, k1, aln->bfk_n * sizeof(uint64_t));
	dna_shl_seqs(k2, aln->g->bf_ksize, link);
	memcpy(k3, k2, aln->bfk_n * sizeof(uint64_t));
	dna_rev_seqs(k3, aln->g->bf_ksize);
	if(dna_cmp_seqs(k2, k3, aln->g->bf_ksize) < 0) k3 = k2;
	hidx = dbg_kmer_smear(k3[0]) % aln->g->nh;
	return get_cbf(aln->g->bfs[hidx], k3, aln->bfk_n * sizeof(uint64_t));
}

inline void cpy_bf_kmer_dbgaln(dbg_aligner *aln, uint32_t dpidx){
	uint64_t *k1, *k2;
	encap_u64list(aln->bfks, aln->bfk_n);
	k1 = aln->bfks->buffer + ((size_t)dpidx) * aln->bfk_n;
	k2 = aln->bfks->buffer + (aln->nodes->size) * aln->bfk_n;
	memcpy(k2, k1, aln->bfk_n * sizeof(uint64_t));
}

inline void print_dp_kmers(dbg_aligner *aln, uint32_t dpidx){
	dbg_dp_t *d;
	uint64_t kmer, *k;
	uint32_t i;
	d = ref_dbgdpv(aln->nodes, dpidx);
	kmer = d->dir? dna_rev_seq(d->k->val, aln->g->ksize) : d->k->val;
	for(i=0;i<aln->g->ksize;i++) fputc(bit_base_table[(kmer >> ((aln->g->ksize - 1 - i) << 1)) & 0x03], stdout);
	fputc('\t', stdout);
	k = aln->bfks->buffer + dpidx * aln->bfk_n;
	for(i=0;i<aln->g->bf_ksize;i++) fputc(bit_base_table[bits2bit(k, i)], stdout);
}

static const int hashing_dp_path[2][4] = {{1, 1, 1, 1}, {1, 0, 0, 0}};

#define DBG_DP_CORE_FINISH	0
#define DBG_DP_CORE_EXTEND	1
#define DBG_DP_CORE_DISCARD	2
#define DBG_DP_CORE_DEADEND	3
#define DBG_DP_CORE_FOUND	4

#define REMOVE_HEAP(aln, idx) array_heap_remove((aln)->heap->ptrs->buffer, (aln)->heap->ptrs->size, (aln)->heap->ptrs->cap, uint32_t, idx, num_cmp((aln)->nodes->buffer[b].score, (aln)->nodes->buffer[a].score))
#define PUSH_HEAP(aln, idx) array_heap_push((aln)->heap->ptrs->buffer, (aln)->heap->ptrs->size, (aln)->heap->ptrs->cap, uint32_t, idx, num_cmp((aln)->nodes->buffer[b].score, (aln)->nodes->buffer[a].score))

// TODO: Supports we know some or all of variants in population or especially in multiple copies of repetitive sequences,
// TODO: need to lower the score penalty for known variants. Finding bubbles on DBG? will help reduce false alignment step
// TODO: No, if the right path exists on DBG, we only need to care about sequencing errors
inline int dbg_aln_core(dbg_aligner *aln){
	dbg_dp_t *dp, *dp2, *dp1;
	kmer_cache_t *kc;
	kmer_t *k, KMER;
	uint64_t kmer, knew, *bk, BK[DBG_MAX_BF_KSIZE / 32], BR[DBG_MAX_BF_KSIZE / 32];
	uint32_t idx, i, d, q, hidx, cnt, *ptr, didx, pidx, path, _qpos;
	uint8_t lst_base;
	int score, exists, kc_exists, found, qpos;
	if((idx = peer_heap(aln->heap)) == 0xFFFFFFFFU) return DBG_DP_CORE_FINISH;
	encap_dbgdpv(aln->nodes, 9);
	dp = ref_dbgdpv(aln->nodes, idx);
	if(verbose > 3){
		printf("DP: %u<-%u\tpos=%d\tpath=%d\tscore=%d\theap=%d\tsearch=%d\tqmax=%d", idx, dp->bt_idx, dp->qpos, dp->path, (int)dp->score, (int)aln->heap->ptrs->size, (int)aln->nodes->size, aln->qmaxs->buffer[dp->qpos]);
		if(verbose > 4){
			if(verbose > 5){
				dbg_hit_t R = stat_aln_paths(aln, idx);
				printf("\t[%d,%d,%d,%d,%d]", R.score, R.alns[0], R.alns[3], R.alns[1], R.alns[2]);
			}
			printf("\t");
			print_dp_kmers(aln, idx);
		}
		printf("\n");
	}
	if(dp->score > aln->local_max){
		aln->local_max = dp->score;
		aln->local_idx = idx;
	}
	if(dp->qpos >= aln->qlen){
		return DBG_DP_CORE_FOUND;
	}
	{
		if(dp->score < aln->smin){
			REMOVE_HEAP(aln, 0); return DBG_DP_CORE_DISCARD;
		}
		//if(dp->mat < (dp->qpos - aln->qbeg) * aln->min_mat_ratio){
			//REMOVE_HEAP(aln, 0); return DBG_DP_CORE_DISCARD;
		//}
		if(dp->qpos + aln->W < aln->qtop){ // W bandwidth
			REMOVE_HEAP(aln, 0); return DBG_DP_CORE_DISCARD;
		} else if(dp->qpos > aln->qtop){
			aln->qtop = dp->qpos;
			// Tidy heap to improve the speed of heap_push
			for(i=aln->heap->ptrs->size-1;i>0;i--){
				if(aln->nodes->buffer[aln->heap->ptrs->buffer[i]].qpos + aln->W < aln->qtop){
					REMOVE_HEAP(aln, i);
				}
			}
		}
		if(dp->score < aln->qmaxs->buffer[dp->qpos]){ // Z bandwidth
			REMOVE_HEAP(aln, 0); return DBG_DP_CORE_DISCARD;
		} else if(dp->score + aln->Z * (aln->has_5q? aln->QMIS : aln->X) > aln->qmaxs->buffer[dp->qpos]){
			aln->qmaxs->buffer[dp->qpos] = dp->score + aln->Z * (aln->has_5q? aln->QMIS : aln->X);
		}
	}
	if(1){
		ptr = prepare_dphash(aln->hash, idx, &exists);
		if(exists){ // Found cached search
			dp2 = ref_dbgdpv(aln->nodes, *ptr);
			if(dp2->cached){
				return process_cached_dps_dbgaln(aln, idx, *ptr);
			} else if(dp->score > dp2->score){
				*ptr = idx;
			} else{
				REMOVE_HEAP(aln, 0); return DBG_DP_CORE_DISCARD;
			}
		} else {
			*ptr = idx;
		}
	}
	// extending
	REMOVE_HEAP(aln, 0);
	cnt = 0;
	for(i=0;i<4;i++){
		if(dp->k->links[dp->dir * 4 + i] < aln->link_cov_cutoff) continue;
		cnt ++;
	}
	if(cnt == 0) return DBG_DP_CORE_DEADEND;
	dp->fw_idx = aln->nodes->size;
	kmer = dp->dir? dna_rev_seq(dp->k->val, aln->g->ksize) : dp->k->val;
	lst_base = kmer & 0x03;
	if(dp->qdir){
		_qpos = aln->qlen - 1 - dp->qpos;
		q = (~aln->qry[_qpos]) & 0x03;
	} else {
		_qpos = dp->qpos;
		q = aln->qry[_qpos];
	}
	memset(&KMER, 0, sizeof(kmer_t));
	cnt = 0;
	dc = prepare_dbgcachehash(aln->gcache, (dbg_cache_t){dp->k, (kmer_t*[8]){0, 0, 0, 0, 0, 0, 0, 0}}, &dc_exists);
	if(!dc_exists) memset(dc->links, 0, sizeof(kmer_t*) * 8);
	for(i=0;i<4;i++){
		if(dp->k->links[dp->dir * 4 + i] < aln->link_cov_cutoff) continue;
		if(aln->bfk_n && (prepare_bf_kmer_dbgaln(aln, idx, i) < aln->kmer_cov_cutoff && dp->qpos - aln->qbeg >= (int)aln->g->bf_ksize)) continue; // Check whether this move exists in big-kmer
		// match or mismatch
		if(aln->has_5q){
			if(i == q) score = dp->score + aln->QMAT;
			else if(dp->qdir){
				score = dp->score + ((aln->qvs[5]->buffer[_qpos] < 4 && i == ((~aln->qvs[5]->buffer[_qpos]) & 0x03U))? - aln->qvs[1]->buffer[_qpos] : aln->QMIS);
			} else {
				score = dp->score + ((i == aln->qvs[5]->buffer[_qpos])? - aln->qvs[1]->buffer[_qpos] : aln->QMIS);
			}
			if(dp->path == 2){
				if(i == q && i == lst_base){
					if(aln->qvs[4]->buffer[_qpos] < - dp->inc){
						score = score - (dp->inc) + (- aln->qvs[4]->buffer[_qpos]);
					}
				}
			}
		} else score = dp->score + (i == q? aln->M : aln->X) + dp->k->links[dp->dir * 4 + i] * aln->cov_bonus;
		if(dc_exists) k = dc->links[dp->dir * 4 + i];
		if(k == NULL){
			knew = ((kmer << 2) | i) & aln->g->kmask;
			KMER.val = dna_rev_seq(knew, aln->g->ksize);
			if(KMER.val < knew){ d = 1; } else { KMER.val = knew; d = 0; }
			hidx = (dbg_kmer_smear(KMER.val) % aln->g->nh);
			k = get_kmerhash(aln->g->hashs[hidx], KMER);
			dc->links[dp->dir * 4 + i] = k;
		}
		if(k == NULL){
			fprintf(stderr, " -- Unexpected error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
		}
		dp2 = next_ref_dbgdpv(aln->nodes);
		dp2->aux1 = 0;
		aln->bfks->size += aln->bfk_n;
		dp2->score  = score;// + (i == q? FLT_MIN : 0); // favor on match
		if(dp->k->links[dp->dir * 4 + i] < aln->min_link_cov) dp2->score += aln->lowcov_bonus;
		dp2->inc    = dp2->score - dp->score;
		//dp2->mat    = dp->mat + (i == q? 1 : 0);// + (dp->k->links[dp->dir * 4 + i] / 256.0f); // favor on high coverage
		dp2->qpos   = dp->qpos + 1;
		dp2->k      = k;
		dp2->bt_idx = idx;
		dp2->path   = (i == q? 0 : 3);
		dp2->qdir   = dp->qdir;
		dp2->dir    = d;
		PUSH_HEAP(aln, aln->nodes->size - 1);
		dp->link |= 0x1 << (i * 2 + 0);
		cnt ++;
		// deletion
		if(dp->path != 1){ // affine gap alignment
			if(aln->bfk_n) cpy_bf_kmer_dbgaln(aln, aln->nodes->size - 1);
			aln->bfks->size += aln->bfk_n;
			dp2 = next_ref_dbgdpv(aln->nodes);
			if(aln->has_5q){ // complicated code for deletions may be more than 1bp, caution of the deletionTag matching
				if(dp->qdir){
					score = dp->score + ((aln->qvs[6]->buffer[_qpos] < 4 && i == ((~aln->qvs[6]->buffer[_qpos]) & 0x03))? -aln->qvs[3]->buffer[_qpos] : (dp->path == 2? aln->QEXT : aln->QDEL));
				} else {
					score = dp->score + (i == aln->qvs[6]->buffer[_qpos + 1]? -aln->qvs[3]->buffer[_qpos + 1] : (dp->path == 2? aln->QEXT : aln->QDEL));
				}
			} else score = dp->score + aln->E + (dp->path == 2? 0 : aln->D);
			dp2->aux1   = 0;
			dp2->score  = score;
			if(dp->k->links[dp->dir * 4 + i] < aln->min_link_cov) dp2->score += aln->lowcov_bonus;
			dp2->inc    = dp2->score - dp->score;
			dp2->qpos   = dp->qpos;
			dp2->k      = k;
			dp2->bt_idx = idx;
			dp2->path   = 2;
			dp2->qdir   = dp->qdir;
			dp2->dir    = d;
			PUSH_HEAP(aln, aln->nodes->size - 1);
			dp->link |= 0x1 << (i * 2 + 1);
		}
	}
	if(cnt == 0) return DBG_DP_CORE_DEADEND;
	// inserion
	if(dp->path != 2){ // affine gap alignment
		if(aln->has_5q){
			score = dp->score + (- aln->qvs[2]->buffer[_qpos]);
		} else score = dp->score + aln->E + (dp->path == 1? 0 : aln->I);
		if(aln->bfk_n) cpy_bf_kmer_dbgaln(aln, idx);
		aln->bfks->size += aln->bfk_n;
		dp2 = next_ref_dbgdpv(aln->nodes);
		dp2->aux1   = 0;
		dp2->score  = score;// + FLT_MIN; // just a bit favor on insertion
		dp2->inc    = dp2->score - dp->score;
		dp2->qpos   = dp->qpos + 1;
		dp2->k      = dp->k;
		dp2->bt_idx = idx;
		dp2->path   = 1;
		dp2->qdir   = dp->qdir;
		dp2->dir    = dp->dir;
		PUSH_HEAP(aln, aln->nodes->size - 1);
		dp->link |= 0x1 << 8;
	}
	return DBG_DP_CORE_EXTEND;
}

void dbg_aln(dbg_aligner *aln, uint8_t *qry, uint32_t qlen, u8list *qvs[7]){
	dbg_dp_t *d, *d1, *d2, *d3;
	kmer_t *k, KMER;
	uint64_t kmer;
	uint32_t ii, i, j, idx, sidx, hidx, *ptr;
	int score, mat, dir, smax, ret, goon, found;
	uint8_t lst;
	clear_dbgcachehash(aln->gcache);
	clear_dbgdpv(aln->qmers);
	clear_dphash(aln->qhash);
	clear_dbghspv(aln->hsps);
	clear_dbgdpv(aln->paths);
	next_ref_dbgdpv(aln->paths);
	clear_dphash(aln->phash);
	clear_b32list(aln->qmaxs);
	for(i=0;i<=qlen;i++) push_b32list(aln->qmaxs, -0x0FFFFFFF);
	aln->best_idx = 0xFFFFFFFFU;
	aln->max_score = - 0x0FFFFFFF;
	aln->qry  = qry;
	aln->qlen = qlen;
	aln->qvs  = qvs;
	aln->has_5q = (qvs != NULL);
	memset(&KMER, 0, sizeof(kmer_t));
	kmer = 0;
	for(i=0;i<qlen;i++){
		kmer = ((kmer << 2) | qry[i]) & aln->g->kmask;
		if(i + 1 < aln->g->ksize) continue;
		KMER.val = dna_rev_seq(kmer, aln->g->ksize);
		lst = i >= aln->g->ksize? 3 - (qry[i - aln->g->ksize] & 0x03) : 8;
		if(KMER.val < kmer){ dir = 1; }
		else { dir = 0; lst = lst == 8? 8 : 4 + lst; KMER.val = kmer; }
		hidx = dbg_kmer_smear(KMER.val) % aln->g->nh;
		k = get_kmerhash(aln->g->hashs[hidx], KMER);
		if(k == NULL) continue;
		//if(!aln->all_kmer && lst < 8){
			//if(k->links[lst] > aln->link_cov_cutoff) continue;
		//}
		d = next_ref_dbgdpv(aln->qmers);
		memset(d, 0, sizeof(dbg_dp_t));
		d->qpos   = i + 1;
		d->k      = k;
		d->dir    = dir;
		push_dbghspv(aln->hsps, (dbg_hsp_t){aln->qmers->size - 1, -0x0FFFFFFF, i + 1, i + 1});
		put_dphash(aln->qhash, aln->qmers->size - 1);
	}
	encap_dbgdpv(aln->qmers, 1);
	if(verbose){
		fprintf(stdout, "Found %d seed-kmers\n", (int)aln->qmers->size); fflush(stdout);
	}
	for(dir=0;dir<2;dir++){
		for(ii=0;ii<aln->qmers->size;ii++){
			i = dir? ii : aln->qmers->size - 1 - ii;
			if(debug_opt1 && (int)i != atoi(debug_opt1)) continue;
			d = ref_dbgdpv(aln->qmers, i);
			if(d->cached) continue;
			score = 0;
			mat   = 0;
			goon = 1;
			clear_dbgdpv(aln->nodes);
			clear_u64list(aln->bfks);
			next_ref_dbgdpv(aln->nodes);
			encap_and_inc_u64list(aln->bfks, aln->bfk_n);
			clear_heap(aln->heap);
			clear_dphash(aln->hash);
			d1 = next_ref_dbgdpv(aln->nodes);
			d1->aux1   = 0;
			d1->aux2   = 0;
			d1->score  = aln->g->ksize * (aln->has_5q? aln->QMAT : aln->M);
			d1->k      = d->k;
			d1->path   = 0;
			d1->qdir   = dir;
			d1->dir    = d->dir ^ dir;
			d1->qpos   = dir? (int)(aln->qlen - (d->qpos - aln->g->ksize)) : d->qpos;
			if(aln->bfk_n){ init_bf_kmer_dbgaln(aln, d1); aln->bfks->size += aln->bfk_n; }
			for(j=d1->qpos;(int)j<=aln->qlen;j++) set_b32list(aln->qmaxs, j, -0x0FFFFFFF);
			aln->local_max = 0;
			aln->local_idx = 1;
			aln->qbeg = d1->qpos - aln->g->ksize;
			aln->qtop = d1->qpos;
			aln->last_cached_pos = 0;
			push_heap(aln->heap, 1);
			smax = - 0x1FFFFFFF;
			sidx = 0;
			found = 0;
			while(1){
				ret = dbg_aln_core(aln);
				if(ret == DBG_DP_CORE_FINISH) break;
				switch(ret){
					case DBG_DP_CORE_EXTEND:
					case DBG_DP_CORE_DISCARD: break; // Do nothing
					//TODO: if only partial of qry can be aligned against DBG, how to do?
					//TODO: Now, skip partial alignment. If no full-length, give up
					case DBG_DP_CORE_DEADEND: break; // Do nothing
					case DBG_DP_CORE_FOUND:
						idx = pop_heap(aln->heap);
						d2  = ref_dbgdpv(aln->nodes, idx);
						if((int)d2->score > smax){
							smax = d2->score;
							sidx = idx;
						}
						found = 1;
						break;
				}
				if(aln->nodes->size > aln->max_nodes){
					if(verbose > 1){
						fprintf(stdout, "memory limit, break\n"); fflush(stdout);
					}
					break;
				}
			}
			if(verbose > 1){
				fprintf(stdout, "[%s] pos=%d %d/%d\tdir=%d\tfound=%d\tsearch=%u\tglobal=%d local=%d(%d-%d)", date(), d->qpos, i, (uint32_t)aln->qmers->size, dir, found, (unsigned)aln->nodes->size, smax, aln->local_max,
				dir? aln->qlen - aln->nodes->buffer[1].qpos : aln->nodes->buffer[1].qpos, dir? aln->qlen - aln->nodes->buffer[aln->local_idx].qpos : aln->nodes->buffer[aln->local_idx].qpos);
				if(found){
					dbg_hit_t R = stat_aln_paths(aln, sidx);
					fprintf(stdout, "\t[%d,%d,%d,%d,%d,cover=%d]", smax, R.alns[0], R.alns[3], R.alns[1], R.alns[2], count_covered_qmers(aln, sidx));
				}
				if(aln->aln_mode){
					dbg_hit_t R = stat_aln_paths(aln, aln->local_idx);
					fprintf(stdout, "\t[%d,%d,%d,%d,%d,cover=%d]", aln->local_max, R.alns[0], R.alns[3], R.alns[1], R.alns[2], count_covered_qmers(aln, aln->local_idx));
				}
				fprintf(stdout, "\n"); fflush(stdout);
			}
			if(aln->aln_mode){
				smax = aln->local_max;
				sidx = aln->local_idx;
				found = 1;
			}
			score += smax;
			if(!found){ goon = 0; break; }
			mat += ref_dbgdpv(aln->nodes, sidx)->mat;
			//Backtrace sidx: add the best path to aln->paths and aln->phash
			idx = 0;
			while(1){
				d1 = ref_dbgdpv(aln->nodes, sidx);
				sidx = d1->bt_idx;
				d2 = next_ref_dbgdpv(aln->paths);
				*d2 = *d1;
				d2->bt_idx = idx; // Here, bt_idx is not backtrace, but forwardly pointer to next dp_t
				d2->score = 0;
				idx = aln->paths->size - 1;
				if(d2->path == 0 && aln->all_kmer == 0){
					d3 = aln->qmers->buffer + aln->qmers->size;
					*d3 = *d2;
					if(d3->qdir){
						d3->qdir = 0;
						d3->dir  = !d3->dir;
						d3->qpos = aln->qlen - (d3->qpos - aln->g->ksize);
					}
					if((ptr = get_dphash(aln->qhash, aln->qmers->size))){
						aln->qmers->buffer[*ptr].cached = 1;
					}
				}
				if(sidx == 0){
					put_dphash(aln->phash, idx);
					break;
				}
			}
		}
			if(goon){
				dbg_hit_t rs;
				rs = call_correct_seq(aln, i, NULL, NULL);
				d->score = score;
				ref_dbghspv(aln->hsps, i)->score = score;
				ref_dbghspv(aln->hsps, i)->x     = rs.qb;
				ref_dbghspv(aln->hsps, i)->y     = rs.qe;
				if(score > aln->max_score){
					aln->max_score = score; aln->best_idx = i;
				}
				if(verbose > 1){
					fprintf(stdout, "[%s] Correction = %d-%d\t%d\t%d\t%d\t%d\t%d\n", date(), rs.qb, rs.qe, score, rs.alns[0], rs.alns[3], rs.alns[1], rs.alns[2]); fflush(stdout);
				}
			}
	}
}

// Get best alignment(s), if two hsp doesn't overlap, output both, and three, so on
int best_hsps_dbgaln(dbg_aligner *aln, u32list *rs, int min_hsp_len, int max_hsp_ovl){
	dbg_hsp_t *p1, *p2;
	uint32_t i, j, f;
	int b, e, ret;
	ret = 0;
	clear_u32list(rs);
	sort_array(aln->hsps->buffer, aln->hsps->size, dbg_hsp_t, b.score > a.score);
	for(i=0;i<aln->hsps->size;i++){
		p1 = ref_dbghspv(aln->hsps, i);
		if(p1->score <= - 0x0FFFFFFF) continue;
		if(p1->x + min_hsp_len > p1->y) continue;
		f = 0;
		for(j=0;j<rs->size;j++){
			p2 = ref_dbghspv(aln->hsps, rs->buffer[j]);
			b = num_max(p1->x, p2->x);
			e = num_min(p1->y, p2->y);
			if(b + max_hsp_ovl < e){ f = 1; break; }
		}
		if(f) continue;
		push_u32list(rs, i);
		ret ++;
	}
	sort_array(rs->buffer, rs->size, uint32_t, aln->hsps->buffer[a].x > aln->hsps->buffer[b].x);
	return ret;
}

thread_beg_def(maln);
dbg_aligner *aln;
String *tag;
u8list *qry;
u8list *qvs[7];
int has_5q;
int min_hsp_len;
int max_hsp_ovl;
u8list *cns;
u32list *cigars;
u32list *corrs;
int is_corrected;
thread_end_def(maln);

thread_beg_func(maln);
uint32_t i, idx;
thread_beg_loop(maln);
if(maln->qry->size){
	dbg_aln(maln->aln, maln->qry->buffer, maln->qry->size, maln->has_5q? (u8list**)maln->qvs : NULL);
	clear_u8list(maln->cns);
	clear_u32list(maln->corrs);
	maln->is_corrected = 0;
	if(debug_opt2){
		idx = atoi(debug_opt2);
		for(i=0;i<maln->aln->hsps->size;i++){
			if(i == idx) continue;
			maln->aln->hsps->buffer[i].score = - 0x0FFFFFFF;
		}
	}
	maln->is_corrected = best_hsps_dbgaln(maln->aln, maln->corrs, maln->min_hsp_len, maln->max_hsp_ovl);
		//maln->rs = call_correct_seq(maln->aln, maln->aln->best_idx, maln->cns, maln->cigars);
}
thread_end_loop(maln);
thread_end_func(maln);

int usage(){
	printf(
	"Program: wtcorr\n"
	" Long sequence corrector using DBG from accurate short reads\n"
	"Author: Jue Ruan\n"
	"Version: 1.0 20150715\n"
	"Usage: wtcorr [options] [long_seq_file(fa/fq/f5q)]\n"
	"Options:\n"
	" -h          Show this document\n"
	" -t <int>    Number of threads, [1]\n"
	" -i <string> Short reads file(fa/fq), +, *\n"
	" -r <string> Load DBG from file, will ignore -i\n"
	" -R <string> Similar with -r, will share DBG with other process\n"
	"             It is faster to use -R, but will generate a description file for DBG file\n"
	" -w <string> Dump DBG into file, will exit without correcting, ignore -o\n"
	" -k <int>    kmer-size, odd, 15 <= k <= 31, [25]\n"
	" -K <int>    Big-kmer-size, odd, k < K <= 127, [45]\n"
	"             kmer-size is smaller, used to build DBG\n"
	"             big-kmer is stored in bloomfilter, used to verify the kmer moving\n"
	" -c <int>    Low edge coverage in DBG, <= 255 [5]\n"
	" -C <int>    Penalty score for move along edge with cov less than low cov [-1]\n"
	" -0 <int>    Discard all edges with coverage < <value>, as 0(zero) [3]\n"
	"             NOTE: -c is somehow smaller than expected cov, but still a big number,\n"
	"                   -0 cuts the edges caused by sequencing errors\n"
	"                   If the DBG is built using error-free sequences, e.g. contigs, please set -0 0 -c 1\n"
	" -1 <int>    Disacard big-kmers with coverage < <value>, value = 1/2/3, [2]\n"
	" -T <string> Test whether sequences from file exist in dbg, print message on stdout\n"
	" -o <string> Output corrected sequences, *\n"
	" -f          Force overwrite\n"
	" -L          Local alignment mode, default: yes\n"
	" -G          Global alignment mode\n"
	" -a          Try all matched kmers, default: skip matched kmers in HSP\n"
	"             Will be very sloooooow\n"
	" -M <int>    Alignment score: match, [1]\n"
	" -X <int>    Alignment score: mismatch, [-5]\n"
	" -I <int>    Alignment score: insertion, [-2]\n"
	" -D <int>    Alignment score: deletion, [-3]\n"
	" -E <int>    Alignment score: gap extension, [-1]\n"
	" -u <int>    Max memory usage for each thread in correction, in M byte, [4096]\n"
	" -Z <int>    Something like alignment score bandwidth, [3]\n"
	"             When coming to one position of query, discarding alignment moves have score less than <max> + <-Z> * (<-X> or <-x>(in quality mode))\n"
	"             Try increase -Z to get more accurate result, but will take much much more time\n"
	" -W <int>    Something like alignment query move bandwidth, [64]\n"
	"             Discard all moves are slower more than <-W> bp on query than the fastest\n"
	" -m <int>    Match bonus, used in quality guided alignment, [5]\n"
	" -x <int>    Negative phredQV for anonymous mismatch, used in quality guided alignment, [-20]\n"
	" -d <int>    Negative phredQV for anonymous default deletion, used in quality guided alignment, [-15]\n"
	" -e <int>    Negative phredQV for gap extension, used in quality guided alignment, [-5]\n"
	" -Q          Disable pacbio quality even provided\n"
	" -v          Verbose, +\n"
	" -8 <string> Option for debug, flexible meaning\n"
	" -9 <string> Option for debug, flexible meaning\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	DBG *g;
	FileReader *fr;
	Sequence *seq;
	cplist *srfs;
	String *alns[2];
	char *dbg_load, *dbg_dump, *chkf, *outf, *leftf;
	FILE *dump, *out, *left;
	size_t max_mem, cur_mem;
	dbg_hsp_t *hsp;
	uint32_t num;
	int ncpu, ksize, bf_ksize, cut_cov[2], min_cov, lowcov_bonus, M, X, I, D, E, Z, W, T, m, x, d, e, mode, all_kmer, dbg_share, overwrite, use_5q, min_hsp, ovl_hsp, c;
	int i, j, end, eof_fr;
	thread_preprocess(maln);
	ncpu = 1;
	max_mem = 4096;
	srfs = init_cplist(4);
	dbg_load = NULL;
	dbg_dump = NULL;
	ksize = 0;
	bf_ksize = 0;
	mode = 1; // local alignment
	min_hsp = 500;
	ovl_hsp = 100;
	min_cov = 5;
	cut_cov[0] = 3;
	cut_cov[1] = 2;
	lowcov_bonus = -1;
	chkf = NULL;
	outf = NULL;
	all_kmer = 0;
	use_5q = 1;
	overwrite = 0;
	dbg_share = 0;
	M = 1; X = -5; I = -2; D = -3; E = -1; Z = 3; W = 64; m = 5; x = -20; d = -15; e = -5; T = -10;
	while((c = getopt(argc, argv, "vht:u:i:r:R:w:k:K:c:C:0:1:T:o:fGLaAM:X:I:D:E:Z:W:m:x:d:e:Q8:9:")) != -1){
		switch(c){
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(srfs, optarg); break;
			case 'u': max_mem = atoi(optarg); break;
			case 'r': dbg_load = optarg; dbg_share = 0; break;
			case 'R': dbg_load = optarg; dbg_share = 1; break;
			case 'w': dbg_dump = optarg; break;
			case 'k': ksize = atoi(optarg); break;
			case 'K': bf_ksize = atoi(optarg); break;
			case 'c': min_cov = atoi(optarg); break;
			case 'C': lowcov_bonus = atoi(optarg); break;
			case '0': cut_cov[0] = atoi(optarg); break;
			case '1': cut_cov[1] = atoi(optarg); break;
			case 'T': chkf = optarg; break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'G': mode = 0; break;
			case 'L': mode = 1; break;
			case 'a': all_kmer = 1; break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'Z': Z = atoi(optarg); break;
			case 'W': W = atoi(optarg); break;
			case 'm': m = atoi(optarg); break;
			case 'x': x = atoi(optarg); break;
			case 'd': d = atoi(optarg); break;
			case 'e': e = atoi(optarg); break;
			case 'Q': use_5q = 0; break;
			case 'v': verbose ++; break;
			case '8': debug_opt1 = optarg; break;
			case '9': debug_opt2 = optarg; break;
			default: return usage();
		}
	}
	if(!overwrite && outf && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, " --  File '%s' already exists in %s -- %s:%d --\n", outf, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(!overwrite && dbg_dump && strcmp(dbg_dump, "-") && file_exists(dbg_dump)){
		fprintf(stderr, " --  File '%s' already exists in %s -- %s:%d --\n", dbg_dump, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(dbg_load && !file_exists(dbg_load)){
		fprintf(stderr, " --  File '%s' not exists in %s -- %s:%d --\n", dbg_load, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(dbg_dump == NULL && outf == NULL && chkf == NULL){
		fprintf(stderr, " -- Miss output file (-o your_output_file) in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		return usage();
	}
	if(dbg_load){
		if(verbose) fprintf(stdout, "[%s] Loading DBG from %s\n", date(), dbg_load);
		if(srfs->size){
			fprintf(stderr, " --  When provided DBG file to be loaded, ignores data from short reads files in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		}
		if(dbg_share) g = mem_load_obj_file(&DBG_obj_desc, dbg_load, NULL, NULL, NULL);
		else g = mem_read_obj_file(&DBG_obj_desc, dbg_load, NULL, NULL, NULL);
		if(ksize){
			if(ksize != (int)g->ksize){
				fprintf(stderr, " --  parameter ksize(%d) is set to the same in loaded DBG(%d) in %s -- %s:%d --\n", ksize, g->ksize, __FUNCTION__, __FILE__, __LINE__);
			}
		} else ksize = g->ksize;
		if(bf_ksize){
			if(bf_ksize < 0) bf_ksize = 0;
			else if(bf_ksize != (int)g->bf_ksize){
				fprintf(stderr, " --  parameter ksize(%d) is set to the same in loaded DBG(%d) in %s -- %s:%d --\n", ksize, g->ksize, __FUNCTION__, __FILE__, __LINE__);
			}
		} else bf_ksize = g->bf_ksize;
		if(verbose) fprintf(stdout, "[%s] Done\n", date());
		if(g->ksize == g->bf_ksize){// it is only used in debug
			kmer_t *hk;
			uint64_t bk[2], *kb;
			uint32_t i, hidx;
			reset_iter_kmerhash(g->hashs[0]);
			while((hk = ref_iter_kmerhash(g->hashs[0]))){
				bk[0] = 0;
				for(i=0;i<g->ksize;i++) dna_shl_seqs(bk, g->bf_ksize, (hk->val >> ((g->ksize - 1 - i) << 1)) & 0x03U);
				bk[1] = bk[0];
				dna_rev_seqs(bk + 1, g->bf_ksize);
				if(dna_cmp_seqs(bk + 0, bk + 1, g->bf_ksize) < 0) kb = bk + 0;
				else kb = bk + 1;
				hidx = dbg_kmer_smear(kb[0]) % g->nh;
				if(!get_cbf(g->bfs[hidx], kb, 8)){
					fprintf(stdout, " -- Failture in checking bloomfilter in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
				}
			}
			fprintf(stdout, " -- Success in checking bloomfilter in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		}
	} else {
		if(ksize == 0) ksize = 25;
		if(bf_ksize == 0) bf_ksize = 45;
		if(srfs->size == 0){
			fprintf(stderr, " --  There is neither short reads files nor DBG file to be loaded  in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
			return usage();
		}
		if(verbose) fprintf(stdout, "[%s] Building DBG, %d threads\n", date(), ncpu);
		fr = fopen_m_filereader(srfs->size, srfs->buffer);
		g = init_dbg(ksize, bf_ksize, ncpu, min_cov);
		build_dbg(g, 64 * 1024 * 1024, fr);
		fclose_filereader(fr);
		if(verbose) fprintf(stdout, "[%s] Done\n", date());
	}
	free_cplist(srfs);
	if(bf_ksize < 0) bf_ksize = 0;
	if(dbg_dump){
		if(verbose) fprintf(stdout, "[%s] Dumping DBG into %s\n", date(), dbg_dump);
		if((dump = fopen(dbg_dump, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s for write in %s -- %s:%d --\n", dbg_dump, __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		if(dbg_load){
			mem_dump_obj_file(g, 1, &DBG_obj_desc, 1, 0, dump);
		} else {
			mem_dump_free_obj_file(g, 1, &DBG_obj_desc, 1, 0, dump);
		}
		fclose(dump);
		if(verbose) fprintf(stdout, "[%s] Done\n", date());
		return 0;
	}
	if(chkf){
			fr = fopen_filereader(chkf);
			seq = NULL;
			while(fread_seq(&seq, fr)){
				chk_sequences_dbg(g, seq->tag.string, seq->seq.string, seq->seq.size, cut_cov[0], cut_cov[1], stdout);
			}
			fclose_filereader(fr);
	} else {
		max_mem = max_mem * 1024 * 1024;
		cur_mem = mem_size_obj(g, 1, &DBG_obj_desc, 0, 1);
		if(verbose) fprintf(stdout, "DBG memory: %llu M\n", (unsigned long long)((cur_mem + 1024 * 1024 - 1) / (1024 * 1024)));
		if(verbose) fprintf(stdout, "Max correction memory: %llu M\n", (unsigned long long)((max_mem + 1024 * 1024 - 1) / (1024 * 1024)));
		max_mem = max_mem / (sizeof(dbg_dp_t) + sizeof(uint32_t) + sizeof(uint32_t) * 2 + (bf_ksize + 31) / 32 * 8); // nodes + heap + hash
		if(max_mem > DBG_MAX_BT_IDX - 1) max_mem = DBG_MAX_BT_IDX - 1;
		if(verbose){ fprintf(stdout, "Maximum search: %llu\n", (unsigned long long)(max_mem)); fflush(stdout); }
		if(optind == argc) fr = stdin_filereader();
		else fr = fopen_m_filereader(argc - optind, argv + optind);
		if(outf == 0 || strcmp(outf, "-") == 0){ out = stdout; left = NULL; }
		else if((out = fopen(outf, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s for write in %s -- %s:%d --\n", outf, __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		} else {
			leftf = malloc(strlen(outf) + 10);
			sprintf(leftf, "%s.intact", outf);
			left = fopen(leftf, "w");
			free(leftf);
		}
		thread_beg_init(maln, ncpu);
		maln->aln = init_dbgaln(g, Z, W, T, M, X, I, D, E, m, x, d, e);
		maln->aln->aln_mode = mode;
		maln->aln->min_link_cov = min_cov;
		maln->aln->link_cov_cutoff = cut_cov[0];
		maln->aln->kmer_cov_cutoff = cut_cov[1];
		maln->aln->lowcov_bonus = lowcov_bonus;
		maln->aln->all_kmer = all_kmer;
		maln->aln->max_nodes = max_mem;
		maln->aln->bfk_n = (bf_ksize + 31) / 32;
		maln->aln->smin = mode? 0 : - 0x0FFFFFFF;
		maln->tag = init_string(64);
		maln->qry = init_u8list(1024);
		for(i=0;i<7;i++) maln->qvs[i] = init_u8list(1024);
		maln->has_5q = 0;
		maln->cns = init_u8list(1024);
		maln->cigars = init_u32list(1024);
		maln->is_corrected = 0;
		maln->corrs = init_u32list(8);
		maln->min_hsp_len = min_hsp;
		maln->max_hsp_ovl = ovl_hsp;
		thread_end_init(maln);
		if(verbose) fprintf(stdout, "[%s] Correcting, %d threads\n", date(), ncpu);
		seq = NULL;
		alns[0] = init_string(1024);
		alns[1] = init_string(1024);
		eof_fr = 0;
		while((eof_fr? 0 : fread_seq(&seq, fr)) || !thread_test_all(maln, maln->qry->size == 0)){
			eof_fr = (seq == NULL);
			thread_wait_one(maln);
			if(maln->qry->size){
				if(maln->is_corrected){
					for(num=0;num<maln->corrs->size;num++){
						hsp = maln->aln->hsps->buffer + maln->corrs->buffer[num];
						dbg_hit_t rs;
						rs = call_correct_seq(maln->aln, hsp->idx, maln->cns, maln->cigars);
						fprintf(out, ">%s/%d_%d corr=Y ori_len=%d len=%d seed=%d score=%d mat=%d mis=%d ins=%d del=%d\n", maln->tag->string, rs.qb, rs.qe, (int)maln->qry->size, (int)maln->cns->size,
							maln->aln->qmers->buffer[hsp->idx].qpos, hsp->score, rs.alns[0], rs.alns[3], rs.alns[1], rs.alns[2]);
						for(i=0;i<(int)maln->cns->size;i++) fputc("ACGTacgtN"[maln->cns->buffer[i]], out);
						fputc('\n', out);
						fflush(out);
						if(verbose){
							fprintf(stdout, "[%s] %s %d-%d corr=Y ori_len=%d len=%d seed=%d score=%d mat=%d mis=%d ins=%d del=%d\n", date(), maln->tag->string, rs.qb, rs.qe, (int)maln->qry->size, (int)maln->cns->size,
								maln->aln->qmers->buffer[hsp->idx].qpos, hsp->score, rs.alns[0], rs.alns[3], rs.alns[1], rs.alns[2]);
							if(verbose > 2){
								clear_string(alns[0]);
								clear_string(alns[1]);
								for(i=0;i<(int)maln->cns->size;i++) maln->cns->buffer[i] &= 0x03;
								kswx_cigar2pairwise(alns, (uint8_t*[2]){maln->cns->buffer, maln->qry->buffer + rs.qb}, maln->cigars->size, maln->cigars->buffer);
								for(i=0;i<(int)alns[0]->size;i+=100){
									char ct;
									end = i + 100;
									if(end > (int)alns[0]->size) end = alns[0]->size;
									ct = alns[0]->string[end]; alns[0]->string[end] = 0;
									fprintf(stdout, "%s\n", alns[0]->string + i);
									alns[0]->string[end] = ct;
									for(j=i;j<end;j++){
										if(alns[0]->string[j] == '-') fputc('-', stdout);
										else if(alns[1]->string[j] == '-') fputc('-', stdout);
										else if(alns[0]->string[j] == alns[1]->string[j]) fputc('|', stdout);
										else fputc('*', stdout);
									}
									fputc('\n', stdout);
									ct = alns[1]->string[end]; alns[1]->string[end] = 0;
									fprintf(stdout, "%s\n\n", alns[1]->string + i);
									alns[1]->string[end] = ct;
								}
							}
							fflush(stdout);
						}
					}
				} else {
					if(left){
						fprintf(left, ">%s corr=N ori_len=%d\n", maln->tag->string, (int)maln->qry->size);
						for(i=0;i<(int)maln->qry->size;i++) fputc("ACGTacgtN"[maln->qry->buffer[i]], left);
						fputc('\n', left);
						fflush(left);
					}
					if(verbose){
						fprintf(stdout, "[%s] %s corr=N ori_len=%d\n", date(), maln->tag->string, (int)maln->qry->size);
					}
				}
			}
			clear_u8list(maln->qry);
			if(seq && seq->seq.size){
				if(verbose){
					fprintf(stdout, "[%s] %s ori_len=%d\n", date(), seq->name.string, seq->seq.size);
					fflush(stdout);
				}
				clear_string(maln->tag); append_string(maln->tag, seq->name.string, seq->name.size);
				for(i=0;i<seq->seq.size;i++) push_u8list(maln->qry, base_bit_table[(int)seq->seq.string[i]]);
				maln->has_5q = 0;
				if(use_5q && seq->qual.size == 7 * seq->seq.size){
					for(i=0;i<5*seq->seq.size;i++) seq->qual.string[i] -= 33;
					for(;i<seq->qual.size;i++) seq->qual.string[i] = base_bit_table[(int)seq->qual.string[i]];
					for(i=0;i<7;i++) clear_u8list(maln->qvs[i]);
					for(i=0;i<7;i++) append_array_u8list(maln->qvs[i], (uint8_t*)(seq->qual.string + i * seq->seq.size), seq->seq.size);
					// curate deletion QV for convenience
					push_u8list(maln->qvs[6], 4);
					push_u8list(maln->qvs[3], maln->aln->QDEL);
					set_u8list(maln->qvs[3], 0, maln->aln->QDEL);
					for(i=0;i<seq->seq.size;i++) maln->qvs[6]->buffer[i] = 4;
					maln->has_5q = 1;
				}
				thread_wake(maln);
			}
		}
		fclose_filereader(fr);
		fclose(out);
		free_string(alns[0]);
		free_string(alns[1]);
		thread_beg_close(maln);
		free_dbgaln(maln->aln);
		free_string(maln->tag);
		free_u8list(maln->qry);
		for(i=0;i<7;i++) free_u8list(maln->qvs[i]);
		free_u8list(maln->cns);
		free_u32list(maln->cigars);
		free_u32list(maln->corrs);
		thread_end_close(maln);
	}
	if(verbose) fprintf(stdout, "[%s] Done\n", date());
	if(!dbg_load) free_dbg(g);
	else if(!dbg_share) free(g);
	return 0;
}

