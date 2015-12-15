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
 
#include <sys/stat.h>
#include <sys/mman.h>
#include "timer.h"
#include "pbcorr_dbg.h"
#include "file_reader.h"
#include "thread.h"

#define PBM_MAX_N_SEQ	512

thread_beg_def(midx);
kmerv *job;
khash *dbg;
String *seqs;
kmerv **kmers;
uint8_t ksize;
uint64_t kmask;
thread_end_def(midx);

thread_beg_func(midx);
int off, len, i, exists;
uint32_t m;
kmer_t *k, *kmer, KMER;
uint64_t K, R;
char *seq;
thread_beg_loop(midx);
for(m=0;m<midx->job->size;m++){
	kmer = ref_kmerv(midx->job, m);
	k = prepare_khash(midx->dbg, *kmer, &exists);
	if(exists){
		k->lnk |= kmer->lnk;
		if(k->cnt < MAX_KCNT) k->cnt ++;
	} else {
		k->kmer = kmer->kmer;
		k->lnk  = kmer->lnk;
		k->cnt = 1;
	}
}
clear_kmerv(midx->job);
KMER.cnt = 1;
for(off=0;off<midx->seqs->size;){
	seq = midx->seqs->string + off;
	len = strlen(seq);
	off += len + 1;
	beg_seq2kmers(seq, len, midx->ksize, midx->kmask, K, i);
	R = dna_rev_seq(K, midx->ksize);
	KMER.lnk = 0;
	if(R > K){
		if(i){ KMER.lnk |= (1U << 4) << ((~base_bit_table[(int)seq[i - 1]]) & 0x03); }
		if(i + midx->ksize < len){ KMER.lnk |= (1U) << base_bit_table[(int)seq[i + midx->ksize]]; }
		KMER.kmer = K;
	} else {
		if(i){ KMER.lnk |= (1U) << ((~base_bit_table[(int)seq[i - 1]]) & 0x03); }
		if(i + midx->ksize < len){ KMER.lnk |= (1U << 4) << base_bit_table[(int)seq[i + midx->ksize]]; }
		KMER.kmer = R;
	}
	push_kmerv(midx->kmers[kmer_idx(KMER.kmer) % midx->n_cpu], KMER);
	end_seq2kmers;
}
thread_end_loop(midx);
thread_end_func(midx);

khash** build_seqdbg(uint8_t ksize, FileReader *fr, int ncpu){
	Sequence *seq;
	khash **dbgs;
	kmerv **kmers;
	uint64_t n_rd;
	int i, idle;
	thread_preprocess(midx);
	kmers = malloc(sizeof(kmerv*) * ncpu);
	for(i=0;i<ncpu;i++) kmers[i] = init_kmerv(1024);
	dbgs = malloc(sizeof(khash*) * ncpu);
	for(i=0;i<ncpu;i++) dbgs[i] = init_khash(1024);
	thread_beg_init(midx, ncpu);
	midx->job = init_kmerv(1024);
	midx->dbg = dbgs[midx->t_idx];
	midx->seqs = init_string(1024);
	midx->kmers = malloc(sizeof(kmerv*) * ncpu);
	for(i=0;i<ncpu;i++) midx->kmers[i] = init_kmerv(1024);
	midx->ksize = ksize;
	midx->kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
	thread_end_init(midx);
	seq = NULL;
	n_rd = 0;
	idle = 0;
	while(idle < 2 * ncpu){
		thread_waitfor_one_idle(midx);
		for(i=0;i<ncpu;i++){ append_kmerv(kmers[i], midx->kmers[i]); clear_kmerv(midx->kmers[i]); }
		append_kmerv(midx->job, kmers[midx->t_idx]); clear_kmerv(kmers[midx->t_idx]);
		clear_string(midx->seqs);
		for(i=0;i<PBM_MAX_N_SEQ;i++){
			if(fread_seq_adv(&seq, fr, SEQ_FLAG_NO_NAME | SEQ_FLAG_NO_QUAL) == 0) break;
			append_string(midx->seqs, seq->seq.string, seq->seq.size + 1);
			n_rd ++;
			if((n_rd & 0xFFFFFU) == 0){
				fprintf(stderr, "[%s] indexed %llu reads\n", date(), (unsigned long long)n_rd); fflush(stderr);
			}
		}
		thread_wake(midx);
		if(i == 0) idle ++;
	}
	thread_waitfor_all_idle(midx);
	fprintf(stderr, "[%s] finished %llu\n", date(), (unsigned long long)n_rd); fflush(stderr);
	for(i=0;i<ncpu;i++) free_kmerv(kmers[i]);
	free(kmers);
	thread_beg_close(midx);
	free_kmerv(midx->job);
	for(i=0;i<ncpu;i++) free_kmerv(midx->kmers[i]);
	free(midx->kmers);
	free_string(midx->seqs);
	thread_end_close(midx);
	return dbgs;
}

int exists_seqdbg(khash **dbg, uint8_t n_dbg, uint8_t ksize, uint64_t kmask, char *seq, uint32_t len){
	kmer_t K, *k;
	uint64_t kmer, krev;
	uint32_t i, ret;
	ret = 0;
	K.kmer = 0;
	K.lnk = 0;
	K.cnt = 0;
	fprintf(stderr, "%s\n", seq);
	beg_seq2kmers(seq, len, ksize, kmask, kmer, i);
	krev = dna_rev_seq(kmer, ksize);
	if(krev > kmer) K.kmer = kmer;
	else K.kmer = krev;
	k = get_khash(dbg[kmer_idx(K.kmer) % n_dbg], K);
	if(k){
		if(k->cnt < 10) fprintf(stderr, "%d", k->cnt);
		else fprintf(stderr, "+");
		ret ++;
	} else {
		fprintf(stderr, "^");
	}
	end_seq2kmers;
	for(i=0;i+1<ksize;i++) fprintf(stderr, "0");
	fprintf(stderr, "\n");
	return ret + ksize - 1;
}

define_list(kmerpv, kmer_t*);

void walk_along_seqdbg(khash **dbg, uint8_t n_dbg, uint8_t ksize, uint64_t kmask, kmerpv *kmers, u8list *dirs, uint32_t max_len){
	kmer_t *k, *n, K;
	uint64_t f, r;
	uint8_t dir, ndir, lnk, v;
	K.cnt = K.lnk = 0;
	while(kmers->size < max_len){
		k = get_kmerpv(kmers, kmers->size - 1);
		dir = get_u8list(dirs, dirs->size - 1);
		if(dir){
			lnk = k->lnk >> 4;
			f = dna_rev_seq(k->kmer, ksize);
		} else {
			lnk = k->lnk & 0x0FU;
			f = k->kmer;
		}
		if(byte_ones_table[lnk] != 1) break;
		v = __builtin_ctz(lnk);
		f = ((f << 2) | v) & kmask;
		r = dna_rev_seq(f, ksize);
		if(r < f){
			K.kmer = r;
			ndir = 1;
		} else {
			K.kmer = f;
			ndir = 0;
		}
		n = get_khash(dbg[kmer_idx(K.kmer) % n_dbg], K);
		if(n == NULL) break;
		lnk = ndir? (n->lnk & 0x0FU) : (n->lnk >> 4) ;
		if(byte_ones_table[lnk] != 1) break;
		push_kmerpv(kmers, n);
		push_u8list(dirs, ndir);
	}
}

size_t remove_tips_seqdbg(khash **dbg, uint8_t n_dbg, uint8_t ksize){
	kmer_t *k;
	kmerpv *kmers;
	u8list *dirs;
	size_t ret;
	uint32_t i;
	uint8_t idx;
	kmers = init_kmerpv(ksize * 2 + 1);
	dirs  = init_u8list(ksize * 2 + 1);
	ret = 0;
	for(idx=0;idx<n_dbg;idx++){
		reset_iter_khash(dbg[idx]);
		while((k = ref_iter_khash(dbg[idx]))){
			if(k->cnt == 0) continue;
			if((k->lnk & 0xF0U) && (k->lnk & 0x0FU)) continue;
			clear_kmerpv(kmers);
			if(k->lnk & 0x0FU){
				push_kmerpv(kmers, k);
				push_u8list(dirs, 0);
			} else if(k->lnk & 0xF0){
				push_kmerpv(kmers, k);
				push_u8list(dirs, 1);
			} else continue;
			walk_along_seqdbg(dbg, n_dbg, ksize, 0xFFFFFFFFFFFFFFFFLLU >> (64 - 2 * ksize), kmers, dirs, 2 * ksize + 1);
			if(kmers->size > 2 * ksize) continue;
			ret ++;
			for(i=0;i<kmers->size;i++) get_kmerpv(kmers, i)->cnt = 0;
		}
	}
	free_kmerpv(kmers);
	free_u8list(dirs);
	return ret;
}

DBGAligner* init_dbgaln(khash **dbg, uint8_t n_dbg, uint8_t ksize, uint8_t kfreq, OnlineSWParam *param, int min_score){
	DBGAligner *aln;
	aln = malloc(sizeof(DBGAligner));
	aln->dbg = dbg;
	aln->n_dbg = n_dbg;
	aln->ksize = ksize;
	aln->kmask = 0xFFFFFFFFFFFFFFFFLLU >> (64 - ksize * 2);
	aln->kfreq = kfreq;
	aln->sw = init_olsw(param, min_score, param->mat * 2000);
	aln->hash = init_u64hash(1023);
	aln->kmers = init_u64hash(1023);
	aln->min_score = min_score;
	aln->s_best = 0;
	aln->z_best = 50;
	aln->fast_mode = 0;
	aln->verbose = 0;
	aln->query = NULL;
	aln->qlen  = 0;
	aln->qdir  = 0;
	aln->qvals = init_u8list(1024);
	aln->q[0] = init_string(1024);
	aln->q[1] = init_string(1024);
	aln->t[0] = init_string(1024);
	aln->t[1] = init_string(1024);
	aln->c[0] = init_string(1024);
	aln->c[1] = init_string(1024);
	aln->e = init_string(1024);
	return aln;
}

void free_dbgaln(DBGAligner *aln){
	free_olsw(aln->sw);
	free_u64hash(aln->hash);
	free_u64hash(aln->kmers);
	free_u8list(aln->qvals);
	free_string(aln->q[0]);
	free_string(aln->q[1]);
	free_string(aln->c[0]);
	free_string(aln->c[1]);
	free_string(aln->t[0]);
	free_string(aln->t[1]);
	free_string(aln->e);
	free(aln);
}

uint32_t align_dbgaln_core(DBGAligner *aln, uint64_t kmer, uint8_t lnk, uint32_t qpos, int qdir, int *score){
	kmer_t MER, *mer;
	ol_sw_t *cur, *nxt;
	uint64_t r, k, i, *e;
	uint32_t cur_idx, nxt_idx, max_idx;
	int j, max_len, s, S, exist, hit, n_hit, rev, step;
	if(qdir){
		if(qpos == 0) return 0xFFFFFFFFU;
		if(beg_olsw(aln->sw, aln->qvals->buffer, qpos, 1, aln->ksize * aln->sw->param->mat) == 0) return 0xFFFFFFFFU;
	} else {
		if(qpos + aln->ksize >= aln->qlen) return 0xFFFFFFFFU;
		if(beg_olsw(aln->sw, aln->qvals->buffer + qpos + aln->ksize, aln->qlen - (qpos + aln->ksize), 0, aln->ksize * aln->sw->param->mat) == 0) return 0xFFFFFFFFU;
	}
	cur = ref_ol_swv(aln->sw->steps, 0);
	cur->aux = kmer << 4 | lnk;
	max_idx = 0xFFFFFFFFU;
	max_len = 0;
	n_hit = aln->sw->param->bw * 2 * 1000;
	hit = 0;
	MER.lnk = 0;
	MER.cnt = 0;
	step = 0;
	while(beg_next_olsw(aln->sw) && aln->sw->q_off < (int)aln->sw->tlen + 1){
		j = 0;
		S = 0;
		if((step % 1) == 0) clear_u64hash(aln->kmers);
		step ++;
		while((cur_idx = iter_olsw(aln->sw))!=0xFFFFFFFFU){
			cur = ref_ol_swv(aln->sw->steps, cur_idx);
			e = prepare_u64hash(aln->kmers, cur->aux, &exist);
			if(exist) continue;
			*e = cur->aux;
			s = max_score_olsw(aln->sw, cur_idx);
			if(s < aln->min_score) break;
			j ++;
			if(s > S) S = s;
			if(s + aln->s_best <= S && j > aln->z_best) break;
			if(cur->t_off + cur->max + 1 > max_len){
				max_len = cur->t_off + cur->max + 1;
				max_idx = cur_idx;
				*score  = s;
				if(max_len >= (int)aln->sw->tlen) hit = 1;
			} else if(cur->t_off + cur->max + 1 == max_len && s > *score){
				max_idx = cur_idx;
				*score  = s;
			}
			kmer = cur->aux >> 4;
			lnk  = cur->aux & 0x0FU;
			for(i=0;i<4;i++){
				if(((lnk >> i) & 0x01) == 0) continue;
				k = ((kmer << 2) | i) & aln->kmask;
				r = dna_rev_seq(k, aln->ksize);
				if(r > k){ MER.kmer = k; rev = 0; }
				else { MER.kmer = r; rev = 1; }
				mer = get_khash(aln->dbg[kmer_idx(MER.kmer) % aln->n_dbg], MER);
				if(mer == NULL){
					continue;
					fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
				} else if(mer->cnt < aln->kfreq) continue;
				nxt_idx = push_olsw(aln->sw, cur_idx, bit_base_table[qdir? ((~i) & 0x03) : i], (mer->cnt > 2 * aln->kfreq)? 1 : 0);
				nxt = ref_ol_swv(aln->sw->steps, nxt_idx);
				nxt->aux = (k << 4) | (rev? (mer->lnk >> 4) : (mer->lnk & 0x0F));
			}
		}
		if(hit){
			n_hit --;
			if(n_hit <= 0) break;
		}
	}
	return max_idx;
}

void cache_seq2kmerhash(DBGAligner *aln, char *seq, int len){
	uint64_t kmer, krev;
	int i;
	if(len <= 0) return;
	beg_seq2kmers(seq, len, aln->ksize, aln->kmask, kmer, i);
	krev = dna_rev_seq(kmer, aln->ksize);
	if(krev > kmer) krev =  kmer;
	put_u64hash(aln->hash, krev);
	end_seq2kmers;
}

int align_dbgaln(DBGAligner *aln, char *query, uint32_t qlen, pbhitv *hits){
	kmer_t K, *k;
	PBHit *hit;
	uint64_t kmer, krev;
	uint32_t i, idx1, idx2, j;
	int score[2], off, tlen, len[2], rev;
	char seq[32];
	aln->query = query;
	aln->qlen = qlen;
	clear_and_encap_u8list(aln->qvals, qlen);
	for(i=0;i<qlen;i++) aln->qvals->buffer[i] = base_bit_table[(int)query[i]];
	clear_u64hash(aln->hash);
	K.lnk = 0;
	K.cnt = 0;
	if(aln->verbose) exists_seqdbg(aln->dbg, aln->n_dbg, aln->ksize, aln->kmask, query, qlen);
	beg_seq2kmers(query, qlen, aln->ksize, aln->kmask, kmer, i);
	krev = dna_rev_seq(kmer, aln->ksize);
	if(krev > kmer){ K.kmer = kmer; rev = 0; }
	else { K.kmer = krev; rev = 1; }
	k = get_khash(aln->dbg[kmer_idx(K.kmer) % aln->n_dbg], K);
	if(k == NULL || k->cnt < aln->kfreq) continue;
	if(exists_u64hash(aln->hash, K.kmer)) continue;
	if(aln->verbose) fprintf(stderr, "Try kmer %d\n", i);

	score[0] = 0;
	idx1 = align_dbgaln_core(aln, krev, rev? (k->lnk & 0x0F) : (k->lnk >> 4), i, 1, score + 0);
	if(idx1 == 0xFFFFFFFFU) continue;
	clear_string(aln->q[0]);
	clear_string(aln->t[0]);
	clear_string(aln->c[0]);
	tostring_olsw(aln->sw, idx1, aln->q[0], aln->t[0], aln->c[0]);
	for(j=0;j<aln->ksize;j++){
		add_char_string(aln->q[0], query[i + j]);
		add_char_string(aln->t[0], query[i + j]);
		add_char_string(aln->c[0], '|');
	}
	off = ((int)i) - (aln->sw->steps->buffer[idx1].t_off + aln->sw->steps->buffer[idx1].max + 1);
	tlen = (aln->sw->steps->buffer[idx1].t_off + aln->sw->steps->buffer[idx1].max + 1) + aln->ksize;
	len[0] = (aln->sw->steps->buffer[idx1].t_off + aln->sw->steps->buffer[idx1].max + 1);

	idx2 = align_dbgaln_core(aln, kmer, rev? (k->lnk >> 4) : (k->lnk & 0x0F), i, 0, score + 1);
	if(idx2 != 0xFFFFFFFFU){
		clear_string(aln->q[1]);
		clear_string(aln->t[1]);
		clear_string(aln->c[1]);
		tostring_olsw(aln->sw, idx2, aln->q[1], aln->t[1], aln->c[1]);
		reverse_string(aln->q[1]);
		reverse_string(aln->t[1]);
		reverse_string(aln->c[1]);
		append_string(aln->q[0], aln->q[1]->string, aln->q[1]->size);
		append_string(aln->c[0], aln->c[1]->string, aln->c[1]->size);
		append_string(aln->t[0], aln->t[1]->string, aln->t[1]->size);
		tlen += (aln->sw->steps->buffer[idx2].t_off + aln->sw->steps->buffer[idx2].max + 1);
		len[1] = (aln->sw->steps->buffer[idx2].t_off + aln->sw->steps->buffer[idx2].max + 1);
	} else {
		score[1] = aln->ksize * aln->sw->param->mat;
		len[1] = 0;
	}
	clear_string(aln->e);
	tidy_string(aln->q[0], aln->e, '-');
	cache_seq2kmerhash(aln, aln->e->string, aln->e->size);
	hit = next_ref_pbhitv(hits);
	hit->off  = off;
	hit->qlen = aln->e->size;
	hit->tlen = tlen;
	hit->score = score[0] + score[1] - aln->ksize * aln->sw->param->mat;
	hit->q = clone_string(aln->q[0]);
	hit->t = clone_string(aln->t[0]);
	hit->c = clone_string(aln->c[0]);
	hit->closed = 0;
	if(aln->verbose){
		kmer2seq(seq, kmer, aln->ksize);
		fprintf(stderr, "%d %d {%d %d %d} kmer[%d %s] - trace[%u+%u] score=%d [%d+%d+%d] qlen=%d tlen=%d [%d+%d+%d]\n", off, off + tlen, 
				occ_str(aln->c[0]->string, aln->c[0]->size, 'x'),
				occ_str(aln->q[0]->string, aln->q[0]->size, '-'),
				occ_str(aln->t[0]->string, aln->t[0]->size, '-'),
				i, seq, idx1, idx2, score[0] + score[1] - aln->ksize * aln->sw->param->mat,
				score[0] - aln->ksize * aln->sw->param->mat, aln->ksize * aln->sw->param->mat, score[1] - aln->ksize * aln->sw->param->mat,
				aln->e->size, tlen, len[0], aln->ksize, len[1]);
		fprintf(stderr, "exists %d\n", exists_seqdbg(aln->dbg, aln->n_dbg, aln->ksize, aln->kmask, aln->e->string, aln->e->size));
		fprintf(stderr, "Q\t%d\t%s\n", aln->e->size, aln->e->string);
		clear_string(aln->e);
		tidy_string(aln->t[0], aln->e, '-');
		//fprintf(stderr, "exists %d\n", exists_seqdbg(aln->dbg, aln->ksize, aln->kmask, aln->e->string, aln->e->size));
		fprintf(stderr, "T\t%d\t%s\n", aln->e->size, aln->e->string);
		fprintf(stderr, "%s\n", aln->q[0]->string);
		fprintf(stderr, "%s\n", aln->c[0]->string);
		fprintf(stderr, "%s\n", aln->t[0]->string);
		fflush(stderr);
	}
	if(aln->fast_mode && (int)tlen >= (int)qlen * 0.99) break;
	end_seq2kmers;
	return 0;
}

PBCorr* init_pbcorr(khash **dbg, uint8_t n_dbg, uint8_t ksize, uint8_t kfreq, OnlineSWParam *param, int min_score){
	PBCorr *corr;
	corr = malloc(sizeof(PBCorr));
	corr->aln = init_dbgaln(dbg, n_dbg, ksize, kfreq, param, min_score);
	corr->hits = init_pbhitv(32);
	corr->seqs = init_string(1024);
	corr->q = init_string(1024);
	corr->t = init_string(1024);
	corr->c = init_string(1024);
	corr->covs = init_bitvec(1024);
	return corr;
}

void free_pbcorr(PBCorr *corr){
	free_dbgaln(corr->aln);
	free_pbhitv(corr->hits);
	free_bitvec(corr->covs);
	free_string(corr->seqs);
	free_string(corr->q);
	free_string(corr->t);
	free_string(corr->c);
	free(corr);
}

void pbcorr(PBCorr *corr, char *header, char *query, uint32_t qlen){
	PBHit *hit;
	uint32_t i;
	int j, f, off, off2;
	if(corr->aln->verbose) fprintf(stderr, "#%s\t%d\t%s\n", header, qlen, query);
	align_dbgaln(corr->aln, query, qlen, corr->hits);
	sort_array(corr->hits->buffer, corr->hits->size, PBHit, (b.tlen > a.tlen + 50)? 1 : ((b.tlen + 50 >= a.tlen)? (b.score - a.score) : -1));
	encap_bitvec(corr->covs, qlen);
	zeros_bitvec(corr->covs);
	for(i=0;i<corr->hits->size;i++){
		hit = ref_pbhitv(corr->hits, i);
		f = 0;
		for(j=0;j<hit->tlen;j++){
			if(get_bitvec(corr->covs, j + hit->off) == 0){
				one_bitvec(corr->covs, j + hit->off);
				f = 1;
			}
		}
		if(f == 0) hit->closed = 1;
	}
	sort_array(corr->hits->buffer, corr->hits->size, PBHit, a.off > b.off);
	clear_string(corr->seqs);
	clear_string(corr->q);
	clear_string(corr->t);
	clear_string(corr->c);
	off = 0;
	for(i=0;i<corr->hits->size;i++){
		hit = ref_pbhitv(corr->hits, i);
		if(hit->closed) continue;
		if(off >= hit->off + hit->tlen) continue;
		while(off < hit->off){
			add_char_string(corr->q, lc(query[off]));
			add_char_string(corr->t, 'N');
			add_char_string(corr->c, 'N');
			off ++;
		}
		if(corr->aln->verbose){
			fprintf(stderr, "@ %d\t%d\n", hit->off, hit->tlen);
		}
		for(off2=j=0;j<hit->t->size;j++){
			if(off2 == off - hit->off) break;
			if(hit->t->string[j] == '-') continue;
			off2 ++;
		}
		append_string(corr->q, hit->q->string + j, hit->q->size - j);
		append_string(corr->t, hit->t->string + j, hit->t->size - j);
		append_string(corr->c, hit->c->string + j, hit->c->size - j);
		off = hit->off + hit->tlen;
	}
	while(off < (int)qlen){
		add_char_string(corr->q, lc(query[off]));
		add_char_string(corr->t, 'N');
		add_char_string(corr->c, 'N');
		off ++;
	}
	for(i=0;i<corr->hits->size;i++){
		hit = ref_pbhitv(corr->hits, i);
		free_string(hit->q);
		free_string(hit->t);
		free_string(hit->c);
	}
	clear_pbhitv(corr->hits);
	tidy_string(corr->q, corr->seqs, '-');
	//fprintf(out, ">%s\n%s\n", header, corr->seqs->string);
	if(corr->aln->verbose){
		fprintf(stderr, "#Q %s\n", corr->q->string);
		fprintf(stderr, "#C %s\n", corr->c->string);
		fprintf(stderr, "#T %s\n", corr->t->string);
	}
}

thread_beg_def(mcorr);
Sequence *seq;
PBCorr *corr;
thread_end_def(mcorr);

thread_beg_func(mcorr);
thread_beg_loop(mcorr);
pbcorr(mcorr->corr, mcorr->seq->header.string, mcorr->seq->seq.string, mcorr->seq->seq.size);
thread_end_loop(mcorr);
thread_end_func(mcorr);

int usage(){
	fprintf(stdout,
"Usage: pbcorr-dbg [command]\n"
" index [options] <files> # Build DBG index\n"
"  -t <int>    Number of threads, [8]\n"
"  -o <string> Name of BF index file, [my.dbg]\n"
"  -k <int>    Kmer size, MUST be odd, <= 31 [21]\n"
" load <bf_index> # load bf_index into memory\n"
" corr [options] <bf_index>\n"
"  -t <int>    Number of threads, [8]\n"
"  -i <string> Input pacbio sequences file, fa/fq, [STDIN]\n"
"  -o <string> Output file, [STDOUT]\n"
"  -f <int>    Minimum frequency of kmer, [2]\n"
"  -v          Verbose mode, print to stderr\n"
"  -F          Fast mode: don't try to find better correction when having one cover >= 99%% bases of query\n"
"  -Z <int>    Z-best [100]\n"
"  -S <int>    distance to the top score, [0]\n"
"              Branches satisfy score distance will ignore z-best limitation\n"
"  -B <int>    Bandwidth in Smith-Waterman alignment, <= 64 [16]\n"
"  -M <int>    Score match, [1]\n"
"  -X <int>    Score mismatch, [-4]\n"
"  -I <int>    Score insertion, [-2]\n"
"  -D <int>    Score deletion, [-3]\n"
);
	return 1;
}

static obj_desc_t khashv_obj_desc = {sizeof(khash*), 1, {0}, {(obj_desc_t*)&khash_obj_desc}, NULL};

int index_invoker(int argc, char **argv){
	khash **dbg;
	FileReader *fr;
	size_t aux, cnt;
	uint8_t ksize;
	char *out_file;
	FILE *out;
	int c, i, ncpu;
	ksize = 21;
	ncpu = 8;
	out_file = NULL;
	while((c = getopt(argc, argv, "ht:k:o:m:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'o': out_file = optarg; break;
			default: printf("Unknown option -%c\n", c); return usage();
		}
	}
	if(optind == argc) return usage();
	if(out_file == NULL) out_file = "my.dbg";
	fr = fopen_m_filereader(argc - optind, argv + optind);
	if((out = fopen(out_file, "w")) == NULL){
		fprintf(stderr, " -- Cannot write %s in %s -- %s:%d --\n", out_file, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	dbg = build_seqdbg(ksize, fr, ncpu);
	cnt = remove_tips_seqdbg(dbg, ncpu, ksize);
	fprintf(stderr, "[%s] remove %lu tips\n", date(), (unsigned long)cnt);
	aux = ksize;
	cnt = ncpu;
	mem_dump_obj_file(dbg, &khashv_obj_desc, cnt, aux, out);
	fclose(out);
	fprintf(stderr, "[%s] wrote %s\n", date(), out_file);
	for(i=0;i<ncpu;i++) free_khash(dbg[i]);
	free(dbg);
	return 0;
	
}

int load_invoker(int argc, char **argv){
	void *addr;
	size_t cnt, aux;
	if(argc < 2) return usage();
	addr = mem_load_obj_file(&khashv_obj_desc, argv[1], &cnt, &aux);
	fprintf(stdout, "[%s] shared dbg_index \"%s\"[ksize=%d, n_dbg=%d] at %p\n", date(), argv[1], (int)aux, (int)cnt, addr);
	fprintf(stdout, "kill this process to unshare\n");
	while(1){
		microsleep(100);
	}
	return 0;
}

int corr_invoker(int argc, char **argv){
	PBCorr *corr;
	khash **dbg;
	OnlineSWParam param;
	FileReader *fr;
	char *inf, *ouf;
	size_t aux, cnt;
	uint8_t ksize, kfreq;
	FILE *out;
	int lo, hi, c, ncpu, s_best, z_best, verbose, fast;
	thread_preprocess(mcorr);
	param = (OnlineSWParam){1, -4, -2, -3, 16};
	inf = ouf = NULL;
	lo = 0;
	hi = 1024;
	ncpu = 8;
	kfreq = 1;
	s_best = 0;
	z_best = 100;
	verbose = 0;
	fast = 0;
	while((c = getopt(argc, argv, "ht:i:o:f:FvB:Z:S:M:X:I:D:L:H:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'i': inf = optarg; break;
			case 'o': ouf = optarg; break;
			case 'f': kfreq = atoi(optarg); break;
			case 'F': fast = 1; break;
			case 'v': verbose = 1; break;
			case 'B': param.bw = atoi(optarg); break;
			case 'Z': z_best = atoi(optarg); break;
			case 'M': param.mat = atoi(optarg); break;
			case 'X': param.mis = atoi(optarg); break;
			case 'I': param.ins = atoi(optarg); break;
			case 'D': param.del = atoi(optarg); break;
			case 'L': lo = atoi(optarg); break;
			case 'H': hi = atoi(optarg); break;
			case 'S': s_best = atoi(optarg); break;
			default: printf("Unknown option -%c\n", c); return usage();
		}
	}
	if(optind + 1 != argc) return usage();
	dbg = (khash**)mem_load_obj_file(&khashv_obj_desc, argv[optind], &cnt, &aux);
	ksize = aux;
	if(inf == NULL) inf = "-";
	if((fr = fopen_filereader(inf)) == NULL){
		fprintf(stderr, " -- Cannot read %s in %s -- %s:%d --\n", inf, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	if((out = ouf? fopen(ouf, "w") : stdout) == NULL){
		fprintf(stderr, " -- Cannot write %s in %s -- %s:%d --\n", ouf, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	thread_beg_init(mcorr, ncpu);
	corr = init_pbcorr(dbg, cnt, ksize, kfreq, &param, lo);
	corr->aln->z_best = z_best;
	corr->aln->s_best = s_best;
	corr->aln->verbose = verbose;
	corr->aln->fast_mode = fast;
	mcorr->corr = corr;
	mcorr->seq  = NULL;
	thread_end_init(mcorr);
	while(1){
		thread_waitfor_one_idle(mcorr);
		if(mcorr->seq){ fprintf(out, ">%s\n%s\n", mcorr->seq->header.string, mcorr->corr->seqs->string); fflush(out); }
		if((fread_seq(&mcorr->seq, fr))){
			thread_wake(mcorr);
		} else {
			break;
		}
	}
	thread_waitfor_all_idle(mcorr);
	thread_beg_close(mcorr);
	if(mcorr->seq){
		fprintf(out, ">%s\n%s\n", mcorr->seq->header.string, mcorr->corr->seqs->string); fflush(out);
		free_sequence(mcorr->seq);
	}
	free_pbcorr(mcorr->corr);
	thread_end_close(mcorr);
	fclose_filereader(fr);
	if(ouf) fclose(out);
	return 0;
}

int main(int argc, char **argv){
	if(argc < 2) return usage();
	if(strcasecmp(argv[1], "index") == 0){
		return index_invoker(argc -1 , argv + 1);
	} else if(strcasecmp(argv[1], "load") == 0){
		return load_invoker(argc -1, argv + 1);
	} else if(strcasecmp(argv[1], "corr") == 0){
		return corr_invoker(argc -1, argv + 1);
	}
	return 1;
}
