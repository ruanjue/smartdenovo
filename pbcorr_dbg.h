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
 
#ifndef __PACBIO_ERROR_CORRECTION_DBG_RJ_H
#define __PACBIO_ERROR_CORRECTION_DBG_RJ_H

#include "OnlineSW.h"
#include "hashset.h"

#define MAX_KSIZE	25
#define MAX_KCNT	63

typedef struct {
	uint64_t kmer:50, lnk:8, cnt:6;
} kmer_t;
#define kmer_hashcode(K) u64hashcode((K).kmer)
#define kmer_equals(K1, K2) ((K1).kmer == (K2).kmer)
define_hashset(khash, kmer_t, kmer_hashcode, kmer_equals);
#define kmer_idx(v) u32hashcode(v)
define_list(kmerv, kmer_t);

typedef struct {
	khash    **dbg;
	uint8_t  ksize, kfreq, n_dbg;
	uint64_t kmask;
	int min_score;
	int s_best, z_best, fast_mode;
	int verbose;
	OnlineSW *sw;
	u64hash  *hash;
	u64hash  *kmers;

	char     *query;
	uint32_t qlen, qdir;
	u8list *qvals;

	String *q[2], *t[2], *c[2], *e;
} DBGAligner;

typedef struct {
	int off, qlen, tlen;
	int score, closed;
	String *q, *t, *c;
} PBHit;
define_list(pbhitv, PBHit);

typedef struct {
	DBGAligner *aln;
	pbhitv *hits;
	String *seqs, *q, *t, *c;
	BitVec *covs;
} PBCorr;

#endif

