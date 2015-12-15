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

#ifndef __WATCHTOWER_SEQ_DB_RJ_H
#define __WATCHTOWER_SEQ_DB_RJ_H

#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "file_reader.h"
#include "bitvec.h"

#define wt_roundup8x(n) (((n) + 0x7LLU) & 0xFFFFFFFFFFFFFFF8LLU)

#define WT_MAX_RD	0x7FFFFFFFU
#define WT_MAX_RD_len	0xFFFFFU

typedef struct {
	uint32_t n_rd;
	BaseBank *rdseqs;
	u32list  *rdords;
	u32list  *rdlens;
	u64list  *rdoffs;
	cplist   *rdnames;
	cuhash   *rdname2id;
	BitVec   *masked;
} WTSEQDB;

static inline WTSEQDB* init_wtseqdb(){
	WTSEQDB* wt;
	wt = malloc(sizeof(WTSEQDB));
	wt->n_rd = 0;
	wt->rdseqs = init_basebank();
	wt->rdords = init_u32list(1024);
	wt->rdlens = init_u32list(1024);
	wt->rdoffs = init_u64list(1024);
	wt->rdnames = init_cplist(1024);
	wt->rdname2id = init_cuhash(1023);
	wt->masked = init_bitvec(1024);
	return wt;
}

static inline void free_wtseqdb(WTSEQDB *wt){
	uint32_t i;
	free_basebank(wt->rdseqs);
	free_u32list(wt->rdords);
	free_u32list(wt->rdlens);
	free_u64list(wt->rdoffs);
	for(i=0;i<wt->rdnames->size;i++) free((char*)get_cplist(wt->rdnames, i));
	free_cplist(wt->rdnames);
	free_cuhash(wt->rdname2id);
	free_bitvec(wt->masked);
	free(wt);
}

static inline void push_wtseqdb(WTSEQDB *wt, char *name, int name_len, char *seq, int seq_len){
	char *ptr;
	push_u32list(wt->rdlens, seq_len);
	push_u64list(wt->rdoffs, wt->rdseqs->size);
	seq2basebank(wt->rdseqs, seq, seq_len);
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_cplist(wt->rdnames, ptr);
	kv_put_cuhash(wt->rdname2id, ptr, wt->n_rd);
	zero2bitvec(wt->masked);
	wt->n_rd ++;
}

static inline void set_read_clip_wtseqdb(WTSEQDB *wt, char *name, int coff, int clen){
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->rdname2id, name)) == 0xFFFFFFFFU) return;
	if(coff < 0 || coff + clen > (int)wt->rdlens->buffer[pbid]) return;
	wt->rdoffs->buffer[pbid] += coff;
	wt->rdlens->buffer[pbid]  = clen;
}

static inline void ready_wtseqdb(WTSEQDB *wt){
	uint32_t i;
	clear_u32list(wt->rdords);
	for(i=0;i<wt->n_rd;i++) push_u32list(wt->rdords, i);
	sort_array(wt->rdords->buffer, wt->rdords->size, uint32_t, wt->rdlens->buffer[b] > wt->rdlens->buffer[a]);
}

typedef struct {
	uint32_t pb1:31, dir1:1, pb2:31, dir2:1;
	int qb, qe, tb, te;
	int score, mat, mis, ins, del, aln;
	char *cigar;
} wt_ovl_t;
define_list(wtovlv, wt_ovl_t);

static inline size_t print_wtovl_hits(WTSEQDB *db1, WTSEQDB *db2, wtovlv *hits, FILE *out){
	wt_ovl_t *hit;
	size_t i, ret;
	ret = 0;
	{
		for(i=0;i<hits->size;i++){
			hit = ref_wtovlv(hits, i);
			if(hit->score == -19830203){ if(hit->cigar){ free(hit->cigar); hit->cigar = NULL; } continue; }
			ret ++;
			if(hit->aln == 0) hit->aln = 1;
			if(out){
				fprintf(out, "%s\t%c\t%d\t%d\t%d", get_cplist(db1->rdnames, hit->pb1), "+-"[hit->dir1], db1->rdlens->buffer[hit->pb1], hit->tb, hit->te);
				fprintf(out, "\t%s\t%c\t%d\t%d\t%d", get_cplist(db2->rdnames, hit->pb2), "+-"[hit->dir2], db2->rdlens->buffer[hit->pb2], hit->qb, hit->qe);
				fprintf(out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d", hit->score, 1.0 * hit->mat / hit->aln, hit->mat, hit->mis, hit->ins, hit->del);
			}
			if(hit->cigar){
				if(out) fprintf(out, "\t%s\n", hit->cigar);
				free(hit->cigar);
				hit->cigar = NULL;
			} else {
				if(out) fprintf(out, "\t0M\n");
			}
		}
	}
	return ret;
}

typedef struct {
	int W, M, X, O, E, T;
	int min_local_block;
	int min_score;
	float min_identity;
} WTALNPARAM;

static const WTALNPARAM default_alnparam = (WTALNPARAM){800, 2, -5, -2, -1, -100, 100, 200, 0.6};

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

#endif
