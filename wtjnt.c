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

#include "string.h"
#include "list.h"
#include "sort.h"
#include "hashset.h"
#include "dna.h"
#include "hzm_aln.h"
#include "dagcns.h"
#include "thread.h"
#include "file_reader.h"

#define WT_MAX_RD_CNT	0x0FFFFFFFU
#define WT_MAX_RD_LEN	0x3FFFF
#define WT_MAX_OL_CNT	0xFFFF

static int verbose = 0;

typedef struct {
	uint64_t seqoff:46, seqlen:18;
	uint64_t hitoff:46, hitcnt:16;
	uint64_t pbi:28, pbx:18, pby:18;
	char *tag;
} wt_rd_t;
define_list(wtrdv, wt_rd_t);

typedef struct {
	uint32_t rid:31, dir:1;
	int beg, end;
} wt_reg_t;

typedef struct {
	wt_reg_t pair[2];
} wt_ol_t;
define_list(wtolv, wt_ol_t);

typedef struct {
	uint64_t seqoff:46, seqlen:18;
	uint32_t rdoff, rdcnt;
	char *tag;
} wt_pb_t;
define_list(wtpbv, wt_pb_t);

typedef struct {
	BaseBank *rdseqs;
	cuhash   *rd2ids;
	BaseBank *pbseqs;
	cuhash   *pb2ids;
	wtrdv    *rds;
	wtpbv    *pbs;
	wtolv    *ols;
} WTJNT;

WTJNT* init_wtjnt(){
	WTJNT *wt;
	wt = calloc(1, sizeof(WTJNT));
	wt->rdseqs = init_basebank();
	wt->rd2ids = init_cuhash(1023);
	wt->pbseqs = init_basebank();
	wt->pb2ids = init_cuhash(1023);
	wt->rds = init_wtrdv(1024);
	wt->pbs = init_wtpbv(1024);
	wt->ols = init_wtolv(1024);
	return wt;
}

void free_wtjnt(WTJNT *wt){
	uint32_t i;
	free_basebank(wt->rdseqs);
	free_cuhash(wt->rd2ids);
	free_basebank(wt->pbseqs);
	free_cuhash(wt->pb2ids);
	for(i=0;i<wt->rds->size;i++) free(wt->rds->buffer[i].tag);
	free_wtrdv(wt->rds);
	for(i=0;i<wt->pbs->size;i++) free(wt->pbs->buffer[i].tag);
	free_wtpbv(wt->pbs);
	free_wtolv(wt->ols);
	free(wt);
}

void push_pb_wtjnt(WTJNT *wt, char *tag, int taglen, char *seq, int seqlen){
	wt_pb_t *pb;
	pb = next_ref_wtpbv(wt->pbs);
	pb->seqlen = seqlen;
	pb->seqoff = wt->pbseqs->size;
	seq2basebank(wt->pbseqs, seq, seqlen);
	pb->tag = malloc(taglen + 1);
	memcpy(pb->tag, tag, taglen);
	pb->tag[taglen] = 0;
	kv_put_cuhash(wt->pb2ids, pb->tag, wt->pbs->size - 1);
	pb->rdoff = 0;
	pb->rdcnt = 0;
}

void finish_pb_wtjnt(WTJNT *wt){
	if(wt) return;
}

void push_rd_wtjnt(WTJNT *wt, char *tag, int taglen, char *seq, int seqlen){
	wt_rd_t *rd;
	uint32_t idx;
	int i, f, x, y;
	rd = next_ref_wtrdv(wt->rds);
	rd->seqlen = seqlen;
	rd->seqoff = wt->rdseqs->size;
	seq2basebank(wt->rdseqs, seq, seqlen);
	rd->tag = malloc(taglen + 1);
	memcpy(rd->tag, tag, taglen);
	rd->tag[taglen] = 0;
	rd->hitoff = 0;
	rd->hitcnt = 0;
	i = taglen - 1;
	x = y = 0;
	while(1){
		f = 1;
		while(i > 0 && tag[i] >= '0' && tag[i] <= '9'){ y += (tag[i] - '0') * f; f *= 10; i --; }
		if(i <= 0 || f == 1 || tag[i] != '_'){ x = y = 0; break; }
		i --;
		f = 1;
		while(i > 0 && tag[i] >= '0' && tag[i] <= '9'){ x += (tag[i] - '0') * f; f *= 10; i --; }
		if(i <= 0 || f == 1 || tag[i] != '/'){ x = y = 0; break; }
		taglen = i;
		tag[i] = '\0';
		break;
	}
	if((idx = kv_get_cuhash(wt->pb2ids, tag)) == 0xFFFFFFFU){
		fprintf(stderr, " -- Cannot find %s in pacbio sequences, in %s -- %s:%d --\n", tag, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	rd->pbi = idx;
	rd->pbx = x;
	rd->pby = y;
	wt->pbs->buffer[idx].rdcnt ++;
}

void finish_rd_wtjnt(WTJNT *wt){
	uint32_t i, lst;
	// Build hashmap for rdtag=>idx
	sort_array(wt->rds->buffer, wt->rds->size, wt_rd_t, num_cmpgtx(a.pbi, b.pbi, a.pbx, b.pbx));
	for(i=0;i<wt->rds->size;i++) kv_put_cuhash(wt->rd2ids, wt->rds->buffer[i].tag, i);
	lst = 0xFFFFFFFFU;
	for(i=0;i<wt->rds->size;i++){
		if(lst != wt->rds->buffer[i].pbi){
			lst = wt->rds->buffer[i].pbi;
			wt->pbs->buffer[lst].rdoff = i;
		}
	}
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

void load_overlaps_wtjnt(WTJNT *wt, FileReader *fr, float min_sm, int min_mat, int max_margin, int ncpu){
	wt_ol_t *ol1, *ol2;
	wt_rd_t *rd;
	uint64_t i, lst;
	uint32_t rids[2];
	float sm;
	int n, mat, tmp;
	if(ncpu < 1) ncpu = 1;
	while((n = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(n < 16) continue;
		if((rids[0] = kv_get_cuhash(wt->rd2ids, get_col_str(fr, 0))) == 0xFFFFFFFFU){
			fprintf(stderr, " -- Cannot find %s INPUT:%llu line, in %s -- %s:%d --\n", get_col_str(fr, 0), fr->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); continue;
		}
		if((rids[1] = kv_get_cuhash(wt->rd2ids, get_col_str(fr, 5))) == 0xFFFFFFFFU){
			fprintf(stderr, " -- Cannot find %s INPUT:%llu line, in %s -- %s:%d --\n", get_col_str(fr, 5), fr->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); continue;
		}
		if((sm = atof(get_col_str(fr, 11))) < min_sm) continue;
		if((mat = atoi(get_col_str(fr, 12))) < min_mat) continue;
		encap_wtolv(wt->ols, 2);
		ol1 = next_ref_wtolv(wt->ols);
		ol1->pair[0].rid = rids[0];
		ol1->pair[0].dir = (get_col_str(fr, 1)[0] == '-');
		ol1->pair[0].beg = atoi(get_col_str(fr, 3));
		ol1->pair[0].end = atoi(get_col_str(fr, 4));
		ol1->pair[1].rid = rids[1];
		ol1->pair[1].dir = (get_col_str(fr, 6)[0] == '-');
		ol1->pair[1].beg = atoi(get_col_str(fr, 8));
		ol1->pair[1].end = atoi(get_col_str(fr, 9));
		if(num_min(ol1->pair[0].beg, ol1->pair[1].beg) > max_margin || num_min(wt->rds->buffer[rids[0]].seqlen - ol1->pair[0].end, wt->rds->buffer[rids[1]].seqlen - ol1->pair[1].end) > max_margin){
			wt->ols->size --; continue;
		}
		if(ol1->pair[0].dir){
			swap_tmp(ol1->pair[0].beg, ol1->pair[0].end, tmp);
			ol1->pair[0].beg = wt->rds->buffer[rids[0]].seqlen - ol1->pair[0].beg;
			ol1->pair[0].end = wt->rds->buffer[rids[0]].seqlen - ol1->pair[0].end;
		}
		if(ol1->pair[1].dir){
			swap_tmp(ol1->pair[1].beg, ol1->pair[1].end, tmp);
			ol1->pair[1].beg = wt->rds->buffer[rids[1]].seqlen - ol1->pair[1].beg;
			ol1->pair[1].end = wt->rds->buffer[rids[1]].seqlen - ol1->pair[1].end;
		}
		ol2 = next_ref_wtolv(wt->ols);
		ol2->pair[0] = ol1->pair[1];
		ol2->pair[1] = ol1->pair[0];
	}
	if(wt->ols->size == 0) return;
	psort_array(wt->ols->buffer, wt->ols->size, wt_ol_t, ncpu, num_cmpgtx(a.pair[0].rid, b.pair[0].rid, a.pair[1].rid, b.pair[1].rid));
	encap_wtolv(wt->ols, 1);
	ol1 = ref_wtolv(wt->ols, wt->ols->size);
	ol1->pair[0].rid = WT_MAX_RD_CNT;
	ol1 = ref_wtolv(wt->ols, (lst = 0));
	for(i=1;i<=wt->ols->size;i++){
		ol2 = ref_wtolv(wt->ols, i);
		if(ol2->pair[0].rid == ol1->pair[0].rid) continue;
		rd = ref_wtrdv(wt->rds, ol1->pair[0].rid);
		rd->hitoff = lst;
		rd->hitcnt = i - lst;
		lst = i;
		ol1 = ol2;
	}
}

typedef struct {
	uint32_t rdid; // rdid --> rdid + 1
	int type, idx, dir, x, y;
} syn_lnk_t;
define_list(synlnkv, syn_lnk_t);

typedef struct {
	u1v *types; // 0: retained; 1: overlaped; 2: accurate filled; 3: cns filled, please reminder reg_inc;
	u1v *gapseqs;
	u4v *seqoffs;
	b4v *seqlens;
} SynGap;

static const int reg_inc = 100;

SynGap* init_syngap(){
	SynGap *sg;
	sg = calloc(1, sizeof(SynGap));
	sg->types = init_u1v(8);
	sg->gapseqs = init_u1v(1024);
	sg->seqoffs = init_u4v(8);
	sg->seqlens = init_b4v(8);
	return sg;
}

void free_syngap(SynGap *sg){
	free_u1v(sg->types);
	free_u1v(sg->gapseqs);
	free_u4v(sg->seqoffs);
	free_b4v(sg->seqlens);
	free(sg);
}

int make_cns_fillgap(WTJNT *wt, uint32_t pbid, uint32_t ridx, syn_lnk_t *lnks, uint32_t lnksize, SynGap *sg, HZMAux *aux, DAGCNS *dag, u1v *rdseq, float min_sm){
	wt_pb_t *pb, *pbb;
	wt_rd_t *r1, *r2;
	syn_lnk_t *lnk;
	uint32_t i, r1b, r2e;
	int j, x, y;
	pb = ref_wtpbv(wt->pbs, pbid);
	r1 = ref_wtrdv(wt->rds, pb->rdoff + ridx);
	r2 = ref_wtrdv(wt->rds, pb->rdoff + ridx + 1);
	x = r1->pby; y = r2->pbx;
	// init backbone
	reset_hzmaux(aux);
	r1b = num_max((uint32_t)reg_inc, r1->seqlen) - reg_inc;
	for(i=r1b;i<r1->seqlen;i++) add_tseq_hzmaux(aux, bits2bit(wt->rdseqs->bits, r1->seqoff + i));
	for(j=x;j<y;j++) add_tseq_hzmaux(aux, bits2bit(wt->pbseqs->bits, pb->seqoff + j));
	r2e = num_min((uint32_t)reg_inc, r2->seqlen);
	for(i=0;i<r2e;i++) add_tseq_hzmaux(aux, bits2bit(wt->rdseqs->bits, r2->seqoff + i));
	ready_hzmaux(aux);
	reset_dagcns(dag);
	append_u8list(dag->cns, aux->tseq);
	gen_pregraph_dagcns(dag);
	// Alignment
	for(i=0;i<lnksize;i++){
		lnk = lnks + i;
		if(lnk->type == 0) continue;
		pbb = ref_wtpbv(wt->pbs, lnk->idx);
		x = num_max(0, lnk->x - reg_inc);
		y = num_min((int)pbb->seqlen, lnk->y + reg_inc);
		if(x > y) continue;
		clear_and_encap_u1v(rdseq, y - x);
		if(lnk->dir) revbitseq_basebank(wt->pbseqs, pbb->seqoff + x, y - x, rdseq->buffer);
		else            bitseq_basebank(wt->pbseqs, pbb->seqoff + x, y - x, rdseq->buffer);
		rdseq->size = y - x;
		if(align_hzmaux(aux, 0, rdseq->buffer, NULL, rdseq->size, -1, -1, 1, min_sm) == 0) continue;
		if(verbose){
			fprintf(stdout, "%d\t%d\n%s\n%s\n%s\n", aux->hit.tb, aux->hit.te, aux->alns[0]->string, aux->alns[2]->string, aux->alns[1]->string); fflush(stdout);
		}
		alignment2dagcns(dag, aux->hit.tb, aux->hit.te, (char*[2]){aux->alns[0]->string, aux->alns[1]->string}, aux->alns[0]->size);
	}
	merge_nodes_dagcns(dag);
	gen_consensus_dagcns(dag, NULL);
	sg->types->buffer[ridx] = 3;
	sg->seqoffs->buffer[ridx] = sg->gapseqs->size;
	sg->seqlens->buffer[ridx] = dag->cns->size;
	append_array_u1v(sg->gapseqs, dag->cns->buffer, dag->cns->size);
	return 1;
}

// First, judge the synteny of reads within group
// Then, fill in the gaps
int join_read_group_wtjnt(WTJNT *wt, uint32_t pbid, wtolv *hits, synlnkv *lnks, SynGap *sg, HZMAux *aux, DAGCNS *cns, u1v *rdseq, uint32_t min_sup, uint32_t max_sup, float min_sm){
	wt_pb_t *pb;
	wt_rd_t *rd, *rds[4];
	wt_ol_t *ol, *ols[2];
	syn_lnk_t *lnk;
	uint32_t i, j, k, beg, end, idx;
	int dir;
	pb = ref_wtpbv(wt->pbs, pbid);
	clear_and_encap_u1v(sg->types, pb->rdcnt);
	clear_u1v(sg->gapseqs);
	clear_and_encap_u4v(sg->seqoffs, pb->rdcnt);
	clear_and_encap_b4v(sg->seqlens, pb->rdcnt);
	for(i=0;i<pb->rdcnt;i++){
		push_u1v(sg->types, 0);
		push_u4v(sg->seqoffs, 0);
		push_b4v(sg->seqlens, 0);
	}
	if(pb->rdcnt < 2) return 0;
	// collecting aligned reads
	clear_wtolv(hits);
	for(i=0;i<pb->rdcnt;i++){
		rd = ref_wtrdv(wt->rds, pb->rdoff + i);
		for(j=0;j<rd->hitcnt;j++){
			push_wtolv(hits, wt->ols->buffer[rd->hitoff + j]);
		}
	}
	if(hits->size == 0) return 0;
	// grouping matched reads
	sort_array(hits->buffer, hits->size, wt_ol_t, num_cmpgtx(wt->rds->buffer[a.pair[1].rid].pbi, wt->rds->buffer[b.pair[1].rid].pbi, a.pair[0].rid, b.pair[0].rid));
	// counting synteny
	clear_synlnkv(lnks);
	if(verbose) fprintf(stdout, ">%s\tn_local=%d\n", pb->tag, pb->rdcnt);
	for(end=beg=0;end<=hits->size;end++){
		if(end < hits->size && wt->rds->buffer[hits->buffer[end].pair[1].rid].pbi == wt->rds->buffer[hits->buffer[beg].pair[1].rid].pbi) continue;
		if(end - beg > 1){
			for(j=beg;j<end;j++){
				ol = ref_wtolv(hits, j);
				ol->pair[1].dir ^= ol->pair[0].dir;
				ol->pair[0].dir = 0;
				rds[0] = ref_wtrdv(wt->rds, ol->pair[0].rid);
				rds[1] = ref_wtrdv(wt->rds, ol->pair[1].rid);
				if(verbose) fprintf(stdout, "%s\t%d\t%c\t%d\t%d\t%s\t%d\t%c\t%d\t%d\n",
							rds[1]->tag, rds[1]->seqlen, "+-"[ol->pair[1].dir], ol->pair[1].beg, ol->pair[1].end,
							rds[0]->tag, rds[0]->seqlen, "+-"[ol->pair[0].dir], ol->pair[0].beg, ol->pair[0].end
							);
			}
			if(verbose) fprintf(stdout, "\n");
			for(j=beg;j<end;j++){
				ols[0] = ref_wtolv(hits, j);
				rds[0]  = ref_wtrdv(wt->rds, ols[0]->pair[0].rid);
				rds[1]  = ref_wtrdv(wt->rds, ols[0]->pair[1].rid);
				dir = ols[0]->pair[1].dir;
				for(k=j+1;k<end;k++){
					ols[1] = ref_wtolv(hits, k);
					rds[2] = ref_wtrdv(wt->rds, ols[1]->pair[0].rid);
					rds[3] = ref_wtrdv(wt->rds, ols[1]->pair[1].rid);
					if(rds[0] + 1 != rds[2]) continue; // only adjacent reads
					if(ols[0]->pair[1].dir ^ ols[1]->pair[1].dir) continue; // read strand checking
					if(rds[1] != rds[3] && dir ^ (rds[1] > rds[3]? 1 : 0)) continue; // read order checking
					lnk = next_ref_synlnkv(lnks);
					lnk->rdid = ols[0]->pair[0].rid;
					lnk->dir = dir;
					if(rds[1] == rds[3]){
						lnk->type = 0; // direct
						lnk->idx  = ols[0]->pair[1].rid;
						lnk->x = ols[ dir]->pair[1].end;
						lnk->y = ols[!dir]->pair[1].beg;
					} else {
						lnk->type = 1; // indirect
						lnk->idx  = wt->rds->buffer[ols[0]->pair[1].rid].pbi;
						lnk->x = rds[  dir  * 2 + 1]->pbx + ols[ dir]->pair[1].end * (rds[  dir  * 2 + 1]->pby - rds[  dir  * 2 + 1]->pbx) / rds[  dir  * 2 + 1]->seqlen;
						lnk->y = rds[(!dir) * 2 + 1]->pbx + ols[!dir]->pair[1].beg * (rds[(!dir) * 2 + 1]->pby - rds[(!dir) * 2 + 1]->pbx) / rds[(!dir) * 2 + 1]->seqlen;
					}
					if(verbose) fprintf(stdout, "%s\t%c\t%d\t%d\t%d\n", lnk->type? wt->pbs->buffer[lnk->idx].tag : wt->rds->buffer[lnk->idx].tag, "+-"[lnk->dir], lnk->y - lnk->x, lnk->x, lnk->y);
				}
			}
			if(verbose){ fprintf(stdout, "//\n"); fflush(stdout); }
		}
		beg = end;
	}
	if(lnks->size == 0) return 0;
	sort_array(lnks->buffer, lnks->size, syn_lnk_t, num_cmpgtx(a.rdid, b.rdid, a.type, b.type));
	for(end=beg=0;end<=lnks->size;end++){
		if(end < lnks->size && lnks->buffer[end].rdid == lnks->buffer[beg].rdid) continue;
		if(verbose){
			for(i=beg;i<end;i++){
				lnk = ref_synlnkv(lnks, i);
				fprintf(stdout, "%s\t%c\t%d\t%d\t%d\n", lnk->type? wt->pbs->buffer[lnk->idx].tag : wt->rds->buffer[lnk->idx].tag, "+-"[lnk->dir], lnk->y - lnk->x, lnk->x, lnk->y);
			}
			fprintf(stdout, "\n");
			fflush(stdout);
		}
		if(end - beg >= min_sup){
			rds[0] = ref_wtrdv(wt->rds, lnks->buffer[beg].rdid);
			rds[1] = ref_wtrdv(wt->rds, lnks->buffer[beg].rdid + 1);
			idx = lnks->buffer[beg].rdid - pb->rdoff;
			if(rds[0]->pby >= rds[1]->pbx){ // r1 and r2 have overlaped bases
				sg->types->buffer[idx] = 1;
				sg->seqlens->buffer[idx] = ((int)rds[1]->pbx) - ((int)rds[0]->pby);
			} else if((lnk = ref_synlnkv(lnks, beg))->type == 0){ // r1 and r2 can be directly connected by another corrected read
				sg->types->buffer[idx] = 2;
				sg->seqoffs->buffer[idx] = sg->gapseqs->size;
				for(i=lnk->x;(int)i<lnk->y;i++) push_u1v(sg->gapseqs, bits2bit(wt->rdseqs->bits, wt->rds->buffer[lnk->idx].seqoff + i));
				sg->seqlens->buffer[idx] = lnk->y - lnk->x;
				if(lnk->dir && lnk->y > lnk->x){
					for(i=0;(int)i<sg->seqlens->buffer[idx];i++) sg->gapseqs->buffer[sg->seqoffs->buffer[idx] + i] = (~sg->gapseqs->buffer[sg->seqoffs->buffer[idx] + i]) & 0x03;
					sub_reverse_u1v(sg->gapseqs, sg->seqoffs->buffer[idx], sg->seqoffs->buffer[idx] + sg->seqlens->buffer[idx]);
				}
			} else {
				make_cns_fillgap(wt, pbid, idx, lnks->buffer + beg, num_min(end - beg, max_sup), sg, aux, cns, rdseq, min_sm);
			}
		}
		beg = end;
	}
	if(verbose){
		for(i=0;i+1<pb->rdcnt;i++){
			rd = ref_wtrdv(wt->rds, pb->rdoff + i);
			fprintf(stdout, "@%s\t%s\t%d\t", rd->tag,((char*[]){"OPENED", "OVERLAPED", "LINKED", "FILLED"})[sg->types->buffer[i]], sg->seqlens->buffer[i]);
			for(j=0;(int)j<sg->seqlens->buffer[i];j++) fputc(bit_base_table[sg->gapseqs->buffer[sg->seqoffs->buffer[i] + j]], stdout);
			fprintf(stdout, "\n");
		}
		fflush(stdout);
	}
	return 1;
}

thread_beg_def(mjnt);
WTJNT *wt;
uint32_t pbid;
float min_cns_sm, min_pb_cov;
FILE *out, *err;
thread_end_def(mjnt);

thread_beg_func(mjnt);
WTJNT *wt;
wtolv *hits;
synlnkv *lnks;
SynGap *sg;
HZMAux *aux;
DAGCNS *dag;
u1v *rdseq;
wt_pb_t *pb;
wt_rd_t *rd;
FILE *out, *err;
uint32_t j, k, l, her, lst;
int cov, x, y, delta, len;
float min_cns_sm, min_pb_cov, crate;
hits = init_wtolv(1024);
lnks = init_synlnkv(1024);
sg = init_syngap();
aux = init_hzmaux();
dag = init_dagcns();
rdseq = init_u1v(1024);
wt = mjnt->wt;
out = mjnt->out;
err = mjnt->err;
min_cns_sm = mjnt->min_cns_sm;
min_pb_cov = mjnt->min_pb_cov;
thread_beg_loop(mjnt);
pb = ref_wtpbv(wt->pbs, mjnt->pbid);
join_read_group_wtjnt(wt, mjnt->pbid, hits, lnks, sg, aux, dag, rdseq, 2, 5, min_cns_sm);
thread_beg_syn(mjnt);
her = 0;
if(pb->rdcnt == 0){
	her = 1;
	crate = 0;
} else {
	cov = 0;
	lst = 0;
	for(j=0;j<pb->rdcnt;j++){
		rd = ref_wtrdv(wt->rds, pb->rdoff + j);
		if(sg->types->buffer[j] == 0){
			clear_u1v(rdseq);
			delta = 0;
			for(k=lst;k<j;k++){
				if(sg->types->buffer[k] == 1 || sg->seqlens->buffer[k] <= 0){
					len = wt->rds->buffer[k + pb->rdoff].seqlen + sg->seqlens->buffer[k] - delta;
					if(len > 0){
						encap_and_inc_u1v(rdseq, len);
						bitseq_basebank(wt->rdseqs, wt->rds->buffer[k + pb->rdoff].seqoff + delta, len, rdseq->buffer + rdseq->size - len);
					}
					delta = 0;
				} else if(sg->types->buffer[k] == 2){
					len = wt->rds->buffer[k + pb->rdoff].seqlen - delta;
					if(len > 0){
						encap_and_inc_u1v(rdseq, len);
						bitseq_basebank(wt->rdseqs, wt->rds->buffer[k + pb->rdoff].seqoff + delta, len, rdseq->buffer + rdseq->size - len);
					}
					append_array_u1v(rdseq, sg->gapseqs->buffer + sg->seqoffs->buffer[k], sg->seqlens->buffer[k]);
					delta = 0;
				} else if(sg->types->buffer[k] == 3){
					len = wt->rds->buffer[k + pb->rdoff].seqlen - delta - reg_inc;
					if(len > 0){
						encap_and_inc_u1v(rdseq, len);
						bitseq_basebank(wt->rdseqs, wt->rds->buffer[k + pb->rdoff].seqoff + delta, len, rdseq->buffer + rdseq->size - len);
					}
					append_array_u1v(rdseq, sg->gapseqs->buffer + sg->seqoffs->buffer[k], sg->seqlens->buffer[k]);
					delta = reg_inc;
				}
			}
			len = wt->rds->buffer[k + pb->rdoff].seqlen - delta;
			if(len > 0){
				encap_and_inc_u1v(rdseq, len);
				bitseq_basebank(wt->rdseqs, wt->rds->buffer[k + pb->rdoff].seqoff + delta, len, rdseq->buffer + rdseq->size - len);
			}
			x = wt->rds->buffer[lst + pb->rdoff].pbx;
			y = wt->rds->buffer[j   + pb->rdoff].pby;
			fprintf(out, ">%s/%d_%d corr=Y len=%d join=%d\n", pb->tag, x, y, (int)rdseq->size, j + 1 - lst);
			for(l=0;l<rdseq->size;l++) fputc(bit_base_table[(int)rdseq->buffer[l]], out);
			fprintf(out, "\n");
			cov += y - x;
			lst = j + 1;
		}
	}
	if((crate = ((double)cov) / pb->seqlen) < min_pb_cov) her = 1;
}
if(her && err){
	fprintf(err, ">%s corr=N len=%d cov=%0.3f\n", pb->tag, pb->seqlen, crate);
	print_seq_basebank(wt->pbseqs, pb->seqoff, pb->seqlen, err);
	fprintf(err, "\n");
}
thread_end_syn(mjnt);
thread_end_loop(mjnt);
free_wtolv(hits);
free_synlnkv(lnks);
free_syngap(sg);
free_hzmaux(aux);
free_dagcns(dag);
free_u1v(rdseq);
thread_end_func(mjnt);

int usage(){
	printf(
	"Program: wtjnt\n"
	" Join local-alignment corrections of one read from wtcorr\n"
	"Author: Jue Ruan\n"
	"Version: 1.0 20150824\n"
	"Usage: wtjnt [options]\n"
	"Options:\n"
	" -h          Show this document\n"
	" -t <int>    Number of threads, [1]\n"
	" -1 <string> Original reads file, +, *\n"
	" -2 <string> Corrected reads file, +, *\n"
	" -3 <string> Alignment of corrected reads, usually from wtzmo, +, *\n"
	" -o <string> Output of joint reads, [stdout]\n"
	" -x <string> Output uncorrected high-error-rate (HER) seuences, [NULL]\n"
	"             HER might come from complex genomic regions, e.g. high GC\n"
	"             It is necessary to bring them back to assembly\n"
	" -m <int>    Min matches in overlap, [1000]\n"
	" -s <float>  Min similarity of overlap, [0.995]\n"
	//" -M <int>    Minimum of other reads which are syntenic, [2]\n"
	" -S <float>  Min similarity of alignment in consensus calling, [0.60]\n"
	" -X <float>  Min coverage of corrected/total for one original read to be not HER, [0.6]\n"
	"             Otherwise, will output whole original read to <-x> as HER\n"
	" -v          Verbose\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	WTJNT *wt;
	FileReader *fr;
	Sequence *seq;
	cplist *pbfs, *rdfs, *olfs;
	char *outf, *errf;
	FILE *out, *err;
	int c, ncpu, min_mat, overwrite;
	float min_sim, min_cns_sm, min_pb_cov;
	uint32_t i;
	thread_preprocess(mjnt);
	outf = NULL;
	errf = NULL;
	pbfs = init_cplist(4);
	rdfs = init_cplist(4);
	olfs = init_cplist(4);
	min_mat = 1000;
	min_sim = 0.995;
	min_cns_sm = 0.60;
	min_pb_cov = 0.60;
	ncpu = 1;
	overwrite = 0;
	while((c = getopt(argc, argv, "hvft:1:2:3:o:x:m:s:M:S:X:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'v': verbose ++; break;
			case 'f': overwrite = 1; break;
			case 't': ncpu = atoi(optarg); break;
			case '1': push_cplist(pbfs, optarg); break;
			case '2': push_cplist(rdfs, optarg); break;
			case '3': push_cplist(olfs, optarg); break;
			case 'o': outf = optarg; break;
			case 'x': errf = optarg; break;
			case 'm': min_mat = atoi(optarg); break;
			case 's': min_sim = atof(optarg); break;
			case 'S': min_cns_sm = atof(optarg); break;
			case 'X': min_pb_cov = atof(optarg); break;
			default: return usage();
		}
	}
	if(!overwrite && outf && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, " --  File '%s' already exists in %s -- %s:%d --\n", outf, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(!overwrite && errf && strcmp(errf, "-") && file_exists(errf)){
		fprintf(stderr, " --  File '%s' already exists in %s -- %s:%d --\n", errf, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	if(pbfs->size == 0){
		fprintf(stderr, " --  There is none original reads file  in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		return usage();
	}
	if(rdfs->size == 0){
		fprintf(stderr, " --  There is none corrected reads file  in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		return usage();
	}
	if(olfs->size == 0){
		fprintf(stderr, " --  There is none alignments file for corrected reads file  in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		return usage();
	}
	wt = init_wtjnt();
	fr = fopen_m_filereader(pbfs->size, pbfs->buffer);
	seq = NULL;
	while(fread_seq(&seq, fr)) push_pb_wtjnt(wt, seq->tag.string, seq->tag.size, seq->seq.string, seq->seq.size);
	fclose_filereader(fr);
	finish_pb_wtjnt(wt);
	fr = fopen_m_filereader(rdfs->size, rdfs->buffer);
	seq = NULL;
	while(fread_seq(&seq, fr)) if(seq->seq.size >= min_mat)  push_rd_wtjnt(wt, seq->tag.string, seq->tag.size, seq->seq.string, seq->seq.size);
	fclose_filereader(fr);
	finish_rd_wtjnt(wt);
	fr = fopen_m_filereader(olfs->size, olfs->buffer);
	load_overlaps_wtjnt(wt, fr, min_sim, min_mat, 10, ncpu);
	fclose_filereader(fr);
	free_cplist(pbfs);
	free_cplist(rdfs);
	free_cplist(olfs);
	out = (outf && strcmp(outf, "-"))? fopen(outf, "w") : stdout;
	err = errf? (strcmp(errf, outf)? (strcmp(errf, "-")? fopen(errf, "w") : stdout) : out) : NULL;
	thread_beg_init(mjnt, ncpu);
	mjnt->wt = wt;
	mjnt->pbid = 0;
	mjnt->out = out;
	mjnt->err = err;
	mjnt->min_cns_sm = min_cns_sm;
	mjnt->min_pb_cov = min_pb_cov;
	thread_end_init(mjnt);
	for(i=0;i<wt->pbs->size;i++){
		thread_wait_one(mjnt);
		mjnt->pbid = i;
		thread_wake(mjnt);
	}
	thread_wait_all(mjnt);
	thread_beg_close(mjnt);
	thread_end_close(mjnt);
	if(err && err != out && err != stdout) fclose(err);
	if(out && out != stdout) fclose(out);
	free_wtjnt(wt);
	return 0;
}

