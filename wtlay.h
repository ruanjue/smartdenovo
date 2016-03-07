
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
 
#ifndef __PACBIO_ASSEMBLY_WATCH_TOWER_LAYOUT_RJ_H
#define __PACBIO_ASSEMBLY_WATCH_TOWER_LAYOUT_RJ_H

#include <stdint.h>
#include <math.h>
#include "string.h"
#include "list.h"
#include "sort.h"
#include "hashset.h"
#include "dna.h"
#include "file_reader.h"
#include "heap.h"

#define SG_EDGE_BITS	10
#define SG_MAX_EDGE	0x3FF
#define SG_NODE_BTID_MAX	0x7FFFF
#define SG_EDGE_MAX_COV	63

typedef struct {
	uint64_t edge_offs[2];
	uint32_t edge_cnts[2];
	struct { uint64_t bt_visit:19, bt_dir:1, bt_idx:10, bt_off:34; };
	struct { uint8_t bogs[2][2][2]; };
	struct { uint32_t lay_id; uint32_t dir:1, end:1, off:30; };
} sg_node_t;
define_list(sgnodev, sg_node_t);

typedef struct {
	uint32_t node_id[2];
	int off[2];
	int ol[2];
	int8_t dir[2];
	int32_t score;
} sg_biedge_t;
define_list(sgbiedgev, sg_biedge_t);

typedef struct {
	uint32_t node_id;
	uint32_t dir:1, closed:2, mark:1, att:1, tta:1, off:20, cov:6;
	int32_t score;
	int16_t ol_var;
	uint16_t rev_idx;
} sg_edge_t;
define_list(sgedgev, sg_edge_t);

typedef struct {
	uint32_t node_id;
	uint32_t dir:1, contained:1;
	int off:30;
	uint16_t eidx[2];
} sg_layout_t;
define_list(sglayv, sg_layout_t);

typedef struct {
	uint32_t n_rd;
	BaseBank *rdseqs;
	u8list   *rdqvs;
	u32list  *rdlens;
	u64list  *rdoffs;
	u32list  *rdqlens;
	u64list  *rdqoffs;
	cplist   *rdnames;
	cuhash   *rdname2id;

	int min_score;
	float min_id;
	int max_ovl_margin;
	int mat_score;
	sgnodev *nodes;
	BitVec  *node_status, *node_flags, *node_atts;
	sgedgev *edges;
	vplist  *lays;
} StringGraph;

#ifdef _CPLUSPLUS
extern "C" {
#endif

StringGraph* init_strgraph(){
	StringGraph *g;
	g = malloc(sizeof(StringGraph));
	g->n_rd = 0;
	g->rdseqs = init_basebank();
	g->rdqvs  = init_u8list(1024);
	g->rdlens = init_u32list(1024);
	g->rdoffs = init_u64list(1024);
	g->rdqlens = init_u32list(1024);
	g->rdqoffs = init_u64list(1024);
	g->rdnames = init_cplist(1024);
	g->rdname2id = init_cuhash(1023);
	g->min_score = 200;
	g->min_id = 0.6;
	g->max_ovl_margin = 500;
	g->mat_score = 0;
	g->nodes = init_sgnodev(1024);
	g->node_status = init_bitvec(1024);
	g->node_flags = init_bitvec(1024);
	g->node_atts = init_bitvec(1024);
	g->edges = init_sgedgev(1024);
	g->lays = init_vplist(1024);
	return g;
}

void free_strgraph(StringGraph *g){
	uint32_t i;
	free_basebank(g->rdseqs);
	free_u8list(g->rdqvs);
	free_u32list(g->rdlens);
	free_u64list(g->rdoffs);
	free_u32list(g->rdqlens);
	free_u64list(g->rdqoffs);
	for(i=0;i<g->rdnames->size;i++) free(get_cplist(g->rdnames, i));
	free_cplist(g->rdnames);
	free_cuhash(g->rdname2id);
	free_sgnodev(g->nodes);
	free_sgedgev(g->edges);
	free_bitvec(g->node_status);
	free_bitvec(g->node_flags);
	free_bitvec(g->node_atts);
	for(i=0;i<g->lays->size;i++) free_sglayv((sglayv*)get_vplist(g->lays, i));
	free_vplist(g->lays);
	free(g);
}

void push_read_strgraph(StringGraph *g, char *name, int name_len, char *seq, int seq_len){
	char *ptr;
	push_u32list(g->rdlens, seq_len);
	push_u64list(g->rdoffs, g->rdseqs->size);
	seq2basebank(g->rdseqs, seq, seq_len);
	push_u32list(g->rdqlens, 0);
	push_u64list(g->rdqoffs, g->rdqvs->size);
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_cplist(g->rdnames, ptr);
	kv_put_cuhash(g->rdname2id, ptr, g->n_rd);
	g->n_rd ++;
}

void push_read5q_strgraph(StringGraph *g, char *name, int name_len, char *seq, int seq_len, char *qvs){
	char *ptr;
	int i, j;
	push_u32list(g->rdlens, seq_len);
	push_u64list(g->rdoffs, g->rdseqs->size);
	seq2basebank(g->rdseqs, seq, seq_len);
	push_u32list(g->rdqlens, seq_len * 7);
	push_u64list(g->rdqoffs, g->rdqvs->size);
	for(i=0;i<seq_len;i++){
		for(j=0;j<7;j++){
			push_u8list(g->rdqvs, qvs[j * seq_len + i]);
		}
	}
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_cplist(g->rdnames, ptr);
	kv_put_cuhash(g->rdname2id, ptr, g->n_rd);
	g->n_rd ++;
}

void set_read_clip_strgraph(StringGraph *wt, char *name, int coff, int clen){
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->rdname2id, name)) == 0xFFFFFFFFU) return;
	if(coff < 0 || coff + clen > (int)wt->rdlens->buffer[pbid]) return;
	wt->rdoffs->buffer[pbid] += coff;
	wt->rdlens->buffer[pbid]  = clen;
	if(wt->rdqlens->buffer[pbid]){
		wt->rdqoffs->buffer[pbid] += coff * 7;
		wt->rdqlens->buffer[pbid]  = clen * 7;
	}
}

void generate_nodes_strgraph(StringGraph *g){
	clear_sgnodev(g->nodes);
	encap_sgnodev(g->nodes, g->n_rd);
	g->nodes->size = g->n_rd;
	zeros_sgnodev(g->nodes);
	encap_bitvec(g->node_status, g->n_rd);
	zeros_bitvec(g->node_status);
	encap_bitvec(g->node_flags, g->n_rd);
	zeros_bitvec(g->node_flags);
	encap_bitvec(g->node_atts, g->n_rd);
	zeros_bitvec(g->node_atts);
}

#define is_dead_node_strgraph(g, node_id) (get_bitvec((g)->node_status, node_id))

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
typedef struct {
	uint32_t node_id[2];
	int dir[2];
	int beg[2];
	int end[2];
	int score;
	int identity;
} OverlapData;

int parse_overlap_item_strgraph(StringGraph *g, FileReader *fr, OverlapData *dat){
	int n, pb1, pb2;
	while((n = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(n < 16) continue;
		dat->node_id[0] = kv_get_cuhash(g->rdname2id, get_col_str(fr, 0));
		if(dat->node_id[0] == 0xFFFFFFFFU) continue;
		if(g->rdlens->buffer[dat->node_id[0]] == 0) continue;
		dat->dir[0] = (get_col_str(fr, 1)[0] == '-');
		pb1 = atoi(get_col_str(fr, 2));
		if(pb1 != (int)g->rdlens->buffer[dat->node_id[0]]){
			fprintf(stderr, " -- Inconsistent read (%s) length %d != %d in %s -- %s:%d --\n", get_col_str(fr, 0), g->rdlens->buffer[dat->node_id[0]], pb1, __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		dat->beg[0] = atoi(get_col_str(fr, 3));
		dat->end[0] = atoi(get_col_str(fr, 4));
		dat->node_id[1] = kv_get_cuhash(g->rdname2id, get_col_str(fr, 5));
		if(dat->node_id[1] == 0xFFFFFFFFU) continue;
		if(dat->node_id[0] == dat->node_id[1]) continue;
		if(g->rdlens->buffer[dat->node_id[1]] == 0) continue;
		dat->dir[1] = (get_col_str(fr, 6)[0] == '-');
		pb2 = atoi(get_col_str(fr, 7));
		if(pb2 != (int)g->rdlens->buffer[dat->node_id[1]]){
			fprintf(stderr, " -- Inconsistent read (%s) length %d != %d in %s -- %s:%d --\n", get_col_str(fr, 5), g->rdlens->buffer[dat->node_id[1]], pb2, __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		dat->beg[1] = atoi(get_col_str(fr, 8));
		dat->end[1] = atoi(get_col_str(fr, 9));
		if(g->mat_score) dat->score = atoi(get_col_str(fr, 12));
		else dat->score = atoi(get_col_str(fr, 10));
		dat->identity = atof(get_col_str(fr, 11)) * 1000;
		if(dat->score < g->min_score) continue;
		if(dat->identity < 1000 * g->min_id) continue;
		return 1;
	}
	return 0;
}

int overlap_item2biedge_strgraph(StringGraph *g, OverlapData *dat, sg_biedge_t *e){
	int len1, len2, l[3], r[3];
	if(is_dead_node_strgraph(g, dat->node_id[0])) return 0;
	if(is_dead_node_strgraph(g, dat->node_id[1])) return 0;
	len1 = g->rdlens->buffer[dat->node_id[0]];
	len2 = g->rdlens->buffer[dat->node_id[1]];
	l[0] = dat->beg[0];
	l[1] = dat->beg[1];
	r[0] = len1 - dat->end[0];
	r[1] = len2 - dat->end[1];
	l[2] = l[0] < l[1]? l[0] : l[1];
	r[2] = r[0] < r[1]? r[0] : r[1];
	if(l[2] + r[2] > g->max_ovl_margin) return 0;
	if(l[0] >= l[1]){
		e->node_id[0] = dat->node_id[0];
		e->node_id[1] = dat->node_id[1];
		e->dir[0] = dat->dir[0];
		e->dir[1] = dat->dir[1];
		e->off[0] = l[0] - l[2];
		e->off[1] = r[1] - r[2];
		e->ol[0]  = dat->end[0] - dat->beg[0];
		e->ol[1]  = dat->end[1] - dat->beg[1];
	} else  {
		e->node_id[1] = dat->node_id[0];
		e->node_id[0] = dat->node_id[1];
		e->dir[1] = dat->dir[0];
		e->dir[0] = dat->dir[1];
		e->off[1] = r[0] - r[2];
		e->off[0] = l[1] - l[2];
		e->ol[1]  = dat->end[0] - dat->beg[0];
		e->ol[0]  = dat->end[1] - dat->beg[1];
	}
	e->score = dat->score;
	return 1;
}

int overlap_item2biedge_v2_strgraph(StringGraph *g, OverlapData *dat, sg_biedge_t *e){
	int len1, len2, l[3], r[3];
	len1 = g->rdlens->buffer[dat->node_id[0]];
	len2 = g->rdlens->buffer[dat->node_id[1]];
	l[0] = dat->beg[0];
	l[1] = dat->beg[1];
	r[0] = len1 - dat->end[0];
	r[1] = len2 - dat->end[1];
	l[2] = l[0] < l[1]? l[0] : l[1];
	r[2] = r[0] < r[1]? r[0] : r[1];
	if(l[2] + r[2] > g->max_ovl_margin) return 0;
	if(l[0] >= l[1]){
		e->node_id[0] = dat->node_id[0];
		e->node_id[1] = dat->node_id[1];
		e->dir[0] = dat->dir[0];
		e->dir[1] = dat->dir[1];
		e->off[0] = l[0] - l[2];
		e->off[1] = r[1] - r[2];
		e->ol[0]  = dat->end[0] - dat->beg[0];
		e->ol[1]  = dat->end[1] - dat->beg[1];
	} else  {
		e->node_id[1] = dat->node_id[0];
		e->node_id[0] = dat->node_id[1];
		e->dir[1] = dat->dir[0];
		e->dir[0] = dat->dir[1];
		e->off[1] = r[0] - r[2];
		e->off[0] = l[1] - l[2];
		e->ol[1]  = dat->end[0] - dat->beg[0];
		e->ol[0]  = dat->end[1] - dat->beg[1];
	}
	e->score = dat->score;
	return 1;
}

#define _ovl_uniq_long_id(id1, id2, dir) ((((uint64_t)(id1)) << 33) | (((uint64_t)(id2)) << 1) | (dir))
#define ovl_uniq_long_id(id1, id2, dir) (((id1) < (id2))? _ovl_uniq_long_id(id1, id2, dir) : _ovl_uniq_long_id(id2, id1, dir))

void load_overlaps_core_strgraph(StringGraph *g, sgbiedgev *biedges){
	sg_biedge_t *b;
	sg_edge_t *e1, *e2;
	sg_node_t *n;
	uint64_t ret, off, i;
	uint32_t k, idx[2];
	int len, len1, len2;
	ret = 0;
	fprintf(stderr, "building edges\n");
	off = 0;
	for(k=0;k<2;k++){
		for(i=0;i<g->nodes->size;i++){
			n = ref_sgnodev(g->nodes, i);
			n->edge_offs[k] = off;
			off += n->edge_cnts[k];
			n->edge_cnts[k] = 0;
		}
	}
	clear_sgedgev(g->edges);
	encap_sgedgev(g->edges, off);
	g->edges->size = off;
	for(i=0;i<biedges->size;i++){
		b = ref_sgbiedgev(biedges, i);
		if((i % 10000) == 0){ fprintf(stderr, "\r%llu overlaps", (unsigned long long)i); fflush(stderr); }

		len1 = g->rdlens->buffer[b->node_id[0]];
		len2 = g->rdlens->buffer[b->node_id[1]];

		n = ref_sgnodev(g->nodes, b->node_id[0]);
		k = b->dir[0];
		idx[0] = n->edge_cnts[k];
		e1 = ref_sgedgev(g->edges, n->edge_offs[k] + n->edge_cnts[k]);
		n->edge_cnts[k] ++;
		e1->node_id = b->node_id[1];
		e1->dir = b->dir[1];
		e1->closed = 0;
		e1->mark   = 0;
		e1->off = b->off[0];
		len = (e1->off + len2 > len1)? len1 - e1->off : len2;
		e1->ol_var = b->ol[0] - len;
		e1->score = b->score;
		e1->att = 0;
		e1->tta = 0;

		n = ref_sgnodev(g->nodes, b->node_id[1]);
		k = !b->dir[1];
		idx[1] = n->edge_cnts[k];
		e2 = ref_sgedgev(g->edges, n->edge_offs[k] + n->edge_cnts[k]);
		n->edge_cnts[k] ++;
		e2->node_id = b->node_id[0];
		e2->dir = !b->dir[0];
		e2->closed = 0;
		e2->mark   = 0;
		e2->off = b->off[1];
		len = (e2->off + len1 > len2)? len2 - e2->off : len1;
		e2->ol_var = b->ol[1] - len;
		e2->score = b->score;
		e2->att = 0;
		e2->tta = 0;

		e1->rev_idx = idx[1];
		e2->rev_idx = idx[0];

		if(b->off[0] == 0){
			if(b->off[1] == 0){
				if(len1 < len2){
					e1->att = 1;
					e2->tta = 1;
				} else if(len1 > len2){
					e2->att = 1;
					e1->tta = 1;
				} else if(b->node_id[0] < b->node_id[1]){
					e2->att = 1;
					e1->tta = 1;
				} else {
					e1->att = 1;
					e2->tta = 1;
				}
			} else {
				e1->att = 1;
				e2->tta = 1;
			}
		} else if(b->off[1] == 0){
			e2->att = 1;
			e1->tta = 1;
		}
	}
	fprintf(stderr, "\r%llu fine overlaps\n", (unsigned long long)i);
}

sgbiedgev* load_overlaps_strgraph(StringGraph *g, FileReader *fr, u64hash *closed_alns){
	sgbiedgev *biedges;
	sg_biedge_t B, *b;
	OverlapData O;
	uint64_t ret;
	b = &B;
	ret = 0;
	biedges = init_sgbiedgev(1024);
	while(parse_overlap_item_strgraph(g, fr, &O)){
		ret ++;
		if((ret % 10000) == 0){ fprintf(stderr, "\r%llu overlaps", (unsigned long long)ret); fflush(stdout); }
		if(closed_alns){
			put_u64hash(closed_alns, ovl_uniq_long_id(O.node_id[0], O.node_id[1], 0));
			put_u64hash(closed_alns, ovl_uniq_long_id(O.node_id[0], O.node_id[1], 1));
		}
		if(overlap_item2biedge_strgraph(g, &O, b) == 0) continue;
		if(ref_sgnodev(g->nodes, b->node_id[0])->edge_cnts[b->dir[0]] >= SG_MAX_EDGE) continue;
		if(ref_sgnodev(g->nodes, b->node_id[1])->edge_cnts[!b->dir[1]] >= SG_MAX_EDGE) continue;
		ref_sgnodev(g->nodes, b->node_id[0])->edge_cnts[b->dir[0]] ++;
		ref_sgnodev(g->nodes, b->node_id[1])->edge_cnts[!b->dir[1]] ++;
		push_sgbiedgev(biedges, B);
	}
	fprintf(stderr, "\rloaded %llu overlaps\n", (unsigned long long)ret);
	load_overlaps_core_strgraph(g, biedges);
	return biedges;
}

inline sg_edge_t* edge_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx){
	sg_node_t *n;
	sg_edge_t *e;
	n = ref_sgnodev(g->nodes, node_id);
	e = ref_sgedgev(g->edges, n->edge_offs[dir] + eidx);
	return e;
}

sg_edge_t* node2node_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t node_id2){
	sg_node_t *n;
	sg_edge_t *e;
	uint32_t i;
	n = ref_sgnodev(g->nodes, node_id);
	for(i=0;i<n->edge_cnts[dir];i++){
		e = ref_sgedgev(g->edges, n->edge_offs[dir] + i);
		if(e->closed) continue;
		if(e->node_id == node_id2) return e;
	}
	return NULL;
}

inline int cut_edge_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx){
	sg_edge_t *e;
	e = edge_strgraph(g, node_id, dir, eidx);
	if(e->closed) return 0;
	e->closed = 1;
	return 1;
}

inline sg_edge_t* edge_partner_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx){
	sg_edge_t *e1, *e2;
	e1 = edge_strgraph(g, node_id, dir, eidx);
	e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
	return e2;
}

inline int cut_biedge_strgraph_core(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx, int closed){
	edge_strgraph(g, node_id, dir, eidx)->closed = closed;
	edge_partner_strgraph(g, node_id, dir, eidx)->closed = closed;
	return 1;
}

inline int cut_biedge_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx){
	return cut_biedge_strgraph_core(g, node_id, dir, eidx, 1);
}

inline int cut_biedge_strgraph2(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx){
	return cut_biedge_strgraph_core(g, node_id, dir, eidx, 2);
}

uint32_t count_living_edge_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e;
	uint32_t i, cnt;
	n = ref_sgnodev(g->nodes, node_id);
	for(i=cnt=0;i<n->edge_cnts[dir];i++){
		e = ref_sgedgev(g->edges, n->edge_offs[dir] + i);
		if(e->closed == 0) cnt ++;
	}
	return cnt;
}

sg_edge_t* first_living_edge_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e;
	uint32_t i;
	n = ref_sgnodev(g->nodes, node_id);
	for(i=0;i<n->edge_cnts[dir];i++){
		e = ref_sgedgev(g->edges, n->edge_offs[dir] + i);
		if(e->closed == 0) return e;
	}
	return NULL;
}

sg_edge_t* single_living_edge_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e, *ret;
	uint32_t i;
	n = ref_sgnodev(g->nodes, node_id);
	ret = NULL;
	for(i=0;i<n->edge_cnts[dir];i++){
		e = ref_sgedgev(g->edges, n->edge_offs[dir] + i);
		if(e->closed) continue;
		if(ret) return NULL;
		else ret = e;
	}
	return ret;
}

uint32_t edge_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx){
	sg_edge_t *e;
	int len1, len2, len;
	e = edge_strgraph(g, node_id, dir, eidx);
	len1 = g->rdlens->buffer[node_id];
	len2 = g->rdlens->buffer[e->node_id];
	len = (e->off + len2 > len1)? len1 - e->off : len2;
	return len + e->ol_var;
}

void print_edge_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t eidx, FILE *out){
	sg_node_t *n;
	sg_edge_t *e;
	n = ref_sgnodev(g->nodes, node_id);
	e = ref_sgedgev(g->edges, n->edge_offs[dir] + eidx);
	fprintf(out, "EDGE %s\t%c\t%s\t%c\t%d\t%d\t%d\n", (char*)get_cplist(g->rdnames, node_id), "+-"[dir], (char*)get_cplist(g->rdnames, e->node_id), "+-"[e->dir], e->off, edge_overlap_strgraph(g, node_id, dir, eidx), e->score);
}

void mask_node_strgraph(StringGraph *g, uint32_t node_id){
	sg_node_t *n;
	uint32_t k, j;
	n = ref_sgnodev(g->nodes, node_id);
	for(k=0;k<2;k++){
		for(j=0;j<n->edge_cnts[k];j++){
			cut_biedge_strgraph(g, node_id, k, j);
		}
	}
	one_bitvec(g->node_status, node_id);
}

void mask_node_strgraph2(StringGraph *g, uint32_t node_id){
	sg_node_t *n;
	uint32_t k, j;
	n = ref_sgnodev(g->nodes, node_id);
	for(k=0;k<2;k++){
		for(j=0;j<n->edge_cnts[k];j++){
			cut_biedge_strgraph2(g, node_id, k, j);
		}
	}
	one_bitvec(g->node_status, node_id);
}

uint64_t remove_duplicate_edges_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_edge_t *e, *e2;
	uuhash *hash;
	uint64_t ret;
	uint32_t i, j, k;
	uuhash_t *c;
	int exists;
	ret = 0;
	hash = init_uuhash(31);
	for(i=0;i<g->nodes->size;i++){
		n = ref_sgnodev(g->nodes, i);
		if(is_dead_node_strgraph(g, i)) continue;
		for(k=0;k<2;k++){
			clear_uuhash(hash);
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				c = kv_prepare_uuhash(hash, e->node_id, &exists);
				if(exists){
					ret ++;
					e2 = ref_sgedgev(g->edges, n->edge_offs[k] + c->val);
					if(e->score < e2->score){
						cut_biedge_strgraph(g, i, k, j);
					} else {
						cut_biedge_strgraph(g, i, k, c->val);
						c->val = j;
					}
				} else {
					c->key = e->node_id;
					c->val = j;
				}
			}
		}
	}
	free_uuhash(hash);
	return ret;
}

void cal_edge_coverage_strgraph(StringGraph *g){
	sg_node_t *n, *n2;
	sg_edge_t *e, *e2;
	u32hash *hash;
	uint64_t I;
	uint32_t i, k, j, k2, j2, cov;
	for(I=0;I<g->edges->size;I++) ref_sgedgev(g->edges, I)->cov = SG_EDGE_MAX_COV;
	hash = init_u32hash(1023);
	for(i=0;i<g->nodes->size;i++){
		n = ref_sgnodev(g->nodes, i);
		clear_u32hash(hash);
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed == 1) continue;
				put_u32hash(hash, e->node_id);
			}
		}
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed == 1) continue;
				if(e->cov != SG_EDGE_MAX_COV) continue;
				cov = 0;
				n2 = ref_sgnodev(g->nodes, e->node_id);
				for(k2=0;k2<2;k2++){
					for(j2=0;j2<n2->edge_cnts[k2];j2++){
						e2 = ref_sgedgev(g->edges, n2->edge_offs[k2] + j2);
						if(e2->closed == 1) continue;
						if(exists_u32hash(hash, e2->node_id)) cov ++;
					}
				}
				if(cov + 1 >= SG_EDGE_MAX_COV) cov = SG_EDGE_MAX_COV - 1;
				e->cov = cov;
				e2 = ref_sgedgev(g->edges, n2->edge_offs[!e->dir] + e->rev_idx);
				e2->cov = cov;
			}
		}
	}
	free_u32hash(hash);
}

uint64_t mask_low_cov_edge_strgraph(StringGraph *g, uint32_t cutoff){
	sg_edge_t *e;
	uint64_t i, ret;
	ret = 0;
	if(cutoff == 0) return ret;
	for(i=0;i<g->edges->size;i++){
		e = ref_sgedgev(g->edges, i);
		if(e->closed == 1) continue;
		if(e->cov >= cutoff) continue;
		e->closed = 1;
		ret ++;
	}
	return ret;
}

uint32_t mask_contained_reads_strgraph(StringGraph *g, FILE *log){
	sg_node_t *n;
	sg_edge_t *e;
	BitVec *flags;
	uint32_t *map;
	uint32_t i, j, k, c, len, ret;
	int max_score;
	ret = 0;
	flags = init_bitvec(g->nodes->size);
	map = malloc(sizeof(uint32_t) * g->nodes->size);
	memset(map, 0xff, sizeof(uint32_t) * g->nodes->size);
	for(i=0;i<g->nodes->size;i++){
		n = ref_sgnodev(g->nodes, i);
		if(is_dead_node_strgraph(g, i)) continue;
		len = g->rdlens->buffer[i];
		c = 0;
		for(k=0;c==0&&k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed == 1) continue;
				if(e->att){ c = 1; break; }
			}
		}
		if(c == 0) continue;
		one_bitvec(flags, i);
		ret ++;
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_sgnodev(g->nodes, i);
		if(is_dead_node_strgraph(g, i)) continue;
		if(!get_bitvec(flags, i)) continue;
		c = 0xFFFFFFFFU;
		max_score = 0;
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed == 1) continue;
				if(!e->att) continue;
				if(get_bitvec(flags, e->node_id)){
					if(c == 0xFFFFFFFFU) c = e->node_id;
					continue;
				}
				if(e->score > max_score){
					c = e->node_id;
					max_score = e->score;
				}
			}
		}
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->node_id != c) e->att = 0;
			}
		}
		map[i] = c;
	}
	for(i=0;i<g->nodes->size;i++){
		if(get_bitvec(flags, i)){
			mask_node_strgraph(g, i);
			if(log) fprintf(log, "%s\t%s\n", get_cplist(g->rdnames, i), (map[i] == 0xFFFFFFFFU)? "NULL" : get_cplist(g->rdnames, map[i]));
		}
	}
	free_bitvec(flags);
	free(map);
	return ret;
}

uint64_t best_overlap_strgraph(StringGraph *g, float best_score_cutoff){
	sg_node_t *n;
	sg_node_t *m;
	sg_edge_t *e, *p;
	uint64_t ret;
	uint32_t i, k, j, b;
	int best;
	float bestS;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			bestS = 0;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed || e->att || e->tta) continue;
				if(e->score > bestS){ bestS = e->score; }
			}
			bestS = bestS * best_score_cutoff;
			best = g->rdlens->buffer[i]; b = 0xFFFFFFFFU;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed || e->att || e->tta) continue;
				if(e->score < bestS) continue;
				if(e->off < best){ best = e->off; b = j; }
			}
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(j != b){
					e->closed = 1;
					ret ++;
				}
			}
		}
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_sgnodev(g->nodes, i);
		memset(n->bogs, 0, 8);
	}
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				m = ref_sgnodev(g->nodes, e->node_id);
				p = edge_partner_strgraph(g, i, k, j);
				if(p->closed){
					e->mark = 1;
					if(n->bogs[1][k][1] < 0xFFU) n->bogs[1][k][1] ++;
					if(m->bogs[0][e->dir][1] < 0xFFU) m->bogs[0][e->dir][1] ++;
				} else {
					e->mark = 0;
					if(n->bogs[1][k][0] < 0xFFU) n->bogs[1][k][0] ++;
					if(m->bogs[0][e->dir][0] < 0xFFU) m->bogs[0][e->dir][0] ++;
				}
			}
		}
	}
	return ret;
}

#ifdef _CPLUSPLUS
}
#endif

#endif
