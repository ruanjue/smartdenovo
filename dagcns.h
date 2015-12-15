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

#include "dna.h"
#include "string.h"
#include "list.h"
#include "queue.h"
#include "hashset.h"
#include <float.h>

#ifndef PACBIO_DAGCNS_RJ_H
#define PACBIO_DAGCNS_RJ_H

static const int debug = 0;

typedef struct {
	uint32_t nodes[2];
	uint32_t cov:28, visit:1, closed:1, cns:1, alt:1;
	uint32_t links[2];
} dagedge_t;
define_list(dagedgev, dagedge_t);

#define NODE_MAX_FW_EDGE	0xFFFFFFFFU
typedef struct dagnode_t {
	uint32_t pos:29, base:2, cns:1;
	float    aux;
	uint32_t edges[2];
	uint32_t fw_edge;
} dagnode_t;
define_list(dagnodev, dagnode_t);

typedef struct {
	uint32_t pos;
	uint32_t bases[4];
} dagsnp_t;
define_list(dagsnpv, dagsnp_t);

typedef struct {
	u8list   *cns;
	u32list  *deps;
	dagnodev *nodes;
	dagedgev *edges;
	u32list  *trash;
	String   *alns[2];
	float    ref_penalty, alt_penalty; // 0.5 and 0.2
	double   cns_score;
	uint32_t cns_head, backbone_size;
} DAGCNS;

static inline DAGCNS* init_dagcns(){
	DAGCNS *g;
	g = malloc(sizeof(DAGCNS));
	g->cns = init_u8list(1024);
	g->deps = init_u32list(1024);
	g->nodes = init_dagnodev(1024);
	g->edges = init_dagedgev(1024);
	g->trash = init_u32list(1024);
	g->alns[0] = init_string(1024);
	g->alns[1] = init_string(1024);
	g->cns_score = 0;
	g->cns_head = 0xFFFFFFFFU;
	g->backbone_size = 0;
	g->ref_penalty = 0.5;
	g->alt_penalty = 0.2;
	return g;
}

static inline void free_dagcns(DAGCNS *g){
	free_dagnodev(g->nodes);
	free_dagedgev(g->edges);
	free_u32list(g->trash);
	free_u8list(g->cns);
	free_u32list(g->deps);
	free_string(g->alns[0]);
	free_string(g->alns[1]);
	free(g);
}

static inline void reset_dagcns(DAGCNS *g){
	clear_dagnodev(g->nodes);
	clear_dagedgev(g->edges);
	clear_u32list(g->trash);
	clear_u8list(g->cns);
	clear_u32list(g->deps);
	g->cns_score = 0;
	g->cns_head = 0xFFFFFFFFU;
	g->backbone_size = 0;
}

static uint32_t prepare_node_dagcns(DAGCNS *g, uint32_t pos, uint8_t base){
	dagnode_t *n;
	n = next_ref_dagnodev(g->nodes);
	n->pos = pos;
	n->base = base;
	n->cns = 0;
	n->aux  = 0;
	n->edges[0] = 0xFFFFFFFFU;
	n->edges[1] = 0xFFFFFFFFU;
	n->fw_edge = NODE_MAX_FW_EDGE;
	return g->nodes->size - 1;
}

static inline dagedge_t* find_edge_by_node_dagcns(DAGCNS *g, uint32_t nid1, uint32_t nid2, int dir){
	dagnode_t *n;
	dagedge_t *e;
	n = ref_dagnodev(g->nodes, nid1);
	if(n->edges[dir] != 0xFFFFFFFFU){
		e = ref_dagedgev(g->edges, n->edges[dir]);
		while(1){
			if(e->nodes[dir] == nid2) return e;
			if(e->links[dir] == 0xFFFFFFFFU) break;
			e = ref_dagedgev(g->edges, e->links[dir]);
		}
	}
	return NULL;
}

static inline dagedge_t* add_edge_dagcns(DAGCNS *g, uint32_t nid1, uint32_t nid2, int dir){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t eid;
	n = ref_dagnodev(g->nodes, nid1);
	if(pop_u32list(g->trash, &eid)){
		e = ref_dagedgev(g->edges, eid);
	} else {
		eid = g->edges->size;
		e = next_ref_dagedgev(g->edges);
	}
	e->nodes[!dir] = nid1;
	e->nodes[dir] = nid2;
	e->cov = 1;
	e->visit = 0;
	e->closed = 0;
	e->cns = 0;
	e->alt = 0;
	e->links[dir] = n->edges[dir];
	n->edges[dir] = eid;
	n = ref_dagnodev(g->nodes, nid2);
	e->links[!dir] = n->edges[!dir];
	n->edges[!dir] = eid;
	return e;
}

static inline dagedge_t* prepare_edge_dagcns(DAGCNS *g, uint32_t nid1, uint32_t nid2, int dir){
	dagedge_t *e;
	e = find_edge_by_node_dagcns(g, nid1, nid2, dir);
	if(e){ e->cov ++; return e; }
	return add_edge_dagcns(g, nid1, nid2, dir);
}

static inline void gen_pregraph_dagcns(DAGCNS *g){
	dagedge_t *e;
	uint32_t i;
	clear_dagnodev(g->nodes);
	clear_dagedgev(g->edges);
	clear_u32list(g->trash);
	clear_u32list(g->deps);
	g->backbone_size = g->cns->size;
	for(i=0;i<g->cns->size;i++){
		push_u32list(g->deps, 0);
		prepare_node_dagcns(g, i, g->cns->buffer[i]);
		if(i){ // make sure the graph is conntective even the alignment is partial
			e = add_edge_dagcns(g, i - 1, i, 0);
			e->cov = 0;
		}
	}
}

static inline int remove_edge_dagcns(DAGCNS *g, uint32_t eid){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t i, lst;
	for(i=0;i<2;i++){
		e = ref_dagedgev(g->edges, eid);
		lst = e->links[i];
		n = ref_dagnodev(g->nodes, e->nodes[!i]);
		if(n->edges[i] == eid){
			n->edges[i] = lst;
		} else if(n->edges[i] == 0xFFFFFFFFU){
			return 0;
		} else {
			e = ref_dagedgev(g->edges, n->edges[i]);
			while(1){
				if(e->links[i] == eid){
					e->links[i] = lst; break;
				} else if(e->links[i] == 0xFFFFFFFFU){
					return 0;
				} else {
					e = ref_dagedgev(g->edges, e->links[i]);
				}
			}
		}
	}
	push_u32list(g->trash, eid);
	return 1;
}

static inline void _polish_pairwise_aln_dagcns(char *alns[2], int len){
	int i, j, m, c, gaps[2];
	while(1){
		c = 0;
		gaps[0] = gaps[1] = 0;
		for(i=0;i<len;i++){
			for(j=0;j<2;j++){
				if(alns[j][i] == '-'){ gaps[j] ++; continue; }
				if(gaps[j] == 0) continue;
				for(m=i-gaps[j];m<i;m++){
					if(alns[!j][m] == alns[j][i]){
						alns[j][m] = alns[j][i];
						alns[j][i] = '-';
						c ++;
						break;
					}
				}
				gaps[j] = i - m;
			}
		}
		if(c == 0) break;
	}
	if(debug) fprintf(stdout, "%s\n%s\n", alns[0], alns[1]);
}

static inline void polish_pairwise_aln_dagcns(char *alns[2], int len, String *rets[2]){
	int i, size;
	clear_string(rets[0]);
	clear_string(rets[1]);
	encap_string(rets[0], len * 2);
	encap_string(rets[1], len * 2);
	for(i=size=0;i<len;i++){
		if(alns[0][i] != alns[1][i] && alns[0][i] != '-' && alns[1][i] != '-'){
			rets[0]->string[size] = '-';
			rets[1]->string[size] = alns[1][i];
			size ++;
			rets[0]->string[size] = alns[0][i];
			rets[1]->string[size] = '-';
			size ++;
		} else {
			rets[0]->string[size] = alns[0][i];
			rets[1]->string[size] = alns[1][i];
			size ++;
		}
	}
	rets[0]->size = rets[1]->size = size;
	rets[0]->string[size] = rets[1]->string[size] = 0;
	_polish_pairwise_aln_dagcns((char*[2]){rets[0]->string, rets[1]->string}, size);
}

static inline void alignment2dagcns(DAGCNS *g, int beg, int end, char *alns[2], int len){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t lst, cur, eid, base;
	int i, j, x1;
	polish_pairwise_aln_dagcns(alns, len, g->alns);
	len = g->alns[0]->size;
	alns[0] = g->alns[0]->string;
	alns[1] = g->alns[1]->string;
	x1 = beg;
	lst = 0xFFFFFFFFU;
	//while(len && (alns[0][len-1] == '-' || alns[1][len-1] == '-')) len --;
	while(len && alns[0][len-1] == '-') len --;
	for(i=0;i<len;i++){
		if(alns[0][i] == alns[1][i]){
			if(alns[0][i] == '-') continue;
			cur = x1 ++;
			if(lst == 0xFFFFFFFFU){ lst = cur; continue; }
			e = prepare_edge_dagcns(g, lst, cur, 0);
			lst = cur;
		} else if(alns[0][i] == '-'){
			if(lst == 0xFFFFFFFFU) continue;
			base = base_bit_table[(int)alns[1][i]];
			n = ref_dagnodev(g->nodes, lst);
			eid = n->edges[0];
			cur = 0xFFFFFFFFU;
			while(eid != 0xFFFFFFFFU){
				e = ref_dagedgev(g->edges, eid);
				if(e->nodes[0] >= g->cns->size && ref_dagnodev(g->nodes, e->nodes[0])->base == base){
					e->cov ++; cur = e->nodes[0]; break;
				}
				eid = e->links[0];
			}
			if(cur == 0xFFFFFFFFU){
				cur = prepare_node_dagcns(g, x1, base);
				//add_edge_dagcns(g, lst, cur, 0);
				prepare_edge_dagcns(g, lst, cur, 0);
			}
			lst = cur;
		} else {
			x1 ++;
		}
	}
	for(j=beg;j<end;j++) g->deps->buffer[j] ++;
}

static inline void merge_nodes_core_dagcns(DAGCNS *g, uint32_t nid, u32list *stack, u32list *cache[4], int dir){
	dagnode_t *n0, *n2, *n;
	dagedge_t *e, *e2, *e1;
	uint32_t base, eid, nid1, i, ret;
	clear_u32list(stack);
	push_u32list(stack, nid);
	ret = 0;
	while(pop_u32list(stack, &nid)){
		n0 = ref_dagnodev(g->nodes, nid);
		if((eid = n0->edges[dir]) == 0xFFFFFFFFU) continue;
		clear_u32list(cache[0]);
		clear_u32list(cache[1]);
		clear_u32list(cache[2]);
		clear_u32list(cache[3]);
		while(1){
			e = ref_dagedgev(g->edges, eid);
			n = ref_dagnodev(g->nodes, e->nodes[dir]);
			e2 = ref_dagedgev(g->edges, n->edges[!dir]);
			if(e2->links[!dir] == 0xFFFFFFFFU){ // check whether there is only one edge from n -(!dir)-> n0
				push_u32list(cache[n->base], eid);
			}
			if((eid = e->links[dir]) == 0xFFFFFFFFU) break;
		}
		for(base=0;base<4;base++){
			for(i=0;i<cache[base]->size;i++) ref_dagedgev(g->edges, cache[base]->buffer[i])->visit = 1;
			if(cache[base]->size < 2) continue;
			e1 = ref_dagedgev(g->edges, cache[base]->buffer[0]);
			nid1 = e1->nodes[dir];
			for(i=1;i<cache[base]->size;i++){
				e2 = ref_dagedgev(g->edges, cache[base]->buffer[i]);
				n2 = ref_dagnodev(g->nodes, e2->nodes[dir]);
				e1->cov += e2->cov;
				remove_edge_dagcns(g, cache[base]->buffer[i]);
				eid = n2->edges[dir];
				while(eid != 0xFFFFFFFFU){
					e2 = ref_dagedgev(g->edges, eid);
					e  = prepare_edge_dagcns(g, nid1, e2->nodes[dir], dir);
					{
						e1 = ref_dagedgev(g->edges, cache[base]->buffer[0]); // memory referred by e1 may be freed in prepare_edge_dagcns
					}
					e->cov = e->cov - 1 + e2->cov;
					e->visit = 1;
					eid = e2->links[dir];
				}
				eid = n2->edges[dir];
				while(eid != 0xFFFFFFFFU){
					e2 = ref_dagedgev(g->edges, eid);
					remove_edge_dagcns(g, eid); // e2->links retain the same values after removing
					eid = e2->links[dir];
				}
			}
			push_u32list(stack, nid1);
		}
	}
}

static inline int has_non_visited_edge_dagcns(DAGCNS *g, uint32_t nid, int dir){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t eid;
	n = ref_dagnodev(g->nodes, nid);
	eid = n->edges[dir];
	while(eid != 0xFFFFFFFFU){
		e = ref_dagedgev(g->edges, eid);
		if(e->visit == 0) return 1;
		eid = e->links[dir];
	}
	return 0;
}

static inline void print_local_dot_dagcns(DAGCNS *g, uint32_t nid, int distance, FILE *out){
	u32list *stack;
	u32hash *hash;
	dagnode_t *n, *n1, *n2;
	dagedge_t *e;
	uint32_t id1, id2, eid, *u;
	int lo, hi, dir, exists;
	n = ref_dagnodev(g->nodes, nid);
	stack = init_u32list(32);
	hash = init_u32hash(1023);
	lo = n->pos - distance;
	hi = n->pos + distance;
	push_u32list(stack, nid);
	put_u32hash(hash, nid);
	fprintf(out, "digraph {\n");
	while(stack->size){
		id1 = stack->buffer[--stack->size];
		n1 = ref_dagnodev(g->nodes, id1);
		for(dir=0;dir<1;dir++){
			eid = n1->edges[dir];
			while(eid != 0xFFFFFFFFU){
				e = ref_dagedgev(g->edges, eid);
				id2 = e->nodes[dir];
				n2 = ref_dagnodev(g->nodes, id2);
				fprintf(out, "N%d_%d_%c -> N%d_%d_%c [label=\"%d:%d\"]\n", id1, n1->pos, "ACGT"[n1->base], id2, n2->pos, "ACGT"[n2->base], e->cov, e->visit);
				if(n2->pos >= lo && n2->pos <= hi){
					u = prepare_u32hash(hash, id2, &exists);
					if(exists){
					} else {
						*u = id2;
						push_u32list(stack, id2);
					}
				}
				eid = e->links[dir];
			}
		}
	}
	fprintf(out, "}\n");
	fflush(out);
	free_u32list(stack);
	free_u32hash(hash);
}

static inline void merge_nodes_dagcns(DAGCNS *g){
	dagnode_t *n;
	dagedge_t *e;
	u32list *stack, *cache[4];
	u32fifo *queue;
	uint32_t i, nid, eid;
	stack = init_u32list(1024);
	cache[0] = init_u32list(4);
	cache[1] = init_u32list(4);
	cache[2] = init_u32list(4);
	cache[3] = init_u32list(4);
	for(i=0;i<g->edges->size;i++) g->edges->buffer[i].visit = 0;
	queue = init_u32fifo();
	for(i=0;i<g->nodes->size;i++){
		n = ref_dagnodev(g->nodes, i);
		if(n->edges[1] != 0xFFFFFFFFU) continue;
		push_u32fifo(queue, i);
	}
	if(i == 0xFFFFFFFFU){ // un-reachable, but is used to call print_local_dot_dagcns in gdb debug
		print_local_dot_dagcns(g, 0, 10, stdout);
	}
	while(pop_u32fifo(queue, &nid)){
		if(debug > 1) fprintf(stdout, "\npop %u\n", nid);
		merge_nodes_core_dagcns(g, nid, stack, cache, 1);
		merge_nodes_core_dagcns(g, nid, stack, cache, 0);
		n = ref_dagnodev(g->nodes, nid);
		eid = n->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			e->visit = 1;
			eid = e->links[0];
		}
		eid = n->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			if(!has_non_visited_edge_dagcns(g, e->nodes[0], 1)){
				if(debug > 1) fprintf(stdout, "push %u\n", e->nodes[0]);
				push_u32fifo(queue, e->nodes[0]);
			}
			eid = e->links[0];
		}
	}
	free_u32fifo(queue);
	free_u32list(stack);
	free_u32list(cache[0]);
	free_u32list(cache[1]);
	free_u32list(cache[2]);
	free_u32list(cache[3]);
}

static inline void print_seq_dagcns(DAGCNS *g, FILE *out){
	uint32_t i;
	for(i=0;i<g->cns->size;i++){
		fputc(bit_base_table[g->cns->buffer[i]], out);
		if((i % 100) == 99) fputc('\n', out);
	}
	if((i % 100)) fputc('\n', out);
}

static inline void gen_consensus_dagcns(DAGCNS *g, u32list *map){
	dagnode_t *n1, *n2;
	dagedge_t *e;
	u32fifo *queue;
	uint32_t i, lst, nid, eid, best_e;
	float best_s, score;
	queue = init_u32fifo();
	for(i=0;i<g->nodes->size;i++){
		n1 = ref_dagnodev(g->nodes, i);
		if(n1->edges[0] == 0xFFFFFFFFU && n1->edges[1] != 0xFFFFFFFFU){
			push_u32fifo(queue, i);
			n1->fw_edge = NODE_MAX_FW_EDGE;
			n1->aux = 0;
		}
	}
	for(i=0;i<g->edges->size;i++) g->edges->buffer[i].visit = 0;
	while(pop_u32fifo(queue, &nid)){
		best_s = - FLT_MAX;
		best_e = 0xFFFFFFFFU;
		n1 = ref_dagnodev(g->nodes, nid);
		eid = n1->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			n2 = ref_dagnodev(g->nodes, e->nodes[0]);
			if(e->nodes[0] < g->backbone_size){
				score = n2->aux + e->cov - g->ref_penalty * g->deps->buffer[n1->pos];
			} else {
				score = n2->aux + e->cov - g->alt_penalty * g->deps->buffer[n1->pos];
			}
			if(score > best_s){
				best_s = score;
				best_e = eid;
			}
			eid = e->links[0];
		}
		if(best_s > - FLT_MAX) n1->aux = best_s;
		n1->fw_edge = best_e;
		eid = n1->edges[1];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			e->visit = 1;
			if(!has_non_visited_edge_dagcns(g, e->nodes[1], 0)){
				push_u32fifo(queue, e->nodes[1]);
			}
			eid = e->links[1];
		}
	}
	free_u32fifo(queue);
	clear_u8list(g->cns);
	clear_u32list(g->deps);
	if(map) clear_u32list(map);
	g->cns_head = 0;
	n1 = ref_dagnodev(g->nodes, g->cns_head);
	g->cns_score = n1->aux;
	n1->cns = 1;
	lst = 0;
	if(map && g->cns_head < g->backbone_size){
		while(lst < g->cns_head){ push_u32list(map, g->cns->size); lst ++; }
	}
	push_u8list(g->cns, n1->base);
	push_u32list(g->deps, 0);
	while(n1->fw_edge != NODE_MAX_FW_EDGE){
		e = ref_dagedgev(g->edges, n1->fw_edge);
		e->cns = 1;
		if(map && e->nodes[0] < g->backbone_size){
			while(lst < e->nodes[0]){ push_u32list(map, g->cns->size); lst ++; }
		}
		n1 = ref_dagnodev(g->nodes, e->nodes[0]);
		n1->cns = 1;
		push_u8list(g->cns, n1->base);
		push_u32list(g->deps, 0);
	}
	if(map) while(lst <= g->backbone_size){ push_u32list(map, g->cns->size); lst ++; }
}


// min_freq: freq = alternative allel count / consensus allel count
static inline uint32_t call_snv_dagcns(DAGCNS *g, dagsnpv *snps, int min_cnt, float min_freq){
	dagnode_t *n, *n2;
	dagedge_t *e, *e2;
	uint32_t eid, nid[3];
	uint32_t c, cns, ret;
	dagsnp_t P;
	ret = 0;
	P.pos = 0;
	nid[0] = g->cns_head;
	while(P.pos + 2 < g->cns->size){
		P.pos ++;
		P.bases[0] = P.bases[1] = P.bases[2] = P.bases[3] = cns = 0;
		n = ref_dagnodev(g->nodes, nid[0]);
		nid[1] = ref_dagedgev(g->edges, n->fw_edge)->nodes[0];
		n2 = ref_dagnodev(g->nodes, nid[1]);
		cns = n2->base;
		P.bases[cns] = num_min(ref_dagedgev(g->edges, n->fw_edge)->cov, ref_dagedgev(g->edges, n2->fw_edge)->cov);
		nid[2] = ref_dagedgev(g->edges, n2->fw_edge)->nodes[0];
		eid = n->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			n2 = ref_dagnodev(g->nodes, e->nodes[0]);
			if(!e->cns){
				e2 = find_edge_by_node_dagcns(g, e->nodes[0], nid[2], 0);
				if(e2){
					c = num_min(e->cov, e2->cov);
					if(c > P.bases[n2->base]) P.bases[n2->base] = c;
				}
			}
			eid = e->links[0];
		}
		for(c=0;c<4;c++){
			if(c == cns) continue;
			if(P.bases[c] >= (uint32_t)min_cnt && P.bases[c] >= P.bases[cns] * min_freq){
				push_dagsnpv(snps, P);
				ret ++;
				break;
			}
		}
		nid[0] = nid[1];
	}
	return ret;
}

#endif
