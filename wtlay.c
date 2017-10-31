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
 
#include "wtlay.h"
#include "timer.h"

static uint32_t MERGE_BUBBLE_MAX_STEP = 20;
static uint32_t CUT_LOOP_MAX_STEP     = 5;
static uint32_t MIN_LAY_NODES = 4;

define_list(f32list, float);

int _test_clustering(float *dat, uint32_t size){
	double sum, var, mean, std;
	uint32_t i;
	sum = 0;
	var = 0;
	for(i=0;i<size;i++){
		sum += dat[i];
		var += dat[i] * dat[i];
	}
	mean = sum / size;
	std = sqrt((sum * sum - var) / size);
	if(dat[size - 1] < mean - 3 * std) return 0;
	return 1;
}

typedef struct {
	uint32_t node_id;
	int x, y;
} pb_reg_t;
define_list(pbregv, pb_reg_t);

uint32_t detect_chimeric_reads_strgraph(StringGraph *g, FILE *log){
	sg_node_t *n, *n1;
	sg_edge_t *e;
	uuhash *hash;
	pbregv *regs;
	pb_reg_t *r1, *r2;
	uint32_t ret;
	uint32_t nidx, ndir, eidx, f, i, k, j, ol, rid, len, max;
	ret = 0;
	hash = init_uuhash(1023);
	regs = init_pbregv(1024);
	for(nidx=0;nidx<g->nodes->size;nidx++){
		if(is_dead_node_strgraph(g, nidx)) continue;
		n = ref_sgnodev(g->nodes, nidx);
		if(n->edge_cnts[0] == 0 || n->edge_cnts[1] == 0) continue;
		clear_pbregv(regs);
		for(ndir=0,f=1;f&&ndir<2;ndir++){
			for(eidx=0;f&&eidx<n->edge_cnts[ndir];eidx++){
				e = ref_sgedgev(g->edges, n->edge_offs[ndir] + eidx);
				ol = edge_overlap_strgraph(g, nidx, ndir, eidx);
				if(ol >= (uint32_t)(0.95 * g->rdlens->buffer[nidx])){ f = 0; break; }
				r1 = next_ref_pbregv(regs);
				r1->node_id = e->node_id;
				r1->x = e->off;
				r1->y = e->off + ol;
				if(ndir){
					r1->x = r1->x ^ r1->y; r1->y = r1->x ^ r1->y; r1->x = r1->x ^ r1->y;
					r1->x = ((int)g->rdlens->buffer[nidx]) - r1->x;
					r1->y = ((int)g->rdlens->buffer[nidx]) - r1->y;
				}
			}
		}
		if(!f) continue;
		sort_array(regs->buffer, regs->size, pb_reg_t, a.x > b.x);
		clear_uuhash(hash);
		for(i=0;i<regs->size;i++){
			r1 = ref_pbregv(regs, i);
			kv_put_uuhash(hash, r1->node_id, i);
		}
		for(i=0;i<regs->size;i++){
			r1 = ref_pbregv(regs, i);
			n1 = ref_sgnodev(g->nodes, r1->node_id);
			for(k=0;k<2;k++){
				for(j=0;j<n1->edge_cnts[k];j++){
					e = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
					if((rid = kv_get_uuhash(hash, e->node_id)) == 0xFFFFFFFFU) continue;
					r2 = ref_pbregv(regs, rid);
					r1->x = r2->x = (r1->x < r2->x)? r1->x : r2->x;
					r1->y = r2->y = (r1->y > r2->y)? r1->y : r2->y;
				}
			}
		}
		max = 0;
		for(i=0;i<regs->size;i++){
			r1 = ref_pbregv(regs, i);
			len = r1->y - r1->x;
			if(len > max) max = len;
		}
		if(max  >= (uint32_t)(0.95 * g->rdlens->buffer[nidx])) continue;
		fprintf(log, ">%s\t\t%u\t%d/%d\n", (char*)get_cplist(g->rdnames, nidx), nidx, max, g->rdlens->buffer[nidx]);
		r1 = next_ref_pbregv(regs);
		for(ndir=0;ndir<2;ndir++){
			for(eidx=0;eidx<n->edge_cnts[ndir];eidx++){
				e = ref_sgedgev(g->edges, n->edge_offs[ndir] + eidx);
				ol = edge_overlap_strgraph(g, nidx, ndir, eidx);
				r1->node_id = e->node_id;
				r1->x = e->off;
				r1->y = e->off + ol;
				if(ndir){
					r1->x = r1->x ^ r1->y; r1->y = r1->x ^ r1->y; r1->x = r1->x ^ r1->y;
					r1->x = ((int)g->rdlens->buffer[nidx]) - r1->x;
					r1->y = ((int)g->rdlens->buffer[nidx]) - r1->y;
				}
				fprintf(log, "%s\t%u\t%c\t%c\t%d\t%d\t%d\n", (char*)get_cplist(g->rdnames, e->node_id), e->node_id, "+-"[ndir], "+-"[e->dir], g->rdlens->buffer[e->node_id], r1->x, r1->y);
			}
		}
		mask_node_strgraph(g, nidx);
		ret ++;
	}
	return ret;
}

uint64_t find_bi_best_overlap_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_edge_t *e, *p;
	uint64_t ret;
	uint32_t i, k, j, b;
	int best, off;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			best = 0; b = SG_MAX_EDGE;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(e->score > best){ best = e->score; }
			}
			if(best == 0) continue;
			off = 1000000;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(e->score >= 0.75 * best){
					if(off > e->off){ off = e->off; b = j; }
				}
			}
			for(j=0;j<n->edge_cnts[k];j++){
				if(j == b) continue;
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				e->closed = 1;
				ret ++;
			}
		}
	}
	if(0){
		for(i=0;i<g->nodes->size;i++){
			if(is_dead_node_strgraph(g, i)) continue;
			n = ref_sgnodev(g->nodes, i);
			for(k=0;k<2;k++){
				for(j=0;j<n->edge_cnts[k];j++){
					e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
					if(e->closed) continue;
					if(!e->mark) continue;
					p = edge_partner_strgraph(g, i, k, j);
					if(!p->mark) continue;
					cut_biedge_strgraph2(g, i, k, j);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint64_t better_overlap_strgraph(StringGraph *g, float var){
	sg_node_t *n;
	sg_edge_t *e, *p;
	uint64_t ret;
	uint32_t i, k, j, b;
	float best, cutoff, score;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			score = b = 0;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				b ++;
				score += e->score;
			}
			if(b < 2) continue;
			cutoff = score / b;
			best = 0; b = 0;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(e->score < cutoff) continue;
				e->mark = 0;
				score = (1.0 * e->score) / edge_overlap_strgraph(g, i, k, j);
				//score = e->score;
				if(score > best){ best = score; b = j; }
			}
			if(best == 0) continue;
			for(j=0;j<n->edge_cnts[k];j++){
				if(j == b) continue;
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				score = (1.0 * e->score) / edge_overlap_strgraph(g, i, k, j);
				//score = e->score;
				if(score >= (1 - var) * best) continue;
				if(0){
					//cut_biedge_strgraph2(g, i, k, j);
					cut_biedge_strgraph(g, i, k, j);
					ret ++;
				} else {
					e->mark = 1;
				}
			}
		}
	}
	if(1){
		for(i=0;i<g->nodes->size;i++){
			if(is_dead_node_strgraph(g, i)) continue;
			n = ref_sgnodev(g->nodes, i);
			for(k=0;k<2;k++){
				for(j=0;j<n->edge_cnts[k];j++){
					e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
					if(e->closed) continue;
					if(!e->mark) continue;
					p = edge_partner_strgraph(g, i, k, j);
					if(!p->mark) continue;
					cut_biedge_strgraph2(g, i, k, j);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint64_t cut_lonely_edges_strgraph(StringGraph *g){
	sg_node_t *n1, *n2;
	sg_edge_t *e1, *e2, *p;
	uihash *hash;
	uint64_t ret;
	uint32_t i, k, j, k2, j2, a, b;
	int off, ol, off2;
	hash = init_uihash(1023);
	for(ret=0;ret<g->edges->size;ret++) g->edges->buffer[ret].mark = 0;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n1 = ref_sgnodev(g->nodes, i);
		clear_uihash(hash);
		for(k=0;k<2;k++){
			for(j=0;j<n1->edge_cnts[k];j++){
				e1 = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e1->closed) continue;
				if(k){
					ol = edge_overlap_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
					if(ol == (int)g->rdlens->buffer[e1->node_id]) continue;
					p = edge_partner_strgraph(g, i, k, j);
					off = - (((int)g->rdlens->buffer[e1->node_id]) - p->off);
				} else {
					off = e1->off;
				}
				kv_put_uihash(hash, e1->node_id, off);
			}
		}
		for(k=0;k<2;k++){
			for(j=0;j<n1->edge_cnts[k];j++){
				e1 = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e1->closed) continue;
				if(k){
					ol = edge_overlap_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
					if(ol == (int)g->rdlens->buffer[e1->node_id]) continue;
					p = edge_partner_strgraph(g, i, k, j);
					off = - (((int)g->rdlens->buffer[e1->node_id]) - p->off);
				} else {
					off = e1->off;
				}
				n2 = ref_sgnodev(g->nodes, e1->node_id);
				a = b = 0;
				for(k2=0;k2<2;k2++){
					for(j2=0;j2<n2->edge_cnts[k2];j2++){
						e2 = ref_sgedgev(g->edges, n2->edge_offs[k2] + j2);
						if(e2->closed) continue;
						if(k2 ^ e1->dir){
							ol = edge_overlap_strgraph(g, e2->node_id, !e2->dir, e2->rev_idx);
							if(ol == (int)g->rdlens->buffer[e2->node_id]) continue;
							p = edge_partner_strgraph(g, e1->node_id, k2, j2);
							off2 = off - (((int)g->rdlens->buffer[e2->node_id]) - p->off);
						} else {
							off2 = off + e2->off;
						}
						if(off2 + 200 > (int)g->rdlens->buffer[i]) continue;
						if(off2 + g->rdlens->buffer[e2->node_id] < 200) continue;
						a ++;
						if(kv_exists_uihash(hash, e2->node_id)) b ++;
					}
				}
				if(b * 3 < a){
					e1->mark = 1;
				}
			}
		}
	}
	free_uihash(hash);
	if(1){
		for(i=0;i<g->nodes->size;i++){
			if(is_dead_node_strgraph(g, i)) continue;
			n1 = ref_sgnodev(g->nodes, i);
			for(k=0;k<2;k++){
				for(j=0;j<n1->edge_cnts[k];j++){
					e1 = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
					if(e1->closed) continue;
					if(!e1->mark) continue;
					p = edge_partner_strgraph(g, i, k, j);
					if(!p->mark) continue;
					cut_biedge_strgraph2(g, i, k, j);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint64_t cut_nasty_edges_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_edge_t *e, *p;
	f32list *dat;
	u32list *idxs;
	uint64_t ret;
	double  v, cutoff;
	uint32_t node_id, i, k, f, r;
	ret = 0;
	dat = init_f32list(1024);
	idxs = init_u32list(16);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(1){
			// cheating
			for(k=0;k<2;k++){
				for(i=0;i<n->edge_cnts[k];i++){
					e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
					if(e->closed) continue;
					int a, b, c, d;
					char *ptr;
					ptr = (char*)g->rdnames->buffer[node_id];
					a = atoi(ptr + 1);
					b = atoi(index(ptr, '_') + 1);
					ptr = (char*)g->rdnames->buffer[e->node_id];
					c = atoi(ptr + 1);
					d = atoi(index(ptr, '_') + 1);
					a = (a > c)? a : c;
					b = (b < d)? b : d;
					if(a > b){
						ret ++;
						cut_biedge_strgraph2(g, node_id, k, i);
						//fprintf(stderr, "%s\t%s\n", (char*)g->rdnames->buffer[node_id], (char*)g->rdnames->buffer[e->node_id]);
					}
				}
			}
			continue;
		}
		if(n->edge_cnts[0] + n->edge_cnts[1] < 3) continue;
		f = count_living_edge_strgraph(g, node_id, 0);
		r = count_living_edge_strgraph(g, node_id, 1);
		if(f < 2 && r < 2) continue;
		clear_f32list(dat);
		for(k=0;k<2;k++){
			for(i=0;i<n->edge_cnts[k];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
				v = ((double)e->score) / edge_overlap_strgraph(g, node_id, k, i);
				push_f32list(dat, v);
			}
		}
		sort_array(dat->buffer, dat->size, float, b > a);
		i = 0;
		while(i < dat->size && dat->buffer[i] >= 0.8) i ++;
		if(i < dat->size / 4) i = dat->size / 4;
		//i = dat->size / 3 + 1;
		for(;i<=dat->size;i++){
			if(_test_clustering(dat->buffer, i)) break;
		}
		cutoff = dat->buffer[i - 1];
		if(f > 1){
			k = 0;
			clear_u32list(idxs);
			for(i=0;i<n->edge_cnts[k];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
				if(e->closed) continue;
				v = ((double)e->score) / edge_overlap_strgraph(g, node_id, k, i);
				if(v < cutoff) push_u32list(idxs, i);
			}
			if(idxs->size < f){
				for(i=0;i<idxs->size;i++){
					e = ref_sgedgev(g->edges, n->edge_offs[k] + get_u32list(idxs, i));
					e->mark = 1;
					//ret ++;
					//cut_biedge_strgraph(g, node_id, k, get_u32list(idxs, i));
				}
			}
		}
		if(r > 1){
			k = 1;
			clear_u32list(idxs);
			for(i=0;i<n->edge_cnts[k];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
				if(e->closed) continue;
				v = ((double)e->score) / edge_overlap_strgraph(g, node_id, k, i);
				if(v < cutoff) push_u32list(idxs, i);
			}
			if(idxs->size < r){
				for(i=0;i<idxs->size;i++){
					e = ref_sgedgev(g->edges, n->edge_offs[k] + get_u32list(idxs, i));
					e->mark = 1;
					//ret ++;
					//cut_biedge_strgraph(g, node_id, k, get_u32list(idxs, i));
				}
			}
		}
	}
	free_u32list(idxs);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		for(k=0;k<2;k++){
			for(i=0;i<n->edge_cnts[k];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
				if(e->closed) continue;
				if(e->mark){
					p = edge_partner_strgraph(g, node_id, k, i);
					if(p->mark){
						e->closed = 3;
						p->closed = 3;
						//fprintf(stderr, "%s\t%s\n", (char*)g->rdnames->buffer[node_id], (char*)g->rdnames->buffer[e->node_id]);
						ret ++;
					}
				}
			}
		}
	}
	return ret;
}

uint32_t mask_self_circle_reads_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_edge_t *e;
	u32hash *hash;
	uint32_t i, j, c, len, ret;
	ret = 0;
	hash = init_u32hash(31);
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		len = g->rdlens->buffer[i];
		clear_u32hash(hash);
		for(j=0;j<n->edge_cnts[0];j++){
			e = ref_sgedgev(g->edges, n->edge_offs[0] + j);
			if(e->closed) continue;
			if(edge_overlap_strgraph(g, i, 0, j) < len / 3) continue;
			put_u32hash(hash, e->node_id);
		}
		c = 0;
		for(j=0;j<n->edge_cnts[1];j++){
			e = ref_sgedgev(g->edges, n->edge_offs[1] + j);
			if(e->closed) continue;
			if(edge_overlap_strgraph(g, i, 1, j) < len / 3) continue;
			if(exists_u32hash(hash, e->node_id)){ c = 1; break; }
		}
		if(c == 0) continue;
		mask_node_strgraph(g, i);
		ret ++;
	}
	free_u32hash(hash);
	return ret;
}

uint64_t reduce_transitive_strgraph(StringGraph *g){
	sg_node_t *n, *n2;
	sg_edge_t *e, *e2;
	uuhash *hash;
	u32list *vec;
	uint64_t ret;
	uint32_t i, j, k, d, eidx, cnt;
	cnt = 0;
	uint32_t wushigang = cnt;
	uint32_t tmp = wushigang;
	wushigang = tmp;
	int dir;
	hash = init_uuhash(0x1fU);
	ret = 0;
	vec = init_u32list(64);
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(d=0;d<2;d++){
			clear_uuhash(hash);
			clear_u32list(vec);
			cnt = 0;
			for(j=0;j<n->edge_cnts[d];j++){
				push_u32list(vec, j);
			}
			sort_array(vec->buffer, vec->size, uint32_t, ref_sgedgev(g->edges, n->edge_offs[d] + b)->off > ref_sgedgev(g->edges, n->edge_offs[d] + a)->off);
			for(j=0;j<vec->size;j++){
				e = ref_sgedgev(g->edges, n->edge_offs[d] + vec->buffer[j]);
				if(e->closed == 1) continue;
				kv_put_uuhash(hash, e->node_id, j);
			}
			for(j=0;j+1<vec->size;j++){
				e = ref_sgedgev(g->edges, n->edge_offs[d] + vec->buffer[j]);
				if(e->closed) continue;
				dir = !e->dir;
				n2 = ref_sgnodev(g->nodes, e->node_id);
				for(k=0;k<n2->edge_cnts[dir];k++){
					e2 = ref_sgedgev(g->edges, n2->edge_offs[dir] + k);
					if(e2->closed == 1) continue;
					eidx = kv_get_uuhash(hash, e2->node_id);
					if(eidx == 0xFFFFFFFFU) continue;
					if(eidx <= j) continue;
					cut_biedge_strgraph2(g, i, d, vec->buffer[j]);
					ret ++;
					break;
				}
			}
		}
	}
	free_uuhash(hash);
	free_u32list(vec);
	return ret;
}

uint64_t enhanced_reduce_transitive_strgraph(StringGraph *g){
	sg_node_t *n, *n2;
	sg_edge_t *e, *e2;
	u64list *eids;
	uint64_t ret;
	uint32_t k, d, i, j, l, m, c;
	ret = 0;
	eids = init_u64list(64);
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			if(n->edge_cnts[k] < 2) continue;
			clear_u64list(eids);
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				push_u64list(eids, n->edge_offs[k] + j);
			}
			if(eids->size < 2) continue;
			sort_array(eids->buffer, eids->size, uint64_t, g->edges->buffer[b].off > g->edges->buffer[a].off);
			for(m=0;m+1<eids->size;m++){
				c = 1;
				for(j=0;c&&j<n->edge_cnts[k];j++){
					if(eids->buffer[m] == n->edge_offs[k] + j) continue;
					e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
					if(e->closed) continue;
					if(is_dead_node_strgraph(g, e->node_id)) continue;
					n2 = ref_sgnodev(g->nodes, e->node_id);
					d = e->dir;
					{
						for(l=0;c&&l<n2->edge_cnts[d];l++){
							e2 = ref_sgedgev(g->edges, n2->edge_offs[d] + l);
							if(e2->closed == 2 || e2->closed == 3) continue;
							if(e2->node_id == ref_sgedgev(g->edges, eids->buffer[m])->node_id){
								cut_biedge_strgraph2(g, i, k, eids->buffer[m] - n->edge_offs[k]);
								ret ++;
								c = 0;
								break;
							}
						}
					}
				}
			}
		}
	}
	free_u64list(eids);
	return ret;
}

typedef struct {
	uint32_t node_id[2];
	uint32_t dir[2], eidx[2];
	uint64_t off;
	uint32_t step;
} sg_backtrace_t;

int traceback_strgraph(StringGraph *g, sg_backtrace_t bts[2], uint32_t max_step){
	sg_node_t *n1, *n2;
	sg_edge_t *e;
	int f;
	f = 1;
	while(1){
		if(bts[0].off == bts[1].off){
			if(bts[0].node_id[0] == bts[1].node_id[0]) return 1;
			else return 0;
		}
		if(bts[0].off < bts[1].off) f = 1; else f = 0;
		if(bts[f].step >= max_step) return 0;
		bts[f].step ++;
		bts[f].node_id[1] = bts[f].node_id[0];
		bts[f].dir[1] = bts[f].dir[0];
		bts[f].eidx[1] = bts[f].eidx[0];
		n1 = ref_sgnodev(g->nodes, bts[f].node_id[0]);
		e = ref_sgedgev(g->edges, n1->edge_offs[bts[f].dir[0]] + bts[f].eidx[0]);
		bts[f].node_id[0] = e->node_id;
		n2 = ref_sgnodev(g->nodes, e->node_id);
		bts[f].dir[0] = n2->bt_dir;
		bts[f].eidx[0] = n2->bt_idx;
		bts[f].off = n2->bt_off;
	}
}

int _cmp_node_backtrace_offset_strgraph(uint32_t id1, uint32_t id2, void *ref){
	sg_node_t *n1, *n2;
	n1 = ref_sgnodev(((StringGraph*)ref)->nodes, id1);
	n2 = ref_sgnodev(((StringGraph*)ref)->nodes, id2);
	if(n1->bt_off > n2->bt_off) return 1;
	if(n1->bt_off < n2->bt_off) return -1;
	return 0;
}

uint64_t merge_bubbles_strgraph(StringGraph *g, int max_step){
	Heap *heap;
	sg_node_t *n, *m;
	sg_edge_t *e;
	sg_backtrace_t bts[2];
	uint64_t bubble, loop;
	uint32_t node_id, idx, i, btid;
	heap = init_heap(_cmp_node_backtrace_offset_strgraph, g);
	for(node_id=0;node_id<g->nodes->size;node_id++) ref_sgnodev(g->nodes, node_id)->bt_visit = 0;
	bubble = loop = 0;
	btid = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		n = ref_sgnodev(g->nodes, node_id);
		if(n->bt_visit) continue;
		btid ++; if(btid > SG_NODE_BTID_MAX) btid = 1;
		clear_heap(heap);
		n->bt_visit = btid; n->bt_dir = 1; n->bt_off = 0; n->bt_idx = 0;
		push_heap(heap, node_id);
		while((idx = pop_heap(heap)) != 0xFFFFFFFFU){
			n = ref_sgnodev(g->nodes, idx);
			for(i=0;i<n->edge_cnts[!n->bt_dir];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[!n->bt_dir] + i);
				if(e->closed) continue;
				m = ref_sgnodev(g->nodes, e->node_id);
				if(m->bt_visit){
					if(m->bt_visit == btid){
						if((!e->dir) != m->bt_dir) continue;
						bts[0].node_id[0] = bts[0].node_id[1] = e->node_id;
						bts[0].dir[0] = bts[0].dir[1] = m->bt_dir;
						bts[0].eidx[0] = bts[0].eidx[1] = m->bt_idx;
						bts[0].off = m->bt_off;
						bts[0].step = 0;
						bts[1].node_id[0] = bts[1].node_id[1] = e->node_id;
						bts[1].dir[0] = bts[1].dir[1] = !e->dir;
						bts[1].eidx[0] = bts[1].eidx[1] = e->rev_idx;
						bts[1].off = n->bt_off + e->off;
						bts[1].step = 0;
						if(traceback_strgraph(g, bts, max_step) == 0) continue;
						if(bts[0].dir[0] != bts[1].dir[0]) continue;
						if(bts[0].step < bts[1].step){
							cut_biedge_strgraph2(g, e->node_id, m->bt_dir, m->bt_idx);
							cut_biedge_strgraph2(g, bts[0].node_id[1], bts[0].dir[1], bts[0].eidx[1]);
						} else {
							cut_biedge_strgraph2(g, e->node_id, !e->dir, e->rev_idx);
							cut_biedge_strgraph2(g, bts[1].node_id[1], bts[1].dir[1], bts[1].eidx[1]);
						}
						bubble ++;
					}
					continue;
				}
				m->bt_visit = btid; m->bt_dir = !e->dir; m->bt_off = n->bt_off + e->off; m->bt_idx = e->rev_idx;
				push_heap(heap, e->node_id);
			}
		}
	}
	free_heap(heap);
	return bubble + loop;
}

uint64_t best_score_overlap_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_edge_t *e, *p;
	uint64_t ret;
	uint32_t i, k, j, b;
	float best, score;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			best = 0; b = 0;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				//score = (1.0 * e->score) / edge_overlap_strgraph(g, i, k, j);
				score = e->score;
				if(score > best){ best = score; b = j; }
			}
			if(best == 0) continue;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(j == b) e->mark = 0;
				else e->mark = 1;
			}
		}
	}
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(!e->mark) continue;
				p = edge_partner_strgraph(g, i, k, j);
				if(!p->mark) continue;
				cut_biedge_strgraph(g, i, k, j);
				ret ++;
			}
		}
	}
	return ret;
}

uint64_t longest_overlap_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_node_t *m;
	sg_edge_t *e, *p;
	uint64_t ret;
	uint32_t i, k, j, b, c;
	int best;
	float score, bestS;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		for(k=0;k<2;k++){
			best = g->rdlens->buffer[i]; b = 0xFFFFFFFFU;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(e->off < best){ best = e->off; b = j; }
			}
			if(b == 0xFFFFFFFFU) continue;
			//best = best + 0.05 * g->rdlens->buffer[i];
			best = best + 50;
			bestS = 0;
			c = 0xFFFFFFFFU;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
				if(e->off > best) continue;
				score = ((float)e->score);
				if(score > bestS){ bestS = score; c = j; }
			}
			if(c != b){
				score = ((float)ref_sgedgev(g->edges, n->edge_offs[k] + b)->score);
				if(score < 0.95 * bestS){ b = c; }
			}
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed) continue;
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

uint64_t pack_bogs_flag_strgraph(sg_node_t *n){
	uint64_t flag;
	flag = 0;
	flag |= ((uint64_t)n->bogs[0][0][0]) << 56;
	flag |= ((uint64_t)n->bogs[0][0][1]) << 48;
	flag |= ((uint64_t)n->bogs[0][1][0]) << 40;
	flag |= ((uint64_t)n->bogs[0][1][1]) << 32;
	flag |= ((uint64_t)n->bogs[1][0][0]) << 24;
	flag |= ((uint64_t)n->bogs[1][0][1]) << 16;
	flag |= ((uint64_t)n->bogs[1][1][0]) << 8;
	flag |= ((uint64_t)n->bogs[1][1][1]) << 0;
	return flag;
}

void check_node_bogs_strgraph(StringGraph *g, uint32_t node_id){
	sg_node_t *n;
	if(is_dead_node_strgraph(g, node_id)) return;
	n = ref_sgnodev(g->nodes, node_id);
	if(n->bogs[1][0][0] + n->bogs[1][0][1] > 1){
		fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	if(n->bogs[1][1][0] + n->bogs[1][1][1] > 1){
		fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
}

void check_all_nodes_bogs_strgraph(StringGraph *g){
	uint32_t i;
	for(i=0;i<g->nodes->size;i++) check_node_bogs_strgraph(g, i);
}

void cut_edge_bog_strgraph(StringGraph *g, sg_edge_t *e){
	sg_node_t *n1, *n2;
	sg_edge_t *p;
	if(e->closed) return;
	p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
	n1 = ref_sgnodev(g->nodes, p->node_id);
	n2 = ref_sgnodev(g->nodes, e->node_id);
	e->closed = 1;
	if(e->mark){
		n1->bogs[1][!p->dir][1] --;
		n2->bogs[0][e->dir][1] --;
	} else {
		p->mark = 1;
		n1->bogs[1][!p->dir][0] --;
		n1->bogs[0][p->dir][0] --;
		n1->bogs[0][p->dir][1] ++;
		n2->bogs[1][!e->dir][0] --;
		n2->bogs[1][!e->dir][1] ++;
		n2->bogs[0][e->dir][0] --;
	}
}

void repair_one_way_edge_bog_strgraph(StringGraph *g, sg_edge_t *e){
	sg_node_t *n1, *n2;
	sg_edge_t *p;
	if(e->closed) return;
	if(e->mark == 0) return;
	p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
	n1 = ref_sgnodev(g->nodes, p->node_id);
	n2 = ref_sgnodev(g->nodes, e->node_id);
	if(n2->bogs[1][!e->dir][0] || n2->bogs[1][!e->dir][1]){
		fprintf(stderr, " *WARNING* -- try to repair the backward edge of %s-[%c][%u]->%s, but it will break BOG, give up in %s -- %s:%d --\n", 
			g->rdnames->buffer[p->node_id], "+-"[!p->dir], p->rev_idx, g->rdnames->buffer[e->node_id], __FUNCTION__, __FILE__, __LINE__); return;
	}
	p->closed = 0;
	e->mark = 0;
	p->mark = 0;
	n1->bogs[1][!p->dir][1] --;
	n1->bogs[1][!p->dir][0] ++;
	n1->bogs[0][p->dir][0] ++;
	n2->bogs[0][e->dir][1] --;
	n2->bogs[0][e->dir][0] ++;
	n2->bogs[1][!e->dir][0] ++;
}

int revive_edge_bog_strgraph(StringGraph *g, sg_edge_t *e){
	sg_node_t *n1, *n2;
	sg_edge_t *p;
	if(e->closed == 0) return 0;
	p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
	n1 = ref_sgnodev(g->nodes, p->node_id);
	n2 = ref_sgnodev(g->nodes, e->node_id);
	if(n1->bogs[1][!p->dir][0] + n1->bogs[1][!p->dir][1]) return 0;
	if(p->closed){
		e->closed = 0;
		e->mark = 1;
		n1->bogs[1][!p->dir][1] ++;
		n2->bogs[0][e->dir][1] ++;
	} else if(n1->bogs[1][!p->dir][0] == 0){
		e->closed = 0;
		e->mark = 0;
		p->mark = 0;
		n1->bogs[0][p->dir][1] --;
		n1->bogs[0][p->dir][0] ++;
		n1->bogs[1][!p->dir][0] ++;

		n2->bogs[0][e->dir][0] ++;
		n2->bogs[1][!e->dir][1] --;
		n2->bogs[1][!e->dir][0] ++;
	} else return 0;
	return 1;
}

uint32_t count_linear_nodes_bog_strgraph(StringGraph *g, uint32_t node_id, int dir, uint32_t max_nodes){
	sg_node_t *n;
	sg_edge_t *e;
	uint32_t cnt;
	for(cnt=0;cnt<max_nodes;cnt++){
		n = ref_sgnodev(g->nodes, node_id);
		if(n->bogs[1][dir][0] == 0) break;
		if(n->bogs[0][!dir][1]) break;
		e = first_living_edge_strgraph(g, node_id, dir);
		node_id = e->node_id;
		dir = e->dir;
	}
	return cnt;
}

sg_edge_t* first_one_way_input_edge_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e1, *e2;
	uint32_t i;
	n = ref_sgnodev(g->nodes, node_id);
	for(i=0;i<n->edge_cnts[!dir];i++){
		e1 = ref_sgedgev(g->edges, n->edge_offs[!dir] + i);
		if(e1->closed != 1) continue;
		e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
		if(e2->closed) continue;
		return e2;
	}
	return NULL;
}

/**
 * Cut tip1
 * n1->n2
 * n1 has only e1 -> n2, no other output or input
 * cut e1, mask n1
 */
int cut_tip1_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n1;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	//sg_node_t *n2;
	sg_edge_t *e1;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, dir);
	cut_edge_bog_strgraph(g, e1);
	return 1;
}

/**
 * Cut tip4
 * n1<-n2
 * n1 has only input e1 from n2, no other output or input
 * cut e1, mask n1
 */
int cut_tip4_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n2;
	sg_edge_t *e1, *e2;
	e1 = first_one_way_input_edge_strgraph(g, node_id, dir);
	cut_edge_bog_strgraph(g, e1);
	e1 = edge_strgraph(g, node_id, !e1->dir, e1->rev_idx);
	n2 = ref_sgnodev(g->nodes, e1->node_id);
	if(n2->bogs[0][e1->dir][1] != 1) return 1;
	e2 = first_one_way_input_edge_strgraph(g, e1->node_id, e1->dir);
	repair_one_way_edge_bog_strgraph(g, e2);
	return 1;
}


/**
 * Cut tip2
 *  ---n1---
 *  $      $
 * =n2=...=n3=
 * mask n1, cut edges from n1
 */
int cut_tip2_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id){
	sg_node_t *n1;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	//sg_node_t *n2, *n3;
	sg_edge_t *e1, *e2;
	//uint64_t flag;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, 0);
	//n2 = ref_sgnodev(g->nodes, e1->node_id);
	//flag = pack_bogs_flag_strgraph(n2);
	//if((flag & 0x0100010001000100LLu) != 0x0100010001000100LLU) return 0;
	e2 = first_living_edge_strgraph(g, node_id, 1);
	//n3 = ref_sgnodev(g->nodes, e2->node_id);
	//flag = pack_bogs_flag_strgraph(n3);
	//if((flag & 0x0100010001000100LLu) != 0x0100010001000100LLU) return 0;
	cut_edge_bog_strgraph(g, e1);
	cut_edge_bog_strgraph(g, e2);
	return 1;
}

/**
 * Cut tip5 
 *  -->n1---
 *  |      $
 * -n2=...=n3=
 * mask n1, cut edges from n1
 */
int cut_tip5_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n1;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	sg_edge_t *e1, *e2;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, dir);
	e2 = first_one_way_input_edge_strgraph(g, node_id, !dir);
	cut_edge_bog_strgraph(g, e1);
	cut_edge_bog_strgraph(g, e2);
	return 1;
}

/**
 * Cut tip3
 * n1=n2-n3=
 *    ||    
 *    n4
 * mask n1, cut edges from n1
 */
int cut_tip3_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e1, *e2, *p;
	uint64_t flag;
	uint32_t step, max_step;
	max_step = 10;
	e1 = first_living_edge_strgraph(g, node_id, dir);
	step = 0;
	while(1){
		if(++step > max_step) return 0;
		n = ref_sgnodev(g->nodes, e1->node_id);
		flag = pack_bogs_flag_strgraph(n);
		if(flag == 0x0100010001000100LLU){
			e1 = first_living_edge_strgraph(g, e1->node_id, e1->dir);
		} else {
			if(e1->dir){
				if(flag != 0x0100010101000100LLU) return 0;
			} else {
				if(flag != 0x0101010001000100LLU) return 0;
			}
			break;
		}
	}
	e2 = first_one_way_input_edge_strgraph(g, e1->node_id, e1->dir);
	if(e2 == NULL) return 0;
	p = edge_strgraph(g, e2->node_id, !e1->dir, e1->rev_idx);
	cut_edge_bog_strgraph(g, e1);
	cut_edge_bog_strgraph(g, p);
	repair_one_way_edge_bog_strgraph(g, e2);
	return 1;
}

/**
 * Cut long tip
 *        n0
 *        || ...
 *        n1
 *        ||
 * =n=n4->n2=n3=
 * OR
 * n1=...x->n2
 */
int cut_tip6_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n, *n2;
	sg_edge_t *e, *p0, *p;
	uint32_t step, max_tip_nodes;
	max_tip_nodes = 10;
	step = 0;
	e = first_living_edge_strgraph(g, node_id, dir);
	while(1){
		if(++step > max_tip_nodes) return 0;
		if(e->mark == 1){ cut_edge_bog_strgraph(g, e); return 1; }
		n = ref_sgnodev(g->nodes, e->node_id);
		if(n->bogs[0][e->dir][1] == 1){
			p0 = first_one_way_input_edge_strgraph(g, e->node_id, e->dir);
			if(p0 == NULL) return 0;
			p = edge_strgraph(g, p0->node_id, !p0->dir, p0->rev_idx);
			if(count_linear_nodes_bog_strgraph(g, p->node_id, p->dir, max_tip_nodes) < max_tip_nodes) return 0;
			cut_edge_bog_strgraph(g, e);
			p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
			cut_edge_bog_strgraph(g, p);
			repair_one_way_edge_bog_strgraph(g, p0);
			return 1;
		}
		if(n->bogs[0][!e->dir][1] == 1){
			if(n->bogs[0][!e->dir][0]) return 0;
			p0 = first_one_way_input_edge_strgraph(g, e->node_id, !e->dir);
			if(p0 == NULL) return 0;
			p = edge_strgraph(g, p0->node_id, !p0->dir, p0->rev_idx);
			n2 = ref_sgnodev(g->nodes, p->node_id);
			if(n2->bogs[0][p->dir][1] != 1) return 0;
			if(count_linear_nodes_bog_strgraph(g, p->node_id, p->dir, max_tip_nodes) < max_tip_nodes) return 0;
			cut_edge_bog_strgraph(g, p0);
			p0 = first_one_way_input_edge_strgraph(g, p->node_id, p->dir);
			repair_one_way_edge_bog_strgraph(g, p0);
			return 1;
		}
		e = first_living_edge_strgraph(g, e->node_id, e->dir);
		if(e == NULL) return 0;
	}
	return 0;
}

/**
 * Cut nail
 *  n1=..=n2
 *  |     |
 *  $     $
 */
int cut_nail_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_edge_t *e, *e1, *e2;
	uint64_t flag;
	uint32_t step, max_tip_nodes;
	max_tip_nodes = 5;
	step = 0;
	e = first_living_edge_strgraph(g, node_id, !dir);
	while(1){
		if(++step > max_tip_nodes) return 0;
		if(e->mark == 1) break;
		flag = pack_bogs_flag_strgraph(ref_sgnodev(g->nodes, e->node_id));
		if(flag & 0x00FF00FF00000000LLU) return 0;
		e = first_living_edge_strgraph(g, e->node_id, e->dir);
		if(e == NULL) return 0;
	}
	e2 = e;
	e1 = first_living_edge_strgraph(g, node_id, dir);
	cut_edge_bog_strgraph(g, e1);
	cut_edge_bog_strgraph(g, e2);
	return 1;
}

/**
 * long jump
 *     ||
 *  ---n1<--
 *  |      |
 *  |      n3=
 *  $
 * =n2=
 *
 */
int repair_jump_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_edge_t *e1, *e2, *p;
	if(count_linear_nodes_bog_strgraph(g, node_id, !dir, 4) < 4) return 0;
	e1 = first_living_edge_strgraph(g, node_id, dir);
	if(count_linear_nodes_bog_strgraph(g, e1->node_id, 0, 4) < 4) return 0;
	if(count_linear_nodes_bog_strgraph(g, e1->node_id, 1, 4) < 4) return 0;
	e2 = first_one_way_input_edge_strgraph(g, node_id, !dir);
	p = edge_strgraph(g, e2->node_id, !e2->dir, e2->rev_idx);
	if(count_linear_nodes_bog_strgraph(g, p->node_id, p->dir, 4) < 4) return 0;
	cut_edge_bog_strgraph(g, e1);
	repair_one_way_edge_bog_strgraph(g, e2);
	return 1;
}

/**
 * long jump
 *  ?
 *  n1
 *  |
 * =n2=
 *
 */
int cut_nasty_jump_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e, *e1;
	float score, s;
	e1 = first_living_edge_strgraph(g, node_id, dir);
	if(count_linear_nodes_bog_strgraph(g, e1->node_id, 0, 4) < 4) return 0;
	if(count_linear_nodes_bog_strgraph(g, e1->node_id, 1, 4) < 4) return 0;
	score = ((float)e1->score) / edge_overlap_strgraph(g, node_id, dir, (e1 - g->edges->buffer) - g->nodes->buffer[node_id].edge_offs[dir]);
	n = ref_sgnodev(g->nodes, e1->node_id);
	e = first_living_edge_strgraph(g, e1->node_id, 0);
	s = ((float)e->score) / edge_overlap_strgraph(g, e1->node_id, 0, (e - g->edges->buffer) - n->edge_offs[0]);
	if(score >= s) return 0;
	e = first_living_edge_strgraph(g, e1->node_id, 1);
	s = ((float)e->score) / edge_overlap_strgraph(g, e1->node_id, 1, (e - g->edges->buffer) - n->edge_offs[1]);
	if(score >= s) return 0;
	cut_edge_bog_strgraph(g, e1);
	return 1;
}

/*
* n1 and n2 is connected by only n0, but n1 and n2 has their other connective paths, n0 is chimeric
*/
int mask_chimeric_node_best_overlap_strgraph(StringGraph *g, uint32_t node_id){
	sg_node_t *n, *n1, *n2;
	sg_edge_t *e, *e1, *e2;
	uint32_t i, k;
	n = ref_sgnodev(g->nodes, node_id);
	if(n->bogs[1][0][0] + n->bogs[1][0][1] != 1) return 0; // n0 must has forward edge to n1
	if(n->bogs[1][1][0] + n->bogs[1][1][1] != 1) return 0; // n0 must has backward edge to n2
	e1 = first_living_edge_strgraph(g, node_id, 0);
	e2 = first_living_edge_strgraph(g, node_id, 1);
	n1 = ref_sgnodev(g->nodes, e1->node_id);
	k = !e1->dir;
	for(i=0;i<n1->edge_cnts[k];i++){
		e = ref_sgedgev(g->edges, n1->edge_offs[k] + i);
		//if(e->node_id == e2->node_id && e->dir == e2->dir) return 0;
		if(e->node_id == e2->node_id) return 0; // n1 and n2 is conntective
	}
	//whether n1 has its alternative path
	if(n1->bogs[0][e1->dir][1] + n1->bogs[1][!e1->dir][1] + n1->bogs[1][!e1->dir][0] <= 1) return 0;
	//whether n2 has its alternative path
	n2 = ref_sgnodev(g->nodes, e2->node_id);
	if(n2->bogs[0][e2->dir][1] + n2->bogs[1][!e2->dir][1] + n2->bogs[1][!e2->dir][0] <= 1) return 0;
	// found chimera, mask it
	for(k=0;k<2;k++){
		for(i=0;i<n->edge_cnts[k];i++){
			e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
			cut_edge_bog_strgraph(g, e);
		}
	}
	one_bitvec(g->node_status, node_id);
	return 1;
}

int repair_lonely_one_way_edge_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n2;
	sg_edge_t *e1;
	e1 = first_living_edge_strgraph(g, node_id, dir);
	n2 = ref_sgnodev(g->nodes, e1->node_id);
	if(n2->bogs[1][!e1->dir][0] > 0 || n2->bogs[1][!e1->dir][1] > 0) return 0;
	repair_one_way_edge_bog_strgraph(g, e1);
	return 1;
}

uint64_t repair_all_lonely_one_way_edges_bog_strgraph(StringGraph *g){
	sg_node_t *n;
	uint64_t ret;
	uint32_t node_id;
	ret = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(n->bogs[1][0][0] == 0 && n->bogs[1][0][1] == 1) ret += repair_lonely_one_way_edge_best_overlap_strgraph(g, node_id, 0);
		if(n->bogs[1][1][0] == 0 && n->bogs[1][1][1] == 1) ret += repair_lonely_one_way_edge_best_overlap_strgraph(g, node_id, 1);
	}
	return ret;
}

/**
 * Bubble pattern:
 *      ||
 *    --n1<-
 *    |    |
 *    $    |
 *   -n2==n3-
 * 1: =n1=n2=n3+
 * 2: =n1=n3=n2+
 *
 */

int merge_bubble1_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n1, *n2, *n3;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	n3 = NULL;
	wushigang = n3;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	sg_edge_t *e1, *e2, *e3;
	uint64_t flag;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, dir);
	n2 = ref_sgnodev(g->nodes, e1->node_id);
	flag = pack_bogs_flag_strgraph(n2);
	if((flag & 0x0100010001000100LLU) == 0x0100010001000100LLU) return 0;
	e2 = first_living_edge_strgraph(g, e1->node_id, e1->dir);
	if(e2 && e2->closed == 0 && e2->mark == 0){
		n3 = ref_sgnodev(g->nodes, e2->node_id);
		e3 = first_living_edge_strgraph(g, e2->node_id, !e2->dir);
		if(e3 && e3->closed == 0 && e3->node_id == node_id && e3->dir != dir){
			cut_edge_bog_strgraph(g, e3);
			repair_one_way_edge_bog_strgraph(g, e1);
			return 1;
		} else return 0;
	}
	e2 = first_living_edge_strgraph(g, e1->node_id, !e1->dir);
	if(e2 && e2->closed == 0 && e2->mark == 0){
		n3 = ref_sgnodev(g->nodes, e2->node_id);
		e3 = first_living_edge_strgraph(g, e2->node_id, e2->dir);
		if(e3 && e3->closed == 0 && e3->node_id == node_id && e3->dir != dir){
			cut_edge_bog_strgraph(g, e1);
			repair_one_way_edge_bog_strgraph(g, e3);
			return 1;
		} else return 0;
	}
	return 0;
}

/**
 * Bubble pattern:
 * =n0=n1----->n3=
 *  $          ||
 *  -----------n2
 *
 * e0 : n1 -> n0
 * e0': n0 -> n1
 * e1 : n1 -> n3
 * e2 : n3 -> n2
 * e2': n2 -> n3
 * e3 : n2 -> n0
 *
 * cut e2+e2'+e3
 * add e1'
 */
int merge_bubble2_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n0, *n1, *n2, *n3;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	sg_edge_t *es[4], *p;
	uint64_t flag;
	n1 = ref_sgnodev(g->nodes, node_id);
	es[1] = first_living_edge_strgraph(g, node_id, dir);
	es[0] = first_living_edge_strgraph(g, node_id, !dir);
	n0 = ref_sgnodev(g->nodes, es[0]->node_id);
	flag = pack_bogs_flag_strgraph(n0);
	if(es[0]->dir){
		if(flag != 0x0100010101000100LLU) return 0;
	} else {
		if(flag != 0x0101010001000100LLU) return 0;
	}
	n3 = ref_sgnodev(g->nodes, es[1]->node_id);
	flag = pack_bogs_flag_strgraph(n3);
	if(es[1]->dir){
		if(flag != 0x0100010101000100LLU) return 0;
	} else {
		if(flag != 0x0101010001000100LLU) return 0;
	}
	es[2] = first_living_edge_strgraph(g, es[1]->node_id, !es[1]->dir);
	n2 = ref_sgnodev(g->nodes, es[2]->node_id);
	flag = pack_bogs_flag_strgraph(n2);
	if(es[2]->dir){
		if(flag != 0x0000010001000001LLU) return 0;
	} else {
		if(flag != 0x0100000000010100LLU) return 0;
	}
	es[3] = first_living_edge_strgraph(g, es[2]->node_id, es[2]->dir);
	if(es[3]->node_id != es[0]->node_id) return 0;
	cut_edge_bog_strgraph(g, es[2]);
	p = edge_strgraph(g, es[2]->node_id, !es[2]->dir, es[2]->rev_idx);
	cut_edge_bog_strgraph(g, p);
	repair_one_way_edge_bog_strgraph(g, es[1]);
	return 1;
}

/**
 * Bubble pattern:
 * =n0=x=n1----->n3=
 *  $           ||
 *  |           xx
 *  |           ||
 *  ------------n2
 */

int merge_bubble3_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n, *n1, *n3;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	n3 = NULL;
	wushigang = n3;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	sg_edge_t *e, *e0, *e1, *e2, *e3, *p;
	uint32_t cnt1, cnt2, step, max_step;
	max_step = MERGE_BUBBLE_MAX_STEP;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, dir);
	n3 = ref_sgnodev(g->nodes, e1->node_id);
	e2 = first_living_edge_strgraph(g, e1->node_id, !e1->dir);
	if(e2 == NULL) return 0;
	cnt1 = 1;
	cnt2 = 0;
	step = 0;
	e = e2;
	e3 = NULL;
	while(1){
		if(++step > max_step) return 0;
		if(e->mark){
			n = ref_sgnodev(g->nodes, e->node_id);
			if(n->bogs[1][!e->dir][0] + n->bogs[1][!e->dir][1]){
				e3 = e; break;
			}
		}
		p = e;
		e = first_living_edge_strgraph(g, e->node_id, e->dir);
		if(e == NULL){
			e = first_one_way_input_edge_strgraph(g, p->node_id, p->dir);
			if(e == NULL) return 0;
			e3 = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
			break;
		}
	}
	if(0 && e3 == e1){ // loop
		cut_edge_bog_strgraph(g, e1);
		e = e2;
		while(1){
			if(e->mark){
				n = ref_sgnodev(g->nodes, e->node_id);
				if(n->bogs[1][!e->dir][0] == 0 && n->bogs[1][!e->dir][1] == 0) repair_one_way_edge_bog_strgraph(g, e);
				else break;
			}
			e = first_living_edge_strgraph(g, e->node_id, e->dir);
			if(e == NULL) break;
		}
		return 1;
	}
	cnt2 += step;
	e = first_living_edge_strgraph(g, node_id, !dir);
	e0 = NULL;
	step = 0;
	while(1){
		if(++step > max_step) return 0;
		if(e->node_id == e3->node_id){ e0 = e;  break; }
		n = ref_sgnodev(g->nodes, e->node_id);
		if(n->bogs[1][e->dir][0] != 1) return 0;
		e = first_living_edge_strgraph(g, e->node_id, e->dir);
	}
	cnt1 += step;
	if(cnt1 >= cnt2){
		cut_edge_bog_strgraph(g, e3);
		if(e2 != e3){
			p = edge_strgraph(g, e2->node_id, !e2->dir, e2->rev_idx);
			cut_edge_bog_strgraph(g, e2);
			cut_edge_bog_strgraph(g, p);
		}
		e = first_living_edge_strgraph(g, node_id, !dir);
		while(0){
			if(e->mark){
				n = ref_sgnodev(g->nodes, e->node_id);
				if(n->bogs[1][!e->dir][0] == 0 && n->bogs[1][!e->dir][1] == 0) repair_one_way_edge_bog_strgraph(g, e);
			}
			if(e->node_id == e3->node_id){ break; }
			e = first_living_edge_strgraph(g, e->node_id, e->dir);
		}
		repair_one_way_edge_bog_strgraph(g, e1);
	} else {
		cut_edge_bog_strgraph(g, e1);
		p = edge_strgraph(g, e0->node_id, !e0->dir, e0->rev_idx);
		cut_edge_bog_strgraph(g, e0);
		cut_edge_bog_strgraph(g, p);
		e = e2;
		while(0){
			if(e->mark){
				n = ref_sgnodev(g->nodes, e->node_id);
				if(n->bogs[1][!e->dir][0] == 0 && n->bogs[1][!e->dir][1] == 0) repair_one_way_edge_bog_strgraph(g, e);
				else break;
			}
			e = first_living_edge_strgraph(g, e->node_id, e->dir);
		}
		repair_one_way_edge_bog_strgraph(g, e3);
	}
	return 1;
}

/**
* Merge bubble4
* -n0-n1-n2=
*  ||    ||
*  =======
*/
int merge_bubble4_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_edge_t *e1, *e2, *e3, *e4;
	uint32_t step, max_step;
	max_step = MERGE_BUBBLE_MAX_STEP;
	e1 = first_living_edge_strgraph(g, node_id, dir);
	e2 = first_one_way_input_edge_strgraph(g, node_id, dir);
	e4 = edge_strgraph(g, node_id, !e2->dir, e2->rev_idx);
	e3 = first_living_edge_strgraph(g, e1->node_id, !e1->dir);
	step = 0;
	while(1){
		if(++step > max_step) return 0;
		if(e3 == NULL) return 0;
		if(e3->node_id == e4->node_id) break;
		e3 = first_living_edge_strgraph(g, e3->node_id, e3->dir);
	}
	cut_edge_bog_strgraph(g, e1);
	cut_edge_bog_strgraph(g, e2);
	e3 = first_living_edge_strgraph(g, e1->node_id, !e1->dir);
	while(1){
		if(e3->mark && g->nodes->buffer[e3->node_id].bogs[1][!e3->dir][0] + g->nodes->buffer[e3->node_id].bogs[1][!e3->dir][1] == 0){
			repair_one_way_edge_bog_strgraph(g, e3);
			break;
		}
		if(e3->node_id == e4->node_id) break;
		e3 = first_living_edge_strgraph(g, e3->node_id, e3->dir);
	}
	return 1;
}

/**
 * Loop pattern:
 * =n0=x=n1----->n3=
 *  $           ||
 *  |           xx
 *  |           ||
 *  ------------n2
 * direction: n1->n2->n3->n0->n1, a loop
 */

int cut_loop1_core_best_overlap_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n, *n1, *n3;
	n1 = NULL;
	sg_node_t *wushigang = n1;
	n3 = NULL;
	wushigang = n3;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	sg_edge_t *e, *e0, *e1, *e2, *e3, *p;
	uint32_t cnt1, cnt2, step, max_step;
	max_step = 10;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, dir);
	n3 = ref_sgnodev(g->nodes, e1->node_id);
	e2 = first_living_edge_strgraph(g, e1->node_id, e1->dir);
	if(e2 == NULL) return 0;
	cnt1 = 1;
	cnt2 = 0;
	step = 0;
	e = e2;
	e3 = NULL;
	while(1){
		if(++step > max_step) return 0;
		if(e->mark == 0){
			e = first_living_edge_strgraph(g, e->node_id, e->dir);
			if(e == NULL) return 0;
		} else {
			e3 = e; break;
		}
	}
	cnt2 += step;
	e = first_living_edge_strgraph(g, node_id, !dir);
	e0 = NULL;
	step = 0;
	while(1){
		if(++step > max_step) return 0;
		if(e->node_id == e3->node_id){ e0 = e;  break; }
		n = ref_sgnodev(g->nodes, e->node_id);
		if(n->bogs[1][e->dir][0] != 1) return 0;
		e = first_living_edge_strgraph(g, e->node_id, e->dir);
	}
	cnt1 += step;
	{
		cut_edge_bog_strgraph(g, e1);
		p = edge_strgraph(g, e0->node_id, !e0->dir, e0->rev_idx);
		cut_edge_bog_strgraph(g, e0);
		cut_edge_bog_strgraph(g, p);
		repair_one_way_edge_bog_strgraph(g, e3);
	}
	return 1;
}

uint64_t cut_tips_bog_strgraph(StringGraph *g, uint32_t max_step){
	sg_node_t *n;
	sg_edge_t *e;
	uint64_t ret;
	uint32_t i, k, dir, step, nid;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		for(k=0;k<2;k++){
			n = ref_sgnodev(g->nodes, i);
			if(n->bogs[1][!k][0] + n->bogs[1][!k][1]) continue;
			if(n->bogs[0][k][0] + n->bogs[0][k][1]) continue;
			if(n->bogs[1][k][1]){
				if(n->bogs[1][k][0]) continue;
				if(n->bogs[0][!k][1]) continue;
				e = first_living_edge_strgraph(g, i, k);
				cut_edge_bog_strgraph(g, e);
				ret ++;
				continue;
			}
			if(n->bogs[1][k][0] == 0) continue;
			if(n->bogs[0][!k][0] + n->bogs[0][!k][1] > 1) continue;
			step = 0;
			nid = i;
			dir = k;
			while(step ++ < max_step){
				e = first_living_edge_strgraph(g, nid, dir);
				if(e == NULL){
					break; // yes, ignore the below
					e = first_one_way_input_edge_strgraph(g, nid, dir);
					if(e == NULL) break;
					ret ++;
					cut_edge_bog_strgraph(g, e);
					break;
				}
				if(e->mark){
					ret ++;
					cut_edge_bog_strgraph(g, e);
					break;
				}
				nid = e->node_id;
				dir = e->dir;
				n = ref_sgnodev(g->nodes, nid);
				if(n->bogs[0][dir][1]) break;
				if(n->bogs[0][!dir][1]) break;
			}
		}
	}
	return ret;
}

typedef struct {
	uint32_t node_id;
	uint32_t dir:1, mark:1, eidx:10, off:20;
} trace_t;
define_list(tracev, trace_t);

int _cmp_heap_trace_t(uint32_t idx1, uint32_t idx2, void *ref){
	tracev *t;
	t = (tracev*)ref;
	cmp_2nums_proc(ref_tracev(t, idx1)->off, ref_tracev(t, idx2)->off);
	return 0;
}

int merge_bubble_core_best_overlap_strgraph(StringGraph *g, uuhash *hash, tracev *bts[2], uint32_t node_id, int dir){
	trace_t *ts[2], *t;
	sg_node_t *n1;
	sg_edge_t *e1, *e2, *e, *p;
	uint32_t k, step, max_step, dead, idx, found;
	max_step = MERGE_BUBBLE_MAX_STEP;
	n1 = ref_sgnodev(g->nodes, node_id);
	e1 = first_living_edge_strgraph(g, node_id, dir);
	e2 = first_one_way_input_edge_strgraph(g, node_id, !dir);
	if(e2 == NULL) return 0; // strange bubble, skip it, because I have no data to reproduce this. Please see https://github.com/ruanjue/smartdenovo/issues/10
	e2 = edge_strgraph(g, e2->node_id, !e2->dir, e2->rev_idx);
	clear_tracev(bts[0]);
	clear_tracev(bts[1]);
	clear_uuhash(hash);
	t = next_ref_tracev(bts[0]);
	t->node_id = node_id;
	t->dir = dir;
	t->eidx = e1 - edge_strgraph(g, node_id, dir, 0);
	t->mark = 0;
	t = next_ref_tracev(bts[0]);
	t->node_id = e1->node_id;
	t->dir = e1->dir;
	t->mark = 0;
	t->eidx = 0;
	kv_put_uuhash(hash, e1->node_id, ((bts[0]->size) << 1) | 0);
	t = next_ref_tracev(bts[1]);
	t->node_id = node_id;
	t->dir = dir;
	t->mark = 0;
	t->eidx = e2 - edge_strgraph(g, node_id, dir, 0);
	t = next_ref_tracev(bts[1]);
	t->node_id = e2->node_id;
	t->dir = e2->dir;
	t->mark = 0;
	t->eidx = 0;
	kv_put_uuhash(hash, e2->node_id, ((bts[1]->size) << 1) | 1);
	dead = 0;
	encap_tracev(bts[0], max_step);
	encap_tracev(bts[1], max_step);
	step = 0;
	found = 0;
	// try find bubble
	// path1: bts[0], forward travese from node_id
	// path2: bts[1], backward travese from node_id
	while(!found){
		step ++;
		if(step >= max_step) return 0;
		for(k=0;k<2;k++){
			ts[0] = peer_tracev(bts[k]);
			if(((dead >> k) & 0x01) == 1) continue;
			n1 = ref_sgnodev(g->nodes, ts[0]->node_id);
			if(n1->bogs[1][ts[0]->dir][0] || n1->bogs[1][ts[0]->dir][1]){
				e = first_living_edge_strgraph(g, ts[0]->node_id, ts[0]->dir);
			} else if(n1->bogs[0][!ts[0]->dir][1] == 1){
				e = first_one_way_input_edge_strgraph(g, ts[0]->node_id, !ts[0]->dir);
				e = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
			} else {
				dead |= 0x01 << k; if(dead == 0x03) return 0; else continue;
			}
			ts[0]->eidx = e - (g->edges->buffer + n1->edge_offs[ts[0]->dir]);
			ts[1] = next_ref_tracev(bts[k]);
			ts[1]->node_id = e->node_id;
			ts[1]->dir = e->dir;
			ts[1]->eidx = 0;
			if((idx = kv_get_uuhash(hash, e->node_id)) != 0xFFFFFFFFU){
				if((idx&0x01) == k) return 0;
				bts[idx&0x01]->size = idx >> 1;
				found = 1;
				break;
			}
			kv_put_uuhash(hash, e->node_id, ((bts[k]->size) << 1) | k);
		}
	}
	k = (bts[0]->size >= bts[1]->size);
	t = ref_tracev(bts[k], 0);
	e = edge_strgraph(g, t->node_id, t->dir, t->eidx);
	p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
	if(e->closed == 0) cut_edge_bog_strgraph(g, e);
	if(p->closed == 0) cut_edge_bog_strgraph(g, p);
	t = ref_tracev(bts[k], bts[k]->size - 2);
	e = edge_strgraph(g, t->node_id, t->dir, t->eidx);
	p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
	if(e->closed == 0) cut_edge_bog_strgraph(g, e);
	if(p->closed == 0) cut_edge_bog_strgraph(g, p);
	return 1;
}

uint64_t merge_bubbles_bog_strgraph(StringGraph *g){
	uuhash *hash;
	tracev *bts[2];
	sg_node_t *n;
	uint64_t ret;
	uint32_t node_id, k;
	hash = init_uuhash(1023);
	bts[0] = init_tracev(1024);
	bts[1] = init_tracev(1024);
	ret = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		for(k=0;k<2;k++){
			if(n->bogs[0][!k][1] == 0) continue;
			if(n->bogs[1][k][0] + n->bogs[1][k][1] != 1) continue;
			ret += merge_bubble_core_best_overlap_strgraph(g, hash, bts, node_id, k);
		}
	}
	free_uuhash(hash);
	free_tracev(bts[0]);
	free_tracev(bts[1]);
	return ret;
}

int cut_loop_core_bog_strgraph(StringGraph *g, uint32_t node_id, int dir, int max_step){
	sg_node_t *n;
	n = NULL;
	sg_node_t *wushigang = n;
	sg_node_t *tmp = wushigang;
	wushigang = tmp;
	sg_edge_t *e, *p;
	uint32_t nid;
	int step, k;
	nid = node_id;
	k = dir;
	for(step=0;step<max_step;step++){
		n = ref_sgnodev(g->nodes, nid);
		if((e = first_living_edge_strgraph(g, nid, k)) == NULL) return 0;
		if(e->node_id == node_id){
			cut_edge_bog_strgraph(g, e);
			if((p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx))) cut_edge_bog_strgraph(g, p);
			return 1;
		}
		nid = e->node_id;
		k = e->dir;
	}
	return 0;
}

uint64_t cut_loops_bog_strgraph(StringGraph *g){
	sg_node_t *n;
	uint64_t ret;
	uint32_t node_id;
	ret = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(n->bogs[0][0][0] + n->bogs[0][0][1] > 1) ret += cut_loop_core_bog_strgraph(g, node_id, 0, CUT_LOOP_MAX_STEP);
		if(n->bogs[0][1][0] + n->bogs[0][1][1] > 1) ret += cut_loop_core_bog_strgraph(g, node_id, 1, CUT_LOOP_MAX_STEP);
	}
	return ret;
}

uint64_t recover_paired_dead_ends_bog_strgraph(StringGraph *g){
	sg_node_t *n, *n2;
	sg_edge_t *e, *p;
	uuhash *hash;
	uuhash_t *h1, *h2;
	uint64_t flag, ret;
	uint32_t i, k, j, dir, step, min_nodes;
	min_nodes = 10;
	hash = init_uuhash(1023);
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		n = ref_sgnodev(g->nodes, i);
		flag = pack_bogs_flag_strgraph(n);
		if(flag == 0x0000010001000000LLU){
			if(count_linear_nodes_bog_strgraph(g, i, 0, min_nodes) < min_nodes) continue;
		} else if(flag == 0x0100000000000100LLU){
			if(count_linear_nodes_bog_strgraph(g, i, 1, min_nodes) < min_nodes) continue;
		} else continue;
		kv_put_uuhash(hash, i, 0);
	}
	reset_iter_uuhash(hash);
	while((h1 = ref_iter_uuhash(hash))){
		n = ref_sgnodev(g->nodes, h1->key);
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed != 1) continue;
				h2 = get_uuhash(hash, (uuhash_t){e->node_id, 0});
				if(h2 == NULL) continue;
				h1->val ++;
			}
		}
	}
	reset_iter_uuhash(hash);
	while((h1 = ref_iter_uuhash(hash))){
		if(h1->val == 1) continue;
		delete_uuhash(hash, h1);
	}
	reset_iter_uuhash(hash);
	while((h1 = ref_iter_uuhash(hash))){
		n = ref_sgnodev(g->nodes, h1->key);
		h1->val = 0xFFFFFFFFU;
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed != 1) continue;
				if(e->node_id < h1->key) continue;
				h2 = get_uuhash(hash, (uuhash_t){e->node_id, 0});
				if(h2 == NULL) continue;
				h1->val = e->node_id;
				break;
			}
		}
	}
	ret = 0;
	reset_iter_uuhash(hash);
	while((h1 = ref_iter_uuhash(hash))){
		if(h1->val == 0xFFFFFFFFU) continue;
		n = ref_sgnodev(g->nodes, h1->key);
		k = n->bogs[1][0][0];
		n2 = ref_sgnodev(g->nodes, h1->val);
		dir = !n2->bogs[1][0][0];
		step = 0;
		while(1){
			if(++step > min_nodes) break;
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(e->closed != 1) continue;
				if(e->node_id != h1->val) continue;
				if(e->dir != dir) break;
				n2 = ref_sgnodev(g->nodes, h1->val);
				if(n2->bogs[0][dir][0]){
					p = first_living_edge_strgraph(g, h1->val, !dir);
					cut_edge_bog_strgraph(g, p);
					p = edge_strgraph(g, p->node_id, !p->dir, p->rev_idx);
					cut_edge_bog_strgraph(g, p);
				}
				/*
				if(n->bogs[1][k][0]){
					fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
				}
				if(n->bogs[0][!k][0]){
					fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
				}
				if(n2->bogs[1][!e->dir][0]){
					fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
				}
				if(n2->bogs[0][e->dir][0]){
					fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
				}
				*/
				p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
				e->closed = 0;
				p->closed = 0;
				e->mark = 0;
				p->mark = 0;
				n->bogs[1][k][0] ++;
				n->bogs[0][!k][0] ++;
				n2->bogs[1][!e->dir][0] ++;
				n2->bogs[0][e->dir][0] ++;
				ret ++;
				h1->val = 0xFFFFFFFFU;
				break;
			}
			if(h1->val == 0xFFFFFFFFU) break;
			e = first_living_edge_strgraph(g, h1->val, dir);
			h1->val = e->node_id;
			dir = e->dir;
		}
	}
	free_uuhash(hash);
	return ret;
}

uint64_t repair_best_overlap_strgraph(StringGraph *g){
	uuhash *hash;
	tracev *bts[2];
	sg_node_t *n;
	unsigned long long tip, bub, single, rec, chi, flag;
	uint32_t node_id;
	tip = bub = single = rec = chi = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(n->bogs[1][0][1] && n->bogs[0][0][0] + n->bogs[0][0][1] == 0){
			cut_edge_bog_strgraph(g, first_living_edge_strgraph(g, node_id, 0));
			tip ++;
		} else if(n->bogs[1][1][1] && n->bogs[0][1][0] + n->bogs[0][1][1] == 0){
			cut_edge_bog_strgraph(g, first_living_edge_strgraph(g, node_id, 1));
			tip ++;
		}
	}
	if(1){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0000000000010000LLU: tip += cut_tip1_core_best_overlap_strgraph(g, node_id, 0); break;
				case 0x0000000000000001LLU: tip += cut_tip1_core_best_overlap_strgraph(g, node_id, 1); break;
			}
		}
	}
	if(1){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0001000000000000LLU: tip += cut_tip4_core_best_overlap_strgraph(g, node_id, 0); break;
				case 0x0000000100000000LLU: tip += cut_tip4_core_best_overlap_strgraph(g, node_id, 1); break;
			}
		}
	}
	if(1){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0000000000010001LLU: tip += cut_tip2_core_best_overlap_strgraph(g, node_id); break;
			}
		}
	}
	if(1){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0100000000010100LLU: tip += cut_nail_core_best_overlap_strgraph(g, node_id, 0); break;
				case 0x0000010001000001LLU: tip += cut_nail_core_best_overlap_strgraph(g, node_id, 1); break;
			}
		}
	}
	if(1){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0000010001000000LLU: tip += cut_tip6_core_best_overlap_strgraph(g, node_id, 0); break;
				case 0x0100000000000100LLU: tip += cut_tip6_core_best_overlap_strgraph(g, node_id, 1); break;
			}
		}
	}
	if(0){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0100000000010100LLU: bub += merge_bubble3_core_best_overlap_strgraph(g, node_id, 0); break;
				case 0x0000010001000001LLU: bub += merge_bubble3_core_best_overlap_strgraph(g, node_id, 1); break;
				case 0x0001000000010000LLU: bub += merge_bubble4_core_best_overlap_strgraph(g, node_id, 0); break;
				case 0x0000000100000001LLU: bub += merge_bubble4_core_best_overlap_strgraph(g, node_id, 1); break;
			}
		}
	}
	if(0){
		hash = init_uuhash(1023);
		bts[0] = init_tracev(1024);
		bts[1] = init_tracev(1024);
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			n = ref_sgnodev(g->nodes, node_id);
			flag = pack_bogs_flag_strgraph(n);
			switch(flag){
				case 0x0100000100010100LLU: bub += merge_bubble_core_best_overlap_strgraph(g, hash, bts, node_id, 0); break;
				case 0x0001010001000001LLU: bub += merge_bubble_core_best_overlap_strgraph(g, hash, bts, node_id, 1); break;
			}
		}
		free_uuhash(hash);
		free_tracev(bts[0]);
		free_tracev(bts[1]);
	}
	bub += merge_bubbles_bog_strgraph(g);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		flag = pack_bogs_flag_strgraph(n);
		switch(flag){
			case 0x0000010001000000LLU: if(cut_tip3_core_best_overlap_strgraph(g, node_id, 0)) tip ++; break;
			case 0x0100000000000100LLU: if(cut_tip3_core_best_overlap_strgraph(g, node_id, 1)) tip ++; break;
		}
	}
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		flag = pack_bogs_flag_strgraph(n);
		switch(flag){
			case 0x0100000100010100LLU: chi += repair_jump_core_best_overlap_strgraph(g, node_id, 0); break;
			case 0x0001010001000001LLU: chi += repair_jump_core_best_overlap_strgraph(g, node_id, 1); break;
		}
	}
	if(1){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			chi += mask_chimeric_node_best_overlap_strgraph(g, node_id);
		}
	}
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		flag = pack_bogs_flag_strgraph(n);
		switch(flag){
			case 0x0100000000010100LLU: chi += cut_nasty_jump_best_overlap_strgraph(g, node_id, 0); break;
			case 0x0000010001000001LLU: chi += cut_nasty_jump_best_overlap_strgraph(g, node_id, 1); break;
		}
	}
	bub += cut_loops_bog_strgraph(g);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		flag = pack_bogs_flag_strgraph(n);
		switch(flag){
			case 0x0100000000010100LLU: single += repair_lonely_one_way_edge_best_overlap_strgraph(g, node_id, 0); break;
			case 0x0000010001000001LLU: single += repair_lonely_one_way_edge_best_overlap_strgraph(g, node_id, 1); break;
		}
	}
	rec += recover_paired_dead_ends_bog_strgraph(g);
	fprintf(stdout, "%llu tips, %llu bubbles, %llu chimera, %llu non-bog, %llu recoveries\n", tip, bub, chi, single, rec);
	return tip + bub + single + rec;
}

uint32_t cut_read_tips_strgraph(StringGraph *g){
	sg_edge_t *e;
	uint32_t i, f, r, w, ret;
	int dir;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		if(is_dead_node_strgraph(g, i)) continue;
		// find read tip
		f = count_living_edge_strgraph(g, i, 0);
		if(f > 1) continue;
		r = count_living_edge_strgraph(g, i, 1);
		if(r + f != 1) continue;
		dir = (r == 1);
		e = first_living_edge_strgraph(g, i, dir);
		w = count_living_edge_strgraph(g, e->node_id, !e->dir);
		if(w < 2) continue; // not a branch
		// found read tip
		// mask node
		mask_node_strgraph(g, i);
		ret ++;
	}
	return ret;
}

sg_edge_t* bog_boldly_walk_strgraph(StringGraph *g, uint32_t node_id, int dir){
	sg_node_t *n;
	sg_edge_t *e, *p;
	uint32_t i;
	n = ref_sgnodev(g->nodes, node_id);
	if(n->bogs[1][dir][0] + (n->bogs[1][dir][1] + n->bogs[0][!dir][1]) != 1) return NULL;
	for(i=0;i<n->edge_cnts[dir];i++){
		e = ref_sgedgev(g->edges, n->edge_offs[dir] + i);
		if(e->closed == 0) return e;
		p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
		if(p->closed == 0) return e; // if partner p is living, return dead e
	}
	return NULL;
}

uint64_t bog_cut_tips_strgraph(StringGraph *g, uint32_t _max_step){
	sg_node_t *n, *n2;
	sg_edge_t *e, *p, *t;
	uint32_t node_id, nid, dir, d, ret, step, max_step;
	ret = 0;
	for(max_step=1;max_step<=_max_step;max_step++){
		for(node_id=0;node_id<g->nodes->size;node_id++){
			if(is_dead_node_strgraph(g, node_id)) continue;
			for(dir=0;dir<2;dir++){
				n = ref_sgnodev(g->nodes, node_id);
				if(n->bogs[0][dir][0] + n->bogs[0][dir][1] + n->bogs[1][!dir][1]) continue;
				nid = node_id;
				d = dir;
				for(step=0;step<max_step;step++){
					e = bog_boldly_walk_strgraph(g, nid, d);
					if(e == NULL) break;
					n2 = ref_sgnodev(g->nodes, e->node_id);
					if(n2->bogs[0][e->dir][0] + n2->bogs[0][e->dir][1] + n2->bogs[1][!e->dir][1] != 1){
						ret ++;
						p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
						if(e->closed == 0) cut_edge_bog_strgraph(g, e);
						if(p->closed == 0){
							cut_edge_bog_strgraph(g, p);
							if(n2->bogs[e->dir][0] == 0 && n2->bogs[0][e->dir][1] == 1){
								t = first_one_way_input_edge_strgraph(g, e->node_id, e->dir);
								repair_one_way_edge_bog_strgraph(g, t);
							}
						}
						break;
					}
					nid = e->node_id;
					d   = e->dir;
				}
			}
		}
	}
	repair_all_lonely_one_way_edges_bog_strgraph(g);
	return ret;
}

int bog_traceback_strgraph(StringGraph *g, sg_backtrace_t bts[2], uint32_t max_step){
	sg_node_t *n1, *n2;
	sg_edge_t *e, *p;
	int f;
	f = 1;
	while(1){
		if(bts[0].off == bts[1].off){
			if(bts[0].node_id[0] == bts[1].node_id[0]) return 1;
		}
		//if(bts[0].step >= max_step && bts[1].step >= max_step) return 0;
		if(bts[0].off < bts[1].off) f = 1; else f = 0;
		if(bts[f].step >= max_step) return 0;
		bts[f].step ++;
		bts[f].node_id[1] = bts[f].node_id[0];
		bts[f].dir[1] = bts[f].dir[0];
		bts[f].eidx[1] = bts[f].eidx[0];
		n1 = ref_sgnodev(g->nodes, bts[f].node_id[0]);
		e = ref_sgedgev(g->edges, n1->edge_offs[bts[f].dir[0]] + bts[f].eidx[0]);
		if(e->closed){
			p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
			if(p->closed) return 0;
		}
		bts[f].node_id[0] = e->node_id;
		n2 = ref_sgnodev(g->nodes, e->node_id);
		bts[f].dir[0] = n2->bt_dir;
		bts[f].eidx[0] = n2->bt_idx;
		bts[f].off = n2->bt_off;
	}
}

int _cmp_node_bog_orders(StringGraph *g, uint32_t a, uint32_t b){
	sg_node_t *n1, *n2;
	int s1, s2;
	n1 = g->nodes->buffer + a;
	n2 = g->nodes->buffer + b;
	s1 = n1->bogs[0][0][0] + n1->bogs[0][1][0] + n1->bogs[1][0][0] + n1->bogs[1][1][0];
	s2 = n2->bogs[0][0][0] + n2->bogs[0][1][0] + n2->bogs[1][0][0] + n2->bogs[1][1][0];
	return (s1 > s2);
}

uint64_t bog_merge_bubbles_strgraph(StringGraph *g, int _max_step, int with_loop){
	Heap *heap;
	sg_node_t *n, *m;
	sg_edge_t *e, *p, *e1, *e2;
	sg_backtrace_t bts[2];
	u32list *orders, *ranks;
	uint64_t bubble, loop;
	uint32_t node_id, idx, i, btid, rank;
	int max_step, bub, fwd, rev;
	heap = init_heap(_cmp_node_backtrace_offset_strgraph, g);
	orders = init_u32list(g->nodes->size);
	ranks  = init_u32list(g->nodes->size);
	bubble = loop = 0;
	for(max_step=2;max_step<=_max_step;max_step++){
		while(1){
			clear_u32list(orders);
			for(node_id=0;node_id<g->nodes->size;node_id++){
				if(is_dead_node_strgraph(g, node_id)) continue;
				n = ref_sgnodev(g->nodes, node_id);
				n->bt_visit = 0;
				if(n->bogs[1][0][0] && n->bogs[0][1][1] == 0){
					if(n->bogs[1][1][0] && n->bogs[0][0][1] == 0) continue;
					else rank = (count_linear_nodes_bog_strgraph(g, node_id, 0, _max_step) + 1) * ((n->bogs[1][1][0] + n->bogs[1][1][1] + n->bogs[0][0][1]));
				} else if(n->bogs[1][1][0] && n->bogs[0][0][1] == 0){
					rank = (count_linear_nodes_bog_strgraph(g, node_id, 1, _max_step) + 1) * ((n->bogs[1][0][0] + n->bogs[1][0][1] + n->bogs[0][1][1]));
				} else if(n->bogs[1][0][0] + n->bogs[1][0][1] + n->bogs[0][1][1] && n->bogs[1][1][0] + n->bogs[1][1][1] + n->bogs[0][0][1] == 0){
					rank = 1000;
				} else if(n->bogs[1][0][0] + n->bogs[1][0][1] + n->bogs[0][1][1] == 0 && n->bogs[1][1][0] + n->bogs[1][1][1] + n->bogs[0][0][1]){
					rank = 1000;
				} else continue;
				//rank  = count_linear_nodes_bog_strgraph(g, node_id, 0, _max_step * 2);
				//rank += count_linear_nodes_bog_strgraph(g, node_id, 1, _max_step * 2);
				push_u32list(orders, node_id);
				push_u32list(ranks, rank);
			}
			sort_array(orders->buffer, orders->size, uint32_t, ranks->buffer[a] > ranks->buffer[b]);
			btid = 0;
			bub = 0;
			while(pop_u32list(orders, &node_id)){
			//for(node_id=0;node_id<g->nodes->size;node_id++){
				if(is_dead_node_strgraph(g, node_id)) continue;
				n = ref_sgnodev(g->nodes, node_id);
				if(n->bt_visit) continue;
				if(n->bogs[1][0][0]){
					rev = n->bogs[0][1][1]? 1 : 2;
				} else if(n->bogs[1][0][0] + n->bogs[1][0][1] + n->bogs[0][1][1] == 0){
					rev = 4;
				} else rev = 0;
				if(n->bogs[1][1][0]){
					fwd = n->bogs[0][0][1]? 1 : 2;
				} else if(n->bogs[1][1][0] + n->bogs[1][1][1] + n->bogs[0][0][1] == 0){
					fwd = 4;
				} else fwd = 0;
				if(fwd < 2 && rev < 2) continue;
				n->bt_dir = (fwd < rev);
				btid ++; if(btid > SG_NODE_BTID_MAX) btid = 1;
				clear_heap(heap);
				n->bt_visit = btid; n->bt_off = 0; n->bt_idx = 0;
				push_heap(heap, node_id);
				while((idx = pop_heap(heap)) != 0xFFFFFFFFU){
					n = ref_sgnodev(g->nodes, idx);
					for(i=0;i<n->edge_cnts[!n->bt_dir];i++){
						e = ref_sgedgev(g->edges, n->edge_offs[!n->bt_dir] + i);
						if(e->closed){
							p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
							if(p->closed) continue;
						}
						m = ref_sgnodev(g->nodes, e->node_id);
						if(m->bt_visit){
							if(m->bt_visit == btid){
								//if((!e->dir) != m->bt_dir) continue;
								// check whether the traceback of m is valid
								{
									e1 = edge_strgraph(g, e->node_id, m->bt_dir, m->bt_idx);
									e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
									if(e1->closed && e2->closed) continue;
								}
								bts[0].node_id[0] = bts[0].node_id[1] = e->node_id;
								bts[0].dir[0] = bts[0].dir[1] = m->bt_dir;
								bts[0].eidx[0] = bts[0].eidx[1] = m->bt_idx;
								bts[0].off = m->bt_off;
								bts[0].step = 0;
								bts[1].node_id[0] = bts[1].node_id[1] = e->node_id;
								bts[1].dir[0] = bts[1].dir[1] = !e->dir;
								bts[1].eidx[0] = bts[1].eidx[1] = e->rev_idx;
								bts[1].off = n->bt_off + e->off;
								bts[1].step = 0;
								if(traceback_strgraph(g, bts, max_step) == 0) continue;
								//fprintf(stderr, " -- %s in %s -- %s:%d --\n", g->rdnames->buffer[idx], __FUNCTION__, __FILE__, __LINE__);
								if((!e->dir) == m->bt_dir){
									if(bts[0].step < bts[1].step){
										e1 = edge_strgraph(g, e->node_id, m->bt_dir, m->bt_idx);
										e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
										cut_edge_bog_strgraph(g, e1);
										cut_edge_bog_strgraph(g, e2);
										e1 = edge_strgraph(g, bts[0].node_id[1], bts[0].dir[1], bts[0].eidx[1]);
										e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
										cut_edge_bog_strgraph(g, e1);
										cut_edge_bog_strgraph(g, e2);
									} else {
										e1 = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
										e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
										cut_edge_bog_strgraph(g, e1);
										cut_edge_bog_strgraph(g, e2);
										e1 = edge_strgraph(g, bts[1].node_id[1], bts[1].dir[1], bts[1].eidx[1]);
										e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
										cut_edge_bog_strgraph(g, e1);
										cut_edge_bog_strgraph(g, e2);
									}
									bubble ++;
									bub ++;
								} else if(with_loop){
									{
										e1 = edge_strgraph(g, bts[0].node_id[1], bts[0].dir[1], bts[0].eidx[1]);
										e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
										cut_edge_bog_strgraph(g, e1);
										cut_edge_bog_strgraph(g, e2);
									}
									{
										e1 = edge_strgraph(g, bts[1].node_id[1], bts[1].dir[1], bts[1].eidx[1]);
										e2 = edge_strgraph(g, e1->node_id, !e1->dir, e1->rev_idx);
										cut_edge_bog_strgraph(g, e1);
										cut_edge_bog_strgraph(g, e2);
									}
									loop ++;
									bub ++;
								}
							}
							continue;
						}
						m->bt_visit = btid; m->bt_dir = !e->dir; m->bt_off = n->bt_off + e->off; m->bt_idx = e->rev_idx;
						push_heap(heap, e->node_id);
					}
				}
			}
			//fprintf(stderr, "[%s] bubble[%d] %d\n", date(), max_step, bub);
			if(bub == 0) break;
			repair_all_lonely_one_way_edges_bog_strgraph(g);
		}
	}
	free_heap(heap);
	free_u32list(orders);
	return bubble + loop;
}

int bog_step_once_strgraph(StringGraph *g, sglayv *lays){
	sg_layout_t *l, *l2;
	sg_edge_t *e;
	sg_node_t *n1, *n2;
	encap_sglayv(lays, 1);
	l = peer_sglayv(lays);
	n1 = ref_sgnodev(g->nodes, l->node_id);
	if(n1->bogs[1][l->dir][1]) return 0;
	if(n1->bogs[1][l->dir][0] == 0) return 0;
	e = single_living_edge_strgraph(g, l->node_id, l->dir);
	if(get_bitvec(g->node_flags, e->node_id)) return 0;
	n2 = ref_sgnodev(g->nodes, e->node_id);
	if(n2->bogs[0][e->dir][1]) return 0;
	l->eidx[l->dir] = e - edge_strgraph(g, l->node_id, l->dir, 0);
	l2 = next_ref_sglayv(lays);
	l2->node_id = e->node_id;
	l2->dir = e->dir;
	l2->contained = 0;
	l2->eidx[!l2->dir] = e->rev_idx;
	l2->eidx[l2->dir] = 0;
	l2->off = l->off + e->off;
	return 1;
}

int step_once_strgraph(StringGraph *g, sglayv *lays){
	sg_layout_t *l, *l2;
	sg_edge_t *e;
	encap_sglayv(lays, 1);
	l = peer_sglayv(lays);
	if((e = single_living_edge_strgraph(g, l->node_id, l->dir)) == NULL) return 0;
	if(get_bitvec(g->node_flags, e->node_id)) return 0;
	if(count_living_edge_strgraph(g, e->node_id, !e->dir) != 1) return 0;
	l->eidx[l->dir] = e - edge_strgraph(g, l->node_id, l->dir, 0);
	l2 = next_ref_sglayv(lays);
	l2->node_id = e->node_id;
	l2->dir = e->dir;
	l2->contained = 0;
	l2->eidx[!l2->dir] = e->rev_idx;
	l2->eidx[l2->dir] = 0;
	l2->off = l->off + e->off;
	return 1;
}

int length_sglayv(StringGraph *g, sglayv *lays){
	sg_layout_t *l;
	uint32_t i;
	int len;
	if(lays->size == 0) return 0;
	len = 0;
	for(i=0;i<lays->size;i++){
		l = ref_sglayv(lays, i);
		if(l->off + (int)g->rdlens->buffer[l->node_id] > len) len = l->off + g->rdlens->buffer[l->node_id];
	}
	return len;
}

void reverse_flip_sglayv(StringGraph *g, sglayv *lays){
	sg_layout_t *l;
	sg_edge_t *e;
	uint32_t i;
	int off;
	if(lays->size == 0) return;
	reverse_sglayv(lays);
	off = 0;
	for(i=0;i<lays->size;i++){
		l = ref_sglayv(lays, i);
		l->dir = !l->dir;
		l->off = off;
		e = edge_strgraph(g, l->node_id, l->dir, l->eidx[l->dir]);
		off += e->off;
	}
}

void print_simple_dot_strgraph(StringGraph *g, FILE *out){
	sg_node_t *n;
	sg_edge_t *e, *p;
	char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};
	char *dead = "black";
	uint32_t node_id, k, i;
	fprintf(out, "digraph {\n");
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		for(k=0;k<2;k++){
			for(i=0;i<n->edge_cnts[k];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
				if(e->closed) continue;
				fprintf(out, "%s -> %s [label=\"%c%c:%u:%d:%0.3f:%d\" color=\"%s\"]\n",
					(char*)get_cplist(g->rdnames, node_id), (char*)get_cplist(g->rdnames, e->node_id), "+-"[k], "+-"[e->dir], e->off, e->score, e->score * 1.0 / edge_overlap_strgraph(g, node_id, k, i), e->cov, colors[k][e->dir]);
				if(0){
					p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
					if(p->closed){
						fprintf(out, "%s -> %s [label=\"%c%c:%u:%d:%0.3f:%d\" style=dashed color=\"%s\"]\n",
							(char*)get_cplist(g->rdnames, e->node_id), (char*)get_cplist(g->rdnames, node_id), "+-"[!e->dir], "+-"[p->dir], p->off, p->score, p->score * 1.0 / edge_overlap_strgraph(g, e->node_id, !e->dir, e->rev_idx), p->mark, dead);
					}
				}
			}
		}
	}
	fprintf(out, "}\n");
}

void print_dot_strgraph(StringGraph *g, FILE *out){
	sg_node_t *n;
	sg_edge_t *e;
	u32list *stack;
	char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};
	uint32_t node_id, k, i, idx;
	zeros_bitvec(g->node_flags);
	stack = init_u32list(1024);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		if(get_bitvec(g->node_flags, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(count_living_edge_strgraph(g, node_id, 0) + count_living_edge_strgraph(g, node_id, 1) == 0) continue;
		fprintf(out, "digraph %u {\n", node_id);
		clear_u32list(stack);
		push_u32list(stack, node_id);
		one_bitvec(g->node_flags, node_id);
		while(pop_u32list(stack, &idx)){
			n = ref_sgnodev(g->nodes, idx);
			for(k=0;k<2;k++){
				for(i=0;i<n->edge_cnts[k];i++){
					e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
					if(e->closed) continue;
					fprintf(out, "%s -> %s [label=\"%c%c:%u:%d:%0.3f\" color=\"%s\"]\n", (char*)get_cplist(g->rdnames, idx), (char*)get_cplist(g->rdnames, e->node_id), "+-"[k], "+-"[e->dir], e->off, e->score, e->score * 1.0 / edge_overlap_strgraph(g, idx, k, i), colors[k][e->dir]);
					if(get_bitvec(g->node_flags, e->node_id)) continue;
					push_u32list(stack, e->node_id);
					one_bitvec(g->node_flags, e->node_id);
				}
			}
		}
		fprintf(out, "}\n");
	}
	free_u32list(stack);
}

void recurit_contained_reads_strgraph(StringGraph *g, sglayv *lays, sglayv *tmp){
	sg_layout_t *l, *l2;
	sg_node_t *n;
	sg_edge_t *e, *e2;
	uint32_t i, k, j, size;
	int len1;
	size = lays->size;
	clear_sglayv(tmp);
	for(i=0;i<size;i++){
		l = ref_sglayv(lays, i);
		n = ref_sgnodev(g->nodes, l->node_id);
		len1 = g->rdlens->buffer[l->node_id];
		push_sglayv(tmp, *l);
		for(k=0;k<2;k++){
			for(j=0;j<n->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + j);
				if(!is_dead_node_strgraph(g, e->node_id)) continue;
				e2 = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
				if(!e2->att) continue;
				l2 = next_ref_sglayv(tmp);
				l2->node_id = e->node_id;
				l2->dir = l->dir ^ k ^ e->dir;
				l2->contained = 1;
				l2->off = (l->dir ^ k)? (int)(l->off + len1 - (e->off + edge_overlap_strgraph(g, l->node_id, k, j))) : (int)(l->off + e->off);
			}
		}
	}
	clear_sglayv(lays);
	append_sglayv(lays, tmp);
}

uint64_t cut_all_branches_bog_strgraph(StringGraph *g){
	sg_node_t *n;
	sg_edge_t *e, *p;
	uint64_t ret;
	uint32_t node_id, i, k;
	ret = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		for(k=0;k<2;k++){
			if(n->bogs[0][k][1]){
				for(i=0;i<n->edge_cnts[!k];i++){
					e = ref_sgedgev(g->edges, n->edge_offs[!k] + i);
					if(e->closed != 1) continue;
					p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
					if(p->closed) continue;
					cut_edge_bog_strgraph(g, p);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint32_t gen_unitigs_layout_strgraph(StringGraph *g){
	sglayv *lays;
	sg_layout_t *l;
	sg_node_t *n;
	uint32_t node_id, ret, i, j;
	ret = 0;
	zeros_bitvec(g->node_flags);
	for(i=0;i<g->lays->size;i++){
		lays = (sglayv*)get_vplist(g->lays, i);
		free_sglayv(lays);
	}
	clear_vplist(g->lays);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		n = ref_sgnodev(g->nodes, node_id);
		n->lay_id = 0xFFFFFFFFU;
		n->dir = 0;
		n->off = 0;
		n->end = 0;
	}
	cut_all_branches_bog_strgraph(g);
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		if(get_bitvec(g->node_flags, node_id)) continue;
		if(g->rdlens->buffer[node_id] == 0) continue;
		lays = init_sglayv(32);
		l = next_ref_sglayv(lays);
		l->node_id = node_id;
		l->dir = 0;
		l->off = 0;
		l->eidx[0] = 0;
		l->eidx[1] = 0;
		one_bitvec(g->node_flags, node_id);
		while(bog_step_once_strgraph(g, lays)){ one_bitvec(g->node_flags, peer_sglayv(lays)->node_id); }
		reverse_flip_sglayv(g, lays);
		while(bog_step_once_strgraph(g, lays)){ one_bitvec(g->node_flags, peer_sglayv(lays)->node_id); }
		//if(lays->size < MIN_LAY_NODES) continue;
		push_vplist(g->lays, lays);
		ret ++;
	}
	for(i=0;i<g->lays->size;i++){
		lays = (sglayv*)get_vplist(g->lays, i);
		if(lays->size < MIN_LAY_NODES) continue;
		for(j=0;j<lays->size;j++){
			l = ref_sglayv(lays, j);
			n = ref_sgnodev(g->nodes, l->node_id);
			n->lay_id = i;
			n->dir = l->dir;
			n->off = l->off;
			n->end = (j < 2)? 1 : ((j + 1 + 1 > lays->size)? 1 : 0);
		}
	}
	return ret;
}

uint64_t recover_edges_inter_unitigs_strgraph(StringGraph *g, float best_score_cutoff){
	sg_node_t *n1, *n2;
	sg_edge_t *e;
	uint64_t ret;
	uint32_t node_id, k, j, b;
	int best;
	float bestS;
	ret = 0;
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n1 = ref_sgnodev(g->nodes, node_id);
		if(n1->lay_id == 0xFFFFFFFFU) continue;
		if(n1->end == 0) continue;
		for(k=0;k<2;k++){
			bestS = 0;
			for(j=0;j<n1->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e->closed != 0 && e->closed != 1) continue;
				n2 = ref_sgnodev(g->nodes, e->node_id);
				if(n2->lay_id == 0xFFFFFFFFU) continue;
				if(n2->end == 0) continue;
				if(e->score > bestS){ bestS = e->score; }
			}
			if(bestS == 0) continue;
			bestS = bestS * best_score_cutoff;
			best = g->rdlens->buffer[node_id]; b = 0xFFFFFFFFU;
			for(j=0;j<n1->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e->closed != 0 && e->closed != 1) continue;
				n2 = ref_sgnodev(g->nodes, e->node_id);
				if(n2->lay_id == 0xFFFFFFFFU) continue;
				if(n2->end == 0) continue;
				if(e->score < bestS) continue;
				if(e->off < best){ best = e->off; b = j; }
			}
			if(b == 0xFFFFFFFFU) continue;
			e = ref_sgedgev(g->edges, n1->edge_offs[k] + b);
			if(e->closed == 0) continue;
			for(j=0;j<n1->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e->closed) continue;
				cut_edge_bog_strgraph(g, e);
			}
			e = ref_sgedgev(g->edges, n1->edge_offs[k] + b);
			ret += revive_edge_bog_strgraph(g, e);
		}
	}
	repair_all_lonely_one_way_edges_bog_strgraph(g);
	return ret;
}

void print_recovered_dot_strgraph(StringGraph *g, FILE *out){
	sg_node_t *n;
	sg_edge_t *e;
	char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};
	uint32_t node_id, k, i;
	fprintf(out, "digraph {\n");
	for(node_id=0;node_id<g->nodes->size;node_id++){
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(n->lay_id == 0xFFFFFFFFU) continue;
		for(k=0;k<2;k++){
			for(i=0;i<n->edge_cnts[k];i++){
				e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
				if(e->closed) continue;
				if(g->nodes->buffer[e->node_id].lay_id == 0xFFFFFFFFU) continue;
				//fprintf(out, "%s -> %s [color=\"%s\"%s]\n",
					//(char*)get_cplist(g->rdnames, node_id), (char*)get_cplist(g->rdnames, e->node_id), colors[k][e->dir], (n->lay_id == g->nodes->buffer[e->node_id].lay_id)? "" : "style=dashed");
				//continue;
				fprintf(out, "%s -> %s [label=\"%c%c:%u:%d:%0.3f\" color=\"%s\"%s]\n",
					(char*)get_cplist(g->rdnames, node_id), (char*)get_cplist(g->rdnames, e->node_id), "+-"[k], "+-"[e->dir], e->off, e->score, e->score * 1.0 / edge_overlap_strgraph(g, node_id, k, i),
					colors[k][e->dir], (n->lay_id == g->nodes->buffer[e->node_id].lay_id)? "" : "style=dashed");
			}
		}
	}
	fprintf(out, "}\n");
}

int is_duplicated_layout_strgraph(StringGraph *g, sglayv *lays, float min_cov, u64hash *hash, u64list *vec, uint32_t *utg, float *max_cov){
	sg_layout_t *l;
	sg_node_t *n1, *n2;
	sg_edge_t *e;
	uint64_t val;
	uint32_t i, j, k, max, layid, cnt, tot_len, cov_len, cov, x, y;
	clear_u64hash(hash);
	for(i=0;i<lays->size;i++){
		l = ref_sglayv(lays, i);
		n1 = ref_sgnodev(g->nodes, l->node_id);
		for(k=0;k<2;k++){
			for(j=0;j<n1->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e->closed != 1) continue;
				n2 = ref_sgnodev(g->nodes, e->node_id);
				if(n2->lay_id == 0xFFFFFFFFU || n2->lay_id == n1->lay_id) continue;
				val = (((uint64_t)n2->lay_id) << 32) | i;
				put_u64hash(hash, val);
			}
		}
	}
	clear_u64list(vec);
	reset_iter_u64hash(hash);
	while(iter_u64hash(hash, &val)){
		push_u64list(vec, val >> 32);
	}
	if(vec->size < 1) return 0;
	sort_array(vec->buffer, vec->size, uint64_t, a > b);
	max = cnt = 0; layid = 0xFFFFFFFFU;
	for(i=0;i<=vec->size;i++){
		if(i == 0) cnt ++;
		else if(i == vec->size || vec->buffer[i] != vec->buffer[i-1]){
			if(cnt > max){ max = cnt; layid = vec->buffer[i-1]; }
			cnt = 1;
		} else cnt ++;
	}
	*utg = layid;
	tot_len = cov_len = 0;
	for(i=0;i<lays->size;i++){
		l = ref_sglayv(lays, i);
		n1 = ref_sgnodev(g->nodes, l->node_id);
		tot_len += g->rdlens->buffer[l->node_id];
		clear_u64list(vec);
		for(k=0;k<2;k++){
			for(j=0;j<n1->edge_cnts[k];j++){
				e = ref_sgedgev(g->edges, n1->edge_offs[k] + j);
				if(e->closed != 1) continue;
				n2 = ref_sgnodev(g->nodes, e->node_id);
				if(n2->lay_id != layid) continue;
				if(k){
					y = e->off;
					x = y + edge_overlap_strgraph(g, l->node_id, k, j);
					x = g->rdlens->buffer[l->node_id] - x;
					y = g->rdlens->buffer[l->node_id] - y;
				} else {
					x = e->off;
					y = x + edge_overlap_strgraph(g, l->node_id, k, j);
				}
				push_u64list(vec, (((uint64_t)x) << 32) | y);
			}
		}
		if(vec->size == 0) continue;
		sort_array(vec->buffer, vec->size, uint64_t, (a >> 32) > (b >> 32));
		val = vec->buffer[0];
		x = val >> 32;
		y = val & 0xFFFFFFFFU;
		cov = 0;
		for(j=1;j<vec->size;j++){
			val = vec->buffer[j];
			if((val >> 32) > y){
				cov += y - x;
				x = (val >> 32); y = val & 0xFFFFFFFFU;
			} else if((val & 0xFFFFFFFFU) > y){
				y = val & 0xFFFFFFFFU;
			}
		}
		cov += y - x;
		cov_len += cov;
	}
	//fprintf(stderr, "layout[%u nodes, %u bp] %u / %u = %0.3f\n", (uint32_t)lays->size, length_sglayv(g, lays), cov_len, tot_len, 1.0 * cov_len / tot_len);
	*max_cov = 1.0 * cov_len / tot_len;
	return (cov_len >= (uint32_t)(min_cov * tot_len));
}

void output_unitigs_layout_strgraph(StringGraph *g, FILE *out_lay, FILE *dup_lay, FILE *out_utg, FILE *out_lnk, FILE *out_dup, float similar_unitig_cov){
	sglayv *lays, *tmp;
	sg_layout_t *l;
	sg_edge_t *e, *p;
	sg_node_t *n1, *n2;
	String *seq, *ctg;
	u64hash *hash;
	u64list *vec;
	FILE *out_ptr;
	uint8_t *qvs;
	uint32_t i, j, k, m, len, ret, is_dup, dup_utg;
	float dup_cov;
	seq = init_string(1024);
	ctg = init_string(1024);
	tmp = init_sglayv(1024);
	hash = init_u64hash(1023);
	vec = init_u64list(1024);
	ret = 0;
	for(i=0;i<g->lays->size;i++){
		lays = (sglayv*)get_vplist(g->lays, i);
		if(lays->size < MIN_LAY_NODES){
			is_dup = 1; dup_utg = 19830203; dup_cov = 0.0;
		} else is_dup = is_duplicated_layout_strgraph(g, lays, similar_unitig_cov, hash, vec, &dup_utg, &dup_cov);
		recurit_contained_reads_strgraph(g, lays, tmp);
		len = length_sglayv(g, lays);
		if(is_dup){
			fprintf(out_dup, ">utg%u length=%d nodes=%u dup=utg%u cov=%0.3f\n", i, len, (uint32_t)lays->size, dup_utg, dup_cov);
			fprintf(dup_lay, ">utg%u length=%d nodes=%u dup=utg%u cov=%0.3f\n", i, len, (uint32_t)lays->size, dup_utg, dup_cov);
			out_ptr = dup_lay;
		} else {
			fprintf(out_utg, ">utg%u length=%d nodes=%u\n", i, len, (uint32_t)lays->size);
			fprintf(out_lay, ">utg%u length=%d nodes=%u\n", i, len, (uint32_t)lays->size);
			out_ptr = out_lay;
			ret ++;
		}
		clear_string(ctg); encap_string(ctg, len);
		len = 0;
		for(j=0;j<lays->size;j++){
			l = ref_sglayv(lays, j);
			n1 = ref_sgnodev(g->nodes, l->node_id);
			for(k=0;k<2;k++){
				for(m=0;m<n1->edge_cnts[k];m++){
					e = ref_sgedgev(g->edges, n1->edge_offs[k] + m);
					if(e->closed == 2) continue;
					n2 = ref_sgnodev(g->nodes, e->node_id);
					if(n2->lay_id == i) continue;
					if(n2->lay_id == 0xFFFFFFFFU) continue;
					p = edge_strgraph(g, e->node_id, !e->dir, e->rev_idx);
					fprintf(out_lnk, "utg%u\t%s\t%c\t%d\tutg%u\t%s\t%c\t%d",
						n1->lay_id, get_cplist(g->rdnames, l->node_id), "+-"[n1->dir], n1->off,
						n2->lay_id, get_cplist(g->rdnames, e->node_id), "+-"[n2->dir], n2->off);
					fprintf(out_lnk, "\t%c\t%u\t%d\t%d", "+-"[k], g->rdlens->buffer[l->node_id], e->off, e->off + edge_overlap_strgraph(g, l->node_id, k, m));
					fprintf(out_lnk, "\t%c\t%u\t%d\t%d", "+-"[e->dir], g->rdlens->buffer[e->node_id], p->off, p->off + edge_overlap_strgraph(g, e->node_id, !e->dir, e->rev_idx));
					fprintf(out_lnk, "\t%d\n", e->score);
				}
			}
			clear_string(seq); encap_string(seq, g->rdlens->buffer[l->node_id]);
			if(l->dir) revseq_basebank(g->rdseqs, g->rdoffs->buffer[l->node_id], g->rdlens->buffer[l->node_id], seq->string);
			else seq_basebank(g->rdseqs, g->rdoffs->buffer[l->node_id], g->rdlens->buffer[l->node_id], seq->string);
			seq->size = g->rdlens->buffer[l->node_id];
			fprintf(out_ptr, "%c\t%s\t%c\t%d\t%d\t%s", "YN"[l->contained], (char*)get_cplist(g->rdnames, l->node_id), "+-"[l->dir], l->off, g->rdlens->buffer[l->node_id], seq->string);
			if(g->rdqlens->buffer[l->node_id] == 7 * g->rdlens->buffer[l->node_id]){ // f5q format
				fputc('\t', out_ptr);
				qvs = g->rdqvs->buffer + g->rdqoffs->buffer[l->node_id];
				if(l->dir){
					for(k=0;k<5;k++){
						for(m=0;m<g->rdlens->buffer[l->node_id];m++){
							fputc(qvs[(g->rdlens->buffer[l->node_id] - m - 1) * 7 + k], out_ptr);
						}
					}
					for(k=5;k<7;k++){
						for(m=0;m<g->rdlens->buffer[l->node_id];m++){
							fputc(reverse_dna_base(qvs[(g->rdlens->buffer[l->node_id] - m - 1) * 7 + k]), out_ptr);
						}
					}
				} else {
					for(k=0;k<7;k++){
						for(m=0;m<g->rdlens->buffer[l->node_id];m++){
							fputc(qvs[m * 7 + k], out_ptr);
						}
					}
				}
			}
			fputc('\n', out_ptr);
			if(l->contained || l->off + seq->size <= (int)len) continue;
			memcpy(ctg->string + l->off, seq->string, seq->size);
			len = l->off + seq->size;
		}
		ctg->string[len] = 0;
		ctg->size = len;
		print_pretty_seq(is_dup? out_dup : out_utg, ctg, 100);
	}
	fprintf(stdout, "[%s] output %u independent unitigs\n", date(), ret);
	free_string(seq);
	free_string(ctg);
	free_sglayv(tmp);
	free_u64hash(hash);
	free_u64list(vec);
}

int usage(){
	printf(
	"WTLay: Layout of long reads using Best Overlap Graph\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtlay [options]\n"
	"Options:\n"
	" -i <string> Long reads sequences file(s), + *\n"
	" -b <string> Long reads retained region, often from wtobt, +\n"
	"             Format: read_name\\toffset\\tlength\\toriginal_len\n"
	//" -C <string> Read name per line from this file will be marked as contained, +\n"
	" -j <string> Overlap file(s), + *\n"
	"             Format: reads1\\t+/-\\tlen1\\tbeg1\\tend1\\treads2\\t+/-\\tlen2\\tbeg2\\tend2\\tscore\n"
	" -s <int>    Minimum alignment matched bases, [500]\n"
	" -m <float>  Minimum alignment identity, [0.6]\n"
	" -w <int>    Maximum margin of alignment, [100]\n"
	" -o <string> Output file of layouts, *\n"
	" -f          Force overwrite output file\n"
	" -c <int>    Minimum estimated coverage of edge to be trusted, [1]\n"
	"             edge coverage is calculated by counting overlaps that can replace this edge\n"
	" -R          Use number of matches as alignment score\n"
	" -r <float>  Best score cutoff, say best overlap MUST have alignment score >= <-r> * read's best score [0.95]\n"
	" -q <float>  Minimum coverage for similar unitig detection, say the aligned length of unitig A by unitig B, divided by total length of unitig A, [0.4]\n"
	" -u <int>    Min nodes of a layout to be output as independent unitig, [4]\n"
	//" -v <float>  Length variance between two long reads in the aligned region, [0.05]\n"
	" -B <int>    Maximum step in tracing bubbles, [20]\n"
	" -S <float>  Variance threshold of (alignment score / overlap) between ture and false overlaps, [0.50]\n"
	"             Used in better_overlap_strgraph\n"
	"---- DEBUG option ---\n"
	" -Q <string> Commands, [gCwgBgRURg]\n"
	"             G: print_dot_strgraph\n"
	"             g: print_simple_dot_strgraph\n"
	"             w: mask_low_cov_edge_strgraph\n"
	"             C: mask_contained_reads_strgraph\n"
	"             B: best_overlap_strgraph\n"
	"             t: bog_cut_tips_strgraph\n"
	"             M: bog_merge_bubbles_strgraph\n"
	"             U: recover_edges_inter_unitigs_strgraph\n"
	"             R: repair_best_overlap_strgraph\n"
	"             b: better_overlap_strgraph\n"
	//"             D: remove_duplicate_edges_strgraph\n"
	"             O: mask_self_circle_reads_strgraph\n"
	"             T: reduce_transitive_strgraph\n"
	//"             E: enhanced_reduce_transitive_strgraph\n"
	"             L: longest_overlap_strgraph\n"
	//"             S: best_score_overlap_strgraph\n"
	//"             X: detect_chimeric_reads_strgraph\n"
	"             .: exit program\n"
	"\n"
	"Example: \n"
	"$> wtlay -i wt.fa -b wt.zmo.obt -j wt.zmo.ext -o wt.lay -m 0.6 -c 1 -r 0.95\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	obj_desc_t wsg = f32list_obj_desc;
	wsg = pbregv_obj_desc;
	wsg = tracev_obj_desc;
	obj_desc_t ttt = wsg;
	wsg = ttt;
	StringGraph *g;
	FileReader *fr;
	Sequence *seq;
	cplist *pbs, *ovls, *obts, *ctns;
	char *prefix, *dot_file, *lay_file, *day_file, *utg_file, *lnk_file, *dup_file, *log_file;
	char *commands;
	FILE *out_dot, *out_lay, *dup_lay, *out_utg, *out_lnk, *out_dup, *log;
	unsigned long long n;
	int c, edgecov_cutoff, force_overwrite, dot_idx, min_score, margin, mat_score;
	float len_var, score_var, best_score_cutoff, min_id, utg_sm;
	len_var = 0;
	float wushigang = len_var;
	float tmp = wushigang;
	wushigang = tmp;
	pbs = init_cplist(4);
	ovls = init_cplist(4);
	obts = init_cplist(4);
	ctns = init_cplist(4);
	prefix = NULL;
	min_score = 500;
	min_id = 0.6;
	margin = 100;
	len_var = 0.05;
	score_var = 0.50;
	utg_sm = 0.4;
	MERGE_BUBBLE_MAX_STEP = 20;
	edgecov_cutoff = 1;
	mat_score = 0;
	best_score_cutoff = 0.95;
	force_overwrite = 0;
	dot_idx = 0;
	commands = "gCwgBgRURg";
	while((c = getopt(argc, argv, "hfi:b:C:j:s:m:w:o:Rr:q:u:v:c:B:S:Q:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'f': force_overwrite = 1; break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			//case 'C': push_cplist(ctns, optarg); break;
			case 'j': push_cplist(ovls, optarg); break;
			case 's': min_score = atoi(optarg); break;
			case 'm': min_id = atof(optarg); break;
			case 'w': margin = atoi(optarg); break;
			case 'o': prefix = optarg; break;
			case 'R': mat_score = 1; break;
			case 'r': best_score_cutoff = atof(optarg); break;
			case 'q': utg_sm = atof(optarg); break;
			case 'u': MIN_LAY_NODES = atof(optarg); break;
			case 'v': len_var = atof(optarg); break;
			case 'c': edgecov_cutoff = atoi(optarg); break;
			case 'S': score_var = atof(optarg); break;
			case 'B': MERGE_BUBBLE_MAX_STEP = atoi(optarg); break;
			case 'Q': commands = optarg; break;
			default: return usage();
		}
	}
	if(prefix == NULL) return usage();
	if(pbs->size == 0) return usage();
	if(ovls->size == 0) return usage();
	if(!force_overwrite && file_exists(prefix)){
		fprintf(stdout, "File exists! '%s'\n\n", prefix);
		return usage();
	}
	g = init_strgraph();
	g->min_score = min_score;
	g->min_id = min_id;
	g->max_ovl_margin = margin;
	g->mat_score = mat_score;
	if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fprintf(stdout, "[%s] loading reads\n", date());
	seq = NULL;
	while(fread_seq(&seq, fr)){
		if(seq->qual.size == seq->seq.size * 7){ // f5q format
			push_read5q_strgraph(g, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size, seq->qual.string);
		} else push_read_strgraph(g, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
		if((g->n_rd % 1000) == 0){
			fprintf(stdout, "\r%u", g->n_rd); fflush(stdout);
		}
	}
	fclose_filereader(fr);
	fprintf(stdout, "\r[%s] Done, %u reads\n", date(), (unsigned)g->n_rd); fflush(stdout);
	if(obts->size){
		fprintf(stdout, "[%s] loading reads obt information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(fr->line->string[0] == '#') continue;
			if(c < 3) continue;
			set_read_clip_strgraph(g, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stdout, "[%s] Done\n", date()); fflush(stdout);
	} else {
		fprintf(stdout, "[%s] No obt information\n", date()); fflush(stdout);
	}
	generate_nodes_strgraph(g);
	if((fr = fopen_m_filereader(ovls->size, ovls->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", ovls->buffer[0], __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fprintf(stdout, "[%s] loading alignments\n", date()); fflush(stdout);
	free_sgbiedgev(load_overlaps_strgraph(g, fr, NULL));
	fclose_filereader(fr);
	fprintf(stdout, "[%s] Done\n", date()); fflush(stdout);
	fprintf(stdout, "[%s] calculating edge coverage ...\n", date());
	cal_edge_coverage_strgraph(g);
	n = remove_duplicate_edges_strgraph(g);
	fprintf(stdout, "[%s] removed %llu duplicate edges\n", date(), n);
	fprintf(stdout, "[%s] Done\n", date());
	for(c=0;c<(int)strlen(commands);c++){
		fflush(stdout);
		switch(commands[c]){
			case '.':
				fprintf(stderr, "[%s] exit program immediately\n", date());
				return 0;
				break;
			case 'G':
				dot_file = malloc(strlen(prefix) + 30);
				sprintf(dot_file, "%s.%d.dot", prefix, ++dot_idx);
				out_dot = fopen(dot_file, "w");
				print_dot_strgraph(g, out_dot);
				fclose(out_dot);
				free(dot_file);
				break;
			case 'g':
				dot_file = malloc(strlen(prefix) + 30);
				sprintf(dot_file, "%s.%d.dot", prefix, ++dot_idx);
				out_dot = fopen(dot_file, "w");
				print_simple_dot_strgraph(g, out_dot);
				fclose(out_dot);
				free(dot_file);
				break;
			case 'w':
				//fprintf(stdout, "[%s] calculating edge coverage ...\n", date());
				//cal_edge_coverage_strgraph(g);
				n = mask_low_cov_edge_strgraph(g, edgecov_cutoff);
				fprintf(stdout, "[%s] masked %llu low coverage (<%u) edges\n", date(), n, edgecov_cutoff);
				break;
			case 'B':
				n = best_overlap_strgraph(g, best_score_cutoff);
				fprintf(stdout, "[%s] 'best_overlap' cut %llu non-best edges\n", date(), n);
				//check_all_nodes_bogs_strgraph(g);
				break;
			case 'b':
				n = better_overlap_strgraph(g, score_var);
				fprintf(stdout, "[%s] 'better_overlap' cut %llu bad edges\n", date(), n);
				break;
			case 'C':
				log_file = catstr(2, prefix, ".contained_reads");
				log = fopen(log_file, "w");
				n = mask_contained_reads_strgraph(g, log);
				fclose(log);
				free(log_file);
				fprintf(stdout, "[%s] masked %llu contained reads\n", date(), n);
				break;
			case 'O':
				n = mask_self_circle_reads_strgraph(g);
				fprintf(stdout, "[%s] masked %llu self circle reads\n", date(), n);
				break;
			case 'T':
				n = reduce_transitive_strgraph(g);
				fprintf(stdout, "[%s] reduced %llu transitive edges\n", date(), n);
				break;
			case 'M':
				n = bog_cut_tips_strgraph(g, 10);
				fprintf(stdout, "[%s] cut %llu read tips\n", date(), n);
				while(1){
					n = bog_merge_bubbles_strgraph(g, 20, 0);
					fprintf(stdout, "[%s] merge %llu bubbles\n", date(), n);
					if(n == 0) break;
					n = bog_cut_tips_strgraph(g, 10);
					fprintf(stdout, "[%s] cut %llu read tips\n", date(), n);
					n = recover_paired_dead_ends_bog_strgraph(g);
					fprintf(stdout, "[%s] recover %llu paired dead ends\n", date(), n);
				}
				while(1){
					n = bog_merge_bubbles_strgraph(g, 20, 1);
					fprintf(stdout, "[%s] merge %llu loops/bubbles\n", date(), n);
					if(n == 0) break;
					n = bog_cut_tips_strgraph(g, 10);
					fprintf(stdout, "[%s] cut %llu read tips\n", date(), n);
					n = recover_paired_dead_ends_bog_strgraph(g);
					fprintf(stdout, "[%s] recover %llu paired dead ends\n", date(), n);
				}
				break;
			case 't':
				n = bog_cut_tips_strgraph(g, 10);
				fprintf(stdout, "[%s] cut %llu read tips\n", date(), n);
				break;
			case 'U':
				n = gen_unitigs_layout_strgraph(g);
				fprintf(stdout, "[%s] generated %llu unitigs\n", date(), n);
				n = recover_edges_inter_unitigs_strgraph(g, best_score_cutoff);
				fprintf(stdout, "[%s] recovered %llu edges inter unitigs\n", date(), n);
				break;
			case 'L':
				n = longest_overlap_strgraph(g);
				fprintf(stdout, "[%s] 'longest_overlap' cut %llu edges\n", date(), n);
				break;
			case 'R':
				while(1){
					n = repair_best_overlap_strgraph(g);
					//check_all_nodes_bogs_strgraph(g);
					if(n == 0) break;
					fprintf(stdout, "[%s] repair %llu bog elements\n", date(), n);
				}
				break;
			case 'S':
				n = best_score_overlap_strgraph(g);
				fprintf(stdout, "[%s] 'best_score' cut %llu edges\n", date(), n);
				break;
			case 'X':
				log_file = catstr(2, prefix, ".chimeric");
				log = fopen(log_file, "w");
				n = detect_chimeric_reads_strgraph(g, log);
				fclose(log);
				free(log_file);
				fprintf(stdout, "[%s] masked %llu chimeric reads\n", date(), n);
				break;
			default: fprintf(stdout, "Unknown command '%c'\n", commands[c]);
		}
	}
	fflush(stdout);
	lay_file = catstr(2, prefix, "");
	if((out_lay = fopen(lay_file, "w")) == NULL){
		fprintf(stderr, " -- Cannot open(write) %s in %s -- %s:%d --\n", lay_file, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	day_file = catstr(2, prefix, ".dup");
	if((dup_lay = fopen(day_file, "w")) == NULL){
		fprintf(stderr, " -- Cannot open(write) %s in %s -- %s:%d --\n", day_file, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	utg_file = catstr(2, prefix, ".utg");
	if((out_utg = fopen(utg_file, "w")) == NULL){
		fprintf(stderr, " -- Cannot open(write) %s in %s -- %s:%d --\n", utg_file, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	dup_file = catstr(2, prefix, ".utg.dup");
	if((out_dup = fopen(dup_file, "w")) == NULL){
		fprintf(stderr, " -- Cannot open(write) %s in %s -- %s:%d --\n", dup_file, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	lnk_file = catstr(2, prefix, ".lnk");
	if((out_lnk = fopen(lnk_file, "w")) == NULL){
		fprintf(stderr, " -- Cannot open(write) %s in %s -- %s:%d --\n", lnk_file, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	n = gen_unitigs_layout_strgraph(g);
	fprintf(stdout, "[%s] generated %llu unitigs\n", date(), n);
	n = recover_edges_inter_unitigs_strgraph(g, best_score_cutoff);
	fprintf(stdout, "[%s] recover %llu edges inter unitigs\n", date(), n);
	{
		dot_file = malloc(strlen(prefix) + 30);
		sprintf(dot_file, "%s.%d.dot", prefix, ++dot_idx);
		out_dot = fopen(dot_file, "w");
		print_recovered_dot_strgraph(g, out_dot);
		fclose(out_dot);
		free(dot_file);
	}
	output_unitigs_layout_strgraph(g, out_lay, dup_lay, out_utg, out_lnk, out_dup, utg_sm);
	fclose(out_lay);
	fclose(dup_lay);
	fclose(out_utg);
	fclose(out_lnk);
	fclose(out_dup);
	free(lay_file);
	free(day_file);
	free(utg_file);
	free(lnk_file);
	free(dup_file);
	free_cplist(pbs);
	free_cplist(ovls);
	free_cplist(obts);
	free_cplist(ctns);
	free_strgraph(g);
	fprintf(stdout, "[%s] Done\n", date());
	return 0;
}
