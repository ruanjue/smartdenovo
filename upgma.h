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

#ifndef __HEAP_RJ_H
#define __HEAP_RJ_H

#include "list.h"

define_list(f32list, float);

typedef struct {
	uint32_t child, sibling, parent;
} upgma_tre_t;
define_list(upgmatrev, upgma_tre_t);

typedef struct {
	uint32_t n;
	f32list *mat;
	upgmatrev *tre;
} UPGMA;

static inline UPGMA* init_upgma(){
	UPGMA *g;
	g = malloc(UPGMA);
	g->n = 0;
	g->mat = init_f32list(32);
	g->tre = init_upgmatrev(32);
	return g;
}

static inline void reset_upgma(uint32_t n){
	uint32_t i;
	if(n < 2) n = 2;
	g->n = n;
	clear_and_encap_f32list(g->mat, n * (n - 1) / 2);
	g->mat->size = n * (n - 1) / 2;
	for(i=0;i<g->mat->size;i++) g->mat->buffer[i] = 0;
	clear_and_encap_upgmatrev(g->tre, n);
	for(i=0;i<n;i++){
		push_u32list(g->tre, (upgma_tre_t){0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU});
	}
}

// i < j
static inline void set_upgma(UPGMA *g, uint32_t i, uint32_t j, float d){
	uint32_t idx;
	idx = (j * (j - 1)) / 2 + i;
	g->mat->buffer[idx] = d;
}

// i < j
static inline float get_upgma(UPGMA *g, uint32_t i, uint32_t j){
	uint32_t idx;
	idx = (j * (j - 1)) / 2 + i;
	return g->mat->buffer[idx];
}

static inline float get2_upgma(UPGMA *g, uint32_t i, uint32_t j){
	uint32_t idx;
	idx = (i < j)? (j * (j - 1)) / 2 + i : (i * (i - 1)) / 2 + j;
	return g->mat->buffer[idx];
}

static inline void do_upgma(UPGMA *g){
	uint32_t i, j, iter, mi, mj;
	float max, dis;
	for(iter=0;iter<g->n-1;iter++){
		max = 0; mi = mj = 0;
		for(i=0;i<g->tre->size;i++){
			if(g->tre->buffer[i].parent != 0xFFFFFFFFU) continue;
			for(j=i+1;j<g->tre->size;j++){
				if(g->tre->buffer[j].parent != 0xFFFFFFFFU) continue;
				dis = get_upgma(g, i, j);
				if(dis > max){ max = dis; mi = i; mj = j; }
			}
		}
		g->tre->buffer[mi].sibling = mj;
		g->tre->buffer[mi].parent  = g->tre->size;
		g->tre->buffer[mj].sibling = 0xFFFFFFFFU;
		g->tre->buffer[mj].parent  = g->tre->size;
		push_upgmatrev(g->tre, (upgma_tre_t){mi, 0xFFFFFFFFU, 0xFFFFFFFFU});
		// update mat
		for(i=0;i<g->tre->size-1;i++){
			if(g->tre->buffer[i].parent == 0xFFFFFFFFU){
				dis = (get2_upgma(g, i, mi) + get2_upgma(g, i, mj)) / 2;
			} else {
				dis = 0;
			}
			push_f32list(g->mat, dis);
		}
	}
}

// Order: Left, Right, Top
// @return: 0, finish; 1, '('; 2, node; 3, ','; 4, ')'
static inline int iter_upgma(UPGMA *g, uint32_t *_idx, int flag){
	uint32_t idx;
	int flag;
	idx = *_idx;
	if(idx == 0xFFFFFFFFU) idx = g->tre->size - 1;
	if(flag == 1){
		if(g->tre->buffer[idx].child == 0xFFFFFFFFU){
			flag = 2;
		} else {
			idx = g->tre->buffer[idx].child;
			flag = 1;
		}
	} else if(flag == 2){
		if(g->tre->buffer[idx].sibling == 0xFFFFFFFFU){
			flag = 4;
		} else {
			idx = g->tre->buffer[idx].sibling;
			flag = 3;
		}
	} else if(flag == 3){
		flag = 1;
	} else {
		if(g->tre->buffer[idx].parent == 0xFFFFFFFFU){
			flag = 0;
		} else {
			idx = g->tre->buffer[idx].parent;
			flag = 2;
		}
	}
	*_idx = idx;
	return flag;
}

static inline void print_upgma(UPGMA *g, FILE *out){
	uint32_t idx;
	int flag;
	idx = 0xFFFFFFFFU; flag = 1;
	while((flag = iter_upgma(g, &idx, flag))){
		if(flag == 2){
			fprintf(out, "%u", idx);
		} else {
			fputc("*(*,)*"[flag], out);
		}
	}
	fprintf(out, "\n");
}

#endif
