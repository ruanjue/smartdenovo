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

typedef int (*heap_comp_func)(uint32_t idx1, uint32_t idx2, void *ref);

static inline int heap_cmp_u32_func(uint32_t idx1, uint32_t idx2, void *ref){ return idx1 == idx2? 0 : (idx1 < idx2? -1 : 1); if(ref) return 0; }
static inline int heap_cmp_b32_func(uint32_t idx1, uint32_t idx2, void *ref){ int v1, v2; v1 = idx1; v2 = idx2; return v1 == v2? 0 : (v1 < v2? -1 : 1); if(ref) return 0; }

typedef struct {
	u32list *ptrs;
	void *ref;
	heap_comp_func cmp;
} Heap;

static inline Heap* init_heap(heap_comp_func cmp, void *ref){
	Heap *heap;
	heap = malloc(sizeof(Heap));
	heap->ptrs = init_u32list(16);
	heap->cmp  = cmp;
	heap->ref  = ref;
	return heap;
}

static inline void free_heap(Heap *heap){ free_u32list(heap->ptrs); free(heap); }

static inline void clear_heap(Heap *heap){ clear_u32list(heap->ptrs); }

static inline void push_heap(Heap *heap, uint32_t p){
	uint32_t pp;
	size_t i;
	i = count_u32list(heap->ptrs);
	push_u32list(heap->ptrs, p);
	while(i && heap->cmp(get_u32list(heap->ptrs, i), get_u32list(heap->ptrs, (i - 1) >> 1), heap->ref) < 0){
		pp = get_u32list(heap->ptrs, i);
		set_u32list(heap->ptrs, i, get_u32list(heap->ptrs, (i - 1) >> 1));
		set_u32list(heap->ptrs, (i - 1) >> 1, pp);
		i = (i - 1) >> 1;
	}
}

static inline size_t count_heap(Heap *heap){ return count_u32list(heap->ptrs); }

static inline uint32_t peer_heap(Heap *heap){ return (count_u32list(heap->ptrs)? get_u32list(heap->ptrs, 0) : 0xFFFFFFFFU);}

static inline void remove_heap(Heap *heap, size_t idx){
	uint32_t pp;
	size_t swap;
	set_u32list(heap->ptrs, idx, get_u32list(heap->ptrs, count_u32list(heap->ptrs) - 1));
	trunc_u32list(heap->ptrs, 1);
	while((idx << 1) + 1 < count_u32list(heap->ptrs)){
		swap = idx;
		if(heap->cmp(get_u32list(heap->ptrs, swap), get_u32list(heap->ptrs, (idx << 1) + 1), heap->ref) > 0){
			swap = (idx << 1) + 1;
		}
		if((idx << 1) + 2 < count_u32list(heap->ptrs) && heap->cmp(get_u32list(heap->ptrs, swap), get_u32list(heap->ptrs, (idx << 1) + 2), heap->ref) > 0){
			swap = (idx << 1) + 2;
		}
		if(swap == idx) break;
		pp = get_u32list(heap->ptrs, idx);
		set_u32list(heap->ptrs,  idx, get_u32list(heap->ptrs, swap));
		set_u32list(heap->ptrs, swap, pp);
		idx = swap;
	}
}

static inline uint32_t pop_heap(Heap *heap){
	uint32_t p;
	if(count_u32list(heap->ptrs)){
		p = get_u32list(heap->ptrs, 0);
		remove_heap(heap, 0);
		return p;
	} else return 0xFFFFFFFFU;
}

static inline uint32_t replace_heap(Heap *heap, size_t idx, uint32_t p){
	uint32_t ret, pp;
	size_t swap;
	ret = heap->ptrs->buffer[idx];
	heap->ptrs->buffer[idx] = p;
	while((idx << 1) + 1 < count_u32list(heap->ptrs)){
		swap = idx;
		if(heap->cmp(get_u32list(heap->ptrs, swap), get_u32list(heap->ptrs, (idx << 1) + 1), heap->ref) > 0){
			swap = (idx << 1) + 1;
		}
		if((idx << 1) + 2 < count_u32list(heap->ptrs) && heap->cmp(get_u32list(heap->ptrs, swap), get_u32list(heap->ptrs, (idx << 1) + 2), heap->ref) > 0){
			swap = (idx << 1) + 2;
		}
		if(swap == idx) break;
		pp = get_u32list(heap->ptrs, idx);
		set_u32list(heap->ptrs,  idx, get_u32list(heap->ptrs, swap));
		set_u32list(heap->ptrs, swap, pp);
		idx = swap;
	}
	return ret;
}

#endif
