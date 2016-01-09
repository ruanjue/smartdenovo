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
 
#ifndef __LIST_RJ_H
#define __LIST_RJ_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "sort.h"
#include "mem_share.h"

/**
 * * Common staic functions
 * */

#define num_min(n1, n2) (((n1) < (n2))? (n1) : (n2))
#define num_max(n1, n2) (((n1) > (n2))? (n1) : (n2))
#define num_diff(n1, n2) (((n1) < (n2))? ((n2) - (n1)) : ((n1) - (n2)))
#define num_cmp(a, b) (((a) > (b))? 1 : (((a) < (b))? -1 : -0))
#define num_cmpgt(a, b) ((a) > (b))
#define num_cmpx(a, b, c, d) (((a) > (b))? 1 : (((a) < (b))? -1 : (((c) > (d))? 1 : (((c) < (d))? -1 : 0))))
#define num_cmpgtx(a, b, c, d) (((a) > (b))? 1 : (((a) < (b))? 0 : (((c) > (d)))))
#define num_abs(n) ((n) < 0? -(n) : (n))

static inline size_t roundup_power2(size_t v){
	if(v == 0) return 0;
	v --;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

static inline size_t encap_list(void **buffer, size_t e_size, size_t size, size_t cur_cap, size_t inc, int mem_zeros){
	void *ptr;
	size_t cap;
	if(size + inc <= cur_cap) return cur_cap;
	if(size + inc < 0xFFFFFFF){
		cap = roundup_power2(size + inc);
	} else {
		cap = ((size + inc + 0xFFFFFFFLLU - 1LLU) / 0xFFFFFFFLLU) * 0xFFFFFFFLLU;
	}
	ptr = realloc(*buffer, e_size * cap);
	if(ptr == NULL){
		fprintf(stderr, " -- Out of memory, try to allocate %llu bytes, old size %llu, old addr %p in %s -- %s:%d --\n", (unsigned long long)(e_size * cap), (unsigned long long)(e_size * cur_cap), *buffer, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	*buffer = ptr;
	if(mem_zeros) memset(ptr + cur_cap, 0, e_size * (cap - cur_cap));
	return cap;
}

/**
 * Heap macros
 */

//ary, size and cap must be explict variable, not expression
#define array_heap_push(ary, size, cap, e_type, id, cmp_expr)\
do {	\
	e_type pp, p;	\
	p = (e_type)(id);	\
	size_t i, j, a, b;	\
	i = (size);	\
	if((size_t)((size) + 1) > (size_t)(cap)){	\
		if((size_t)((size) + 1) < 0xFFFFFFFFU){	\
			(cap) = roundup_power2((size) + 1);	\
		} else {	\
			(cap) = (((size) + 1 + 0xFFFFFFFLLU - 1LLU) / 0xFFFFFFFLLU) * 0xFFFFFFFLLU;	\
		}	\
		(ary) = realloc((ary), sizeof(e_type) * (cap));	\
		if((ary) == NULL){	\
			fprintf(stderr, " -- Out of memory, try to allocate %llu bytes in %s, -- %s:%d --\n", (unsigned long long)(sizeof(e_type) * (cap)), __FUNCTION__, __FILE__, __LINE__);	\
			exit(1);	\
		}	\
	}	\
	(ary)[(size)++] = p;	\
	while(i){	\
		j = (i - 1) >> 1;	\
		a = (ary)[i]; b = (ary)[j];	\
		if((cmp_expr) >= 0) break;	\
		pp = (ary)[i]; (ary)[i] = (ary)[j]; (ary)[j] = pp;	\
		i = j;	\
	}	\
} while(0)

#define array_heap_remove(ary, len, cap, e_type, _idx, cmp_expr)\
do {	\
	e_type pp;	\
	size_t swap, idx, a, b;	\
	idx = (size_t)(_idx);	\
	(ary)[idx] = (ary)[--(len)];	\
	while((size_t)((idx << 1) + 1) < (size_t)(len)){	\
		swap = idx;	\
		a = (ary)[swap]; b = (ary)[(idx << 1) + 1];	\
		if((cmp_expr) > 0) swap = (idx << 1) + 1;	\
		if(((idx << 1) + 2) < (size_t)(len)){	\
			a = (ary)[swap]; b = (ary)[(idx << 1) + 2];	\
			if((cmp_expr) > 0) swap = (idx << 1) + 2;	\
		}	\
		if(swap == idx) break;	\
		pp = (ary)[idx]; (ary)[idx] = (ary)[swap]; (ary)[swap] = pp;	\
		idx = swap;	\
	}	\
} while(0)

#define array_heap_replace(ary, len, cap, e_type, _idx, _val, cmp_expr)\
do {	\
	e_type pp;	\
	size_t swap, idx, a, b;	\
	idx = (size_t)(_idx);	\
	(ary)[idx] = _val;	\
	while((size_t)((idx << 1) + 1) < (size_t)(len)){	\
		swap = idx;	\
		a = (ary)[swap]; b = (ary)[(idx << 1) + 1];	\
		if((cmp_expr) > 0) swap = (idx << 1) + 1;	\
		if(((idx << 1) + 2) < (size_t)(len)){	\
			a = (ary)[swap]; b = (ary)[(idx << 1) + 2];	\
			if((cmp_expr) > 0) swap = (idx << 1) + 2;	\
		}	\
		if(swap == idx) break;	\
		pp = (ary)[idx]; (ary)[idx] = (ary)[swap]; (ary)[swap] = pp;	\
		idx = swap;	\
	}	\
} while(0)

#define array_heap_pop(ary, len, cap, e_type, cmp_expr)	\
({	\
	size_t ret;	\
	if(len){ ret = (size_t)ary[0]; array_heap_remove(ary, len, cap, e_type, 0, cmp_expr); }	\
	else ret = (size_t)0xFFFFFFFFFFFFFFFFLLU;	\
	ret;	\
})

/**
 * List
 */

#define define_list_core(list_type, e_type, size_type, inc_size)	\
	\
typedef struct { e_type* buffer; size_type size; size_type cap; int mem_zero; } list_type;	\
	\
static inline size_t list_type##_obj_desc_cnt(void *list, int idx){	\
	if(idx == 0) return ((list_type*)list)->size * sizeof(e_type);	\
	else return 1;	\
}	\
	\
static const obj_desc_t list_type##_obj_desc = {.size = sizeof(list_type), .n_child = 1, .mem_type = {1}, .addr = {offsetof(list_type, buffer)}, .desc = {(struct obj_desc_t*)&OBJ_DESC_DATA}, .cnt = list_type##_obj_desc_cnt};	\
	\
static inline list_type* init_##list_type(size_type init_size){	\
	if(init_size == 0) init_size = 2;	\
	list_type *list = (list_type*)malloc(sizeof(list_type));	\
	list->size = 0;	\
	list->cap  = init_size;	\
	list->buffer = (e_type*)calloc(list->cap, sizeof(e_type));	\
	list->mem_zero = 0;	\
	return list;	\
}	\
	\
static inline void list_type##_init(list_type *list, size_type init_size){	\
	if(init_size == 0) init_size = 2;	\
	list->size = 0;	\
	list->cap  = init_size;	\
	list->buffer = (e_type*)calloc(list->cap, sizeof(e_type));	\
	list->mem_zero = 0;	\
}	\
	\
static inline size_type count_##list_type(list_type *list){ return list->size; }	\
	\
static inline void clear_##list_type(list_type *list){ list->size = 0; }	\
	\
static inline void zeros_##list_type(list_type *list){ memset(list->buffer, 0, list->cap * sizeof(e_type)); }	\
	\
static inline void encap_##list_type(list_type *list, size_type n){	\
	list->cap = encap_list((void**)&list->buffer, sizeof(e_type), list->size, list->cap, n, list->mem_zero);	\
}	\
	\
static inline void recap_##list_type(list_type *list, size_type n){	\
	if((size_t)n == (size_t)list->cap) return;	\
	list->cap = n;	\
	if(list->size > n) list->size = n;	\
	list->buffer = realloc(list->buffer, list->cap * sizeof(e_type));	\
}	\
	\
static inline void encap_and_zeros_##list_type(list_type *list, size_type n){	\
	size_t old;	\
	if(list->size + ((size_t)n) <= (size_t)list->cap) return;	\
	if((size_type)(list->size + n) < list->size + ((size_t)n)){	\
		fprintf(stderr, " -- elements size exceed %s's data type %s in %s -- %s:%d --\n", #list_type, #size_type, __FUNCTION__, __FILE__, __LINE__);	\
		fflush(stderr);	\
		exit(1);	\
	}	\
	old = list->cap;	\
	while(list->size + n > list->cap){	\
		if(list->cap < inc_size){	\
			list->cap <<= 1;	\
		} else {	\
			list->cap += inc_size;	\
		}	\
	}	\
	list->buffer = realloc(list->buffer, list->cap * sizeof(e_type));	\
	memset(list->buffer + old, 0, (list->cap - old) * sizeof(e_type));	\
}	\
	\
static inline void clear_and_encap_##list_type(list_type *list, size_type n){	\
	list->size = 0;	\
	encap_##list_type(list, n);	\
}	\
	\
static inline void encap_and_inc_##list_type(list_type *list, size_type n){	\
	encap_##list_type(list, n);	\
	list->size += n;	\
}	\
	\
static inline void trunc_##list_type(list_type *list, size_type size){	\
	if(size > count_##list_type(list)) size = count_##list_type(list);	\
	list->size -= size;	\
}	\
	\
static inline void set_##list_type##_size(list_type *list, size_type size){ list->size = size; }	\
	\
static inline void inc_##list_type(list_type *list, size_type size){	\
	if(size + list->size > list->cap) list->size = list->cap;	\
	else list->size += size;	\
}	\
	\
static inline void lazy_push_##list_type(list_type *list, e_type e){ list->buffer[list->size ++] = e; }	\
	\
static inline void push_##list_type(list_type *list, e_type e){	\
	encap_##list_type(list, 1);	\
	list->buffer[list->size++] = e;	\
}	\
	\
static inline int fpush_##list_type(list_type *list, FILE *inp){	\
	int ret;	\
	encap_##list_type(list, 1);	\
	ret = fread(list->buffer + list->size, sizeof(e_type), 1, inp);	\
	if(ret == 1) list->size ++;	\
	return ret;	\
}	\
	\
static inline e_type* peer_##list_type(list_type *list){	\
	if(count_##list_type(list)){	\
		return list->buffer + list->size - 1;	\
	} else return NULL;	\
}	\
	\
static inline int pop_##list_type(list_type *list, e_type*e){	\
	if(count_##list_type(list)){	\
		list->size --;	\
		*e = list->buffer[list->size];	\
		return 1;	\
	} else return 0;	\
}	\
	\
static inline int fpop_##list_type(list_type *list, FILE *oup){	\
	if(list->size){	\
		fwrite(list->buffer + list->size - 1, sizeof(e_type), 1, oup);	\
		list->size --;	\
		return 1;	\
	} else {	\
		return 0;	\
	}	\
}	\
	\
static inline void insert_##list_type(list_type *list, size_type idx, e_type e){	\
	if(idx > list->size) return;	\
	encap_##list_type(list, 1);	\
	if(idx == list->size){	\
		list->buffer[list->size] = e;	\
	} else {	\
		memmove(list->buffer + idx + 1, list->buffer + idx, (list->size - idx) * sizeof(e_type));	\
		list->buffer[idx] = e;	\
	}	\
	list->size ++;	\
}	\
	\
static inline void insert_array_##list_type(list_type *list, size_type idx, e_type *es, size_type size){	\
	if(idx > list->size) return;	\
	encap_##list_type(list, size);	\
	if(idx == list->size){	\
	} else {	\
		memmove(list->buffer + idx + size, list->buffer + idx, (list->size - idx) * sizeof(e_type));	\
	}	\
	memcpy(list->buffer + idx, es, size * sizeof(e_type));	\
	list->size += size;	\
}	\
	\
static inline void remove_##list_type(list_type *list, size_type idx){	\
	if(idx >= list->size) return;	\
	if(idx + 1 < list->size){	\
		memmove(list->buffer + idx, list->buffer + idx + 1, (list->size - idx - 1) * sizeof(e_type));	\
	}	\
	list->size --;	\
}	\
	\
static inline void remove_array_##list_type(list_type *list, size_type off, size_type len){	\
	if(off >= list->size) return;	\
	if(off + len < list->size){	\
		memmove(list->buffer + off, list->buffer + off + len, (list->size - off - len) * sizeof(e_type));	\
		list->size -= len;	\
	} else { \
		list->size = off;	\
	}	\
}	\
	\
static inline void set_##list_type(list_type *list, size_type idx, e_type e){ list->buffer[idx] = e; }	\
	\
static inline e_type get_##list_type(list_type *list, size_type idx){ return list->buffer[idx]; }	\
	\
static inline e_type* ref_##list_type(list_type *list, size_type idx){ return list->buffer + idx; }	\
	\
static inline e_type* next_ref_##list_type(list_type *list){ encap_##list_type(list, 1); list->size ++; return list->buffer + list->size - 1; }	\
	\
static inline e_type* ref_next_##list_type(list_type *list){ list->size ++; return list->buffer + list->size - 1; }	\
	\
static inline e_type* as_array_##list_type(list_type *list){ return list->buffer; }	\
	\
static inline void reverse_##list_type(list_type *list){	\
	size_type i, j;	\
	e_type t;	\
	if(count_##list_type(list) == 0) return;	\
	i = 0;	\
	j = count_##list_type(list) - 1;	\
	while(i < j){	\
		t = get_##list_type(list, i);	\
		set_##list_type(list, i, get_##list_type(list, j));	\
		set_##list_type(list, j, t);	\
		i ++;	\
		j --;	\
	}	\
}	\
	\
static inline void sub_reverse_##list_type(list_type *list, size_type beg, size_type end){	\
	size_type i, j;	\
	e_type t;	\
	if(end == 0) return;	\
	i = beg;	\
	j = end - 1;	\
	while(i < j){	\
		t = get_##list_type(list, i);	\
		set_##list_type(list, i, get_##list_type(list, j));	\
		set_##list_type(list, j, t);	\
		i ++;	\
		j --;	\
	}	\
}	\
	\
static inline void append_##list_type(list_type *list1, list_type *list2){	\
	encap_##list_type(list1, count_##list_type(list2));	\
	memcpy(list1->buffer + list1->size, list2->buffer, sizeof(e_type) * list2->size);	\
	list1->size += list2->size;	\
}	\
	\
static inline void append_array_##list_type(list_type *list1, e_type *ary, size_type size){	\
	if(size == 0) return;	\
	encap_##list_type(list1, size);	\
	memcpy(list1->buffer + list1->size, ary, sizeof(e_type) * size);	\
	list1->size += size;	\
}	\
	\
static inline size_type write_##list_type(list_type *list, FILE *out){	\
	return fwrite(list->buffer, sizeof(e_type), count_##list_type(list), out);	\
}	\
	\
static inline size_type dump_##list_type(list_type *list, FILE *out){	\
	fwrite(&list->size, sizeof(size_type), 1, out);	\
	fwrite(list->buffer, sizeof(e_type), list->size, out);	\
	return sizeof(size_type) + sizeof(e_type) * list->size;	\
}	\
	\
static inline size_t mem_size_##list_type(list_type *list){	\
	return ((sizeof(list_type) + 7) / 8) * 8 + (list->size * sizeof(e_type) + 7) / 8 * 8;	\
}	\
	\
static inline size_t mem_dump_##list_type(list_type *list, FILE *out, void *addr){	\
	list_type clone;	\
	size_t size;	\
	int8_t i, v;	\
	size = mem_size_##list_type(list);	\
	clone.size = list->size;	\
	clone.cap  = list->size;	\
	clone.buffer = addr + (sizeof(list_type) + 7) / 8 * 8;	\
	fwrite(&clone, sizeof(list_type), 1, out);	\
	v = 0;	\
	for(i=(sizeof(list_type) + 7) / 8 * 8 - sizeof(list_type);i>0;i--) fwrite(&v, 1, 1, out);	\
	fwrite(list->buffer, sizeof(e_type), list->size, out);	\
	for(i=(sizeof(e_type) * list->size + 7) / 8 * 8 - sizeof(e_type) * list->size;i>0;i--) fwrite(&v, 1, 1, out);	\
	return size;	\
}	\
	\
static inline list_type* load_##list_type(FILE *inp){	\
	list_type *list;	\
	size_type n;	\
	list = (list_type*)malloc(sizeof(list_type));	\
	if((n = fread(&list->size, sizeof(size_type), 1, inp)) != 1){	\
		free(list); return NULL;	\
	}	\
	list->cap = list->size;	\
	list->buffer = (e_type*)malloc(sizeof(e_type) * list->cap);	\
	if(list->buffer == NULL){	\
		fprintf(stderr, " Out of memory in load_%s \n", #list_type); fflush(stderr); exit(1);	\
	}	\
	if((n = fread(list->buffer, sizeof(e_type), list->size, inp)) != list->size){	\
		free(list->buffer); free(list); return NULL;	\
	}	\
	return list;	\
}	\
	\
static inline void free_##list_type(list_type *list){ free(list->buffer); free(list); }	\
	\
static inline void list_type##_free(list_type *list){ free(list->buffer); list->buffer = NULL; }	\

#define define_list_ext(list_type, e_type, size_type, cmp_func)	\
static inline size_type delete_##list_type(list_type *list, e_type e){	\
	size_type i, ret;	\
	ret = 0;	\
	for(i=list->size;i>0;i--){	\
		if(cmp_func(list->buffer[i-1], e, NULL) == 0){	\
			if(i < list->size){	\
				memmove(list->buffer + i - 1, list->buffer + i, (list->size - i) * sizeof(e_type));	\
			}	\
			list->size --;	\
			ret ++;	\
		}	\
	}	\
	return ret;	\
}	\
	\
static inline size_type occ_##list_type(list_type *list, e_type e){	\
	size_type i, n;	\
	for(i=0,n=0;i<list->size;i++){	\
		if(cmp_func(list->buffer[i], e, NULL) == 0) n++;	\
	}	\
	return n;	\
}	\
	\
static inline size_type replace_##list_type(list_type *list, e_type from, e_type to){	\
	size_type i, ret;	\
	ret = 0;	\
	for(i=0;i<list->size;i++){	\
		if(cmp_func(list->buffer[i], from, NULL) == 0){	\
			list->buffer[i] = to;	\
			ret ++;	\
		}	\
	}	\
	return ret;	\
}	\
	\
static inline size_type locate_##list_type(list_type *list, e_type e, size_type start){	\
	size_type i;	\
	for(i=start;i<list->size;i++){	\
		if(cmp_func(list->buffer[i], e, NULL) == 0) return i;	\
	}	\
	return i;	\
}	\
	\
define_quick_sort(sort_##list_type##_core, e_type, cmp_func);	\
	\
static inline void sort_##list_type(list_type *list){ sort_##list_type##_core(ref_##list_type(list, 0), count_##list_type(list), NULL); }

#define define_list(name, e_type) define_list_core(name, e_type, size_t, 0xFFFFFU)

#define native_number_cmp(e1, e2, obj) (((e1) == (e2))? 0 : (((e1) < (e2))? -1 : 1))

#define define_native_list(name, e_type)	\
define_list_core(name, e_type, size_t, 0xFFFFFU);	\
define_list_ext(name, e_type, size_t, native_number_cmp);

define_native_list(u8list,  uint8_t);
define_native_list(u1v, uint8_t);
define_native_list(u16list, uint16_t);
define_native_list(u2v, uint16_t);
define_native_list(u32list, uint32_t);
define_native_list(u4v, uint32_t);
define_native_list(u64list, uint64_t);
define_native_list(u8v, uint64_t);

define_native_list(b8list,  int8_t);
define_native_list(b1v, int8_t);
define_native_list(b16list, int16_t);
define_native_list(b2v, int16_t);
define_native_list(b32list, int32_t);
define_native_list(b4v, int32_t);
define_native_list(b64list, int64_t);
define_native_list(b8v, int64_t);

define_native_list(f4v, float);
define_native_list(f8v, long double);

define_list(vplist, void*);
define_list(cplist, char*);
// mem_share for deep dump of cplist
static inline size_t cplist_deep_obj_desc_cnt(void *list, int idx){
	if(idx == 0) return ((cplist*)list)->size;
	else return 1;
}
static const obj_desc_t cplist_deep_obj_desc = {.size = sizeof(cplist), .n_child = 1, .mem_type = {3}, .addr = {offsetof(cplist, buffer)}, .desc = {(struct obj_desc_t*)&OBJ_DESC_CHAR_ARRAY}, .cnt = cplist_deep_obj_desc_cnt};

#define define_recycle_list_array(name, list_type)	\
typedef struct {	\
	vplist *array;	\
	vplist *dustbin;	\
} name;	\
	\
static inline name* init_##name(){	\
	name *g;	\
	g = (name*)malloc(sizeof(name));	\
	g->buffer = init_vplist(4);	\
	g->dustbin = init_vplist(4);	\
	return g;	\
}	\
	\
static inline void free_##name(name *g){	\
	list_type *v;	\
	size_t i;	\
	for(i=0;i<g->array->size;i++){	\
		v = (list_type*)get_vplist(g->array, i);	\
		if(v) free_##list_type(v);	\
	}	\
	for(i=0;i<g->dustbin->size;i++){	\
		v = (list_type*)get_vplist(g->dustbin, i);	\
		if(v) free_##list_type(v);	\
	}	\
	free_vplist(g->array);	\
	free_vplist(g->dustbin);	\
	free(g);	\
}	\
	\
static inline list_type* fetch_##name(name *g){	\
	list_type *v;	\
	if(g->dustbin->size) v = (list_type*)g->dustbin->buffer[--g->dustbin->size];	\
	else v = init_##list_type(4);	\
	return v;	\
}	\
	\
static inline void recyc_##name(name *g, list_type *v){	\
	push_vplist(g->dustbin, v);	\
}	\
	\
static inline void recyc_all_##name(name *g, vplist *vs){	\
	append_vplist(g->dustbin, vs);	\
	vs->size = 0;	\
}

#endif
