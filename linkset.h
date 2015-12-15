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
 
#ifndef __LINKED_SET_RJ_H
#define __LINKED_SET_RJ_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "sort.h"
#include "mem_share.h"
#include "list.h"

#define define_linkset_core(list_type, e_type, size_type, inc_size)	\
	\
typedef struct { e_type* buffer; size_type* fwd_links; size_type* rev_links; size_type head, tail, size; size_type* trash; size_type trash_size, trash_cap; size_type off; size_type cap; } list_type;	\
	\
static inline size_t list_type##_obj_desc_cnt(void *list, int idx){	\
	if(idx == 0) return ((list_type*)list)->size * sizeof(e_type);	\
	else if(idx == 1) return ((list_type*)list)->size * sizeof(size_type);	\
	else if(idx == 2) return ((list_type*)list)->size * sizeof(size_type);	\
	else if(idx == 3) return ((list_type*)list)->trash_size * sizeof(size_type);	\
	else return 1;	\
}	\
	\
static const obj_desc_t list_type##_obj_desc = {.size = sizeof(list_type), .n_child = 4, .mem_type = {1, 1, 1, 1}, .addr = {offsetof(list_type, buffer), offsetof(list_type, fwd_links), offsetof(list_type, rev_links), offsetof(list_type, trash)}, .desc = {(struct obj_desc_t*)&OBJ_DESC_DATA, (struct obj_desc_t*)&OBJ_DESC_DATA, (struct obj_desc_t*)&OBJ_DESC_DATA, (struct obj_desc_t*)&OBJ_DESC_DATA}, .cnt = list_type##_obj_desc_cnt};	\
	\
static inline list_type* init_##list_type(){	\
	list_type *list = (list_type*)malloc(sizeof(list_type));	\
	list->off = 0;	\
	list->cap = 0;	\
	list->buffer = NULL;	\
	list->fwd_links = NULL;	\
	list->rev_links = NULL;	\
	list->trash  = NULL;	\
	list->trash_size = 0;	\
	list->trash_cap = 0;	\
	list->head = 0;	\
	list->tail = 0;	\
	list->size = 0;	\
	return list;	\
}	\
	\
static inline size_type count_##list_type(list_type *list){ return list->size; }	\
	\
static inline void clear_##list_type(list_type *list){ list->trash_size = 0; list->off = 0; list->head = 0; list->tail = 0; list->size = 0; }	\
	\
static inline void push_##list_type(list_type *list, e_type e){	\
	size_type idx;	\
	if(list->trash_size){	\
		idx = list->trash[--list->trash_size];	\
	} else {	\
		encap_list((void**)&(list->buffer), sizeof(e_type), list->off, list->cap, 1, 0);	\
		encap_list((void**)&(list->fwd_links), sizeof(size_type), list->off, list->cap, 1, 0);	\
		list->cap = encap_list((void**)&(list->rev_links), sizeof(size_type), list->off, list->cap, 1, 0);	\
		idx = list->off ++;	\
	}	\
	list->buffer[idx] = e;	\
	list->fwd_links[idx] = 0;	\
	if(list->size == 0){	\
		list->head = idx;	\
		list->rev_links[idx] = 0;	\
	} else {	\
		list->fwd_links[list->tail] = idx;	\
		list->rev_links[idx] = list->tail;	\
	}	\
	list->tail = idx;	\
	list->size ++;	\
}	\
	\
static inline int pop_##list_type(list_type *list, e_type*e){	\
	if(list->size == 0) return 0;	\
	*e = list->buffer[list->tail];	\
	list->trash_cap = encap_list((void**)&list->trash, sizeof(size_type), list->trash_size, list->trash_cap, 1, 0);	\
	list->trash[list->trash_size ++] = list->head;	\
	list->tail = list->rev_links[list->tail];	\
	list->size --;	\
	return 1;	\
}	\
	\
static inline void unshift_##list_type(list_type *list, e_type e){	\
	size_type idx;	\
	if(list->trash_size){	\
		idx = list->trash[--list->trash_size];	\
	} else {	\
		encap_list((void**)&(list->buffer), sizeof(e_type), list->off, list->cap, 1, 0);	\
		list->cap = encap_list((void**)&(list->links), sizeof(size_type), list->off, list->cap, 1, 0);	\
		idx = list->off ++;	\
	}	\
	list->buffer[idx] = e;	\
	list->rev_links[idx] = 0;	\
	if(list->size == 0){	\
		list->fwd_links[idx] = 0;	\
		list->tail = idx;	\
	} else {	\
		list->fwd_links[idx] = list->head;	\
		list->rev_links[list->head] = idx;	\
	}	\
	list->head = idx;	\
	list->size ++;	\
}	\
	\
static inline void shift_##list_type(list_type *list, e_type *e){	\
	if(list->size == 0) return 0;	\
	*e = list->buffer[list->head];	\
	list->trash_cap = encap_list((void**)&list->trash, sizeof(size_type), list->trash_size, list->trash_cap, 1, 0);	\
	list->trash[list->trash_size ++] = list->head;	\
	list->head = list->fwd_links[list->head];	\
	list->size --;	\
	return 1;	\
}	\
	\
	\
static inline void free_##list_type(list_type *list){	\
	if(list->buffer) free(list->buffer);	\
	if(list->fwd_links) free(list->fwd_links);	\
	if(list->rev_links) free(list->rev_links);	\
	if(list->trash) free(list->trash);	\
	free(list);	\
}

#define define_linkset(name, e_type) define_linkset_core(name, e_type, size_t, 0xFFFFFU)

define_linkset(u32linkv, uint32_t);
define_linkset(b32linkv, int32_t);

#endif
