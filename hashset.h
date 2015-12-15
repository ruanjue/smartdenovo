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
 
#ifndef __HASH_SET_RJ
#define __HASH_SET_RJ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "mem_share.h"

static const uint64_t sys_prime_list[61] = {
	0x7LLU, 0xfLLU, 0x1fLLU, 0x43LLU, 0x89LLU,
	0x115LLU, 0x22dLLU, 0x45dLLU, 0x8bdLLU, 0x1181LLU,
	0x2303LLU, 0x4609LLU, 0x8c17LLU, 0x1183dLLU, 0x2307bLLU,
	0x460fdLLU, 0x8c201LLU, 0x118411LLU, 0x230833LLU, 0x461069LLU,
	0x8c20e1LLU, 0x11841cbLLU, 0x2308397LLU, 0x461075bLLU, 0x8c20ecbLLU,
	0x11841da5LLU, 0x23083b61LLU, 0x461076c7LLU, 0x8c20ed91LLU, 0x11841db31LLU,
	0x23083b673LLU, 0x461076d1bLLU, 0x8c20eda41LLU, 0x11841db48dLLU, 0x23083b6937LLU,
	0x461076d27fLLU, 0x8c20eda50dLLU, 0x11841db4a59LLU, 0x23083b694ebLLU, 0x461076d29f1LLU,
	0x8c20eda5441LLU, 0x11841db4a887LLU, 0x23083b69511fLLU, 0x461076d2a2c1LLU, 0x8c20eda54591LLU,
	0x11841db4a8b55LLU, 0x23083b69516c1LLU, 0x461076d2a2da5LLU, 0x8c20eda545b55LLU, 0x11841db4a8b6b5LLU,
	0x23083b69516d91LLU, 0x461076d2a2db3bLLU, 0x8c20eda545b69dLLU, 0x11841db4a8b6d5dLLU, 0x23083b69516daf5LLU,
	0x461076d2a2db5edLLU, 0x8c20eda545b6c5fLLU, 0x11841db4a8b6d8ebLLU, 0x23083b69516db1ffLLU, 0x461076d2a2db643fLLU,
	0x8c20eda545b6c8f3LLU
};

static inline uint64_t _rj_hashset_find_prime(uint64_t n){
	uint32_t i;
	i = 0;
	while(i < 60 && n > sys_prime_list[i]) i ++;
	return sys_prime_list[i];
}

#ifndef HASH_FLAG_MACROS
#define HASH_FLAG_MACROS
#define is_entity_null(flags, idx)    ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x01)
#define is_entity_del(flags, idx)     ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x02)
#define exists_entity(flags, idx)     (!((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x03))
#define set_entity_null(flags, idx)   ((flags)[(idx)>>4] |= (0x01u<<(((idx)&0x0f)<<1)))
#define set_entity_del(flags, idx)    ((flags)[(idx)>>4] |= (0x02u<<(((idx)&0x0f)<<1)))
#define clear_entity_null(flags, idx) ((flags)[(idx)>>4] &= ~(0x01u<<(((idx)&0x0f)<<1)))
#define clear_entity_del(flags, idx)  ((flags)[(idx)>>4] &= ~(0x02u<<(((idx)&0x0f)<<1)))
#endif

#define init_hashset_macro(hash_type, hash_key_type) \
typedef struct { hash_key_type *array;  uint32_t *flags; size_t e_size; size_t ocp; size_t size; size_t count; size_t max; float load_factor; size_t iter_ptr; void *userdata; } hash_type; \
static inline size_t hash_type##_obj_desc_cnt(void *obj, int idx){	\
	switch(idx){	\
		case 0: return ((hash_type*)obj)->size * sizeof(hash_key_type);	\
		case 1: return (((hash_type*)obj)->size + 15) / 16 * 4;	\
		default: return 0;	\
	}	\
}	\
static const obj_desc_t hash_type##_obj_desc = {sizeof(hash_type), 2, {1, 1}, {offsetof(hash_type, array), offsetof(hash_type, flags)}, {(obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&OBJ_DESC_DATA}, hash_type##_obj_desc_cnt};	\
static inline void prepare_mem_dump_##hash_type(hash_type *hash){	\
	size_t i;	\
	for(i=0;i<hash->size;i++){	\
		if(!exists_entity(hash->flags, i)) memset(hash->array + i, 0, sizeof(hash_key_type));	\
	}	\
}	\
static inline int hash_type##_is_prime(uint64_t num){                          \
	uint64_t i, max;                                                           \
	if(num < 4) return 1;                                                      \
	if(num % 2 == 0) return 0;                                                 \
	max = (uint64_t)sqrt((double)num);                                         \
	for(i=3;i<max;i+=2){ if(num % i == 0) return 0; }                          \
	return 1;                                                                  \
}                                                                              \
static inline uint64_t hash_type##_find_next_prime(uint64_t num){              \
	if(num % 2 == 0) num ++;                                                   \
	while(1){ if(hash_type##_is_prime(num)) return num; num += 2; }            \
}                                                                              \
static inline hash_type* init2_##hash_type(uint32_t size, float factor){       \
	hash_type *set;                                                            \
	set = (hash_type*)calloc(1, sizeof(hash_type));                            \
	set->e_size = sizeof(hash_key_type);                                       \
	set->size   = _rj_hashset_find_prime(size);                                \
	set->count  = 0;                                                           \
	set->ocp    = 0;                                                           \
	set->load_factor = factor;                                                 \
	set->max    = set->size * set->load_factor;                                \
	set->iter_ptr    = 0;                                                      \
	set->array       = calloc(set->size, set->e_size);                         \
	set->flags       = malloc((set->size + 15)/16 * 4);                        \
	memset(set->flags, 0x55, (set->size + 15) / 16 * 4);                       \
	set->userdata = NULL;                                                      \
	return set;                                                                \
}                                                                              \
static inline void set_userdata_##hash_type(hash_type *set, void *userdata){ set->userdata = userdata; } \
static inline hash_type* init_##hash_type(uint32_t size){ return init2_##hash_type(size, 0.67f); }

#define get_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline hash_key_type* get_##hash_type(hash_type *set, hash_key_type key){\
	hash_key_type *e;                                                          \
	uint32_t flag;                                                             \
	size_t hc;                                                                 \
	hc = hash_code_macro(key) % set->size;                                     \
	while(1){                                                                  \
		flag = (set->flags[hc >> 4] >> (((hc) & 0x0f) << 1)) & 0x03;           \
		if(flag & 0x01){                                                       \
			return NULL;                                                       \
		} else if(flag & 0x02){                                                \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)) return e;                            \
		}                                                                      \
		if(hc + 1 == set->size) hc = 0; else hc ++;                            \
	}                                                                          \
	return NULL;                                                               \
}                                                                              \
static inline size_t offset_##hash_type(hash_type *set, hash_key_type *ptr){   \
	return ptr - set->array;                                                   \
}

#define prepare_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline void encap_##hash_type(hash_type *set, size_t num);              \
static inline hash_key_type* prepare_##hash_type(hash_type *set, hash_key_type key, int *exists){\
	hash_key_type *e;                                                          \
	size_t hc, d;                                                              \
	encap_##hash_type(set, 1);                                                 \
	hc = hash_code_macro((key)) % set->size;                                   \
	d = set->size;                                                             \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			if(d == set->size){                                                \
				clear_entity_null(set->flags, hc);                             \
				set->ocp ++;                                                   \
			} else {                                                           \
				hc = d;                                                        \
				clear_entity_del(set->flags, hc);                              \
			}                                                                  \
			*exists = 0;                                                       \
			set->count ++;                                                     \
			e = ((hash_key_type*)set->array) + hc;                             \
			return e;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
			if(d == set->size) d = hc;                                         \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro((*e), (key))){                                \
				*exists = 1;                                                   \
				return e;                                                      \
			}                                                                  \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return NULL;                                                               \
}

#define exists_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline int exists_##hash_type(hash_type *set, hash_key_type key){       \
	hash_key_type *e;                                                          \
	size_t hc;                                                                 \
	hc = hash_code_macro(key) % set->size;                                     \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			return 0;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)) return 1;                            \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return 0;                                                                  \
}

#define add_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline hash_key_type* add_##hash_type(hash_type *set, hash_key_type key){       \
	hash_key_type *e;                                                          \
	size_t d, hc;                                                              \
	hc = hash_code_macro(key) % set->size;                                     \
	d  = set->size;                                                            \
	do{                                                                        \
		if(is_entity_null(set->flags, hc)){                                    \
			if(d == set->size){                                                \
				d = hc;                                                        \
				clear_entity_null(set->flags, d);                              \
				set->ocp ++;                                                   \
			} else {                                                           \
				clear_entity_del(set->flags, d);                               \
			}                                                                  \
			e = ((hash_key_type*)set->array) + d;                              \
			*e = key;                                                          \
			set->count ++;                                                     \
			return e;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
			if(d == set->size) d = hc;                                         \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)){                                     \
				*e = key;                                                      \
				return e;                                                      \
			}                                                                  \
		}                                                                      \
		if(hc + 1 == set->size) hc = 0;                                        \
		else hc = hc + 1;                                                      \
	} while(1);                                                                \
	return NULL;                                                                  \
}

#define put_hashset_macro(hash_type, hash_key_type) \
static inline hash_key_type* put_##hash_type(hash_type *set, hash_key_type key){         \
	encap_##hash_type(set, 1);                                                 \
	return add_##hash_type(set, key);                                          \
}

#define remove_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline void delete_##hash_type(hash_type *set, hash_key_type *key){ set_entity_del(set->flags, key - set->array); set->count --; }   \
static inline int remove_##hash_type(hash_type *set, hash_key_type key){       \
	hash_key_type *e;                                                          \
	size_t hc;                                                                 \
	hc = hash_code_macro(key) % set->size;                                     \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			return 0;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)){                                     \
				set->count --;                                                 \
				set_entity_del(set->flags, hc);                                \
				return 1;                                                      \
			}                                                                  \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return 0;                                                                  \
}

#define reset_iter_hashset_macro(hash_type) static inline void reset_iter_##hash_type(hash_type *set){ set->iter_ptr = 0; }

#define iter_hashset_macro(hash_type, hash_key_type) \
static inline int iter_##hash_type(hash_type *set, hash_key_type *ret){        \
	if(set->iter_ptr >= set->size) return 0;                                   \
	while(set->iter_ptr < set->size){                                          \
		if(exists_entity(set->flags, set->iter_ptr)){                          \
			*ret = *(((hash_key_type*)set->array) + set->iter_ptr);            \
			set->iter_ptr ++;                                                  \
			return 1;                                                          \
		}                                                                      \
		set->iter_ptr ++;                                                      \
	}                                                                          \
	return 0;                                                                  \
}

#define ref_iter_hashset_macro(hash_type, hash_key_type) \
static inline hash_key_type* ref_iter_##hash_type(hash_type *set){             \
	if(set->iter_ptr >= set->size) return NULL;                                \
	while(set->iter_ptr < set->size){                                          \
		if(exists_entity(set->flags, set->iter_ptr)){                          \
			return (((hash_key_type*)set->array) + set->iter_ptr++);           \
		}                                                                      \
		set->iter_ptr ++;                                                      \
	}                                                                          \
	return NULL;                                                               \
}

#define ref_rev_iter_hashset_macro(hash_type, hash_key_type) \
static inline hash_key_type* ref_rev_iter_##hash_type(hash_type *set){             \
	if(set->iter_ptr >= set->size) return NULL;                                \
	while(set->iter_ptr < set->size){                                          \
		if(!exists_entity(set->flags, set->iter_ptr)){                          \
			return (((hash_key_type*)set->array) + set->iter_ptr++);           \
		}                                                                      \
		set->iter_ptr ++;                                                      \
	}                                                                          \
	return NULL;                                                               \
}

#define count_hashset_macro(hash_type) static inline int64_t count_##hash_type(hash_type *set){ return set->count; }

#define clear_hashset_macro(hash_type) \
static inline void clear_##hash_type(hash_type *set){                          \
	if(set->ocp == 0) return;                                                  \
	memset(set->flags, 0x55, (set->size + 15) / 16 * 4);                       \
	set->count = 0;                                                            \
	set->ocp   = 0;                                                            \
	set->iter_ptr = 0;                                                         \
}

#define ffwrite(ptr, e_size, size, file) (e_size * fwrite(ptr, e_size, size, file))
#define ffread(ptr, e_size, size, file) (e_size * fread(ptr, e_size, size, file))

#define dump_hashset_macro(hash_type) \
static inline size_t sizeof_##hash_type(hash_type *set){                       \
	return sizeof(size_t) * 3 + sizeof(float) + set->e_size * set->size        \
				+ sizeof(uint32_t) * ((set->size + 15) / 16);                  \
}                                                                              \
static inline size_t dump_##hash_type(hash_type *set, FILE *out){              \
	size_t n;                                                          \
	n =  ffwrite(&set->e_size, sizeof(size_t), 1, out);                        \
	n += ffwrite(&set->size, sizeof(size_t), 1, out);                          \
	n += ffwrite(&set->count, sizeof(size_t), 1, out);                         \
	n += ffwrite(&set->load_factor, sizeof(float), 1, out);                    \
	n += ffwrite(set->array, set->e_size, set->size, out);	\
	n += ffwrite(set->flags, sizeof(uint32_t), (set->size + 15) / 16, out);	\
	return n;                                                                  \
}

#define load_hashset_macro(hash_type) \
static inline hash_type* load_##hash_type(FILE *in){                           \
	hash_type *set;                                                            \
	size_t n;                                                                  \
	set = (hash_type*)malloc(sizeof(hash_type));                               \
	n =  ffread(&set->e_size, sizeof(size_t), 1, in);                          \
	n += ffread(&set->size, sizeof(size_t), 1, in);                            \
	n += ffread(&set->count, sizeof(size_t), 1, in);                           \
	n += ffread(&set->load_factor, sizeof(float), 1, in);                      \
	set->max   = set->size * set->load_factor;                                 \
	set->array = malloc(set->size * set->e_size);                              \
	n += ffread(set->array, set->e_size, set->size, in);                       \
	set->flags = (uint32_t*)malloc((set->size + 15) / 16 * 4);                 \
	n += ffread(set->flags, sizeof(uint32_t), (set->size + 15) / 16, in);      \
	return set;                                                                \
}

#define free_hashset_macro(hash_type) \
static inline void free_##hash_type(hash_type *set){                           \
	free(set->array);                                                          \
	free(set->flags);                                                          \
	free(set);                                                                 \
}

#define encap_hashset_macro(hash_type, hash_key_type, hash_code_macro) \
static inline void encap_##hash_type(hash_type *set, size_t num){             \
	uint32_t *flags, *f;                                                      \
	uint64_t i, n, size, hc;                                                  \
	hash_key_type key;                                                        \
	hash_key_type tmp;                                                        \
	if(set->ocp + num <= set->max) return;                                  \
	n = set->size;                                                            \
	do{ n = _rj_hashset_find_prime(n * 2); } while(n * set->load_factor < set->count + num);    \
	set->array = realloc(set->array, n * set->e_size);                        \
	if(set->array == NULL){                                                   \
		fprintf(stderr, "-- Out of memory --\n");                             \
		exit(1);                                                              \
	}                                                                         \
	memset(set->array + set->size, 0, (n - set->size) * set->e_size);         \
	flags = malloc((n+15)/16 * 4);                                            \
	memset(flags, 0x55, (n+15)/16 * 4);                                       \
	size = set->size;                                                         \
	set->size = n;                                                            \
	set->ocp  = set->count;                                                   \
	set->max = n * set->load_factor;                                          \
	f = set->flags;                                                           \
	set->flags = flags;                                                       \
	flags = f;                                                                \
	for(i=0;i<size;i++){                                                      \
		if(!exists_entity(flags, i)) continue;                                \
		key = ((hash_key_type*)set->array)[i];                                \
		set_entity_del(flags, i);                                             \
		while(1){                                                             \
			hc = hash_code_macro(key) % set->size;                            \
			while(!is_entity_null(set->flags, hc)){ hc = (hc + 1) % set->size; }        \
			clear_entity_null(set->flags, hc);                                \
			if(hc < size && exists_entity(flags, hc)){                        \
				tmp = key;                                                    \
				key = ((hash_key_type*)set->array)[hc];                       \
				((hash_key_type*)set->array)[hc] = tmp;                       \
				set_entity_del(flags, hc);                                    \
			} else {                                                          \
				((hash_key_type*)set->array)[hc] = key;                       \
				break;                                                        \
			}                                                                 \
		}                                                                     \
	}                                                                         \
	free(flags);                                                              \
}                                                                             \
static inline size_t offsetof_##hash_type(hash_type *set, hash_key_type *ptr){ return ptr - set->array; }	\

#define define_key_val_hashset_ext(hash_type, hash_element_type, hash_key_type, hash_val_type) \
static inline int kv_exists_##hash_type(hash_type *set, hash_key_type key){	\
	hash_element_type e;	\
	memset(&e, 0xFFU, sizeof(hash_element_type));	\
	e.key = key;	\
	return exists_##hash_type(set, e);	\
}	\
	\
static inline hash_val_type kv_get_##hash_type(hash_type *set, hash_key_type key){	\
	hash_element_type e, *ee;	\
	memset(&e, 0xFFU, sizeof(hash_element_type));	\
	e.key = key;	\
	ee = get_##hash_type(set, e);	\
	if(ee == NULL) return e.val;	\
	else return ee->val;	\
}	\
	\
static inline void kv_put_##hash_type(hash_type *set, hash_key_type key, hash_val_type val){	\
	hash_element_type e;	\
	e.key = key;	\
	e.val = val;	\
	put_##hash_type(set, e);	\
}	\
	\
static inline hash_element_type* kv_prepare_##hash_type(hash_type *set, hash_key_type key, int *exists){	\
	hash_element_type e;	\
	memset(&e, 0xFFU, sizeof(hash_element_type));	\
	e.key = key;	\
	return prepare_##hash_type(set, e, exists);	\
}

// ---------------------- Define your own hashset ----------------------------------
// Example: 
// typedef struct { int group; int user; } Info;
// #define my_hashcode(val) (val)->group
// #define my_hashequal(v1, v2) (((v1)->group == (v2)->group) && ((v1)->user == (v2)->user))
// define_hashset(myhash, Info, my_hashcode, my_hashequal);

#define define_hashset(hash_type, hash_key_type, hash_code_macro, hash_equal_macro)    \
	init_hashset_macro(hash_type, hash_key_type);                              \
	get_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro);    \
	prepare_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro);\
	exists_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro); \
	add_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro);    \
	put_hashset_macro(hash_type, hash_key_type);                               \
	remove_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro); \
	iter_hashset_macro(hash_type, hash_key_type);                              \
	ref_iter_hashset_macro(hash_type, hash_key_type);                          \
	reset_iter_hashset_macro(hash_type);                                       \
	count_hashset_macro(hash_type);                                            \
	clear_hashset_macro(hash_type);                                            \
	dump_hashset_macro(hash_type);                                             \
	load_hashset_macro(hash_type);                                             \
	free_hashset_macro(hash_type);                                             \
	encap_hashset_macro(hash_type, hash_key_type, hash_code_macro);

/* ------------------ Useful functions ------------------------------------- */

static inline uint32_t __lh3_Jenkins_hash_int(uint32_t key){
	key += (key << 12);
	key ^= (key >> 22);
	key += (key << 4);
	key ^= (key >> 9);
	key += (key << 10);
	key ^= (key >> 2);
	key += (key << 7);
	key ^= (key >> 12);
	return key;
}

static inline uint64_t __lh3_Jenkins_hash_64(uint64_t key){
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

static inline uint32_t jenkins_one_at_a_time_hash(char *key, size_t len){
	uint32_t hash, i;
	for(hash = i = 0; i < len; ++i){
		hash += key[i];
		hash += (hash << 10);
		hash ^= (hash >> 6);
	}
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	return hash;
}

static inline uint64_t hash64shift(uint64_t key){
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}


static inline uint64_t MurmurHash64A(const void * key, int len, uint32_t seed){
	const uint64_t m = 0xc6a4a7935bd1e995LLU;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end){
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7){
	case 7: h ^= ((uint64_t)data2[6]) << 48;
	case 6: h ^= ((uint64_t)data2[5]) << 40;
	case 5: h ^= ((uint64_t)data2[4]) << 32;
	case 4: h ^= ((uint64_t)data2[3]) << 24;
	case 3: h ^= ((uint64_t)data2[2]) << 16;
	case 2: h ^= ((uint64_t)data2[1]) << 8;
	case 1: h ^= ((uint64_t)data2[0]);
	        h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

#define u32hashcode(key) __lh3_Jenkins_hash_int(key)
#define u64hashcode(key) __lh3_Jenkins_hash_64(key)

static inline uint32_t __string_hashcode(const char *s){
	uint32_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

#define u32hash_code(e) u32hashcode(e)
#define u64hash_code(e) u64hashcode(e)
#define uxxhash_equals(e1, e2) ((e1) == (e2))
define_hashset(u32hash, uint32_t, u32hash_code, uxxhash_equals);
define_hashset(u64hash, uint64_t, u64hash_code, uxxhash_equals);

#define i32hash_code(e) u32hashcode((uint32_t)(e))
#define i32hash_equals(e1, e2) ((e1) == (e2))
define_hashset(i32hash, int, i32hash_code, i32hash_equals);

#define chash_code(e) __string_hashcode(e)
#define chash_equals(e1, e2) (strcmp(e1, e2) == 0)
define_hashset(chash, char*, chash_code, chash_equals);

typedef struct { uint32_t key, val; } uuhash_t;
#define uuhash_code(e) u32hashcode((e).key)
#define uuhash_equals(e1, e2) ((e1).key == (e2).key)
define_hashset(uuhash, uuhash_t, uuhash_code, uuhash_equals);
define_key_val_hashset_ext(uuhash, uuhash_t, uint32_t, uint32_t);

typedef struct {uint32_t key; int val; } uihash_t;
#define uihashcode(E) u32hashcode((E).key)
#define uihashequals(E1, E2) (E1).key == (E2).key
define_hashset(uihash, uihash_t, uihashcode, uihashequals);
define_key_val_hashset_ext(uihash, uihash_t, uint32_t, int);

typedef struct {uint64_t key; uint64_t val; } UUhash_t;
#define UUhashcode(E) u64hashcode((E).key)
#define UUhashequals(E1, E2) (E1).key == (E2).key
define_hashset(UUhash, UUhash_t, UUhashcode, UUhashequals);
define_key_val_hashset_ext(UUhash, UUhash_t, uint64_t, uint64_t);

typedef struct { char *key; uint32_t val; } cuhash_t;
#define cuhash_code(e) __string_hashcode((e).key)
#define cuhash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
define_hashset(cuhash, cuhash_t, cuhash_code, cuhash_equals);
define_key_val_hashset_ext(cuhash, cuhash_t, char*, uint32_t);
static const obj_desc_t cuhash_struct_deep_obj_desc = {sizeof(cuhash_t), 1, {1}, {offsetof(cuhash_t, key)}, {(obj_desc_t*)&OBJ_DESC_CHAR_ARRAY}, NULL};
static inline size_t cuhash_deep_obj_desc_cnt(void *obj, int idx){
	switch(idx){
		case 0: return ((cuhash*)obj)->size;
		case 1: return (((cuhash*)obj)->size + 15) / 16 * 4;
		default: return 0;
	}
}
static const obj_desc_t cuhash_deep_obj_desc = {sizeof(cuhash), 2, {1, 1}, {offsetof(cuhash, array), offsetof(cuhash, flags)}, {(obj_desc_t*)&cuhash_struct_deep_obj_desc, (obj_desc_t*)&OBJ_DESC_DATA}, cuhash_deep_obj_desc_cnt};

typedef struct { char *key; int val; } cihash_t;
#define cihash_code(e) __string_hashcode((e).key)
#define cihash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
define_hashset(cihash, cihash_t, cihash_code, cihash_equals);
define_key_val_hashset_ext(cihash, cihash_t, char*, int);

typedef struct { char *key; unsigned long long val; } clhash_t;
#define clhash_code(e) __string_hashcode((e).key)
#define clhash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
define_hashset(clhash, clhash_t, clhash_code, clhash_equals);
define_key_val_hashset_ext(clhash, clhash_t, char*, unsigned long long);

/**
* Example of using userdata in thread-safe mode
* char **strs;
* ... codes init strs
* #define test_hc(E) __string_hashcode(((char**)set->userdata)[E])
* #define test_he(E1, E2) (strcmp(((char**)set->userdata)[E1], ((char**)set->userdata)[E2]) == 0)
* define_hashset(testhash, uint32_t, test_hc, test_he);
* testhash *hash = init_testhash(13);
* set_userdata_testhash(hash, strs);
* ... now, the key of testhash is uint32_t, but refer to strs
*/

#endif
