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

#ifndef __BLOCK_SPARSE_ARRAY_RJ_H
#define __BLOCK_SPARSE_ARRAY_RJ_H

#include "list.h"
#include "hashset.h"

#define define_blocksparsearray(dtype, etype, enull)	\
	\
typedef struct {	\
	UUhash *hash;	\
	etype  *buffer;	\
	size_t buf_size;	\
	size_t buf_cap;	\
	size_t block_size;	\
	u8v    *trash;	\
	int mem_zero;	\
} dtype;	\
	\
static inline dtype* init_##dtype(size_t block_size){	\
	dtype *bsa;	\
	if(block_size == 0) block_size = 1;	\
	bsa = malloc(sizeof(dtype));	\
	bsa->hash = init_UUhash(13);	\
	bsa->buffer = NULL;	\
	bsa->buf_size = 0;	\
	bsa->buf_cap  = 0;	\
	bsa->block_size = block_size;	\
	bsa->trash = init_u8v(16);	\
	bsa->mem_zero = 0;	\
	return bsa;	\
}	\
	\
static inline void free_##dtype(dtype *bsa){	\
	free_UUhash(bsa->hash);	\
	free(bsa->buffer);	\
	free_u8v(bsa->trash);	\
	free(bsa);	\
}	\
	\
static inline void clear_##dtype(dtype *bsa){	\
	free_UUhash(bsa->hash); bsa->hash = init_UUhash(13);	\
	free(bsa->buffer); bsa->buffer = NULL;	\
	free_u8v(bsa->trash); bsa->trash = init_u8v(16);	\
	bsa->buf_size = 0;	\
	bsa->buf_cap  = 0;	\
}	\
	\
static inline size_t count_block_##dtype(dtype *bsa){ return bsa->hash->count; }	\
	\
static inline etype* ref_block_##dtype(dtype *bsa, size_t block_idx){	\
	UUhash_t *u, U;	\
	int exists;	\
	U.key = block_idx;	\
	U.val = 0;	\
	u = prepare_UUhash(bsa->hash, U, &exists);	\
	if(!exists){	\
		u->key = block_idx;	\
		if(bsa->trash->size){	\
			u->val = bsa->trash->buffer[--bsa->trash->size];	\
			if(bsa->mem_zero) memset(bsa->buffer + u->val, 0, sizeof(etype) * bsa->block_size);	\
		} else	\
			u->val = bsa->buf_size;	\
			bsa->buf_cap = encap_list(&bsa->buffer, sizeof(etype), bsa->buf_size, bsa->buf_cap, bsa->block_size, bsa->mem_zero);	\
			bsa->buf_size += bsa->block_size;	\
		}	\
	}	\
	return bsa->buffer + u->val;	\
}	\
	\
static inline etype get_##dtype(dtype *bsa, size_t idx){	\
	etype* blk;	\
	size_t bidx, boff;	\
	bidx = idx / bsa->block_size;	\
	boff = idx % bsa->block_size;	\
	blk = ref_block_##dtype(bsa, bidx);	\
	return blk[boff];	\
}	\
	\
static inline void set_##dtype(dtype *bsa, size_t idx, etype val){	\
	etype* blk;	\
	size_t bidx, boff;	\
	bidx = idx / bsa->block_size;	\
	boff = idx % bsa->block_size;	\
	blk = ref_block_##dtype(bsa, bidx);	\
	blk[boff] = val;	\
}	\
	\
static inline void del_block_##dtype(dtype *bsa, size_t block_idx){	\
	UUhash_t *u, U;	\
	U.key = block_idx;	\
	U.val = 0;	\
	if((u = get_UUhash(bsa->hash, U)) == NULL) return;	\
	push_u8v(bsa->trash, u->val);	\
	delete_UUhash(bsa->hash, u);	\
}

#endif
