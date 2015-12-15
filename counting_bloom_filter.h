
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
 
#ifndef __COUNTING_BLOOM_FILTER_RJ_H
#define __COUNTING_BLOOM_FILTER_RJ_H

#include "bitsvec.h"
#include "hashset.h"
#include "mem_share.h"

static const uint32_t cbf_total_seeds = 20;

static const uint32_t cbf_seeds[20] = 
{
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul
};

typedef struct {
	BitsVec *bits;
	size_t size;
	uint64_t *codes;
	uint8_t *cnts;
	uint32_t n_bit, n_seed, seed_off;
} CBF;

// size MUST be prime number
static inline CBF* init_cbf(size_t size, uint8_t n_bit, uint32_t n_seed){
	CBF *cbf;
	if(n_seed > cbf_total_seeds) n_seed = cbf_total_seeds;
	if(n_seed == 0) n_seed = 1;
	if(n_bit < 2) n_bit = 2;
	else if(n_bit > 8) n_bit = 8;
	size = _rj_hashset_find_prime(size);
	cbf = calloc(1, sizeof(CBF));
	cbf->bits = init_bitsvec(size, n_bit);
	cbf->size = size;
	cbf->n_bit = n_bit;
	cbf->codes = malloc(sizeof(uint64_t) * n_seed);
	cbf->cnts  = malloc(sizeof(uint8_t) * n_seed);
	cbf->n_seed = n_seed;
	cbf->seed_off = 0;
	return cbf;
}

static inline size_t cbf_obj_desc_cnt(void *obj, int idx){
	switch(idx){
		case 0: return 1;
		case 1: return ((CBF*)obj)->n_seed * sizeof(uint64_t);
		case 2: return ((CBF*)obj)->n_seed * sizeof(uint8_t);
		default: return 1;
	}
}

static const obj_desc_t cbf_obj_desc = {sizeof(CBF), 3, {1, 1, 1}, {offsetof(CBF, bits), offsetof(CBF, codes), offsetof(CBF, cnts)}, {(obj_desc_t*)&bitsvec_obj_desc, (obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&OBJ_DESC_DATA}, cbf_obj_desc_cnt};

static inline void clear_cbf(CBF *cbf){ clear_bitsvec(cbf->bits); }

static inline void change_seeds_cbf(CBF *cbf){ cbf->seed_off = (cbf->seed_off + cbf->n_seed) % cbf_total_seeds; }

static inline uint8_t put_cbf(CBF *cbf, const void *key, uint32_t len){
	uint32_t i, flag;
	uint8_t min;
	min = 0xFFU;
	flag = 0;
	for(i=0;i<cbf->n_seed;i++){
		cbf->codes[i] = MurmurHash64A(key, len, cbf_seeds[(i + cbf->seed_off) % cbf_total_seeds]) % cbf->size;
		cbf->cnts[i] = get_bitsvec(cbf->bits, cbf->codes[i]);
		if(cbf->cnts[i] < min){
			min = cbf->cnts[i];
			flag = 1U << i;
		} else if(cbf->cnts[i] == min){
			flag |= 1U << i;
		}
	}
	if(min >= (1U << cbf->n_bit) - 1) return min;
	min ++;
	for(i=0;i<cbf->n_seed;i++){
		if((flag >> i) & 0x01){
			set_bitsvec(cbf->bits, cbf->codes[i], min);
		}
	}
	return min;
}

static inline uint8_t get_cbf(CBF *cbf, const void *key, uint32_t len){
	uint32_t i;
	uint8_t min;
	min = 0xFFU;
	for(i=0;i<cbf->n_seed;i++){
		cbf->codes[i] = MurmurHash64A(key, len, cbf_seeds[(i + cbf->seed_off) % cbf_total_seeds]) % cbf->size;
		cbf->cnts[i] = get_bitsvec(cbf->bits, cbf->codes[i]);
		if(cbf->cnts[i] < min) min = cbf->cnts[i];
	}
	return min;
}

static inline void free_cbf(CBF *cbf){
	free_bitsvec(cbf->bits);
	free(cbf->codes);
	free(cbf->cnts);
	free(cbf);
}

#endif
