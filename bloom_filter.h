
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
 
#ifndef __BLOOM_FILTER_RJ_H
#define __BLOOM_FILTER_RJ_H

#include "bitvec.h"
#include "hashset.h"
#include "mem_share.h"

static const uint32_t bloom_filter_total_seeds = 20;

static const uint32_t bloom_filter_seeds[20] = 
{
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul
};

typedef struct {
	BitVec *bits;
	size_t size;
	uint32_t n_seed, seed_off;
} BloomFilter;

// size MUST be prime number
static inline BloomFilter* init_bloomfilter(size_t size, uint32_t n_seed){
	BloomFilter *bf;
	if(n_seed > bloom_filter_total_seeds) n_seed = bloom_filter_total_seeds;
	if(n_seed == 0) n_seed = 1;
	size = _rj_hashset_find_prime(size);
	bf = malloc(sizeof(BloomFilter));
	bf->bits = init_bitvec(size);
	bf->size = size;
	bf->n_seed = n_seed;
	bf->seed_off = 0;
	return bf;
}

static const obj_desc_t bloomfilter_obj_desc = {sizeof(BloomFilter), 1, {1}, {offsetof(BloomFilter, bits)}, {(obj_desc_t*)&bitvec_obj_desc}, NULL};

static inline void clear_bloomfilter(BloomFilter *bf){ zeros_bitvec(bf->bits); }

static inline void change_seeds_bloomfilter(BloomFilter *bf){ bf->seed_off = (bf->seed_off + bf->n_seed) % bloom_filter_total_seeds; }

static inline void put_bloomfilter(BloomFilter *bf, const void *key, uint32_t len){
	uint32_t i;
	for(i=0;i<bf->n_seed;i++) one_bitvec(bf->bits, MurmurHash64A(key, len, bloom_filter_seeds[(i + bf->seed_off) % bloom_filter_total_seeds]) % bf->size);
}

static inline int  get_bloomfilter(BloomFilter *bf, const void *key, uint32_t len){
	uint32_t i;
	for(i=0;i<bf->n_seed;i++){
		if(get_bitvec(bf->bits, MurmurHash64A(key, len, bloom_filter_seeds[(i + bf->seed_off) % bloom_filter_total_seeds]) % bf->size) == 0) return 0;
	}
	return 1;
}

static inline void free_bloomfilter(BloomFilter *bf){
	free_bitvec(bf->bits);
	free(bf);
}

#endif
