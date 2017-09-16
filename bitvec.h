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
 
#ifndef __BIT_VEC_RJ_H
#define __BIT_VEC_RJ_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "mem_share.h"

#define get_bit8(bits, idx) ((((bits)[(idx) >> 3]) >> ((idx) & 0x07)) & 0x01)
#define get_bit16(bits, idx) ((((bits)[(idx) >> 4]) >> ((idx) & 0x0F)) & 0x01)
#define get_bit32(bits, idx) ((((bits)[(idx) >> 5]) >> ((idx) & 0x1F)) & 0x01)
#define get_bit64(bits, idx) ((((bits)[(idx) >> 6]) >> ((idx) & 0x3F)) & 0x01)

#define get_2bit8(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 1)) & 0x03)
#define get_2bit16(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 1)) & 0x03)
#define get_2bit32(bits, idx) ((((bits)[(idx) >> 4]) >> (((idx) & 0x0F) << 1)) & 0x03)
#define get_2bit64(bits, idx) ((((bits)[(idx) >> 5]) >> (((idx) & 0x1F) << 1)) & 0x03)

#define get_4bit8(bits, idx) ((((bits)[(idx) >> 1]) >> (((idx) & 0x01) << 2)) & 0x0F)
#define get_4bit16(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 2)) & 0x0F)
#define get_4bit32(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 2)) & 0x0F)
#define get_4bit64(bits, idx) ((((bits)[(idx) >> 4]) >> (((idx) & 0x0F) << 2)) & 0x0F)

static const uint8_t byte_ones_table[256] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

static inline unsigned int _bitvec_roundup_power2(unsigned int v){
	if(v == 0) return 0;
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

typedef struct {
	uint64_t *bits;
	uint64_t n_bit;
	uint64_t n_cap;

	uint64_t *sums;
	uint64_t sum_size;
	uint64_t n_ones;
	uint64_t *hash;
	uint64_t hash_size;
	uint64_t hash_mod;
	int64_t iter_idx;
} BitVec;

#if 0

static inline uint32_t count_ones_bit32(uint32_t v){
	v = v - ((v >> 1) & 0x55555555U);                        // reuse input as temporary
	v = (v & 0x33333333U) + ((v >> 2) & 0x33333333U);        // temp
	return (((v + (v >> 4)) & 0xF0F0F0FU) * 0x1010101U) >> 24; // count
}

#define ONES_STEP_4 0x1111111111111111ULL
#define ONES_STEP_8 0x0101010101010101ULL

static inline int count_ones_bit64(const uint64_t x){
	register uint64_t byte_sums = x - ((x & 0xa * ONES_STEP_4) >> 1);
	byte_sums = (byte_sums & 3 * ONES_STEP_4) + ((byte_sums >> 2) & 3 * ONES_STEP_4);
	byte_sums = (byte_sums + (byte_sums >> 4)) & 0x0f * ONES_STEP_8;
	return byte_sums * ONES_STEP_8 >> 56;
}

#else

#define count_ones_bit32(v) __builtin_popcount(v)
#define count_ones_bit64(v) __builtin_popcountll(v)

#endif

static inline size_t bitvec_obj_desc_cnt(void *bitv, int idx){
	switch(idx){
		case 0: return ((BitVec*)bitv)->n_cap / 64 * 8;
		case 1: return ((BitVec*)bitv)->sums? (((BitVec*)bitv)->sum_size * 2 + 1) * 8 : 0;
		case 2: return ((BitVec*)bitv)->hash? (((BitVec*)bitv)->hash_size) * 8 : 0;
		default: return 0;
	}
}

static const obj_desc_t bitvec_obj_desc = {sizeof(BitVec), 3, {1, 1, 1}, {offsetof(BitVec, bits), offsetof(BitVec, sums), offsetof(BitVec, hash)}, {(obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&OBJ_DESC_DATA}, bitvec_obj_desc_cnt};

static inline BitVec* init_bitvec(uint64_t n_bit){
	BitVec *bitv;
	if(n_bit == 0) n_bit = 64 * 8;
	bitv = (BitVec*)malloc(sizeof(BitVec));
	bitv->n_bit = 0;
	bitv->n_cap = (((n_bit + 63) / 64) + 7) / 8 * 64 * 8;
	bitv->bits  = (uint64_t*)calloc(bitv->n_cap / 64 + 1, 8);
	//memset(bitv->bits, 0, bitv->n_cap / 8);
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->sum_size = 0;
	bitv->n_ones = 0;
	bitv->hash_size = 0;
	bitv->hash_mod = 0;
	bitv->iter_idx = 0;
	return bitv;
}

static inline size_t dump_bitvec(BitVec *bitv, FILE *out){
	fwrite(&bitv->n_bit, sizeof(uint64_t), 1, out);
	fwrite(&bitv->n_cap, sizeof(uint64_t), 1, out);
	fwrite(bitv->bits, sizeof(uint64_t), bitv->n_cap / 64, out);
	return sizeof(uint64_t) * (2 + bitv->n_cap / 64);
}

static inline BitVec* load_bitvec(FILE *inp){
	BitVec *bitv;
	size_t n;
	bitv = (BitVec*)malloc(sizeof(BitVec));
	if((n = fread(&bitv->n_bit, sizeof(uint64_t), 1, inp)) != 1){
		free(bitv); return NULL;
	}
	if((n = fread(&bitv->n_cap, sizeof(uint64_t), 1, inp)) != 1){
		free(bitv); return NULL;
	}
	bitv->bits = (uint64_t*)malloc(bitv->n_cap / 8);
	if(bitv->bits == NULL){
		fprintf(stderr, " Out of memeory in load_bitvec\n "); fflush(stderr); exit(1);
	}
	if((n = fread(bitv->bits, sizeof(uint64_t), bitv->n_cap / 64, inp)) != bitv->n_cap / 64){
		free(bitv); free(bitv->bits); return NULL;
	}
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->hash_size = 0;
	return bitv;
}

#if 0
static inline BitVec* mem_load_bitvec(void *mem, FILE *inp){
	BitVec *bitv;
	size_t off, n;
	bitv = mem;
	off = ((sizeof(BitVec) + 7) / 8) * 8;
	if((n = fread(&bitv->n_bit, sizeof(uint64_t), 1, inp)) != 1) return NULL;
	if((n = fread(&bitv->n_cap, sizeof(uint64_t), 1, inp)) != 1) return NULL;
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->hash_size = 0;
	bitv->bits = mem + off;
	off += (bitv->n_cap / 64) * 8;
	if((n = fread(bitv->bits, sizeof(uint64_t), bitv->n_cap / 64, inp)) != bitv->n_cap / 64) return NULL;
	return bitv;
}
#endif

static inline void clear_bitvec(BitVec *bitv){ bitv->n_bit = 0; }

static inline void zeros_bitvec(BitVec *bitv){ memset(bitv->bits, 0, bitv->n_cap / 8); }

static inline void ones_bitvec(BitVec *bitv){ memset(bitv->bits, 0xFFU, bitv->n_cap / 8); }

static inline void flip_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] ^= 1LLU << (idx&0x3FU); }

static inline void one_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] |= 1LLU << (idx&0x3FU); }

static inline void zero_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] &= ~(1LLU << (idx&0x3FU)); }

static inline uint64_t get_bitvec(BitVec *bitv, uint64_t idx){ return (bitv->bits[idx>>6] >> (idx&0x3FU)) & 0x01LLU; }

static inline void encap_bitvec(BitVec *bitv, uint64_t num){
	uint64_t cap;
	if(bitv->n_bit + num < bitv->n_cap) return;
	cap = bitv->n_cap;
	while(bitv->n_bit + num >= bitv->n_cap){
		if(bitv->n_cap < 1024 * 1024 * 8){
			bitv->n_cap <<= 1;
		} else bitv->n_cap += 1024 * 1024 * 8;
	}
	bitv->bits = (uint64_t*)realloc(bitv->bits, bitv->n_cap / 8 + 8);
	memset(((void*)bitv->bits) + cap / 8, 0, (bitv->n_cap - cap) / 8 + 8);
	bitv->bits[cap / 64] = 0x0000000000000001LLU;
}

static inline void recap_bitvec(BitVec *bitv, uint64_t new_cap){
	if(new_cap & 0x3FU) new_cap = (new_cap & 0xFFFFFFFFFFFFFFC0LLU) + 0x40U;
	if(bitv->n_cap == new_cap) return;
	bitv->bits = (uint64_t*)realloc(bitv->bits, new_cap / 8 + 8);
	if(new_cap > bitv->n_cap){
		memset(((void*)bitv->bits) + bitv->n_cap / 8, 0, (new_cap - bitv->n_cap) / 8 + 8);
	}
	bitv->bits[new_cap / 64] = 0x0000000000000001LLU;
	bitv->n_cap = new_cap;
}

static inline void one2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); one_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline void zero2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); zero_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline uint64_t get_2bitvec(BitVec *bitv, uint64_t idx){ return (bitv->bits[idx>>5] >> ((idx&0x1FU) << 1)) & 0x03LLU; }

static inline void set_2bitvec(BitVec *bitv, uint64_t idx, uint64_t v){
	bitv->bits[idx>>5] = (bitv->bits[idx>>5] & (~(0x03LLU << ((idx&0x1FU) << 1)))) | ((v&0x03LLU) << ((idx&0x1FU) << 1));
}

static inline void push_2bitvec(BitVec *bitv, uint64_t v){
	encap_bitvec(bitv, 2);
	set_2bitvec(bitv, bitv->n_bit >> 1, v);
	bitv->n_bit = ((bitv->n_bit >> 1) + 1) << 1;
}

static inline void end_bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); one_bitvec(bitv, bitv->n_bit); }

static inline uint64_t next_one_bitvec(BitVec *bitv, uint64_t idx){
	register uint64_t p, v;
	register uint32_t s;
	p = idx >> 6;
	s = idx & 0x3F;
	while(!(bitv->bits[p] >> s)){ p ++; s = 0; }
	v = bitv->bits[p] >> s;
	s += __builtin_ctzll(v);
	return (p << 6) + s;
}

static const int Mod37BitPosition[] = // map a bit value mod 37 to its position
{
 32,  0,  1, 26,  2, 23, 27,  0,  3, 16,
 24, 30, 28, 11,  0, 13,  4,  7, 17,  0,
 25, 22, 31, 15, 29, 10, 12,  6,  0, 21,
 14,  9,  5, 20,  8, 19, 18
};

static inline uint64_t next_one_bitvec2(BitVec *bitv, uint64_t idx){
	register uint64_t p;
	register uint32_t s, v;
	p = idx >> 6;
	s = idx & 0x3F;
	while(!(bitv->bits[p] >> s)){ p ++; s = 0; }
	if(!((bitv->bits[p] >> s) & 0xFFFFFFFFU)) s += 32;
	v = bitv->bits[p] >> s;
	s += Mod37BitPosition[(-v & v) % 37];
	return (p << 6) + s;
}

static inline uint64_t next_one_bitvec3(BitVec *bitv, uint64_t idx){
	register uint64_t p;
	register uint32_t s;
	p = idx >> 6;
	s = idx & 0x3F;
	while(!(bitv->bits[p] >> s)){ p ++; s = 0; }
	while(!((bitv->bits[p] >> s) & 0xFFU)) s += 8;
	while(!((bitv->bits[p] >> s) & 0x01U)) s ++;
	return (p << 6) + s;
}

static inline void index_bitvec(BitVec *bitv){
	uint64_t i, k, s, t, m;
	m = ((bitv->n_cap + 63) / 64 + 7) / 8;
	bitv->sums = (uint64_t*)calloc((m * 2 + 1), 8);
	t = 0;
	for(i=0;i<bitv->n_cap;i+=64*8){
		k = ((i>>6) >> 3) << 1;
		bitv->sums[k] = t;
		s = 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+0]);
		bitv->sums[k+1] |= s << 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+1]);
		bitv->sums[k+1] |= s << 9;
		s += count_ones_bit64(bitv->bits[(i>>6)+2]);
		bitv->sums[k+1] |= s << 18;
		s += count_ones_bit64(bitv->bits[(i>>6)+3]);
		bitv->sums[k+1] |= s << 27;
		s += count_ones_bit64(bitv->bits[(i>>6)+4]);
		bitv->sums[k+1] |= s << 36;
		s += count_ones_bit64(bitv->bits[(i>>6)+5]);
		bitv->sums[k+1] |= s << 45;
		s += count_ones_bit64(bitv->bits[(i>>6)+6]);
		bitv->sums[k+1] |= s << 54;
		s += count_ones_bit64(bitv->bits[(i>>6)+7]);
		t += s;
	}
	bitv->sums[((i>>6) >> 3) << 1] = t;
	bitv->n_ones = t;
	bitv->sum_size = m;
	bitv->hash_size = (bitv->n_cap / 64 / 8) / 2;
	if(bitv->hash_size == 0) bitv->hash_size = 1;
	bitv->hash_mod = (t + bitv->hash_size) / bitv->hash_size;
	if(bitv->hash_mod == 0) bitv->hash_mod = 1;
	bitv->hash = (uint64_t*)malloc(sizeof(uint64_t) * bitv->hash_size);
	s = 0;
	t = 0;
	for(i=0;i<=m;i++){
		k = bitv->sums[i*2] / bitv->hash_mod;
		if(s < k){
			while(s < k){ bitv->hash[s] = t; s ++; }
			t = i? i - 1 : 0;
		}
	}
	bitv->hash[bitv->sums[m*2] / bitv->hash_mod] = t;
}

static inline uint64_t rank_bitvec(BitVec *bitv, uint64_t idx){
	uint64_t p, s, sum;
	p = (idx>>6)>>3;
	s = (idx >> 6) & 0x07U;
	sum = bitv->sums[p<<1];
	if(s) sum += (bitv->sums[(p<<1)+1] >> (9 * (s - 1))) & 0x1FFU;
	if(idx & 0x3FU) sum += count_ones_bit64(bitv->bits[idx>>6]<<(64-(idx&0x3FU)));
	return sum;
}

static inline uint8_t select_8bytes(uint64_t word, uint8_t n_one){
	uint8_t idx, n, m;
	n = count_ones_bit32((uint32_t)word);
	if(n >= n_one){
		n = 0;
		idx = 0;
		word = word & 0xFFFFFFFFU;
	} else {
		idx = 32;
		word = word >> 32;
	}
	while(1){
		m = byte_ones_table[(uint8_t)word];
		if(n + m >= n_one) break;
		n += m;
		idx += 8;
		word >>= 8;
	}
	m = byte_ones_table[(uint8_t)(word & 0xF)];
	if(n + m < n_one){
		idx += 4;
		word >>= 4;
		n += m;
	}
	while(word){
		idx ++;
		if(word & 0x01){
			n ++;
			if(n == n_one) break;
		}
		word >>= 1;
	}
	return idx;
}

/*
 * To select the 1'st one, use select_bitvec(bitv, 1) - 1
 * */
static inline uint64_t select_bitvec(BitVec *bitv, uint64_t idx){
	uint64_t i, p, s, sum, t;
	p = bitv->hash[idx / bitv->hash_mod];
	while(p + 1 < bitv->sum_size && bitv->sums[(p + 1) << 1] < idx) p ++;
	sum = bitv->sums[p << 1];
	i = 0;
	t = sum;
	while(i < 7){
		s = (bitv->sums[(p << 1) + 1] >> (9 * i)) & 0x1FFU;
		if(sum + s >= idx) break;
		t = sum + s;
		i ++;
	}
	p = p * 8 + i;
	s = idx - t;
	return p * 64 + select_8bytes(bitv->bits[p], s);
}

static inline void begin_iter_bitvec(BitVec *bitv){ bitv->iter_idx = -1; }

static inline uint64_t iter_bitvec(BitVec *bitv){
	if((uint64_t)(bitv->iter_idx + 1) > bitv->n_cap) return 0xFFFFFFFFFFFFFFFFLLU;
	bitv->iter_idx = next_one_bitvec(bitv, bitv->iter_idx + 1);
	return (uint64_t)bitv->iter_idx;
}

static inline void free_bitvec(BitVec *bitv){
	free(bitv->bits);
	if(bitv->sums) free(bitv->sums);
	if(bitv->hash) free(bitv->hash);
	free(bitv);
}

#if 0

static inline size_t mem_size_bitvec(BitVec *bitv){
	size_t m;
	m = (sizeof(BitVec) + 7) / 8 * 8 + ((bitv->n_cap / 64) * 8);
	if(bitv->sums){
		m += (bitv->sum_size * 2 + 1) * 8;
	}
	if(bitv->hash){
		m += bitv->hash_size * 8;
	}
	return m;
}

static inline size_t mem_dump_bitvec(BitVec *bitv, void *mem){
	BitVec *clone;
	size_t off;
	clone = mem;
	memcpy(clone, bitv, sizeof(BitVec));
	off = ((sizeof(BitVec) + 7) / 8) * 8;
	clone->bits = mem + off;
	memcpy(clone->bits, bitv->bits, (bitv->n_cap / 64) * 8);
	off += (bitv->n_cap / 64) * 8;
	if(bitv->sums){
		clone->sums = mem + off;
		memcpy(clone->sums, bitv->sums, (bitv->sum_size * 2 + 1) * 8);
		off += (bitv->sum_size * 2 + 1) * 8;
	}
	if(bitv->hash){
		clone->hash = mem + off;
		memcpy(clone->hash, bitv->hash, bitv->hash_size * 8);
		off += bitv->hash_size * 8;
	}
	return off;
}
#endif

#endif
