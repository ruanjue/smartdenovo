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
 
#ifndef __GOLOMB_CODING_RJ_H
#define __GOLOMB_CODING_RJ_H

#include "bitsvec.h"
#include "sort.h"
#include "list.h"

/**
 * Please use apply_array in sort.h to calculate occs.
 * Example:
 * uint64_t occs[256];
 * memset(occs, 0, 256 * 8);
 * apply_array(my_ary, ary_size, uint8_t, (occs[a] ++));
 */

static inline uint8_t estimate_golomb_param_bit8(uint64_t occs[256], uint64_t *est_bits){
	uint32_t i, j;
	uint8_t m, q, r, b;
	uint64_t bits, min, n;
	m = 0; min = 0xFFFFFFFFFFFFFFFFLLU;
	for(i=1;i<=128;i++){
		bits = 0;
		for(b=0;(1U<<b)<i;b++);
		for(j=0;j<256;j++){
			q = j / i;
			r = j % i;
			if((1U << b) == i){
				n = q + 1 + b;
			} else {
				n = q + 1 + ((r < ((1U << b) - i))? b - 1 : b);
			}
			bits += n * occs[j];
		}
		if(bits < min){ m = i; min = bits; }
	}
	if(est_bits) *est_bits = min;
	return m;
}

static inline uint16_t estimate_golomb_param_bit16(uint64_t *occs, uint64_t *est_bits){
	uint32_t i, j;
	uint16_t m, q, r, b;
	uint64_t bits, min, n;
	m = 0; min = 0xFFFFFFFFFFFFFFFFLLU;
	for(i=1;i<=256;i++){
		bits = 0;
		for(b=0;(1U<<b)<i;b++);
		for(j=0;j<65535;j++){
			q = j / i;
			r = j % i;
			if((1U << b) == i){
				n = q + 1 + b;
			} else {
				n = q + 1 + ((r < ((1U << b) - i))? b - 1 : b);
			}
			bits += n * occs[j];
		}
		if(bits < min){ m = i; min = bits; }
	}
	if(est_bits) *est_bits = min;
	return m;
}

static inline void golomb_encode_bit8(uint8_t *data, uint64_t size, uint8_t param, BitsVec *out){
	uint64_t i;
	uint8_t j, n, m, q, r, b, p;
	m = param;
	for(b=0;(1U<<b)<m;b++);
	p = (1U << b) - m;
	for(i=0;i<size;i++){
		n = data[i];
		q = n / m;
		r = n % m;
		for(j=0;j<q;j++) push_bitsvec(out, 1);
		push_bitsvec(out, 0);
		if(r < p){
			for(j=0;j<b-1;j++) push_bitsvec(out, (r >> j) & 0x01);
		} else{
			r += p;
			for(j=1;j<b;j++) push_bitsvec(out, (r >> j) & 0x01);
			push_bitsvec(out, r);
		}
	}
}

static inline uint64_t golomb_decode_bit8(BitsVec *inp, uint64_t *offset, uint8_t param, uint8_t *data, uint64_t max_size){
	uint64_t i, off;
	uint8_t j, m, u, q, r, b, rice, p;
	m = param;
	for(b=0;(1U<<b)<m;b++);
	rice = ((1U << b) == m);
	p = (1U << b) - m;
	off = *offset;
	u = 0;
	for(i=0;i<max_size&&off<inp->size;i++){
		q = 0;
		r = 0;
		while(1){
			u = get_bitsvec(inp, off);
			off ++;
			if(u == 0) break;
			q ++;
		}
		if(rice){
			r = 0;
			for(j=0;j<b-1;j++) r |= get_bitsvec(inp, off++) << j;
			u = get_bitsvec(inp, off ++);
			r = ((r) << 1) | (u);
		} else {
			r = 0;
			for(j=0;j<b-1;j++) r |= get_bitsvec(inp, off++) << j;
			if(r >= p){
				u = get_bitsvec(inp, off ++);
				r = ((r << 1) | (u & 0x01U)) - p;
			}
		}
		data[i] = q * m + r;
	}
	*offset = off;
	return i;
}

static inline uint8_t compress_golomb_bit8(BitsVec *bits, uint8_t *dat, uint64_t size, int remap){
	uint64_t occs[256], occs2[256];
	uint8_t trans[256], inv_trans[256];
	uint64_t i, j;
	uint8_t m, n;
	memset(occs, 0, 256 * 8);
	apply_array(dat, size, uint8_t, occs[a] ++);
	if(remap){
		for(i=0;i<256;i++) inv_trans[i] = i;
		sort_array(inv_trans, 256, uint8_t, occs[b] > occs[a]);
		for(i=0,n=0;i<256;i++,n++){
			if(occs[inv_trans[i]] == 0) break;
		}
		for(i=0;i<256;i++) trans[inv_trans[i]] = i;
		for(i=0;i<256;i++) occs2[i] = occs[inv_trans[i]];
		for(i=0;i<size;i++) dat[i] = trans[dat[i]];
		m = estimate_golomb_param_bit8(occs2, NULL);
	} else {
		m = estimate_golomb_param_bit8(occs, NULL);
		n = 0;
	}
	for(i=0;i<8;i++) push_bitsvec(bits, m >> i);
	for(i=0;i<8;i++) push_bitsvec(bits, n >> i);
	for(i=0;i<n;i++){
		for(j=0;j<8;j++) push_bitsvec(bits, inv_trans[i] >> j);
	}
	golomb_encode_bit8(dat, size, m, bits);
	return m;
}

static inline void decompress_golomb_bit8(BitsVec *bits, uint64_t *offset, u8list *dat){
	uint8_t inv_trans[256];
	uint64_t i, j, off, size;
	uint8_t m, n;
	off = *offset;
	for(i=0,m=0;i<8;i++) m |= get_bitsvec(bits, off ++) << i;
	for(i=0,n=0;i<8;i++) n |= get_bitsvec(bits, off ++) << i;
	memset(inv_trans, 0, 256);
	for(i=0;i<n;i++){
		for(j=0;j<8;j++) inv_trans[i] |= get_bitsvec(bits, off ++) << j;
	}
	*offset = off;
	size = 1024;
	while(size == 1024){
		encap_u8list(dat, 1024);
		size = golomb_decode_bit8(bits, offset, m, dat->buffer + dat->size, 1024);
		if(n) for(i=0;i<size;i++) dat->buffer[dat->size + i] = inv_trans[dat->buffer[dat->size + i]];
		dat->size += size;
	}
}

#endif
