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

#ifndef WATCHTOWER_INDEX_DB_RJ_H
#define WATCHTOWER_INDEX_DB_RJ_H

#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "file_reader.h"
#include "bitvec.h"
#include "thread.h"
#include <getopt.h>

static const uint64_t WTIDX_MAGIC_NUM = 20150307;

typedef struct {
	uint32_t n_rd;
	uint32_t ksize, max_kmer_freq;
	BaseBank *rdseqs;
	u32list  *rdlens;
	u64list  *rdoffs;
	cplist   *rdnames;
	cuhash   *rdname2id;
	hzmv     *seeds;
	hzmhash  *hash;
} WTIDX;

static const obj_desc_t wtidx_obj_desc = {.size = sizeof(WTIDX), .n_child = 7, .addr = {offsetof(WTIDX, rdseqs), offsetof(WTIDX, rdlens), offsetof(WTIDX, rdoffs), offsetof(WTIDX, rdnames), offsetof(WTIDX, rdname2id), offsetof(WTIDX, seeds), offsetof(WTIDX, hash)}, .desc = {(obj_desc_t*)&basebank_obj_desc, (obj_desc_t*)&u32list_obj_desc, (obj_desc_t*)&u64list_obj_desc, (obj_desc_t*)&cplist_deep_obj_desc, (obj_desc_t*)&cuhash_deep_obj_desc, (obj_desc_t*)&hzmv_obj_desc, (obj_desc_t*)&hzmhash_obj_desc}, .cnt = NULL, .is_array=0};

static inline uint32_t wt_kmer_max(uint32_t ksize, int hk, int ksave){
	uint32_t i, max, save;
	save = ksave? 2 : 1;
	if(hk){
		max = 4;
		for(i=1;i<ksize;i++){
			max *= 3;
		}
	} else {
		max = 0xFFFFFFFFU >> (32 - (ksize << 1));
	}
	return max / save;
}

#endif
