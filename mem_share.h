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
 
#ifndef __MEM_SHARE_RJ_H
#define __MEM_SHARE_RJ_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <alloca.h>

#ifndef OBJ_DESC_MAX_CHILD
#define OBJ_DESC_MAX_CHILD 64
#endif

static inline size_t mem_size_round(size_t size){ return (size + 7) & 0xFFFFFFFFFFFFFFF8LLU; }

static inline uint8_t mem_size_gap(size_t size){ return (size & 0x07U)? 8 - (size & 0x07U) : 0; }

static inline size_t mem_dump(void *mem, size_t len, FILE *out){
	size_t size;
	uint8_t i, v;
	if(mem == NULL) return 0;
	size = mem_size_round(len);
	if(out){
		fwrite(mem, 1, len, out);
		v = 0;
		for(i=0;i<mem_size_gap(len);i++) fwrite(&v, 1, 1, out);
	}
	return size;
}

typedef size_t (*mem_array_count)(void *obj, int idx);

struct obj_desc_t;

#define MEM_PTR_TYPE_DUMP	1
#define MEM_PTR_TYPE_POINTER	2

typedef struct obj_desc_t {
	size_t size; // If size = 0, it is virtual, mem_type will be OR with MEM_PTR_TYPE_DUMP. See OBJ_DESC_CHAR_ARRAY
	int n_child; // <= OBJ_DESC_MAX_CHILD.
	uint8_t mem_type[OBJ_DESC_MAX_CHILD];
	off_t  addr[OBJ_DESC_MAX_CHILD]; // offsetof(type, field)
	const struct obj_desc_t *desc[OBJ_DESC_MAX_CHILD];
	mem_array_count cnt;
} obj_desc_t;

// Basic obj_desc_t, size = 1 byte
static const struct obj_desc_t OBJ_DESC_DATA = {1, 0, {}, {}, {}, NULL};
// Special obj_desc_t for string, set mem_type=0 and addr=0 to call the _char_array_obj_desc_cnt on itself
// so that we know the length of string, then set size=0 to indicate that it is an virtual reference, program should add MEM_PTR_TYPE_DUMP to its mem_type
static inline size_t _char_array_obj_desc_cnt(void *obj, int idx){ if(idx == 0) return strlen((char*)obj) + 1; else return 0; }
static const struct obj_desc_t OBJ_DESC_CHAR_ARRAY = {0, 1, {0}, {0}, {&OBJ_DESC_DATA}, _char_array_obj_desc_cnt};

static inline size_t mem_size_obj(void *obj, uint8_t mem_type, const obj_desc_t *desc, size_t size, size_t cnt){
	size_t m;
	void *ref;
	int i;
	if(desc == NULL) return size;
	if(obj == NULL) return size;
	switch(mem_type){
		case 3: size += mem_size_round(sizeof(void*) * cnt);
		case 2:
			for(m=0;m<cnt;m++) if(((void**)obj)[m]) size += mem_size_round(desc->size);
			break;
		case 1: size += mem_size_round(cnt * desc->size); // TODO: if desc == &OBJ_DESC_DATA, mem_size_round may waste many memory
		case 0: break;
	}
	if(desc->n_child == 0) return size;
	for(m=0;m<cnt;m++){
		if(mem_type & 0x02){
			ref = ((void**)obj)[m];
			if(ref == NULL) continue;
		} else {
			ref = obj + m * desc->size;
		}
		for(i=0;i<desc->n_child;i++){
			if(desc->mem_type[i] & 0x01){
				size += mem_size_obj(*((void**)(ref + desc->addr[i])), desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], 0, desc->cnt? desc->cnt(ref, i) : 1);
			} else {
				size += mem_size_obj(ref + desc->addr[i], desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], 0, desc->cnt? desc->cnt(ref, i) : 1);
			}
		}
	}
	return size;
}

static inline size_t mem_dump_obj(void *obj, uint8_t mem_type, const obj_desc_t *desc, size_t offset, size_t cnt, FILE *out, int mem_free){
	void *ref;
	size_t size, m;
	int i;
	if(obj == NULL) return offset;
	size = offset;
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			size += mem_dump(obj, cnt * sizeof(void*), out);
		} else {
			size += mem_dump(obj, cnt * desc->size, out);
		}
	}
	if((mem_type & 0x02) == 0 && desc->n_child == 0){
	} else {
		for(m=0;m<cnt;m++){
			if(mem_type & 0x02){
				ref = ((void**)obj)[m];
				if(ref == NULL) continue;
				size += mem_dump(ref, desc->size, out);
			} else {
				ref = obj + m * desc->size;
			}
			for(i=0;i<desc->n_child;i++){
				if(desc->mem_type[i] & 0x01){
					size = mem_dump_obj(*((void**)(ref + desc->addr[i])), desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], size, desc->cnt? desc->cnt(ref, i) : 1, out, mem_free);
				} else {
					size = mem_dump_obj(ref + desc->addr[i], desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], size, desc->cnt? desc->cnt(ref, i) : 1, out, mem_free);
				}
			}
		}
	}
	if(mem_free){
		if((mem_type & 0x02) && desc->size) for(m=0;m<cnt;m++) free(((void**)obj)[m]);
		if((mem_type & 0x01) && ((mem_type & 0x02) | desc->size)) free(obj);
	}
	return size;
}

static inline size_t mem_load_obj(void *obj, uint8_t mem_type, const obj_desc_t *desc, size_t addr_beg, size_t cnt){
	size_t addr, m;
	int i;
	void *ref, **ptr;
	if(obj == NULL) return 0;
	addr = addr_beg? : (size_t)obj;
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			addr += mem_size_round(cnt * sizeof(void*));
		} else {
			addr += mem_size_round(cnt * desc->size);
		}
	}
	if(desc->n_child == 0){
		switch(mem_type){
			case 2:
			case 3:
				for(m=0;m<cnt;m++){
					ptr = ((void**)obj) + m;
					if(*ptr == NULL) continue;
					*ptr = (void*)addr;
					addr += mem_size_round(desc->size);
				}
			default: return addr;
		}
	}
	for(m=0;m<cnt;m++){
		if(mem_type & 0x02){
			ptr = ((void**)obj) + m;
			if(*ptr == NULL) continue;
			ref = *ptr = (void*)addr;
			addr += mem_size_round(desc->size);
		} else {
			ref = obj + m * desc->size;
		}
		for(i=0;i<desc->n_child;i++){
			ptr = ref + desc->addr[i];
			if(desc->mem_type[i] & 0x01){
				if(*ptr == NULL) continue;
				*ptr = (void*)addr;
				addr = mem_load_obj(*ptr, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1);
			} else {
				addr = mem_load_obj((void*)ptr, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1);
			}
		}
	}
	return addr;
}

static inline size_t mem_dump_obj_file(void *obj, size_t mem_type, const obj_desc_t *desc, size_t cnt, size_t aux_data, FILE *out){
	size_t size;
	if(desc == NULL) return 0;
	if((mem_type & 0x01) == 0){
		fprintf(stderr, " -- Illegal mem_type (%u) to call mem_dump, object should have standalone memory in %s -- %s:%d --\n", (int)mem_type, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	size = mem_size_obj(obj, mem_type, desc, 0, cnt);
	if(out){
		fwrite(&size, sizeof(size_t), 1, out);
		fwrite(&mem_type, sizeof(size_t), 1, out);
		fwrite(&cnt, sizeof(size_t), 1, out);
		fwrite(&aux_data, sizeof(size_t), 1, out);
	}
	size = 4 * sizeof(size_t);
	size = mem_dump_obj(obj, mem_type, desc, size, cnt, out, 0);
	if(out) fflush(out);
	return size;
}

static inline size_t mem_dump_free_obj_file(void *obj, size_t mem_type, const obj_desc_t *desc, size_t cnt, size_t aux_data, FILE *out){
	size_t size;
	if(desc == NULL) return 0;
	if((mem_type & 0x01) == 0){
		fprintf(stderr, " -- Illegal mem_type (%u) to call mem_dump, object should have standalone memory in %s -- %s:%d --\n", (int)mem_type, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	size = mem_size_obj(obj, mem_type, desc, 0, cnt);
	if(out){
		fwrite(&size, sizeof(size_t), 1, out);
		fwrite(&mem_type, sizeof(size_t), 1, out);
		fwrite(&cnt, sizeof(size_t), 1, out);
		fwrite(&aux_data, sizeof(size_t), 1, out);
	}
	size = 4 * sizeof(size_t);
	size = mem_dump_obj(obj, mem_type, desc, size, cnt, out, 1);
	if(out) fflush(out);
	return size;
}

static char *mem_share_locks = NULL;
static int mem_share_lock_size = 0;

static inline void cleanup_mem_share_file_locks(){
	int off;
	off = 0;
	while(off < mem_share_lock_size){
		unlink(mem_share_locks + off);
		off += strlen(mem_share_locks + off) + 1;
	}
	if(mem_share_locks) free(mem_share_locks);
	mem_share_lock_size = 0;
	mem_share_locks = NULL;
}

#ifndef sighandler_t
typedef void (*sighandler_t)(int sig);
#endif
static sighandler_t sig_term = SIG_IGN;
static sighandler_t sig_int  = SIG_IGN;
static sighandler_t sig_hup  = SIG_IGN;
static sighandler_t sig_kill = SIG_IGN;
static volatile sig_atomic_t cleanup_mem_share_in_progress = 0;

static inline void sig_cleanup_mem_share_file_locks(int sig){
	if(cleanup_mem_share_in_progress) raise(sig);
	cleanup_mem_share_in_progress = 1;
	cleanup_mem_share_file_locks();
	signal(SIGTERM, sig_term);
	signal(SIGINT , sig_int);
	signal(SIGHUP, sig_hup);
	signal(SIGKILL, sig_kill);
	raise(sig);
}

static inline void register_mem_share_file_lock(char *file){
	int len;
	if(mem_share_lock_size == 0){
		if((sig_term = signal(SIGTERM, sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGTERM, SIG_IGN);
		if((sig_int  = signal(SIGINT , sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGINT , SIG_IGN);
		if((sig_hup  = signal(SIGHUP , sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGHUP , SIG_IGN);
		if((sig_kill = signal(SIGKILL, sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGKILL, SIG_IGN);
		atexit(cleanup_mem_share_file_locks);
	}
	len = strlen(file);
	mem_share_locks = realloc(mem_share_locks, mem_share_lock_size + len + 1);
	strcpy(mem_share_locks + mem_share_lock_size, file);
	mem_share_lock_size += len + 1;
}

// Directly read from file, don't share this object
static inline void* mem_read_obj_file(const obj_desc_t *desc, char *path, size_t *mem_type, size_t *cnt, size_t *aux_data){
	void *mem;
	size_t size, nin;
	FILE *file;
	if(desc == NULL) return NULL;
	if((file = fopen(path, "r")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	if(mem_type == NULL) mem_type = alloca(sizeof(size_t));
	if(cnt == NULL) cnt = alloca(sizeof(size_t));
	if(aux_data == NULL) aux_data = alloca(sizeof(size_t));
	fread(&size, sizeof(size_t), 1, file);
	fread(mem_type, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	mem = malloc(size);
	if(mem == NULL){
		fprintf(stderr, " -- Cannot alloc %llu bytes memory for %s in %s -- %s:%d --\n", (unsigned long long)size, path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	if((nin = fread(mem, 1, size, file)) != size){
		fprintf(stderr, " -- Read %llu bytes, not %llu bytes in %s -- %s:%d --\n", (unsigned long long)nin, (unsigned long long)size, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fclose(file);
	mem_load_obj(mem, *mem_type, desc, 0, *cnt);
	return mem;
}

static inline void* mem_load_obj_file_core(const obj_desc_t *desc, char *path, size_t *mem_type, size_t *cnt, size_t *aux_data){
	void *mem;
	size_t size, psize;
	char *lock;
	char hostname[65];
	FILE *file;
	int fd;
	if(desc == NULL) return NULL;
	if((file = fopen(path, "r+")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fread(&size, sizeof(size_t), 1, file);
	fread(mem_type, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	fd = fileno(file);
	psize = getpagesize();
	mem = mmap(0, (size + 4 * sizeof(size_t) + psize - 1) / psize * psize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if(mem == MAP_FAILED){
		perror("Cannot mmap");
		exit(1);
	}
	fclose(file);
	mem_load_obj(mem + 4 * sizeof(size_t), *mem_type, desc, 0, *cnt);
	gethostname(hostname, 64);
	lock = alloca(strlen(path) + strlen(hostname) + 20);
	sprintf(lock, "%s.mem_share.%s.%ld", path, hostname, gethostid());
	if((file = fopen(lock, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", lock, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fprintf(file, "%p\n", mem);
	fclose(file);
	fprintf(stderr, "-- mem_share '%s' at %p --\n", path, mem); fflush(stderr);
	register_mem_share_file_lock(lock);
	return mem + 4 * sizeof(size_t);
}

static inline void* mem_load_obj_file(const obj_desc_t *desc, char *path, size_t *mem_type, size_t *cnt, size_t *aux_data){
	char *lock;
	char hostname[65];
	void *addr, *mem;
	size_t size, psize;
	FILE *file;
	int fd;
	gethostname(hostname, 64);
	lock = alloca(strlen(path) + strlen(hostname) + 32);
	sprintf(lock, "%s.mem_share.%s.%ld", path, hostname, gethostid());
	if(mem_type == NULL) mem_type = alloca(sizeof(size_t));
	if(cnt == NULL) cnt = alloca(sizeof(size_t));
	if(aux_data == NULL) aux_data = alloca(sizeof(size_t));
	if((file = fopen(lock, "r")) == NULL){
		return mem_load_obj_file_core(desc, path, mem_type, cnt, aux_data);
	}
	fscanf(file, "%p", &addr);
	fclose(file);
	if((file = fopen(path, "r+")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fd = fileno(file);
	fread(&size, sizeof(void*), 1, file);
	fread(mem_type, sizeof(void*), 1, file);
	fread(cnt, sizeof(void*), 1, file);
	fread(aux_data, sizeof(void*), 1, file);
	psize = getpagesize();
	mem = mmap(addr, (size + 4 * sizeof(size_t) + psize - 1) / psize * psize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_FIXED, fd, 0);
	fclose(file);
	if(mem == MAP_FAILED){
		perror("Cannot map shared object");
		return mem_load_obj_file_core(desc, path, mem_type, cnt, aux_data);
	}
	mem_load_obj(mem + 4 * sizeof(size_t), *mem_type, desc, 0, *cnt);
	fprintf(stderr, "-- mem_map '%s' at %p --\n", path, addr); fflush(stderr);
	return mem + 4 * sizeof(size_t);
}

#endif

/*
 * An example to use mem_dump
*/

/*

#include "mem_share.h"

typedef struct {
	char *str;
	int val;
} Type1;

size_t type1_count(void *obj, int idx){ if(idx == 0){ return strlen(((Type1*)obj)->str) + 1; } else return 0; }
//const obj_desc_t type1_obj_desc = {sizeof(Type1), 1, {1}, {offsetof(Type1, str)}, {&OBJ_DESC_DATA}, type1_count};
const obj_desc_t type1_obj_desc = {sizeof(Type1), 1, {1}, {offsetof(Type1, str)}, {&OBJ_DESC_CHAR_ARRAY}, NULL};

typedef struct {
	int a, b, c;
	Type1 d1, d2[10], *d3, *d4[10], **d5;
	char **strs;
	int d3len, d5len;
} Type2;

size_t type2_count(void *obj, int idx){
	switch(idx){
		case 0: return 1;
		case 1: return 10;
		case 2: return ((Type2*)obj)->d3len;
		case 3: return 10;
		case 4: return ((Type2*)obj)->d5len;
		default: return 10;
	}
}

const obj_desc_t type2_obj_desc = {sizeof(Type2), 6, {0, 0, 1, 2, 3, 3}, {offsetof(Type2, d1), offsetof(Type2, d2), offsetof(Type2, d3), offsetof(Type2, d4), offsetof(Type2, d5), offsetof(Type2, strs)}, {&type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &OBJ_DESC_CHAR_ARRAY}, type2_count};

int main(){
	Type2 *t2, *t3;
	t2 = calloc(1, sizeof(Type2));
	int idx = 0;
	t2->d1.val = idx ++;
	t2->d1.str = strdup("d1");
	int i;
	for(i=0;i<10;i++){
		t2->d2[i].val = idx ++;
		t2->d2[i].str = strdup("d2");
	}
	t2->d3len = 10;
	t2->d3 = malloc(sizeof(Type1) * t2->d3len);
	for(i=0;i<t2->d3len;i++){
		t2->d3[i].val = idx ++;
		t2->d3[i].str = strdup("d3");
	}
	for(i=0;i<10;i++){
		t2->d4[i] = malloc(sizeof(Type1));
		t2->d4[i]->val = idx ++;
		t2->d4[i]->str = strdup("d4");
	}
	t2->d5len = 10;
	t2->d5 = malloc(sizeof(Type1*) * t2->d5len);
	for(i=0;i<t2->d5len;i++){
		t2->d5[i] = malloc(sizeof(Type1));
		t2->d5[i]->val = idx ++;
		t2->d5[i]->str = strdup("d5");
	}
	t2->strs = malloc(sizeof(char*) * 10);
	for(i=0;i<10;i++){
		t2->strs[i] = malloc(32);
		sprintf(t2->strs[i], "strs[%d,%d]", i, idx ++);
	}
	size_t aux_data, size, cnt, mem_type;
	FILE *file;
	size = mem_size_obj(t2, 1, &type2_obj_desc, 0, 1);
	fprintf(stdout, " -- size = %d in %s -- %s:%d --\n", (int)size, __FUNCTION__, __FILE__, __LINE__);
	aux_data = 1000999900;
	file = fopen("test.mem_share", "w");
	size = mem_dump_free_obj_file(t2, 1, &type2_obj_desc, 1, aux_data, file);
	fclose(file);
	fprintf(stdout, " -- size = %d in %s -- %s:%d --\n", (int)(size - 4 * sizeof(size_t)), __FUNCTION__, __FILE__, __LINE__);
	t3 = mem_read_obj_file(&type2_obj_desc, "test.mem_share", &mem_type, &cnt, &aux_data);
	fprintf(stdout, " -- aux_data = %d in %s -- %s:%d --\n", (int)aux_data, __FUNCTION__, __FILE__, __LINE__);
	printf("%d %s\n", t3->d1.val, t3->d1.str);
	for(i=0;i<10;i++) printf("%d %s\n", t3->d2[i].val, t3->d2[i].str);
	for(i=0;i<t3->d3len;i++) printf("%d %s\n", t3->d3[i].val, t3->d3[i].str);
	for(i=0;i<10;i++) printf("%d %s\n", t3->d4[i]->val, t3->d4[i]->str);
	for(i=0;i<t3->d5len;i++) printf("%d %s\n", t3->d5[i]->val, t3->d5[i]->str);
	for(i=0;i<10;i++) printf("%s\n", t3->strs[i]);
	free(t3);
	return 0;
}

*/
