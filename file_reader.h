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
 
#ifndef __FILE_READER_RJ_H
#define __FILE_READER_RJ_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "string.h"
#include "list.h"

/**
 * Sequence IO
 */

typedef struct {
	union { String name; String tag; };
	String header;
	String seq;
	String qual;
} Sequence;

typedef struct {
	FILE *file;
	char *filename;
	int is_proc;
} fr_file_t;
define_list(fr_filev, fr_file_t);

typedef struct {
	fr_filev *files;
	uint32_t fidx;
	char *buffer;
	int size;
	int capacity;
	int ptr;
	int last_brk;
	char line_breaker;
	char delimiter;
	unsigned long long n_line;
	String *line;
	String *vline;
	VStrv  *tabs;
	int seq_type;
} FileReader;

#define free_sequence(sequence) { if(sequence->name.string) free(sequence->name.string);\
	if(sequence->header.string) free(sequence->header.string);\
	if(sequence->seq.string) free(sequence->seq.string);\
	if(sequence->qual.string) free(sequence->qual.string);\
	free(sequence); }

FileReader* fopen_filereader(char *filename);

FileReader* fopen_filereader2(char *prefix, char *postfix);

FileReader* fopen_m_filereader(int n_file, char **file_names);

FileReader* stdin_filereader();

/**
 * Read characters from a copy of string
 */

FileReader* string_filereader(char *string);

void fclose_filereader(FileReader *fr);

size_t fseek_filereader(FileReader *fr, size_t pos);

int reset_filereader(FileReader *fr);

int fread_line(String *line, FileReader *fr);
int froll_back(FileReader *fr);

int fread_table(FileReader *fr);
#define get_col_str(fr, col) ref_VStrv((fr)->tabs, col)->string
#define get_col_len(fr, col) ref_VStrv((fr)->tabs, col)->size
#define col_length(fr, col) get_col_len(fr, (col) - 1)
#define col_string(fr, col) get_col_str(fr, (col) - 1)

typedef struct {
	int is_fq;
	int avg_seq_len;
	int min_seq_len;
	int max_seq_len;
} SeqFileAttr;

int guess_seq_file_type(FileReader *fr);

void guess_seq_file(FileReader *fr, SeqFileAttr *attr);

#define FASTA_FLAG_NORMAL		0
#define FASTA_FLAG_NO_NAME		1
#define FASTA_FLAG_NO_SEQ		2

int fread_fasta_adv(Sequence **seq, FileReader *fr, int flag);

#define fread_fasta(seq, fr) fread_fasta_adv(seq, fr, FASTA_FLAG_NORMAL)

#define FASTQ_FLAG_NORMAL		0
#define FASTQ_FLAG_NO_NAME		1
#define FASTQ_FLAG_NO_SEQ		2
#define FASTQ_FLAG_NO_QUAL		4

int fread_fastq_adv(Sequence **seq, FileReader *fr, int flag);

#define fread_fastq(seq, fr) fread_fastq_adv(seq, fr, FASTQ_FLAG_NORMAL)

#define SEQ_FLAG_NORMAL	0
#define SEQ_FLAG_NO_NAME	1
#define SEQ_FLAG_NO_SEQ	2
#define SEQ_FLAG_NO_QUAL	4

int fread_seq_adv(Sequence **seq, FileReader *fr, int flag);
#define fread_seq(seq, fr) fread_seq_adv(seq, fr, SEQ_FLAG_NORMAL)

char * fread_all(FileReader *fr);

static inline int file_exists(const char *filename){
	struct stat s;
	if(stat(filename, &s) == -1) return 0;
	switch(s.st_mode & S_IFMT){
		//case S_IFBLK:
		//case S_IFCHR:
		//case S_IFDIR:
		//case S_IFIFO:
		//case S_IFSOCK:
		case S_IFLNK:
		case S_IFREG: return 1;
		default: return 0;

	}
}

static inline void print_pretty_seq(FILE *out, String *seq, int line_width){
	char c;
	int i, j;
	i = 0;
	while(i < seq->size){
		j = i + line_width;
		if(j > seq->size) j = seq->size;
		c  = seq->string[j];
		seq->string[j] = '\0';
		fprintf(out, "%s\n", seq->string + i);
		seq->string[j] = c;
		i = j;
	}
}

static inline void print_pretty_str(FILE *out, char *seq, int size, int line_width){
	char c;
	int i, j;
	i = 0;
	while(i < size){
		j = i + line_width;
		if(j > size) j = size;
		c  = seq[j];
		seq[j] = '\0';
		fprintf(out, "%s\n", seq + i);
		seq[j] = c;
		i = j;
	}
}

static inline FILE* open_file_for_read(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(name == NULL && suffix == NULL){
		full_name = "-";
	} else if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	if(strcmp(full_name, "-") == 0){
		file = stdin;
	} else {
		file = fopen(full_name, "r");
	}
	if(file == NULL){
		fprintf(stderr, "Cannot open file for read: %s\n", full_name);
		perror(NULL);
		exit(1);
	}
	return file;
}

static inline FILE* open_file_for_write(char *name, char *suffix, int overwrite){
	char *full_name;
	FILE *file;
	if(name == NULL && suffix == NULL){
		full_name = "-";
	} else if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	if(strcmp(full_name, "-") == 0){
		file = stdout;
	} else if(!overwrite && file_exists(full_name)){
		fprintf(stderr, "File exists: %s\n", full_name); exit(1);
	} else {
		file = fopen(full_name, "w+");
	}
	if(file == NULL){
		fprintf(stderr, "Cannot open file for write: %s\n", full_name);
		perror(NULL);
		exit(1);
	}
	return file;
}

static inline FILE* open_file_for_append(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(name == NULL && suffix == NULL){
		full_name = "-";
	} else if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	if(strcmp(full_name, "-") == 0){
		file = stdout;
	} else {
		file = fopen(full_name, "a+");
	}
	if(file == NULL){
		fprintf(stderr, "Cannot open file for append: %s\n", full_name);
		perror(NULL);
		exit(1);
	}
	return file;
}

static inline void close_file_safely(FILE *file){
	if(file == NULL) return;
	if(file == stdin || file == stdout || file == stderr) return;
	if(fclose(file)) perror("Error on close file");
}

#endif
