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

#include "wtidx.h"

WTIDX* init_wtidx(uint32_t ksize, uint32_t kcut){
	WTIDX *wt;
	wt = malloc(sizeof(WTIDX));
	wt->rdseqs = init_basebank();
	wt->rdlens = init_u32list(1024);
	wt->rdoffs = init_u64list(1024);
	wt->rdnames = init_cplist(1024);
	wt->rdname2id = init_cuhash(1023);
	wt->seeds = init_hzmv(1024);
	wt->hash  = init_hzmhash(1023);
	wt->n_rd = 0;
	wt->ksize = ksize;
	wt->max_kmer_freq = kcut;
	return wt;
}

void free_wtidx(WTIDX *wt){
	uint32_t i;
	free_basebank(wt->rdseqs);
	free_u32list(wt->rdlens);
	free_u64list(wt->rdoffs);
	for(i=0;i<wt->rdnames->size;i++) free((char*)get_cplist(wt->rdnames, i));
	free_cplist(wt->rdnames);
	free_cuhash(wt->rdname2id);
	free_hzmv(wt->seeds);
	free_hzmhash(wt->hash);
	free(wt);
}

void push_long_read_wtidx(WTIDX *wt, char *name, int name_len, char *seq, int seq_len){
	char *ptr;
	push_u32list(wt->rdlens, seq_len);
	push_u64list(wt->rdoffs, wt->rdseqs->size);
	seq2basebank(wt->rdseqs, seq, seq_len);
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_cplist(wt->rdnames, ptr);
	kv_put_cuhash(wt->rdname2id, ptr, wt->n_rd);
	wt->n_rd ++;
}

void set_read_clip_wtidx(WTIDX *wt, char *name, int coff, int clen){
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->rdname2id, name)) == 0xFFFFFFFFU) return;
	if(coff < 0 || coff + clen > (int)wt->rdlens->buffer[pbid]) return;
	wt->rdoffs->buffer[pbid] += coff;
	wt->rdlens->buffer[pbid]  = clen;
}

#define homo_compress(bs, zs) if((bs) == (zs)){ continue; }

void index_wtidx(WTIDX *wt, uint32_t ncpu){
	u32list *hzoff;
	hzm_t   *m;
	hzmh_t  *u, U;
	uint64_t kmer, krev, kmask, off, idx, nflt, nrem, none;
	uint32_t pbid, pblen, kcnt, i, j;
	uint8_t b, c, dir;
	int exists;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - wt->ksize) << 1);
	hzoff = init_u32list(1024);
	memset(&U, 0, sizeof(hzmh_t));
	clear_hzmv(wt->seeds);
	fprintf(stderr, "[%s] - spliting kmers (%d bp)\n", date(), wt->ksize);
	for(pbid=0;pbid<wt->n_rd;pbid++){
		if((pbid % 1000) == 0){ fprintf(stderr, "\r%u", pbid); fflush(stderr); }
		pblen = wt->rdlens->buffer[pbid];
		b = 4;
		kmer = 0;
		off = wt->rdoffs->buffer[pbid];
		clear_u32list(hzoff);
		for(i=j=0;j<pblen;j++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			homo_compress(c, b);
			b = c;
			i ++;
			push_u32list(hzoff, j);
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->ksize) continue;
			krev = dna_rev_seq(kmer, wt->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			m = next_ref_hzmv(wt->seeds);
			m->rd_id = pbid;
			m->mer   = krev;
			m->dir   = dir;
			m->off   = hzoff->buffer[i - wt->ksize];
			m->len   = j + 1 - m->off;
		}
	}
	fprintf(stderr, "\r%u\n", pbid); fflush(stderr);
	fprintf(stderr, "[%s] - sorting kmers, %u threads\n", date(), ncpu);
	psort_array(wt->seeds->buffer, wt->seeds->size, hzm_t, ncpu, (a.mer > b.mer)? 1 : ((a.mer < b.mer)? 0 : (((((uint64_t)a.rd_id) << 32) | a.off) > ((((uint64_t)b.rd_id) << 32) | b.off))));
	nrem = none = nflt = kcnt = 0; kmer = kmask;
	for(idx=0;idx<wt->seeds->size;idx++){
		if(wt->seeds->buffer[idx].mer != kmer){
			kmer = wt->seeds->buffer[idx].mer;
			if(kcnt == 0);
			else if(kcnt == 1) none ++;
			else if(kcnt > wt->max_kmer_freq) nflt ++;
			else nrem ++;
			kcnt = 1;
		} else kcnt ++;
	}
	fprintf(stderr, "[%s] - filtered %llu single kmers (==1)\n", date(), (unsigned long long)none);
	fprintf(stderr, "[%s] - filtered %llu high frequency kmers (>=%d)\n", date(), (unsigned long long)nflt, wt->max_kmer_freq);
	fprintf(stderr, "[%s] - hashing %llu kmers\n", date(), (unsigned long long)nrem);
	clear_hzmhash(wt->hash);
	encap_hzmhash(wt->hash, nrem / wt->hash->load_factor + 1024);
	kmer = kmask;
	u = &U;
	kcnt = 0;
	for(idx=off=0;idx<wt->seeds->size;idx++){
		if(wt->seeds->buffer[idx].mer != kmer){
			if(kcnt > 1 && kcnt <= wt->max_kmer_freq){
				U.mer = kmer;
				u = prepare_hzmhash(wt->hash, U, &exists);
				if(exists){
					fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); exit(1);
				}
				u->mer = kmer;
				u->off = off;
				u->cnt = kcnt;
			}
			kmer = wt->seeds->buffer[idx].mer;
			kcnt = 1;
			off = idx;
		} else kcnt ++;
	}
	fprintf(stderr, "[%s] - Done\n", date());
	free_u32list(hzoff);
}

int usage(){
	printf(
	"WTIDX: indexing long reads using homopolymer compressed k-mer seeding\n"
	"WatchTower: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtidx [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -i <string> Long reads sequences file, + *\n"
	" -b <string> Long reads retained region, often from wtobt/wtcyc, +\n"
	"             Format: read_name\\toffset\\tlength\\toriginal_len\n"
	" -o <string> Index file of long reads, *\n"
	" -f          Force overwrite\n"
	" -k <int>    Kmer size, 5 <= <-k> <= 32, [16]\n"
	" -K <int>    Filter high frequency kmers, maybe repetitive, [10000]\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	WTIDX *wt;
	cplist *pbs, *obts;
	FileReader *fr;
	Sequence *seq;
	char *output;
	FILE *out;
	int c, ncpu, ksize, kcut, overwrite;
	output = NULL;
	ncpu = 1;
	ksize = 16;
	kcut = 10000;
	overwrite = 0;
	pbs = init_cplist(4);
	obts = init_cplist(4);
	while((c = getopt(argc, argv, "ht:i:b:o:fk:K:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			case 'o': output = optarg; break;
			case 'f': overwrite = 1; break;
			case 'k': ksize = atoi(optarg); break;
			case 'K': kcut = atoi(optarg); break;
			default: return usage();
		}
	}
	if(output == NULL) return usage();
	if(!overwrite && strcmp(output, "-") && file_exists(output)){
		fprintf(stderr, "File exists! '%s'\n\n", output);
		return usage();
	}
	if(pbs->size == 0) return usage();
	if(ksize > 32 || ksize < 5) return usage();
	wt = init_wtidx(ksize, kcut);
	if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	fprintf(stderr, "[%s] loading long reads\n", date());
	seq = NULL;
	while(fread_seq(&seq, fr)){
		push_long_read_wtidx(wt, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
		if((wt->n_rd % 1000) == 0){
			fprintf(stderr, "\r%u", wt->n_rd); fflush(stderr);
		}
	}
	fclose_filereader(fr);
	fprintf(stderr, "\r[%s] Done, %u reads\n", date(), (unsigned)wt->n_rd);
	if(obts->size){
		fprintf(stderr, "[%s] loading reads obt information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(c < 3) continue;
			set_read_clip_wtidx(wt, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date());
	} else {
		fprintf(stderr, "[%s] No obt information\n", date());
	}
	fprintf(stderr, "[%s] indexing\n", date());
	index_wtidx(wt, ncpu);
	fprintf(stderr, "[%s] Done\n", date());
	out = strcmp(output, "-")? fopen(output, "w") : stdout;
	fprintf(stderr, "[%s] Output index to %s\n", date(), strcmp(output, "-")? output : "STDOUT");
	prepare_mem_dump_cuhash(wt->rdname2id);
	mem_dump_obj_file(wt, &wtidx_obj_desc, 1, WTIDX_MAGIC_NUM, out);
	if(strcmp(output, "-")) fclose(out);
	fprintf(stderr, "[%s] Done\n", date());
	free_wtidx(wt);
	free_cplist(pbs);
	free_cplist(obts);
	return 0;
}

