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

#include "dna.h"
#include "timer.h"
#include "file_reader.h"

int usage(){
	printf(
	"WTMER: K-mer analysis of long reads\n"
	"SMARTdenvo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtmer [options]\n"
	"Options:\n"
	" -i <string> Long reads sequences file, + *\n"
	" -o <string> Output file of kmer_frequency, *\n"
	" -f          Force overwrite\n"
	" -H          Disable homopolymer compression\n"
	" -k <int>    Kmer size, 5 <= <-k> <= 16, [16]\n"
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	Sequence *seq;
	uint64_t *cnts, cnt, sum;
	uint16_t *freqs;
	uint64_t kmer, krev, kmask;
	cplist *pbs;
	char *output;
	FILE *out;
	uint32_t pbid, pblen, i, j, z, b, ksize, n_freq;
	int c, hz, overwrite;
	ksize = 16;
	hz = 1;
	n_freq = 256;
	overwrite = 0;
	pbs = init_cplist(4);
	output = NULL;
	while((c = getopt(argc, argv, "hi:o:fHk:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'i': push_cplist(pbs, optarg); break;
			case 'o': output = optarg; break;
			case 'f': overwrite = 1; break;
			case 'H': hz = 0; break;
			case 'k': ksize = atoi(optarg); break;
			default: return usage();
		}
	}
	if(output == NULL) return usage();
	if(!overwrite && strcmp(output, "-") && file_exists(output)){
		fprintf(stderr, "File exists! '%s'\n\n", output);
		return usage();
	}
	if(pbs->size == 0) return usage();
	if(ksize > 16 || ksize < 5) return usage();
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
	freqs = calloc(kmask, sizeof(uint16_t));
	fr = fopen_m_filereader(pbs->size, pbs->buffer);
	fprintf(stderr, "[%s] loading long reads\n", date());
	seq = NULL;
	pbid = 0;
	while(fread_seq(&seq, fr)){
		if((pbid % 1000) == 0){ fprintf(stderr, "\r%u", pbid); fflush(stderr); }
		pbid ++;
		pblen = seq->seq.size;
		b = 4;
		kmer = 0;
		for(i=j=0;j<pblen;j++){
			z = base_bit_table[(int)seq->seq.string[j]];
			if(hz && z == b) continue;
			b = z;
			i ++;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < ksize) continue;
			krev = dna_rev_seq(kmer, ksize);
			krev = krev > kmer? kmer : krev;
			if(freqs[krev] < 0xFFFFU) freqs[krev] ++;
		}
	}
	fclose_filereader(fr);
	fprintf(stderr, "\r%u\n", pbid); fflush(stderr);
	fprintf(stderr, "[%s] Done\n", date());
	cnts = calloc(0xFFFFU, sizeof(uint64_t));
	for(kmer=0;kmer<kmask;kmer++) cnts[freqs[kmer]] ++;
	for(i=1;i<=n_freq;i++){
		fprintf(stderr, "%u\t%llu\n", i, (unsigned long long)cnts[i]);
	}
	cnt = sum = 0;
	while(i < 0xFFFFU){
		cnt += cnts[i];
		sum += cnts[i] * i;
		i ++;
	}
	fprintf(stderr, "+\t%llu\t%llu\n", (unsigned long long)cnt, (unsigned long long)sum);
	free(cnts);
	out = strcmp(output, "-")? fopen(output, "w") : stdout;
	for(kmer=0;kmer<kmask;kmer++){
		if(freqs[kmer] == 0) continue;
		fprintf(out, "0x%X\t%u\n", (uint32_t)kmer, freqs[kmer]);
	}
	if(strcmp(output, "-")) fclose(out);
	fprintf(stderr, "[%s] Done\n", date());
	free_cplist(pbs);
	return 0;
}

