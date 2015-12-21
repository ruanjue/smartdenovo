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
#include "kswx.h"
#include "file_reader.h"

int usage(){
	printf(
	"Program: pairaln\n"
	" pairaln read paired sequences (fasta or fastq) from STDIN, here is an example input\n"
	" >read1\n"
	" AAGGCCTT\n"
	" >read2\n"
	" AAGCCTT\n"
	" and so on, read3, read4, ...\n"
	" pairaln will perform alignment on read1 and read2, read3 and read4, ...\n"
	" The default parameters is used for pacbio long reads\n"
	"Usage: pairaln [options]\n"
	"Author: Jue Ruan\n"
	"Options:\n"
	" -s          Try both strands\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: insertion or deletion, [-3]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping, 0: distable HSP extension, otherwise set to -100 or other [-100]\n"
	" -W <int>    Bandwidth, [800]\n"
	" -a          Output alignment\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	String *alns[3];
	Sequence *seq1, *seq2;
	u8list *s1, *s2;
	kswx_t x;
	int W, M, X, O, E, T, bidir, has_aln, d, c, i;
	int n_cigar;
	uint32_t *cigar;
	int8_t matrix[4 * 4];
	W = 800;
	M = 2;
	X = -5;
	O = -3;
	E = -1;
	T = -100;
	bidir = 1;
	has_aln = 0;
	while((c = getopt(argc, argv, "hsW:M:X:O:E:T:a")) != -1){
		switch(c){
			case 'h': return usage();
			case 's': bidir = 3; break;
			case 'a': has_aln = 1; break;
			case 'W': W = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			default: return usage();
		}
	}
	for(i=0;i<4*4;i++) matrix[i] = ((i % 4) == (i / 4))? M : X;
	if(optind < argc){
		fr = fopen_m_filereader(argc - optind, argv + optind);
	} else fr = stdin_filereader();
	seq1 = seq2 = NULL;
	s1 = init_u8list(1024);
	s2 = init_u8list(1024);
	alns[0] = init_string(1024);
	alns[1] = init_string(1024);
	alns[2] = init_string(1024);
	while(1){
		if(! fread_seq(&seq1, fr) || ! fread_seq(&seq2, fr)) break;
		clear_and_encap_u8list(s1, seq1->seq.size);
		clear_and_encap_u8list(s2, seq2->seq.size);
		for(i=0;i<seq1->seq.size;i++) s1->buffer[i] = base_bit_table[(int)seq1->seq.string[i]];
		for(i=0;i<seq2->seq.size;i++) s2->buffer[i] = base_bit_table[(int)seq2->seq.string[i]];
		s1->size = seq1->seq.size;
		s2->size = seq2->seq.size;
		for(d=0;d<2;d++){
			if(((bidir >> d) & 0x01) == 0)  continue;
			if(d){ reverse_u8list(s2); for(i=0;(size_t)i<s2->size;i++) s2->buffer[i] = (~s2->buffer[i]) & 0x03; }
			n_cigar = 0;
			cigar = NULL;
			x = kswx_align_with_cigar(s1->size, s1->buffer, s2->size, s2->buffer, 4, matrix, W, O, O, E, T, &n_cigar, &cigar);
			fprintf(stdout, "%s\t%c\t%d\t%d\t%d"  , seq1->name.string, '+', seq1->seq.size, x.qb, x.qe);
			fprintf(stdout, "\t%s\t%c\t%d\t%d\t%d", seq2->name.string, "+-"[d], seq2->seq.size, x.tb, x.te);
			fprintf(stdout, "\t%d\t%0.2f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
			if(has_aln){
				clear_string(alns[0]);
				clear_string(alns[1]);
				clear_string(alns[2]);
				kswx_cigar2pairwise((String*[2]){alns[0], alns[1]}, (uint8_t*[2]){s1->buffer + x.qb, s2->buffer + x.tb}, n_cigar, cigar);
				for(i=0;i<alns[0]->size;i++){
					if(alns[0]->string[i] == '-'){
						add_char_string(alns[2], '-');
					} else if(alns[0]->string[i] == alns[1]->string[i]){
						add_char_string(alns[2], '|');
					} else if(alns[1]->string[i] == '-'){
						add_char_string(alns[2], '-');
					} else {
						add_char_string(alns[2], '*');
					}
				}
				fprintf(stdout, "%s\n%s\n%s\n", alns[0]->string, alns[2]->string, alns[1]->string);
			}
			if(cigar) free(cigar);
		}
	}
	fclose_filereader(fr);
	free_u8list(s1);
	free_u8list(s2);
	free_string(alns[0]);
	free_string(alns[1]);
	free_string(alns[2]);
	return 0;
}
