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

#include "file_reader.h"

int usage(){
	printf(
	"WTPRE: Prepare raw reads for assembly\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtpre [options] <raw_reads_file:fq/fa>\n"
	"Options:\n"
	" -o <string> Output of processed reads, [-]\n"
	" -f          Force overwrite output file\n"
	" -L          Keep all subreads in a well, default: the longest one\n"
	" -J <int>    Jack knife of read length, [0]\n"
	" -c <int>    Clip <-c> bases at both ends, [0]\n"
	" -p <string> Change the read name into {\"%%s%%012d\", <-p>}, [pb]\n"
	"\n"
	"Example: \n"
	"$> wtpre -J 5000 -p pb my_raw_reads_1.fq my_raw_reads_2.fq >wt.fa\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	Sequence *seq;
	String *lst_tag, *lst_dsc, *lst_seq;
	char *seqstr;
	int longest, min_len, clp_len, seqlen, c, overwrite, max;
	unsigned long long idx;
	char *prefix, *outf;
	FILE *out;
	longest = 1;
	min_len = 0;
	clp_len = 0;
	overwrite = 0;
	prefix = "pb";
	outf = NULL;
	while((c = getopt(argc, argv, "ho:fLJ:c:p:")) >= 0){
		switch(c){
			case 'h': return usage();
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'L': longest = 0; break;
			case 'J': min_len = atoi(optarg); break;
			case 'c': clp_len = atoi(optarg); break;
			case 'p': prefix = optarg; break;
			default: return usage();
		}
	}
	if(optind == argc) return usage();
	if(!overwrite && outf && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	fr = fopen_m_filereader(argc - optind, argv + optind);
	seq = NULL;
	idx = 0;
	lst_tag = init_string(64);
	lst_dsc = init_string(64);
	lst_seq = init_string(1024);
	max = 0;
	out = outf? ((strcmp(outf, "-") == 0)? stdout : fopen(outf, "w")) : stdout;
	int print;
	print = 0;
	int wushigang = print;
	int tmp = wushigang;
	wushigang = tmp;
	while(fread_seq(&seq, fr)){
		seqstr = seq->seq.string + clp_len;
		seqlen = seq->seq.size - 2 * clp_len;
		if(seqlen < min_len) continue;
		seqstr[seqlen] = '\0';
		if(longest){
			int size, f;
			size = seq->tag.size;
			f = 0;
			while(size){
				if(seq->tag.string[size-1] <= '9' && seq->tag.string[size-1] >= '0'){
					size --;
				} else if(seq->tag.string[size-1] == '_'){
					if(f){
						break;
					} else {
						size --; f = 1;
					}
				} else if(seq->tag.string[size-1] == '/'){
					if(f == 1){
						size --; f = 2; break;
					} else {
						break;
					}
				} else break;
			}
			if(size <= 0 || f < 2) size = seq->tag.size;
			if(lst_tag->size == size && strncmp(lst_tag->string, seq->tag.string, size) == 0){
				if(seqlen > max){
					clear_string(lst_tag); append_string(lst_tag, seq->tag.string, size);
					clear_string(lst_dsc); append_string(lst_dsc, seq->header.string + seq->tag.size, seq->header.size - seq->tag.size);
					clear_string(lst_seq); append_string(lst_seq, seqstr, seqlen);
					max = seq->seq.size;
				}
			} else {
				if(lst_tag->size) fprintf(out, ">%s%012llu%s\n%s\n", prefix, idx ++, lst_dsc->string, lst_seq->string);
				clear_string(lst_tag); append_string(lst_tag, seq->tag.string, size);
				clear_string(lst_dsc); append_string(lst_dsc, seq->header.string + seq->tag.size, seq->header.size - seq->tag.size);
				clear_string(lst_seq); append_string(lst_seq, seqstr, seqlen);
				max = seq->seq.size;
			}
		} else if(seqlen >= min_len){
			fprintf(out, ">%s%012llu%s\n%s\n", prefix, idx ++, seq->header.string + seq->tag.size, seqstr);
		}
	}
	if(lst_tag->size) fprintf(out, ">%s%012llu%s\n%s\n", prefix, idx ++, lst_dsc->string, lst_seq->string);
	free_string(lst_tag);
	free_string(lst_dsc);
	free_string(lst_seq);
	fclose_filereader(fr);
	if(out != stdout) fclose(out);
	return 0;
}

