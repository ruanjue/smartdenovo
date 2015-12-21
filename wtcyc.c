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

#include "kswx.h"
#include "dna.h"
#include "list.h"
#include "file_reader.h"
#include "thread.h"

thread_beg_def(mcyc);
String *n;
u8list *t, *q;
int8_t *matrix;
kswx_t x;
int W, M, X, O, E, T;
thread_end_def(mcyc);

thread_beg_func(mcyc);
int i;
thread_beg_loop(mcyc);
clear_u8list(mcyc->q); append_u8list(mcyc->q, mcyc->t);
reverse_u8list(mcyc->q);
for(i=0;i<(int)mcyc->q->size;i++) mcyc->q->buffer[i] = (~mcyc->q->buffer[i]) & 0x03;
mcyc->x = kswx_align(mcyc->q->size, mcyc->q->buffer, mcyc->t->size, mcyc->t->buffer, 4, mcyc->matrix, mcyc->W, mcyc->O, mcyc->O, mcyc->E, mcyc->T);
thread_end_loop(mcyc);
thread_end_func(mcyc);

int usage(){
	printf(
	"WTCYC: Align long read against its reverse complementary\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtcyc [options] <long_read_file>\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -P <int>    Total parallel jobs, [1]\n"
	" -p <int>    Index of current job (0-based), [0]\n"
	"             Suppose to run it parallelly in 60 nodes. For node1, -P 60 -p 0; node2, -P 60 -p 1, ...\n"
	" -o <string> Output of reads' regions after trimming, [-]\n"
	" -a <string> Output of alignments, [NULL]\n"
	" -f          Force overwrite output file\n"
	" -s <int>    Mininum alignment score, [400]\n"
	" -m <int>    Mininum alignment identity, [0.7]\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: gap open, [-3]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping, 0: distable HSP extension, otherwise set to -30 or other [-100]\n"
	" -W <int>    Bandwidth, [800]\n"
	"\n"
	"Example: \n"
	"$> wtcyc -t 32 wt.fa -fo wt.zmo.cyc -a wt.zmo.cyc.info\n"
	"\n"
	);
	return 1;
}

void output_alignment(String *name, int seqlen, kswx_t x, int ms, float mi, FILE *alno, FILE *cbto){
	int bp;
	float identity;
	identity = (1.0 * x.mat) / (x.aln + 1);
	if(alno) fprintf(alno, "%s\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", name->string, seqlen, x.score, identity, x.tb, x.te, seqlen - x.qe, seqlen - x.qb);
	if(x.score >= ms && identity >= mi){
		if((x.tb == 0 || x.te == seqlen) && (num_diff(x.tb, seqlen - x.qe) < 50 && num_diff(x.te, seqlen - x.qb) < 50)){
			bp = (x.tb + x.te) / 2;
			if(bp < seqlen / 2){
				fprintf(cbto, "%s\t%d\t%d\t%d\n", name->string, bp, seqlen - bp, seqlen);
			} else {
				fprintf(cbto, "%s\t%d\t%d\t%d\n", name->string, 0, bp, seqlen);
			}
		}
	}
}

int main(int argc, char **argv){
	FileReader *fr;
	Sequence *seq;
	char *alnf, *cbtf;
	FILE *alno, *cbto;
	int8_t matrix[4 * 4];
	uint64_t n_rd;
	int i, c, ncpu, W, M, X, O, E, T, ms, overwrite, n_job, i_job;
	float mi;
	thread_preprocess(mcyc);
	ncpu = 1;
	n_job = 1;
	i_job = 0;
	W = 800;
	M = 2;
	X = -5;
	O = -3;
	E = -1;
	T = -100;
	ms = 400;
	mi = 0.7;
	alnf = NULL;
	cbtf = NULL;
	overwrite = 0;
	while((c = getopt(argc, argv, "ht:P:p:o:a:fs:m:M:X:O:E:W:T:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'P': n_job = atoi(optarg); break;
			case 'p': i_job = atoi(optarg); break;
			case 'a': alnf = optarg; break;
			case 'o': cbtf = optarg; break;
			case 'f': overwrite = 1; break;
			case 's': ms = atoi(optarg); break;
			case 'm': mi = atof(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'W': W = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			default: return usage();
		}
	}
	if(optind == argc) return usage();
	if(ms < 1) ms = 1;
	if(alnf && !overwrite && strcmp(alnf, "-") && file_exists(alnf)){
		fprintf(stderr, "File exists! '%s'\n\n", alnf);
		return usage();
	}
	if(cbtf && !overwrite && strcmp(cbtf, "-") && file_exists(cbtf)){
		fprintf(stderr, "File exists! '%s'\n\n", cbtf);
		return usage();
	}
	alno = alnf? ((strcmp(alnf, "-") == 0)? stdout : fopen(alnf, "w")) : NULL;
	cbto = cbtf? ((strcmp(cbtf, "-") == 0)? stdout : fopen(cbtf, "w")) : stdout;
	for(i=0;i<4*4;i++) matrix[i] = ((i / 4) == (i % 4))? M : X;
	fr = fopen_m_filereader(argc - optind, argv + optind);
	thread_beg_init(mcyc, ncpu);
	mcyc->n = init_string(1024);
	mcyc->t = init_u8list(1024);
	mcyc->q = init_u8list(1024);
	mcyc->matrix  = matrix;
	mcyc->W = W;
	mcyc->M = M;
	mcyc->X = X;
	mcyc->O = O;
	mcyc->E = E;
	mcyc->T = T;
	mcyc->x = KSWX_NULL;
	mcyc->x.score = -1;
	thread_end_init(mcyc);
	seq = NULL;
	n_rd = 0;
	while(fread_seq(&seq, fr)){
		if(((int)(n_rd ++) % n_job) != i_job) continue;
		thread_waitfor_one_idle(mcyc);
		if(mcyc->x.score >= 0) output_alignment(mcyc->n, mcyc->t->size, mcyc->x, ms, mi, alno, cbto);
		clear_string(mcyc->n); append_string(mcyc->n, seq->name.string, seq->name.size);
		clear_u8list(mcyc->t);
		for(i=0;i<seq->seq.size;i++){ push_u8list(mcyc->t, base_bit_table[(int)seq->seq.string[i]]); }
		mcyc->x.score = -1;
		thread_wake(mcyc);
	}
	fclose_filereader(fr);
	thread_waitfor_all_idle(mcyc);
	thread_beg_close(mcyc);
	if(mcyc->x.score >= 0) output_alignment(mcyc->n, mcyc->t->size, mcyc->x, ms, mi, alno, cbto);
	free_string(mcyc->n);
	free_u8list(mcyc->t);
	free_u8list(mcyc->q);
	thread_end_close(mcyc);
	if(cbto && cbto != stdout) fclose(cbto);
	if(alno && alno != stdout) fclose(alno);
	return 0;
}

