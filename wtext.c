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
#include "list.h"
#include "hashset.h"
#include "file_reader.h"
#include "kswx.h"
#include "thread.h"

typedef struct {
	uint32_t pb1:31, dir1:1, pb2:31, dir2:1;
	int qb, qe, tb, te;
	int score, aln, mat, mis, ins, del;
	char *cigar;
} wt_ovl_t;
define_list(wtovlv, wt_ovl_t);

typedef struct {
	uint32_t n_pb;
	BaseBank *pbseqs;
	u32list  *prev_clp_offs;
	u32list  *prev_clp_lens;
	u32list  *clp_offs;
	u32list  *clp_lens;
	u32list  *pblens;
	u64list  *pboffs;
	cplist   *pbnames;
	cuhash   *pbname2id;
	int      W, M, X, O, E, T;
	int8_t   matrix[4 * 4];
	int      min_score;
	int      max_ext;
	float    min_id;
} WTEXT;

WTEXT* init_wtext(int W, int M, int X, int O, int E, int T, int min_score, float min_id, int max_ext){
	WTEXT *wt;
	int i;
	wt = malloc(sizeof(WTEXT));
	wt->pbseqs = init_basebank();
	wt->prev_clp_offs = init_u32list(1024);
	wt->prev_clp_lens = init_u32list(1024);
	wt->clp_offs = init_u32list(1024);
	wt->clp_lens = init_u32list(1024);
	wt->pblens = init_u32list(1024);
	wt->pboffs = init_u64list(1024);
	wt->pbnames = init_cplist(1024);
	wt->pbname2id = init_cuhash(1023);
	wt->n_pb = 0;
	wt->W = W;
	wt->M = M;
	wt->X = X;
	wt->O = O;
	wt->E = E;
	wt->T = T;
	for(i=0;i<4*4;i++) wt->matrix[i] = ((i % 4) == (i / 4))? M : X;
	wt->min_score = min_score;
	wt->min_id = min_id;
	wt->max_ext = max_ext;
	return wt;
}

void free_wtext(WTEXT *wt){
	uint32_t i;
	free_basebank(wt->pbseqs);
	free_u32list(wt->prev_clp_offs);
	free_u32list(wt->prev_clp_lens);
	free_u32list(wt->clp_offs);
	free_u32list(wt->clp_lens);
	free_u32list(wt->pblens);
	free_u64list(wt->pboffs);
	free_cuhash(wt->pbname2id);
	for(i=0;i<wt->pbnames->size;i++) free(get_cplist(wt->pbnames, i));
	free_cplist(wt->pbnames);
	free(wt);
}

void push_read_wtext(WTEXT *wt, char *name, int name_len, char *seq, int seq_len){
	char *ptr;
	push_u32list(wt->pblens, seq_len);
	push_u32list(wt->prev_clp_offs, 0);
	push_u32list(wt->prev_clp_lens, seq_len);
	push_u32list(wt->clp_offs, 0);
	push_u32list(wt->clp_lens, seq_len);
	push_u64list(wt->pboffs, wt->pbseqs->size);
	seq2basebank(wt->pbseqs, seq, seq_len);
	ptr = malloc(name_len + 1);
	memcpy(ptr, name, name_len);
	ptr[name_len] = 0;
	push_cplist(wt->pbnames, ptr);
	kv_put_cuhash(wt->pbname2id, ptr, wt->n_pb);
	wt->n_pb ++;
}

void set_read_previous_clip_wtext(WTEXT *wt, char *name, int coff, int clen){
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->pbname2id, name)) == 0xFFFFFFFFU) return;
	if(coff < 0 || coff + clen > (int)wt->pblens->buffer[pbid]) return;
	set_u32list(wt->prev_clp_offs, pbid, coff);
	set_u32list(wt->prev_clp_lens, pbid, clen);
}

void set_read_clip_wtext(WTEXT *wt, char *name, int coff, int clen){
	uint32_t pbid;
	if((pbid = kv_get_cuhash(wt->pbname2id, name)) == 0xFFFFFFFFU) return;
	if(coff < 0 || coff + clen > (int)wt->pblens->buffer[pbid]) return;
	set_u32list(wt->clp_offs, pbid, coff);
	set_u32list(wt->clp_lens, pbid, clen);
}

void extending_overlap_wtext(WTEXT *wt, wt_ovl_t *hit, u8list *pb1, u8list *pb2, u32list *cigar){
	kswx_t x0, x1;
	int seqlens[2], clpoffs[2], clplens[2];
	int clp[2], x[2], y[2], dx[2], dy[2], cx[2], cy[2], nx[2], ny[2], i, j, d;
	int len, op;
	int nc;
	uint32_t *cs;
	String *str;
	seqlens[0] = wt->pblens->buffer[hit->pb1];
	clpoffs[0] = wt->clp_offs->buffer[hit->pb1];
	clplens[0] = wt->clp_lens->buffer[hit->pb1];
	seqlens[1] = wt->pblens->buffer[hit->pb2];
	clpoffs[1] = wt->clp_offs->buffer[hit->pb2];
	clplens[1] = wt->clp_lens->buffer[hit->pb2];
	clear_and_encap_u8list(pb1, clplens[0]);
	if(hit->dir1) revbitseq_basebank(wt->pbseqs, wt->pboffs->buffer[hit->pb1] + clpoffs[0], clplens[0], pb1->buffer);
	else             bitseq_basebank(wt->pbseqs, wt->pboffs->buffer[hit->pb1] + clpoffs[0], clplens[0], pb1->buffer);
	clear_and_encap_u8list(pb2, clplens[1]);
	if(hit->dir2) revbitseq_basebank(wt->pbseqs, wt->pboffs->buffer[hit->pb2] + clpoffs[1], clplens[1], pb2->buffer);
	else             bitseq_basebank(wt->pbseqs, wt->pboffs->buffer[hit->pb2] + clpoffs[1], clplens[1], pb2->buffer);
	clear_u32list(cigar);
	kswx_string2cigar(cigar, hit->cigar);
	free(hit->cigar); hit->cigar = NULL;
	str = init_string(32);
	x0 = KSWX_NULL;
	clp[0] = clp[1] = 0;
	x[0] = hit->tb;
	x[1] = hit->qb;
	y[0] = seqlens[0] - hit->te;
	y[1] = seqlens[1] - hit->qe;
	dx[0] = hit->dir1? seqlens[0] - clpoffs[0] - clplens[0] : clpoffs[0];
	dx[1] = hit->dir2? seqlens[1] - clpoffs[1] - clplens[1] : clpoffs[1];
	dy[0] = hit->dir1? clpoffs[0] : seqlens[0] - clpoffs[0] - clplens[0];
	dy[1] = hit->dir2? clpoffs[1] : seqlens[1] - clpoffs[1] - clplens[1];
	cx[0] = dx[0] > x[0]? dx[0] - x[0] : 0;
	cx[1] = dx[1] > x[1]? dx[1] - x[1] : 0;
	cy[0] = dy[0] > y[0]? dy[0] - y[0] : 0;
	cy[1] = dy[1] > y[1]? dy[1] - y[1] : 0;
	nx[0] = nx[1] = 0;
	ny[0] = ny[1] = 0;
	while(clp[0] < (int)cigar->size){
		op  = cigar->buffer[clp[0]] & 0xFU;
		len = cigar->buffer[clp[0]] >> 4;
		if(op == 1){
			// pb2 ++
			nx[1] += len;
		} else if(op == 2){
			// pb1 ++
			nx[0] += len;
		} else {
			// pb1 ++ and pb2 ++
			if(nx[0] >= cx[0] && nx[1] >= cx[1]) break;
			d = cx[0] - nx[0] > cx[1] - nx[1]? cx[0] - nx[0] : cx[1] - nx[1];
			d = d > len? len : d;
			nx[0] += d;
			nx[1] += d;
			if(d < len){
				len -= d;
				cigar->buffer[clp[0]] = (len << 4) | op; break;
			}
		}
		clp[0] ++;
	}
	if(nx[0] < cx[0] || nx[1] < cx[1]) goto RET;
	while(clp[1] < (int)cigar->size){
		op  = cigar->buffer[cigar->size - 1 - clp[1]] & 0xFU;
		len = cigar->buffer[cigar->size - 1 - clp[1]] >> 4;
		if(op == 1){
			// pb2 ++
			ny[1] += len;
		} else if(op == 2){
			// pb1 ++
			ny[0] += len;
		} else {
			// pb1 ++ and pb2 ++
			if(ny[0] >= cy[0] && ny[1] >= cy[1]) break;
			d = cy[0] - ny[0] > cy[1] - ny[1]? cy[0] - ny[0] : cy[1] - ny[1];
			d = d > len? len : d;
			ny[0] += d;
			ny[1] += d;
			if(d < len){
				len -= d;
				cigar->buffer[cigar->size - 1 - clp[1]] = (len << 4) | op; break;
			}
		}
		clp[1] ++;
	}
	if(ny[0] < cy[0] || ny[1] < cy[1]) goto RET;
	if(clp[0] + clp[1] >= (int)cigar->size) goto RET;
	x0.tb = hit->tb + nx[0] - dx[0];
	x0.qb = hit->qb + nx[1] - dx[1];
	x0.te = hit->te - ny[0] - dx[0];
	x0.qe = hit->qe - ny[1] - dx[1];
	cx[0] = x0.tb;
	cx[1] = x0.qb;
	x0.score = 0;
	for(i=clp[0];i+clp[1]<(int)cigar->size;i++){
		op  = cigar->buffer[i] & 0xFU;
		len = cigar->buffer[i] >> 4;
		x0.aln += len;
		if(op == 1){
			x0.ins += len;
			cx[1] += len;
			x0.score += wt->O + wt->E * len;
		} else if(op == 2){
			x0.del += len;
			cx[0] += len;
			x0.score += wt->O + wt->E * len;
		} else {
			for(j=0;j<len;j++){
				if(pb1->buffer[cx[0] + j] == pb2->buffer[cx[1] + j]) x0.mat ++;
				else x0.mis ++;
			}
			cx[0] += len;
			cx[1] += len;
		}
	}
	x0.score += x0.mat * wt->M;
	x0.score += x0.mis * wt->X;
	// left extension
	if(x0.qb <= wt->max_ext || x0.tb <= wt->max_ext){
		nc = 0; cs = NULL;
		x1 = kswx_extend_align(x0.qb, pb2->buffer + x0.qb - 1, x0.tb, pb1->buffer + x0.tb - 1, -1, x0.score, wt->W, wt->M, wt->X, wt->O, wt->O, wt->E, wt->T, &nc, &cs);
		x0.score  = x1.score;
		x0.aln   += x1.aln;
		x0.mat   += x1.mat;
		x0.mis   += x1.mis;
		x0.ins   += x1.ins;
		x0.del   += x1.del;
		x0.qb    -= x1.qe;
		x0.tb    -= x1.te;
		revseq_4bytes(cs, nc);
		kswx_cigar2string(str, nc, cs);
		if(cs) free(cs);
	}
	// append core
	kswx_cigar2string(str, cigar->size - clp[0] - clp[1], cigar->buffer + clp[0]);
	// right extension
	if(clplens[1] - x0.qe <= wt->max_ext || clplens[0] - x0.te <= wt->max_ext){
		nc = 0; cs = NULL;
		x1 = kswx_extend_align(clplens[1] - x0.qe, pb2->buffer + x0.qe, clplens[0] - x0.te, pb1->buffer + x0.te, 1, x0.score, wt->W, wt->M, wt->X, wt->O, wt->O, wt->E, wt->T, &nc, &cs);
		x0.score  = x1.score;
		x0.aln   += x1.aln;
		x0.mat   += x1.mat;
		x0.mis   += x1.mis;
		x0.ins   += x1.ins;
		x0.del   += x1.del;
		x0.qe    += x1.qe;
		x0.te    += x1.te;
		kswx_cigar2string(str, nc, cs);
		if(cs) free(cs);
	}
	RET:
	hit->score = x0.score;
	hit->qb = x0.qb;
	hit->qe = x0.qe;
	hit->tb = x0.tb;
	hit->te = x0.te;
	hit->aln = x0.aln;
	hit->mat = x0.mat;
	hit->mis = x0.mis;
	hit->ins = x0.ins;
	hit->del = x0.del;
	hit->cigar = str->string;
	free(str); // NB: Cannot free_string
}

thread_beg_def(mext);
WTEXT *wt;
wtovlv *hits;
thread_end_def(mext);

thread_beg_func(mext);
uint32_t i;
u8list *ex1, *ex2;
u32list *cigar;
ex1 = init_u8list(1024);
ex2 = init_u8list(1024);
cigar = init_u32list(1024);
thread_beg_loop(mext);
for(i=0;i<mext->hits->size;i++){
	//wt_ovl_t *hit = ref_wtovlv(mext->hits, i);
	//fprintf(stdout, "%s\t%c\t%d\t%d\t%d", get_cplist(mext->wt->pbnames, hit->pb1), "+-"[hit->dir1], mext->wt->clp_lens->buffer[hit->pb1], hit->tb, hit->te);
	//fprintf(stdout, "\t%s\t%c\t%d\t%d\t%d\n", get_cplist(mext->wt->pbnames, hit->pb2), "+-"[hit->dir2], mext->wt->clp_lens->buffer[hit->pb2], hit->qb, hit->qe);
	//fflush(stdout);
	extending_overlap_wtext(mext->wt, ref_wtovlv(mext->hits, i), ex1, ex2, cigar);
}
thread_end_loop(mext);
free_u8list(ex1);
free_u8list(ex2);
free_u32list(cigar);
thread_end_func(mext);

uint32_t output_alignments_wtext(WTEXT *wt, wtovlv *hits, FILE *out){
	wt_ovl_t *hit;
	uint32_t i, ret;
	ret = hits->size;
	for(i=0;i<hits->size;i++){
		hit = ref_wtovlv(hits, i);
		if(hit->aln > 0){
			fprintf(out, "%s\t%c\t%d\t%d\t%d", get_cplist(wt->pbnames, hit->pb1), "+-"[hit->dir1], wt->clp_lens->buffer[hit->pb1], hit->tb, hit->te);
			fprintf(out, "\t%s\t%c\t%d\t%d\t%d", get_cplist(wt->pbnames, hit->pb2), "+-"[hit->dir2], wt->clp_lens->buffer[hit->pb2], hit->qb, hit->qe);
			fprintf(out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%s\n", hit->score, 1.0 * hit->mat / hit->aln, hit->mat, hit->mis, hit->ins, hit->del, hit->cigar);
		}
		if(hit->cigar){ free(hit->cigar); hit->cigar = NULL; }
	}
	clear_wtovlv(hits);
	return ret;
}


int usage(){
	printf(
	"WTEXT: Extending and clipping overlaps\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtext [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -P <int>    Total parallel jobs, [1]\n"
	" -p <int>    Index of current job (0-based), [0]\n"
	"             Suppose to run it parallelly in 60 nodes. For node1, -P 60 -p 0; node2, -P 60 -p 1, ...\n"
	" -i <string> Long reads sequences file, + *\n"
	" -B <string> Long reads previous retained region, often from wtcyc, +\n"
	"             Format: read_name\\toffset\\tlength\n"
	" -b <string> Long reads retained region, often from wtobt, +\n"
	"             Format: read_name\\toffset\\tlength\n"
	" -j <string> Overlap file(s), + *\n"
	"             Format: reads1\\t+/-\\tlen1\\tbeg1\\tend1\\treads2\\t+/-\\tlen2\\tbeg2\\tend2\\tscore\\tidentity<float>\\tmat\\tmis\\tins\\tdel\\tcigar\n"
	" -o <string> Output file of extended alignments, -:stdout, *\n"
	" -f          Force overwrite\n"
	" -W <float>  Bandwidth, [800]\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: insertion or deletion, [-3]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping [-100]\n"
	" -S <int>    Maximum extension (bp) in each end, [400]\n"
	//" -s <int>    Minimum alignment score, [100]\n"
	//" -m <float>  Minimum alignment identity, [50]\n"
	"\n"
	"Example: \n"
	"$> wtext -t 32 -i wt.fa -b wt.zmo.obt -j wt.zmo.ovl -o wt.zmo.ext\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	obj_desc_t wsg = wtovlv_obj_desc;
	obj_desc_t ttt = wsg;
	wsg = ttt;
	WTEXT *wt;
	FileReader *fr;
	Sequence *seq;
	cplist *pbs, *ovls, *obts, *cycs;
	wt_ovl_t HIT;
	uint32_t i, idx, batch_size;
	unsigned long long n, nb;
	char *outf;
	FILE *out;
	int c, f, W, M, X, O, E, T, force_overwrite, min_score, max_ext, ncpu, n_job, i_job;
	float min_id;
	thread_preprocess(mext);
	pbs = init_cplist(4);
	ovls = init_cplist(4);
	obts = init_cplist(4);
	cycs = init_cplist(4);
	ncpu = 1;
	n_job = 1;
	i_job = 0;
	W = 800; M = 2; X = -5; O = -3; E = -1; T = -100;
	min_score = 200;
	max_ext = 400;
	min_id = 0.5;
	batch_size = 100;
	force_overwrite = 0;
	outf = NULL;
	while((c = getopt(argc, argv, "hft:P:p:i:B:b:j:o:W:M:X:O:E:T:s:m:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'f': force_overwrite = 1; break;
			case 't': ncpu = atoi(optarg); break;
			case 'P': n_job = atoi(optarg); break;
			case 'p': i_job = atoi(optarg); break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'B': push_cplist(cycs, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			case 'j': push_cplist(ovls, optarg); break;
			case 'o': outf = optarg; break;
			case 'W': W = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			case 'S': max_ext = atoi(optarg); break;
			case 's': min_score = atoi(optarg); break;
			case 'm': min_id    = atof(optarg) / 100.0; break;
			default: return usage();
		}
	}
	if(pbs->size == 0 || ovls->size == 0) return usage();
	if(outf == NULL) return usage();
	if(!force_overwrite && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	wt = init_wtext(W, M, X, O, E, T, min_score, min_id, max_ext);
	if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	fprintf(stderr, "[%s] loading long reads\n", date());
	seq = NULL;
	while(fread_seq(&seq, fr)){
		push_read_wtext(wt, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
	}
	fclose_filereader(fr);
	fprintf(stderr, "[%s] Done, %u reads\n", date(), (unsigned)wt->n_pb);
	if(cycs->size){
		fprintf(stderr, "[%s] loading reads previous obt information\n", date());
		if((fr = fopen_m_filereader(cycs->size, cycs->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(c < 3) continue;
			set_read_previous_clip_wtext(wt, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date());
	}
	if(obts->size){
		fprintf(stderr, "[%s] loading reads obt information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(fr->line->string[0] == '#') continue;
			if(c < 3) continue;
			set_read_clip_wtext(wt, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date());
	}
	out = strcmp(outf, "-")? fopen(outf, "w") : stdout;
	if((fr = fopen_m_filereader(ovls->size, ovls->buffer)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", ovls->buffer[0], __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	thread_beg_init(mext, ncpu);
	mext->wt = wt;
	mext->hits = init_wtovlv(1024);
	thread_end_init(mext);
	fprintf(stderr, "[%s] extending alignments\n", date());
	n = nb = 0;
	memset(&HIT, 0, sizeof(wt_ovl_t));
	while(1){
		thread_waitfor_one_idle(mext);
		output_alignments_wtext(wt, mext->hits, out);
		clear_wtovlv(mext->hits);
		f = ((int)((nb ++) % n_job) == i_job);
		for(i=0;i<batch_size;i++){
			if((c = fread_table(fr)) == -1) break;
			if(f == 0) continue;
			if(fr->line->string[0] == '#') continue;
			if(c < 17) continue;
			if((idx = kv_get_cuhash(wt->pbname2id, get_col_str(fr, 0))) == 0xFFFFFFFFU) continue;
			HIT.pb1 = idx;
			HIT.dir1 = (get_col_str(fr, 1)[0] == '-');
			if(HIT.dir1){
				HIT.tb = atoi(get_col_str(fr, 3)) + wt->pblens->buffer[HIT.pb1] - (wt->prev_clp_offs->buffer[HIT.pb1] + wt->prev_clp_lens->buffer[HIT.pb1]);
				HIT.te = atoi(get_col_str(fr, 4)) + wt->pblens->buffer[HIT.pb1] - (wt->prev_clp_offs->buffer[HIT.pb1] + wt->prev_clp_lens->buffer[HIT.pb1]);
			} else {
				HIT.tb = atoi(get_col_str(fr, 3)) + wt->prev_clp_offs->buffer[HIT.pb1];
				HIT.te = atoi(get_col_str(fr, 4)) + wt->prev_clp_offs->buffer[HIT.pb1];
			}
			if((idx = kv_get_cuhash(wt->pbname2id, get_col_str(fr, 5))) == 0xFFFFFFFFU) continue;
			HIT.pb2 = idx;
			HIT.dir2 = (get_col_str(fr, 6)[0] == '-');
			if(HIT.dir2){
				HIT.qb = atoi(get_col_str(fr, 8)) + wt->pblens->buffer[HIT.pb2] - (wt->prev_clp_offs->buffer[HIT.pb2] + wt->prev_clp_lens->buffer[HIT.pb2]);
				HIT.qe = atoi(get_col_str(fr, 9)) + wt->pblens->buffer[HIT.pb2] - (wt->prev_clp_offs->buffer[HIT.pb2] + wt->prev_clp_lens->buffer[HIT.pb2]);
			} else {
				HIT.qb = atoi(get_col_str(fr, 8)) + wt->prev_clp_offs->buffer[HIT.pb2];
				HIT.qe = atoi(get_col_str(fr, 9)) + wt->prev_clp_offs->buffer[HIT.pb2];
			}
			HIT.score = atoi(get_col_str(fr, 10));
			HIT.mat   = atoi(get_col_str(fr, 12));
			HIT.mis   = atoi(get_col_str(fr, 13));
			HIT.ins   = atoi(get_col_str(fr, 14));
			HIT.del   = atoi(get_col_str(fr, 15));
			HIT.aln   = 0;
			HIT.cigar = strdup(get_col_str(fr, 16));
			push_wtovlv(mext->hits, HIT);
		}
		n += mext->hits->size;
		if(f){
			fprintf(stderr, "\r%llu", n); fflush(stderr);
			thread_wake(mext);
		}
		if(i < batch_size) break;
	}
	fprintf(stderr, "\r%llu", n); fflush(stderr);
	fclose_filereader(fr);
	thread_waitfor_all_idle(mext);
	thread_beg_close(mext);
	output_alignments_wtext(wt, mext->hits, out);
	free_wtovlv(mext->hits);
	thread_end_close(mext);
	fprintf(stderr, "\r[%s] Done, %llu\n", date(), n);
	free_wtext(wt);
	free_cplist(pbs);
	free_cplist(ovls);
	free_cplist(obts);
	if(out != stdout) fclose(out);
	return 0;
}
