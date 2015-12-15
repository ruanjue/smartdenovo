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
#include "file_reader.h"

typedef struct {
	int x, y;
} pw_mat_t;

define_list(pwmatv, pw_mat_t);
define_list(kswxv, kswx_t);
define_list(strpv, String*);

typedef struct {
	String *name;
	uint32_t n_rd;
	u8list  *rdseqs;
	u64list *seqoffs;
	u32list *seqlens;
	u8list  *namestr;
	u32list *seqtags;
	u32list *rdranks;
	BitVec  *confids;
	BitVec  *varmaps;
	kswxv   *swxs;
	pwmatv  *mats;
	vplist  *cigs;
	uint32_t hlen;
	u8list  *hmsa;
	uint32_t vlen;
	u8list  *vmsa;

	int W, M, X, O, E;

} WTDIF;

WTDIF* init_wtdif(int W, int M, int X, int O, int E){
	WTDIF *g;
	g = malloc(sizeof(WTDIF));
	g->name = init_string(1024);
	g->n_rd = 0;
	g->rdseqs  = init_u8list(1024);
	g->seqoffs = init_u64list(1024);
	g->seqlens = init_u32list(1024);
	g->namestr = init_u8list(1024);
	g->seqtags = init_u32list(1024);
	g->rdranks = init_u32list(1024);
	g->confids = init_bitvec(1024);
	g->varmaps = init_bitvec(1024);
	g->swxs = init_kswxv(32);
	g->mats = init_pwmatv(1024);
	g->cigs = init_vplist(32);
	g->hlen = 0;
	g->hmsa = init_u8list(1024);
	g->vlen = 0;
	g->vmsa = init_u8list(1024);
	g->W = W;
	g->M = M;
	g->X = X;
	g->O = O;
	g->E = E;
	return g;
}

void free_wtdif(WTDIF *g){
	uint32_t i;
	free_string(g->name);
	free_u8list(g->rdseqs);
	free_u64list(g->seqoffs);
	free_u32list(g->seqlens);
	free_u8list(g->namestr);
	free_u32list(g->seqtags);
	free_u32list(g->rdranks);
	free_bitvec(g->confids);
	free_bitvec(g->varmaps);
	free_kswxv(g->swxs);
	free_pwmatv(g->mats);
	if(g->cigs){ for(i=0;i<g->cigs->size;i++) free_u32list((u32list*)g->cigs->buffer[i]); free_vplist(g->cigs); g->cigs = NULL; }
	free_u8list(g->hmsa);
	free_u8list(g->vmsa);
	free(g);
}

void reset_wtdif(WTDIF *g){
	uint32_t i;
	//clear_string(g->name);
	g->n_rd = 0;
	clear_u8list(g->rdseqs);
	clear_u64list(g->seqoffs);
	clear_u32list(g->seqlens);
	clear_u8list(g->namestr);
	clear_u32list(g->seqtags);
	clear_u32list(g->rdranks);
	clear_bitvec(g->confids);
	clear_bitvec(g->varmaps);
	clear_kswxv(g->swxs);
	clear_pwmatv(g->mats);
	for(i=0;i<g->cigs->size;i++) free_u32list((u32list*)g->cigs->buffer[i]);
	clear_vplist(g->cigs);
	clear_u8list(g->hmsa);
	clear_u8list(g->vmsa);
}

// The first sequence is the reference, x = KSWX_NULL, cigar like '10000M'
void push_wtdif(WTDIF *g, char *tag, int tag_len, char *seq, int len, kswx_t x, u32list *cigar){
	int i;
	u32list *cig;
	push_u32list(g->seqlens, len);
	push_u64list(g->seqoffs, g->rdseqs->size);
	encap_u8list(g->rdseqs, len);
	for(i=0;i<len;i++) lazy_push_u8list(g->rdseqs, base_bit_table[(int)seq[i]]);
	push_u32list(g->seqtags, g->namestr->size);
	append_array_u8list(g->namestr, (uint8_t*)tag, tag_len);
	push_u8list(g->namestr, 0);
	push_kswxv(g->swxs, x);
	if(cigar){
		cig = init_u32list(cigar->size); append_u32list(cig, cigar);
	} else cig = init_u32list(4);
	push_vplist(g->cigs, cig);
	g->n_rd ++;
}

void polish_alignment_wtdif(WTDIF *g, uint32_t rid, u32list *cigar_cache, u8list *mem_cache){
	kswx_t x, y;
	u32list *cig;
	uint8_t *q, *t;
	if(rid == 0) return;
	x = get_kswxv(g->swxs, rid);
	cig = (u32list*)get_vplist(g->cigs, rid);
	clear_u32list(cigar_cache); append_u32list(cigar_cache, cig);
	q = g->rdseqs->buffer + g->seqoffs->buffer[rid];
	t = g->rdseqs->buffer + g->seqoffs->buffer[0];
	y = kswx_refine_alignment(q, x.qb, t, x.tb, g->W, g->M, g->X, g->O, g->O, g->E, cigar_cache, mem_cache, cig);
	set_kswxv(g->swxs, rid, y);
}

void scan_long_matches_wtdif(WTDIF *g, uint32_t min_len, uint32_t min_cov){
	u32list *cig, *regs;
	kswx_t *s;
	uint32_t i, j;
	uint32_t x, c, op, sz, margin;
	margin = 2;
	clear_pwmatv(g->mats);
	regs = init_u32list(1024);
	for(i=1;i<g->n_rd;i++){
		s = ref_kswxv(g->swxs, i);
		cig = (u32list*)get_vplist(g->cigs, i);
		x = s->tb;
		for(j=0;j<cig->size;j++){
			op = cig->buffer[j] & 0x0F;
			sz = cig->buffer[j] >> 4;
			switch(op){
				case 1: break;
				case 2: x += sz; break;
				default:
				if(sz > margin * 2){
					push_u32list(regs, x + margin);
					push_u32list(regs, 0x80000000U | (x + margin + sz - margin));
				}
				x += sz;
			}
		}
	}
	sort_array(regs->buffer, regs->size, uint32_t, (a&0x7FFFFFFFU) > (b&0x7FFFFFFFU));
	if(min_cov < 1) min_cov = 1;
	c = 0;
	x = 0;
	for(i=0;i<regs->size;i++){
		if(regs->buffer[i] >> 31){
			if(c == min_cov && x + min_len <= (regs->buffer[i] & 0x7FFFFFFFU)) push_pwmatv(g->mats, (pw_mat_t){x, regs->buffer[i] & 0x7FFFFFFFU});
			c --;
		} else {
			c ++;
			if(c == min_cov) x = regs->buffer[i] & 0x7FFFFFFFU;
		}
		//fprintf(stdout, "%c\t%d\t%d\n", "+-"[regs->buffer[i]>> 31], regs->buffer[i] & 0x7FFFFFFFU, c);
	}
	if(c >= min_cov && x + min_len <= g->seqlens->buffer[0]) push_pwmatv(g->mats, (pw_mat_t){x, regs->buffer[regs->size-1] & 0x7FFFFFFFU});
	free_u32list(regs);
}

void gen_hard_alignments_wtdif(WTDIF *g){
	u32list *cig;
	kswx_t x;
	uint8_t *seq, *q, *t;
	uint32_t i, j, a, x1, x2, op, sz;
	clear_and_encap_u8list(g->hmsa, g->hlen * g->n_rd);
	for(i=0;i<g->n_rd;i++){
		x = get_kswxv(g->swxs, i);
		cig = (u32list*)get_vplist(g->cigs, i);
		seq = g->rdseqs->buffer + g->seqoffs->buffer[i];
		x1 = x.qb;
		for(x2=0;(int)x2<x.tb;x2++) if(get_bitvec(g->confids, x2)) push_u8list(g->hmsa, '^');
		for(j=0;j<cig->size;j++){
			op = cig->buffer[j] & 0x0F;
			sz = cig->buffer[j] >> 4;
			switch(op){
				case 1: x1 += sz; break;
				case 2:
				for(a=0;a<sz;a++) if(get_bitvec(g->confids, x2 + a)) push_u8list(g->hmsa, '-');
				x2 += sz; break;
				default:
				for(a=0;a<sz;a++) if(get_bitvec(g->confids, x2 + a)) push_u8list(g->hmsa, bit_base_table[seq[x1 + a]]);
				x1 += sz; x2 += sz;
			}
		}
		for(;x2<g->seqlens->buffer[0];x2++) if(get_bitvec(g->confids, x2)) push_u8list(g->hmsa, '$');
		push_u8list(g->hmsa, '\0');
		t = g->hmsa->buffer;
		q = g->hmsa->buffer + i * g->hlen;
		while(*q){
			if(*q == *t){
			} else if(*q < 'A' || *q >= 'Z'){
			} else {
				*q = (*q) + 'a' - 'A';
			}
			q ++; t ++;
		}
	}
}

void gen_var_alignments_wtdif(WTDIF *g, uint32_t min_var, float min_freq){
	uint8_t *q, *t, c;
	uint32_t i, j, cnts[5], max, sum;
	clear_bitvec(g->varmaps);
	clear_u8list(g->vmsa);
	g->vlen = 1;
	for(i=0;i<g->hlen-1;i++){
		cnts[0] = cnts[1] = cnts[2] = cnts[3] = cnts[4] = 0;
		for(j=0;j<g->n_rd;j++){
			c = g->hmsa->buffer[j * g->hlen + i];
			if(c == '^' || c == '$') continue;
			cnts[base_bit_table[c]] ++;
		}
		max = 0; sum = cnts[0];
		for(j=1;j<4;j++){
			if(cnts[j] > cnts[max]) max = j;
			sum += cnts[j];
		}
		if(cnts[4] < sum && sum - cnts[max] >= min_var && sum - cnts[max] >= (uint32_t)(sum * min_freq)){
			one2bitvec(g->varmaps);
			g->vlen ++;
		} else {
			zero2bitvec(g->varmaps);
		}
	}
	t = g->hmsa->buffer;
	for(i=0;i<g->n_rd;i++){
		q = g->hmsa->buffer + i * g->hlen;
		for(j=0;j<g->hlen-1;j++){
			if(get_bitvec(g->varmaps, j)){
				if(q[j] == t[j]) push_u8list(g->vmsa, '.');
				else if(q[j] == '^' || q[j] == '$') push_u8list(g->vmsa, '-');
				else push_u8list(g->vmsa, q[j]);
			}
		}
		push_u8list(g->vmsa, '\0');
	}
}

void run_wtdif(WTDIF *g, uint32_t min_len, uint32_t min_cov, uint32_t min_var, float min_freq, FILE *out){
	u32list *cigar_cache;
	u8list *mem_cache;
	pw_mat_t *m;
	kswx_t x;
	uint32_t i, j, idx;
	cigar_cache = init_u32list(1024);
	mem_cache = init_u8list(1024);
	for(i=0;i<g->n_rd;i++) push_u32list(g->rdranks, i);
	sort_array(g->rdranks->buffer + 1, g->n_rd - 1, uint32_t, g->swxs->buffer[a].tb > g->swxs->buffer[b].tb);
	for(i=1;i<g->n_rd;i++){
		idx = g->rdranks->buffer[i];
		polish_alignment_wtdif(g, idx, cigar_cache, mem_cache);
		x = get_kswxv(g->swxs, idx);
		fprintf(out, "%s\t%d\t%d\t%d"  , g->namestr->buffer + g->seqtags->buffer[idx], g->seqlens->buffer[idx], x.qb, x.qe);
		fprintf(out, "\t%s\t%d\t%d\t%d", g->namestr->buffer + g->seqtags->buffer[0], g->seqlens->buffer[0], x.tb, x.te);
		fprintf(out, "\t%d\t%0.2f\t%d\t%d\t%d\t%d\n", x.score, 1.0 * x.mat / x.aln, x.mat, x.mis, x.ins, x.del);
	}
	free_u32list(cigar_cache);
	free_u8list(mem_cache);
	scan_long_matches_wtdif(g, min_len, min_cov);
	g->hlen = 1;
	fprintf(out, "RULER\t0\t");
	for(i=0;i<g->mats->size;i++){
		m = ref_pwmatv(g->mats, i);
		g->hlen += m->y - m->x;
		fputc('|', out);
		for(j=m->x+1;(int)j<m->y;j++) fputc('-', out);
		//fprintf(out, "%d\t%d\t%d\n", m->x, m->y, m->y - m->x);
	}
	fputc('\n', out);
	encap_bitvec(g->confids, g->hlen);
	zeros_bitvec(g->confids);
	for(i=0;i<g->mats->size;i++){
		m = ref_pwmatv(g->mats, i);
		for(j=m->x;(int)j<m->y;j++) one_bitvec(g->confids, j);
	}
	gen_hard_alignments_wtdif(g);
	for(i=0;i<g->n_rd;i++){
		idx = g->rdranks->buffer[i];
		x = get_kswxv(g->swxs, idx);
		fprintf(out, "%s\t%d\t%s\n"  , g->namestr->buffer + g->seqtags->buffer[idx], x.tb, g->hmsa->buffer + idx * g->hlen);
	}
	gen_var_alignments_wtdif(g, min_var, min_freq);
	for(i=0;i<g->n_rd;i++){
		idx = g->rdranks->buffer[i];
		fprintf(out, "VAR\t%s\t%s\n"  , g->namestr->buffer + g->seqtags->buffer[idx], g->vmsa->buffer + idx * g->vlen);
	}
}

int usage(){
	printf(
	"WTDIF: Diff reference pacbio read and other other pacbio reads\n"
	"WatchTower: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtdif [options] [wt.fid]\n"
	"Options:\n"
	" -i <string> Input file, +\n"
	" -o <string> Output file, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: insertion or deletion, used in first round [-2]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -W <int>    Basic bandwidth in refine-alignment, [8]\n"
	" -l <int>    Min length of non-indel block, [5]\n"
	" -c <int>    Min coverage of non-indel block, [4]\n"
	" -V <float>  Variants calling, -V 2.05 mean: min_allele_count>=2,min_allele_freq>=0.05, [2.2]\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	WTDIF *g;
	FileReader *fr;
	cplist *infs;
	char *outf;
	FILE *out;
	kswx_t x;
	u32list *cig;
	int c, M, X, O, E, W, overwrite, min_len, min_cov, min_var;
	float min_freq;
	min_len = 5;
	min_cov = 4;
	W = 8;
	M = 2;
	X = -5;
	O = -2;
	E = -1;
	min_var = 2;
	min_freq = 0.2;
	overwrite = 0;
	infs = init_cplist(5);
	outf = NULL;
	while((c = getopt(argc, argv, "hi:o:fM:X:O:E:W:l:c:V:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'W': W  = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'c': min_cov = atoi(optarg); break;
			case 'V': min_freq = atof(optarg); min_var = min_freq; min_freq -= min_var; break;
			default: return usage();
		}
	}
	if(outf && !overwrite && strcmp(outf, "-") && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	for(c=optind;c<argc;c++) push_cplist(infs, argv[c]);
	if(infs->size){
		fr = fopen_m_filereader(infs->size, infs->buffer);
	} else fr = stdin_filereader();
	if(outf) out = fopen(outf, "w");
	else out = stdout;
	g = init_wtdif(W, M, X, O, E);
	x = KSWX_NULL;
	cig = init_u32list(1024);
	while(1){
		c = fread_table(fr);
		if(c == -1 || fr->line->string[0] == '>'){
			if(g->n_rd){
				run_wtdif(g, min_len, min_cov, min_var, min_freq, out);
				g->n_rd = 0;
			}
			if(c == -1) break;
			clear_string(g->name);
			append_string(g->name, fr->line->string + 1, fr->line->size - 1);
			trim_string(g->name);
			fprintf(out, ">%s\n", g->name->string);
			continue;
		}
		if(fr->line->string[0] == '#') continue;
		if(c < 5) continue;
		// rname, roff, toff, seq, cigar
		x.qb = atoi(get_col_str(fr, 1));
		x.tb = atoi(get_col_str(fr, 2));
		clear_u32list(cig); kswx_string2cigar(cig, get_col_str(fr, 4));
		push_wtdif(g, get_col_str(fr, 0), get_col_len(fr, 0), get_col_str(fr, 3), get_col_len(fr, 3), x, cig);
	}
	free_u32list(cig);
	fclose_filereader(fr);
	free_cplist(infs);
	if(outf) fclose(out);
	free_wtdif(g);
	return 0;
}

