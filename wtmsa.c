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

#include "pomsa.h"
#include "bit2vec.h"
#include "kswx.h"
#include "hzm_aln.h"
#include "file_reader.h"
#include "timer.h"
#include "thread.h"

static int msa_debug = 0;

define_list(kswxv, kswx_t);
define_list(strpv, String*);

typedef struct {
	String *name;
	uint32_t n_rd;
	int offset_delta;
	u8list  *rdseqs;
	u8list  *rdqvs;
	u64list *seqoffs;
	u32list *seqlens;
	u64list *rdqoffs;
	u32list *rdqlens;
	u8list  *namestr;
	u32list *seqtags;
	Bit2Vec *flags;
	kswxv   *kswxs;
	u4v     *begs;
	float min_sm;
	u1v     *ref;
	POMSA   *msa;
} WTMSA;

WTMSA* init_wtmsa(int W, int gw, int M, int X, int O, int I, int D, int E, int T, float min_sm){
	WTMSA *g;
	g = malloc(sizeof(WTMSA));
	g->name = init_string(1024);
	g->n_rd = 0;
	g->offset_delta = 0;
	g->rdseqs  = init_u8list(1024);
	g->rdqvs   = init_u8list(1024);
	g->seqoffs = init_u64list(1024);
	g->seqlens = init_u32list(1024);
	g->rdqoffs = init_u64list(1024);
	g->rdqlens = init_u32list(1024);
	g->namestr = init_u8list(1024);
	g->seqtags = init_u32list(1024);
	g->flags   = init_bit2vec(1024);
	g->kswxs   = init_kswxv(1024);
	g->begs    = init_u4v(1024);
	g->msa     = init_pomsa();
	g->ref     = init_u1v(1024);
	g->msa->W  = gw;
	g->msa->aux->W = W;
	g->msa->aux->M = M;
	g->msa->aux->X = X;
	g->msa->aux->O = O;
	g->msa->aux->I = I;
	g->msa->aux->D = D;
	g->msa->aux->E = E;
	g->msa->aux->T = T;
	g->min_sm = min_sm;
	return g;
}

void free_wtmsa(WTMSA *g){
	free_string(g->name);
	free_u8list(g->rdseqs);
	free_u8list(g->rdqvs);
	free_u64list(g->seqoffs);
	free_u32list(g->seqlens);
	free_u64list(g->rdqoffs);
	free_u32list(g->rdqlens);
	free_u8list(g->namestr);
	free_u32list(g->seqtags);
	free_bit2vec(g->flags);
	free_kswxv(g->kswxs);
	free_u4v(g->begs);
	free_u1v(g->ref);
	free_pomsa(g->msa);
	free(g);
}

void reset_wtmsa(WTMSA *g){
	clear_string(g->name);
	g->n_rd = 0;
	clear_u8list(g->rdseqs);
	clear_u8list(g->rdqvs);
	clear_u64list(g->seqoffs);
	clear_u32list(g->seqlens);
	clear_u64list(g->rdqoffs);
	clear_u32list(g->rdqlens);
	clear_u8list(g->namestr);
	clear_u32list(g->seqtags);
	clear_bit2vec(g->flags);
	clear_kswxv(g->kswxs);
	clear_u4v(g->begs);
	clear_u1v(g->ref);
}

void push_wtmsa(WTMSA *g, int flag, char *tag, int tag_len, char *seq, int len, int offset){
	int i;
	push_bit2vec(g->flags, flag);
	push_u32list(g->seqlens, len);
	push_u64list(g->seqoffs, g->rdseqs->size);
	push_u32list(g->rdqlens, 0);
	push_u64list(g->rdqoffs, g->rdqvs->size);
	encap_u8list(g->rdseqs, len);
	for(i=0;i<len;i++) lazy_push_u8list(g->rdseqs, base_bit_table[(int)seq[i]]);
	push_kswxv(g->kswxs, (kswx_t){0, offset, offset + len, 0, len, 0, 0, 0, 0, 0});
	push_u32list(g->seqtags, g->namestr->size);
	append_array_u8list(g->namestr, (uint8_t*)tag, tag_len);
	push_u8list(g->namestr, 0);
	g->n_rd ++;
}

void push5q_wtmsa(WTMSA *g, int flag, char *tag, int tag_len, char *seq, int len, int offset, char *qvs){
	int i;
	push_bit2vec(g->flags, flag);
	push_u32list(g->seqlens, len);
	push_u64list(g->seqoffs, g->rdseqs->size);
	push_u32list(g->rdqlens, len * 7);
	push_u64list(g->rdqoffs, g->rdqvs->size);
	encap_u8list(g->rdseqs, len);
	for(i=0;i<len;i++) lazy_push_u8list(g->rdseqs, base_bit_table[(int)seq[i]]);
	encap_u8list(g->rdqvs, len * 7);
	for(i=0;i<5*len;i++) lazy_push_u8list(g->rdqvs, qvs[i] - 33);
	for(;i<7*len;i++) lazy_push_u8list(g->rdqvs, base_bit_table[(int)qvs[i]]);
	push_kswxv(g->kswxs, (kswx_t){0, offset, offset + len, 0, len, 0, 0, 0, 0, 0});
	push_u32list(g->seqtags, g->namestr->size);
	append_array_u8list(g->namestr, (uint8_t*)tag, tag_len);
	push_u8list(g->namestr, 0);
	g->n_rd ++;
}

int gen_backbone_wtmsa(WTMSA *g){
	kswx_t *p, *s;
	uint32_t i;
	uint8_t *t, *q;
	t = NULL;
	uint8_t *wushigang = t;
	uint8_t *tmp = wushigang;
	wushigang = tmp;
	int tlen, qlen, tb, lstx, lsty;
	clear_u1v(g->ref);
	lstx = lsty = 0;
	for(i=0;i<g->n_rd;i++){
		p = ref_kswxv(g->kswxs, i);
		if(get_bit2vec(g->flags, i) == 0){
			t = g->ref->buffer;
			tlen = g->ref->size;
			q = g->rdseqs->buffer + g->seqoffs->buffer[i];
			qlen = g->seqlens->buffer[i];
			if(tlen == 0){
				append_array_u1v(g->ref, q, qlen);
				g->offset_delta = - p->tb;
				p->tb = 0;
				p->te = qlen;
				p->qb = 0;
				p->qe = qlen;
				p->aln = -1;
			} else if(p->tb >= 0){
				tb = p->tb + g->offset_delta;
				p->tb = (tb - lstx) + lsty;
				p->te = p->tb + qlen;
				reset_hzmaux(g->msa->aux);
				app_tseq_hzmaux(g->msa->aux, g->ref->buffer, g->ref->size);
				ready_hzmaux(g->msa->aux);
				if(align_hzmaux(g->msa->aux, 0, q, NULL, qlen, p->tb, p->te, 0, g->min_sm)){
					s = &g->msa->aux->hit;
					if(msa_debug >= 1) fprintf(stderr, "BACKBONE\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", g->namestr->buffer + g->seqtags->buffer[i],  s->qb, s->qe, s->tb, s->te, s->score, s->mat, s->mis, s->ins, s->del);
					if(tlen - s->te >= qlen - s->qe) continue; // maybe inside or partial
					g->ref->size = s->te;
					append_array_u1v(g->ref, q + s->qe, qlen - s->qe);
					p->tb = s->tb;
					p->te = g->ref->size;
					p->qb = s->qb;
					p->qe = qlen;
				} else {
					if(msa_debug >= 1) fprintf(stderr, "BACKBONE\t%s\tfailed\n", g->namestr->buffer + g->seqtags->buffer[i]);
					p->tb = g->ref->size;
					append_array_u1v(g->ref, q, qlen);
					p->te = g->ref->size;
					p->qb = 0;
					p->qe = qlen;
				}
				lstx = tb;
				lsty = p->tb;
				p->aln = -1;
			} else {
				// mark it as non-backbone
				set_bit2vec(g->flags, i, 1);
				p->tb  = -1;
				p->te  = -1;
				p->aln = -1;
			}
		}
	}
	return 1;
}

thread_beg_def(mhzm);
WTMSA *g;
HZMAux *aux;
uint32_t qid, tid;
uint8_t *q, *t;
int qlen, tlen;
int tb, te;
int refine;
float min_sm;
int ret;
thread_end_def(mhzm);

thread_beg_func(mhzm);
thread_beg_loop(mhzm);
if(mhzm->tlen){
	reset_hzmaux(mhzm->aux);
	app_tseq_hzmaux(mhzm->aux, mhzm->t, mhzm->tlen);
	ready_hzmaux(mhzm->aux);
	mhzm->tlen = 0;
}
if(mhzm->te < 0) mhzm->te = mhzm->aux->tseq->size;
//if(strcmp("pb000000003230", (char*)mhzm->g->namestr->buffer + mhzm->g->seqtags->buffer[mhzm->qid]) == 0){
	//if(strcmp("pb000000005637", (char*)mhzm->g->namestr->buffer + mhzm->g->seqtags->buffer[mhzm->tid]) == 0){
		//fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
	//}
//}
mhzm->ret = align_hzmaux(mhzm->aux, 0, mhzm->q, NULL, mhzm->qlen, mhzm->tb, mhzm->te, mhzm->refine, mhzm->min_sm);
thread_end_loop(mhzm);
thread_end_func(mhzm);

int mp_gen_backbone_wtmsa(WTMSA *g, int ncpu){
	kswx_t *p, *s, *k;
	uint32_t i, j;
	uint8_t *t, *q;
	t = NULL;
	uint8_t *wushigang = t;
	uint8_t *tmp = wushigang;
	wushigang = tmp;
	int tlen, qlen, tb;
	thread_preprocess(mhzm);
	thread_beg_init(mhzm, ncpu);
	mhzm->g = g;
	mhzm->aux = init_hzmaux();
	copy_param_hzmaux(mhzm->aux, g->msa->aux);
	mhzm->q = NULL;
	mhzm->t = NULL;
	mhzm->qlen = 0;
	mhzm->tlen = 0;
	mhzm->refine = 0;
	mhzm->min_sm = g->min_sm;
	mhzm->tb = 0;
	mhzm->te = -1;
	mhzm->ret = 0;
	mhzm->qid = 0xFFFFFFFFU;
	mhzm->tid = 0xFFFFFFFFU;
	thread_end_init(mhzm);
	s = NULL;
	for(i=j=0;i<g->n_rd;i++){
		p = ref_kswxv(g->kswxs, i);
		p->score = 0;
		if(get_bit2vec(g->flags, i) != 0) continue;
		if(s){
			tb = p->tb - s->tb;
			thread_wait_one(mhzm);
			if(mhzm->ret){
				k = ref_kswxv(g->kswxs, mhzm->qid);
				k->score = 1;
				k->aln = mhzm->tid;
				k->mat = mhzm->aux->hit.qb;
				k->mis = mhzm->aux->hit.qe;
				k->ins = mhzm->aux->hit.tb;
				k->del = mhzm->aux->hit.te;
				if(msa_debug >= 1){
					s = &mhzm->aux->hit;
					fprintf(stderr, "BACKBONE\t%s\t%d\t%d\t%d", g->namestr->buffer + g->seqtags->buffer[mhzm->qid], g->seqlens->buffer[mhzm->qid], s->qb, s->qe);
					fprintf(stderr, "\t%s\t%d\t%d\t%d", g->namestr->buffer + g->seqtags->buffer[mhzm->tid], g->seqlens->buffer[mhzm->tid], s->tb, s->te);
					fprintf(stderr, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", s->score, 1.0 * s->mat / s->aln, s->mat, s->mis, s->ins, s->del);
				}
			} else if(mhzm->qlen){
				k = ref_kswxv(g->kswxs, mhzm->qid);
				k->score = -1;
				if(msa_debug >= 1){
					fprintf(stderr, "BACKBONE\t%s\t%d\t%s\t%d\tfailed\n",
					g->namestr->buffer + g->seqtags->buffer[mhzm->qid], g->seqlens->buffer[mhzm->qid],
					g->namestr->buffer + g->seqtags->buffer[mhzm->tid], g->seqlens->buffer[mhzm->tid]);
				}
			}
			mhzm->tid = j;
			mhzm->qid = i;
			mhzm->tb  = 0;
			mhzm->te  = -1;
			mhzm->t = g->rdseqs->buffer + g->seqoffs->buffer[j];
			mhzm->q = g->rdseqs->buffer + g->seqoffs->buffer[i];
			mhzm->tlen = g->seqlens->buffer[j];
			mhzm->qlen = g->seqlens->buffer[i];
			thread_wake(mhzm);
		}
		s = p;
		j = i;
	}
	thread_wait_all(mhzm);
	thread_beg_close(mhzm);
	if(mhzm->ret){
		k = ref_kswxv(g->kswxs, mhzm->qid);
		k->score = 1;
		k->aln = mhzm->tid;
		k->mat = mhzm->aux->hit.qb;
		k->mis = mhzm->aux->hit.qe;
		k->ins = mhzm->aux->hit.tb;
		k->del = mhzm->aux->hit.te;
		if(msa_debug >= 1){
			s = &mhzm->aux->hit;
			fprintf(stderr, "BACKBONE\t%s\t%d\t%d\t%d", g->namestr->buffer + g->seqtags->buffer[mhzm->qid], g->seqlens->buffer[mhzm->qid], s->qb, s->qe);
			fprintf(stderr, "\t%s\t%d\t%d\t%d", g->namestr->buffer + g->seqtags->buffer[mhzm->tid], g->seqlens->buffer[mhzm->tid], s->tb, s->te);
			fprintf(stderr, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\n", s->score, 1.0 * s->mat / s->aln, s->mat, s->mis, s->ins, s->del);
		}
	} else if(mhzm->qlen){
		k = ref_kswxv(g->kswxs, mhzm->qid);
		k->score = -1;
		if(msa_debug >= 1){
			fprintf(stderr, "BACKBONE\t%s\t%d\tfailed\n", g->namestr->buffer + g->seqtags->buffer[mhzm->qid], g->seqlens->buffer[mhzm->qid]);
		}
	}
	free_hzmaux(mhzm->aux);
	thread_end_close(mhzm);
	clear_u1v(g->ref);
	s = NULL;
	for(i=0;i<g->n_rd;i++){
		p = ref_kswxv(g->kswxs, i);
		if(get_bit2vec(g->flags, i)) continue;
		t = g->ref->buffer;
		tlen = g->ref->size;
		q = g->rdseqs->buffer + g->seqoffs->buffer[i];
		qlen = g->seqlens->buffer[i];
		if(tlen == 0){
			append_array_u1v(g->ref, q, qlen);
			p->tb = 0;
			p->te = qlen;
			p->qb = 0;
			p->qe = qlen;
			//p->aln = -1;
		} else if(p->score == 1){
			p->qb = p->mat;
			p->qe = qlen;
			p->tb = (s->tb - s->mat) + p->ins; // (s->tb - s->mat) is global offset of s, p->ins is the tb of p to s
			p->te = (s->tb - s->mat) + p->del + qlen - p->mis;
			//p->aln = -1;
			g->ref->size = tlen - (g->seqlens->buffer[p->aln] - p->del); // p->aln record the idx of previous read
			append_array_u1v(g->ref, q + p->mis, qlen - p->mis);
		} else {
			tb = tlen - qlen * 1.5;
			if(tb < 0) tb = 0;
			//call sw-align
			{
				int8_t matrix[16];
				kswr_t r;
				for(j=0;j<16;j++) matrix[j] = ((j / 4) == (j % 4))? g->msa->aux->M : g->msa->aux->X;
				r = kswx_align_no_stat(qlen, q, tlen - tb, g->ref->buffer + tb, 4, matrix, g->msa->aux->ew, g->msa->aux->I, g->msa->aux->D, g->msa->aux->E, g->msa->aux->T);
				if(r.qb < r.qe){
					p->qb = r.qb;
					p->qe = qlen;
					p->tb = tb + r.tb;
					p->te = tb + r.te + qlen - r.qe;
					p->score = r.score;
					g->ref->size = tb + r.te;
					append_array_u1v(g->ref, q + r.qe, qlen - r.qe);
					if(msa_debug >= 1){
						fprintf(stderr, "WARNNING\t%s\t%d\tsuccess to rescue\n", g->namestr->buffer + g->seqtags->buffer[i], g->seqlens->buffer[i]);
					}
				} else {
					// just append it, unhappy
					p->qb = 0;
					p->qe = qlen;
					p->tb = tlen;
					p->te = tlen + qlen;
					p->mat = 0;
					//p->aln = -1;
					append_array_u1v(g->ref, q, qlen);
					if(msa_debug >= 1){
						fprintf(stderr, "WARNNING\t%s\t%d\tfail to find overlap, just append it to backbone\n", g->namestr->buffer + g->seqtags->buffer[i], g->seqlens->buffer[i]);
					}
				}
			}
		}
		s = p;
	}
	return 1;
}

void run_wtmsa(WTMSA *g, int max_iter, int ncpu, FILE *backbone, FILE *dot){
	kswx_t *s, *p;
	uint32_t i, rdlen;
	uint8_t *rdseq;
	int iter, ret;
	thread_preprocess(mhzm);
	if(msa_debug >= 0) fprintf(stderr, "[%s]\"%s\"\n", date(), g->name->string);
	mp_gen_backbone_wtmsa(g, ncpu);
	if(msa_debug >= 0) fprintf(stderr, "[%s] generated backbone length=%d\n", date(), (int)g->ref->size);
	if(backbone){
		fprintf(backbone, ">backbone_%s\n", g->name->string);
		for(i=0;i<g->ref->size;i++){
			fputc(bit_base_table[(int)g->ref->buffer[i]], backbone);
			if((i % 100) == 99) fputc('\n', backbone);
		}
		if((i % 100) < 99) fputc('\n', backbone);
	}
	for(iter=0;iter<max_iter;iter++){
		beg_update_pomsa(g->msa, g->ref->buffer, g->ref->size);
		//if(iter == 0){
		if(1){
			thread_beg_init(mhzm, ncpu);
			mhzm->g = g;
			mhzm->aux = init_hzmaux();
			copy_param_hzmaux(mhzm->aux, g->msa->aux);
			share_index_hzmaux(mhzm->aux, g->msa->aux);
			mhzm->q = NULL;
			mhzm->t = NULL;
			mhzm->qlen = 0;
			mhzm->tlen = 0;
			mhzm->refine = 0;
			mhzm->min_sm = g->min_sm;
			mhzm->tb = 0;
			mhzm->te = -1;
			mhzm->ret = 0;
			mhzm->qid = 0xFFFFFFFFU;
			mhzm->tid = 0xFFFFFFFFU;
			thread_end_init(mhzm);
			for(i=0;i<g->n_rd;i++){
				//if(get_bit2vec(g->flags, i) == 3) continue;
				s = ref_kswxv(g->kswxs, i);
				s->score = s->mat = s->mis = s->ins = s->del = 0;
				thread_wait_one(mhzm);
				if(mhzm->ret){
					g->kswxs->buffer[mhzm->qid] = mhzm->aux->hit;
					if(msa_debug >= 1){
						p = g->kswxs->buffer + mhzm->qid;
						fprintf(stderr, "ALN\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", g->namestr->buffer + g->seqtags->buffer[mhzm->qid], g->seqlens->buffer[mhzm->qid],
							p->qb, p->qe, p->tb, p->te, p->score, p->mat, p->mis, p->ins, p->del);
					}
				} else if(mhzm->qlen){
					set_bit2vec(g->flags, mhzm->qid, 3);
					if(msa_debug >= 1){
						fprintf(stderr, "ALN\t%s\tfailed\n", g->namestr->buffer + g->seqtags->buffer[mhzm->qid]);
					}
				}
				mhzm->qid  = i;
				mhzm->q    = g->rdseqs->buffer + g->seqoffs->buffer[i];
				mhzm->qlen = g->seqlens->buffer[i];
				if(get_bit2vec(g->flags, i) == 0){
					mhzm->tb = s->tb;
					mhzm->te = s->te;
				} else {
					mhzm->tb = 0;
					mhzm->te = -1;
				}
				thread_wake(mhzm);
			}
			thread_wait_all(mhzm);
			thread_beg_close(mhzm);
			if(mhzm->ret){
				g->kswxs->buffer[mhzm->qid] = mhzm->aux->hit;
				if(msa_debug >= 1){
					p = g->kswxs->buffer + mhzm->qid;
					fprintf(stderr, "ALN\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", g->namestr->buffer + g->seqtags->buffer[mhzm->qid], g->seqlens->buffer[mhzm->qid],
						p->qb, p->qe, p->tb, p->te, p->score, p->mat, p->mis, p->ins, p->del);
				}
			} else if(mhzm->qlen){
				set_bit2vec(g->flags, mhzm->qid, 3);
				if(msa_debug >= 1){
					fprintf(stderr, "ALN\t%s\tfailed\n", g->namestr->buffer + g->seqtags->buffer[mhzm->qid]);
				}
			}
			free_hzmaux(mhzm->aux);
			thread_end_close(mhzm);
		}
		clear_u4v(g->begs);
		for(i=0;i<g->n_rd;i++){
			if(get_bit2vec(g->flags, i) == 3){ push_u4v(g->begs, 0); continue; }
			s = ref_kswxv(g->kswxs, i);
			if(msa_debug >= 1) fprintf(stderr, "[%s] alignment %s len=%d %u/%u\n", date(), g->namestr->buffer + g->seqtags->buffer[i], g->seqlens->buffer[i], i + 1, g->n_rd);
			if(msa_debug >= 1) fprintf(stderr, "PRE\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s->qb, s->qe, s->tb, s->te, s->score, s->mat, s->mis, s->ins, s->del);
			rdlen = g->seqlens->buffer[i];
			rdseq = g->rdseqs->buffer + g->seqoffs->buffer[i];
			ret = update_pomsa(g->msa, rdseq, rdlen, g->min_sm, s);
			push_u4v(g->begs, g->msa->beginning);
			if(ret == 0){ set_bit2vec(g->flags, i, 3); continue; }
			*s = g->msa->aux->hit;
			if(msa_debug >= 1) fprintf(stderr, "RST\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s->qb, s->qe, s->tb, s->te, s->score, s->mat, s->mis, s->ins, s->del);
			if(0){
				uint32_t j;
				for(j=0;j<=i;j++){
					if(get_bit2vec(g->flags, j) == 3) continue;
					if((ret = match_pomsa(g->msa, j, g->rdseqs->buffer + g->seqoffs->buffer[j] + g->kswxs->buffer[j].qb, g->kswxs->buffer[j].qe - g->kswxs->buffer[j].qb, g->begs->buffer[j], NULL)) == 0){
						fprintf(stderr, " -- i=%u j=%u in %s -- %s:%d --\n", i, j, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						exit(1);
					}
				}
			}
		}
		if(1){
			beg_match_pomsa(g->msa);
			for(i=0;i<g->n_rd;i++){
				if(get_bit2vec(g->flags, i) == 3){ continue; }
				s = ref_kswxv(g->kswxs, i);
				rdseq = g->rdseqs->buffer + g->seqoffs->buffer[i];
				ret = match_pomsa(g->msa, i, rdseq + s->qb, s->qe - s->qb, g->begs->buffer[i], NULL);
				if(msa_debug >= 1) fprintf(stderr, "[%s] match %s len=%d %u/%u ret=%d\n", date(), g->namestr->buffer + g->seqtags->buffer[i], g->seqlens->buffer[i], i + 1, g->n_rd, ret);
			}
		}
		call_consensus_pomsa(g->msa);
		if(msa_debug >= 0) fprintf(stderr, "[%s] iter%d length=%d aln_score=%d cns_score=%d\n", date(), iter + 1, (int)g->msa->cns->size, g->msa->aln_score, g->msa->cns_score);
		//printf_local_dot_pomsa(g->msa, 0, g->ref->size, "1.dot");
		clear_u1v(g->ref); append_u1v(g->ref, g->msa->cns);
		if(iter + 1 < max_iter){
			// update alignment coordinates
			for(i=0;i<g->n_rd;i++){
				if(get_bit2vec(g->flags, i) == 3) continue;
				s = ref_kswxv(g->kswxs, i);
				s->tb = s->tb? g->msa->map->buffer[s->tb - 1] : 0;
				s->te = s->te < (int)g->msa->ref->size? g->msa->map->buffer[s->te] : g->msa->cns->size;
			}
		}
	}
	if(msa_debug >= 0) fprintf(stderr, "[%s]\"%s\" Done\n\n", date(), g->name->string);
	if(dot){
		print_local_dot_pomsa(g->msa, NULL, 0, 0xFFFFFFFFU, dot);
	}
}

int usage(){
	printf(
	"WTMSA: Consensus caller using POA\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtmsa [options]\n"
	"Options:\n"
	//" -t <int>    Number of threads, [16]\n"
	" -P <int>    Total parallel jobs, [1]\n"
	" -p <int>    Index of current job (0-based), [0]\n"
	"             Suppose to run wtmsa for the same layout file parallelly in 60 cpu. For cpu1, -P 60 -p 0; cpu2, -P 60 -p 1, ...\n"
	" -i <string> Input file, layout from wtlay, +, *\n"
	" -o <string> Output file, consensus sequences, *\n"
	" -B <string> Print backbone sequences on file for debug [NULL]\n"
	" -G <stirng> Print dot graph on file, H U G E, be careful, [NULL]\n"
	" -f          Force overwrite\n"
	" -H          Trun off homopolymer compression\n"
	" -z <int>    Zmer size, 5 <= <-z> <= 16, [10]\n"
	" -y <int>    Zmer window, [800]\n"
	" -R <int>    Minimum size of seeding region within zmer window, [200]\n"
	" -l <int>    Maximum variant of uncompressed sizes between two matched zmer, [2]\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -I <int>    Alignment penalty: insertion, [-2]\n"
	" -D <int>    Alignment penalty: deletion, [-3]\n"
	" -V <int>    turn on homopolymer merge penalty\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping [-10]\n"
	" -F          Disable PhreadQV in refine-alignment\n"
	" -w <int>    Minimum bandwidth of pairwise alignment, iteratively doubled to maximum [50]\n"
	" -W <int>    Maximum bandwidth of pairwise alignment, [3200]\n"
	" -e <int>    Maximum bandwidth at graph alignment and ending extension, [800]\n"
	" -g <int>    Basic bandwidth in graph alignment, [100]\n"
	" -m <float>  Minimum alignment identity, [0.5]\n"
	//" -Y <float>  Penalty of backbone edge in calling consensus, [0.5]\n"
	//" -N <float>  Penalty of alternative edge in calling consensus, [0.2]\n"
	//"             The above two options control whether the consensus look like backbone or alternative\n"
	//"             Default 0.5 and 0.2, will let the consensus don't look like backbone\n"
	" -n <int>    Number of iterations for consensus calling, the more, the accurater, the slower [2]\n"
	//" -a <string> Align reads against final consensus, and output to <-a>\n"
	//" -A          Disable fast zmer align in final aligning (see -a), use standard smith-waterman\n"
	//"             More than once -A, will disable fast zmer align in all process\n"
	//" -V <float> Ouput call variants and print to <-a>, -V 2.05 mean: min_allele_count>=2,min_allele_freq>=0.05\n"
	" -v          Verbose, +\n"
	"\n"
	"Example: \n"
	"$> wtmsa -i wt.lay -o wt.lay.cns.fa 2>log.msa\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	obj_desc_t wsg = kswxv_obj_desc;
	wsg = strpv_obj_desc;
	obj_desc_t ttt = wsg;
	wsg = ttt;
	FileReader *fr;
	cplist *lays;
	char *outf, *print_ref, *print_dot;
	WTMSA *g;
	FILE *out, *refo, *doto;
	long job;
	int M, X, I, D, O, E, T, W, ew, rw, gw, w, zsize, hz, zwin, zovl, zvar, n, c, flag, n_job, i_job, ncpu, work, overwrite;
	int use_qv, has_merge;
	float min_sm;
	HZM_FAST_WINDOW_KMER_CHAINING = 0;
	ncpu = 1;
	zsize = 10;
	zvar = 2;
	hz = 1;
	zwin = 800;
	zovl = 200;
	use_qv = 1;
	w = 50;
	ew = 800;
	rw = 8; // basic bandwidth for refine-alignment
	gw = 100;
	W = 3200;
	M = 2;
	X = -5;
	O = -2;
	I = -2;
	D = -3;
	E = -1;
	T = -10;
	has_merge = 0;
	min_sm = 0.5;
	n = 2;
	n_job = 1;
	i_job = 0;
	outf = NULL;
	print_ref = NULL;
	print_dot = NULL;
	overwrite = 0;
	lays = init_cplist(4);
	while((c = getopt(argc, argv, "ht:P:p:i:o:B:G:fHz:y:R:l:M:X:O:I:D:VE:T:Fw:W:e:g:m:n:v")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'P': n_job = atoi(optarg); break;
			case 'p': i_job = atoi(optarg); break;
			case 'i': push_cplist(lays, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'B': print_ref = optarg; break;
			case 'G': print_dot = optarg; break;
			//case 'A': if(fast_align[1]) fast_align[1] = 0; else fast_align[0] = 0; break;
			case 'H': hz = 0; break;
			case 'z': zsize = atoi(optarg); break;
			case 'y': zwin  = atoi(optarg); break;
			case 'R': zovl  = atoi(optarg); break;
			case 'l': zvar  = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'V': has_merge = 1; break;
			case 'E': E = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			case 'F': use_qv = 0; break;
			case 'w': w  = atoi(optarg); break;
			case 'W': W  = atoi(optarg); break;
			case 'e': ew = atoi(optarg); break;
			//case 'r': rw = atoi(optarg); break;
			case 'g': gw = atoi(optarg); break;
			case 'm': min_sm = atof(optarg); break;
			//case 'Y': ref_penalty = atof(optarg); break;
			//case 'N': alt_penalty = atof(optarg); break;
			case 'n': n = atoi(optarg); break;
			//case 'a': alnf = optarg; break;
			//case 'V': vmsa = 1; min_freq = atof(optarg); min_cnt = min_freq; min_freq -= min_cnt; break;
			case 'v': msa_debug ++; break;
			default: return usage();
		}
	}
	if(msa_debug > 1) hzm_debug = msa_debug - 1;
	if(outf == NULL) return usage();
	for(c=optind;c<argc;c++) push_cplist(lays, argv[c]);
	if(lays->size){
		fr = fopen_m_filereader(lays->size, lays->buffer);
	} else return usage();
	out = open_file_for_write(outf, NULL, overwrite);
	refo = print_ref? open_file_for_write(print_ref, NULL, overwrite) : NULL;
	doto = print_dot? open_file_for_write(print_dot, NULL, overwrite) : NULL;
	g = init_wtmsa(W, gw, M, X, O, I, D, E, T, min_sm);
	g->msa->has_merge = has_merge;
	g->msa->MW = ew;
	g->msa->aux->hz    = hz;
	g->msa->aux->zsize = zsize;
	g->msa->aux->zwin  = zwin;
	g->msa->aux->zstep = zwin / 2;
	g->msa->aux->zovl = zovl;
	g->msa->aux->zvar = zvar;
	g->msa->aux->ew = ew;
	g->msa->aux->rw = rw;
	g->msa->aux->w = w;
	work = 0;
	job = 0;
	while(1){
		c = fread_table(fr);
		if(c == -1 || fr->line->string[0] == '>'){
			if(g->n_rd){
				run_wtmsa(g, n, ncpu, refo, doto);
				fprintf(out, ">%s\n", g->name->string);
				print_cnsseq_pomsa(g->msa, out);
				fflush(out);
				reset_wtmsa(g);
			}
			if((job % n_job) == i_job) work = 1;
			else work = 0;
			job ++;
			if(c == -1) break;
			if(work){
				append_string(g->name, fr->line->string + 1, fr->line->size - 1);
				trim_string(g->name);
			}
			continue;
		}
		if(work == 0) continue;
		if(fr->line->string[0] == '#') continue;
		if(c < 6) continue;
		switch(get_col_str(fr, 0)[0]){
			case 'Y': flag = 0; break;
			case 'N': flag = 1; break;
			case 'n': flag = 2; break;
			default: flag = 3;
		}
		if(use_qv && c > 6 && get_col_len(fr, 6) == 7 * get_col_len(fr, 5)){ // f5q
			push5q_wtmsa(g, flag, get_col_str(fr, 1), get_col_len(fr, 1), get_col_str(fr, 5), get_col_len(fr, 5), atoi(get_col_str(fr, 3)), get_col_str(fr, 6));
		} else {
			push_wtmsa(g, flag, get_col_str(fr, 1), get_col_len(fr, 1), get_col_str(fr, 5), get_col_len(fr, 5), atoi(get_col_str(fr, 3)));
		}
	}
	fclose_filereader(fr);
	free_cplist(lays);
	close_file_safely(out);
	close_file_safely(refo);
	close_file_safely(doto);
	free_wtmsa(g);
	return 0;
}
