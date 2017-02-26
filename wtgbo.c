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

#include "wtlay.h"
#include "hzm_aln.h"
#include "thread.h"

thread_beg_def(mgbo);
StringGraph *g;
HZMAux *aux;
float min_sm;
int refine;
uint32_t obj_id, qry_id;
int dir, ret, contained;
thread_end_def(mgbo);

thread_beg_func(mgbo);
u1v *rdseq;
uint64_t rdoff;
uint32_t i, rdlen;
rdseq = init_u1v(1024);
thread_beg_loop(mgbo);
if(mgbo->aux->tseq->size == 0){
	rdoff = mgbo->g->rdoffs->buffer[mgbo->obj_id];
	rdlen = mgbo->g->rdlens->buffer[mgbo->obj_id];
	for(i=0;i<rdlen;i++) add_tseq_hzmaux(mgbo->aux, get_basebank(mgbo->g->rdseqs, rdoff + i));
	ready_hzmaux(mgbo->aux);
}
rdoff = mgbo->g->rdoffs->buffer[mgbo->qry_id];
rdlen = mgbo->g->rdlens->buffer[mgbo->qry_id];
encap_u1v(rdseq, rdlen);
if(mgbo->dir) revbitseq_basebank(mgbo->g->rdseqs, rdoff, rdlen, rdseq->buffer);
else             bitseq_basebank(mgbo->g->rdseqs, rdoff, rdlen, rdseq->buffer);
mgbo->contained = 0;
if((mgbo->ret = align_hzmaux(mgbo->aux, mgbo->qry_id, rdseq->buffer, NULL, rdlen, 0, 0, mgbo->refine, mgbo->min_sm))){
	if(mgbo->aux->hit.tb == 0 && mgbo->aux->hit.te == (int)rdlen){
		mgbo->contained = 1;
	}
}
thread_end_loop(mgbo);
free_u1v(rdseq);
thread_end_func(mgbo);

#define WTGBO_GRAPH_TRACE_LEVEL	2

void collect_graph_candidates_wtgbo(StringGraph *g, uint32_t node_id, uint32_t max_ext, u64hash *closed_alns, u4v *candidates, u8v *heap, u32hash *hash){
	sg_node_t *n1;
	sg_edge_t *e;
	uint64_t idx, lval, *lptr;
	uint32_t nid, dir, k, i, j, f, val, off1, off2, max, lv, *ptr;
	int exists, att;
	max = max_ext + g->rdlens->buffer[node_id];
	clear_u8v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint64_t, (0x0LLU << 40) | (0x0LLU << 32) | (node_id << 2) | (0 << 1) | 0, num_cmp((a >> 40), (b >> 40)));
	array_heap_push(heap->buffer, heap->size, heap->cap, uint64_t, (0x0LLU << 40) | (0x0LLU << 32) | (node_id << 2) | (1 << 1) | 1, num_cmp((a >> 40), (b >> 40)));
	f = 0;
	while(heap->size){
		idx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint64_t, num_cmp((a >> 40), (b >> 40)));
		nid = (idx & 0xFFFFFFFFU) >> 2;
		dir = (idx >> 1) & 0x01;
		k = idx & 0x01;
		lv = (idx >> 32) & 0xFFU;
		off1 = idx >> 40;
		n1 = ref_sgnodev(g->nodes, nid);
		att = 0;
		for(i=0;i<=n1->edge_cnts[k];i++){
			if(i == n1->edge_cnts[k] && att == 0){
				// special operation on contained reads
				if(!is_dead_node_strgraph(g, nid)) break;
				for(j=0;j<n1->edge_cnts[!k];j++){
					e = ref_sgedgev(g->edges, n1->edge_offs[!k] + j);
					if(e->att) break;
				}
				if(j == n1->edge_cnts[!k]) break;
			} else {
				e = ref_sgedgev(g->edges, n1->edge_offs[k] + i);
			}
			if(e->att) att = 1;
			off2 = off1 + e->off;
			if(off2 > max) continue;
			val = (e->node_id << 1) | (dir ^ e->dir);
			ptr = prepare_u32hash(hash, val, &exists);
			if(exists) continue;
			*ptr = val;
			if(f && !is_dead_node_strgraph(g, e->node_id)){
				lval = ovl_uniq_long_id(node_id, e->node_id, (e->dir ^ dir));
				lptr = prepare_u64hash(closed_alns, lval, &exists);
				if(!exists){
					*lptr = lval;
					push_u4v(candidates, val);
				}
			}
			if(lv < WTGBO_GRAPH_TRACE_LEVEL){
				array_heap_push(heap->buffer, heap->size, heap->cap, uint64_t,
					(((uint64_t)off2) << 40) | (((uint64_t)lv + 1) << 32) | (e->node_id << 2) | (dir << 1) | e->dir, num_cmp((a >> 40), (b >> 40)));
			}
		}
		f = 1;
	}
}

int output_wtgbo(StringGraph *g, uint32_t id1, uint32_t id2, int dir, HZMAux *aux, sgbiedgev *biedges, FILE *out){
	OverlapData O;
	sg_biedge_t B, *b;
	fprintf(out,   "%s\t%c\t%d\t%d\t%d", g->rdnames->buffer[id1], '+',       g->rdlens->buffer[id1], aux->hit.tb, aux->hit.te);
	fprintf(out, "\t%s\t%c\t%d\t%d\t%d", g->rdnames->buffer[id2], "+-"[dir], g->rdlens->buffer[id2], aux->hit.qb, aux->hit.qe);
	fprintf(out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\t", aux->hit.score, 1.0 * aux->hit.mat / aux->hit.aln, aux->hit.mat, aux->hit.mis, aux->hit.ins, aux->hit.del);
	kswx_print_cigars(aux->cigars->buffer, aux->cigars->size, out);
	fprintf(out, "\n");
	b = &B;
	O.node_id[0] = id1;
	O.node_id[1] = id1;
	O.dir[0] = 0;
	O.dir[1] = dir;
	O.beg[0] = aux->hit.tb;
	O.beg[1] = aux->hit.qb;
	O.end[0] = aux->hit.te;
	O.end[1] = aux->hit.qe;
	O.score  = aux->hit.score;
	O.identity = 1000.0 * aux->hit.mat / aux->hit.aln;
	if(overlap_item2biedge_v2_strgraph(g, &O, b)){
		push_sgbiedgev(biedges, B);
		return 1;
	} else return 0;
}

uint64_t gbo_core_wtgbo(StringGraph *g, HZMAux *aux, int max_ext, float min_sm, int refine, int ncpu, u64hash *closed_alns, sgbiedgev *biedges, FILE *out){
	sg_node_t *n;
	u4v *candidates;
	u8v *heap;
	u32hash *hash;
	uint64_t ret, ncand, rdoff;
	uint32_t node_id, rdlen, i;
	def_counter(cnt);
	thread_preprocess(mgbo);
	candidates = init_u4v(1024);
	hash = init_u32hash(1023);
	heap = init_u8v(1024);
	ret = 0;
	thread_beg_init(mgbo, ncpu);
	mgbo->g = g;
	mgbo->aux = init_hzmaux();
	mgbo->min_sm = min_sm;
	mgbo->refine = refine;
	mgbo->obj_id = 0xFFFFFFFFU;
	mgbo->ret = 0;
	thread_end_init(mgbo);
	beg_counter(cnt);
	ncand = 0;
	for(node_id=0;node_id<g->n_rd;node_id++){
		run_counter(cnt, 1000, stderr, 0);
		if(is_dead_node_strgraph(g, node_id)) continue;
		n = ref_sgnodev(g->nodes, node_id);
		if(n->bogs[1][0][0] && n->bogs[1][1][0]) continue;
		clear_u4v(candidates);
		clear_u32hash(hash);
		put_u32hash(hash, (node_id << 1) | 0);
		put_u32hash(hash, (node_id << 1) | 1);
		collect_graph_candidates_wtgbo(g, node_id, max_ext, closed_alns, candidates, heap, hash);
		//fprintf(stderr, " -- %u in %s -- %s:%d --\n", (unsigned)candidates->size, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		if(candidates->size == 0) continue;
		ncand += candidates->size;
		rdoff = g->rdoffs->buffer[node_id];
		rdlen = g->rdlens->buffer[node_id];
		reset_hzmaux(aux);
		for(i=0;i<rdlen;i++) add_tseq_hzmaux(aux, get_basebank(g->rdseqs, rdoff + i));
		ready_hzmaux(aux);
		for(i=0;i<candidates->size;i++){
			thread_wait_one(mgbo);
			if(mgbo->ret) ret += output_wtgbo(g, mgbo->obj_id, mgbo->qry_id, mgbo->dir, mgbo->aux, biedges, out);
			if(mgbo->obj_id != node_id){
				copy_index_hzmaux(mgbo->aux, aux);
				mgbo->obj_id = node_id;
			} else if(mgbo->contained){
				break;
			}
			mgbo->qry_id = candidates->buffer[i] >> 1;
			mgbo->dir    = candidates->buffer[i] & 0x01;
			thread_wake(mgbo);
		}
	}
	thread_wait_all(mgbo);
	end_counter(cnt, stderr);
	fprintf(stderr, "[%s] %llu candidates\n", date(), (unsigned long long)ncand);
	thread_beg_close(mgbo);
	if(mgbo->ret) ret += output_wtgbo(g, mgbo->obj_id, mgbo->qry_id, mgbo->dir, mgbo->aux, biedges, out);
	free_hzmaux(mgbo->aux);
	thread_end_close(mgbo);
	free_u4v(candidates);
	free_u8v(heap);
	free_u32hash(hash);
	return ret;
}

typedef struct {
	uint32_t node_id:30, dir:1, end:1;
	int pos;
} mark_t;
define_list(markv, mark_t);

void collect_anchor_candidates_wtgbo(StringGraph *g, uint32_t node_id, u64hash *closed_alns, u8v *candidates, markv *marks, uuhash *hash){
	sg_node_t *n;
	sg_edge_t *e;
	mark_t *m;
	uuhash_t *u, U;
	uint64_t lval, *lptr;
	uint32_t i, k;
	int exists, len, beg, end, dir;
	n = ref_sgnodev(g->nodes, node_id); // whatever node
	len = g->rdlens->buffer[node_id];
	clear_markv(marks);
	for(k=0;k<2;k++){
		for(i=0;i<n->edge_cnts[k];i++){
			e = ref_sgedgev(g->edges, n->edge_offs[k] + i);
			if(is_dead_node_strgraph(g, e->node_id)) continue; // avoid to calculate overlap on contained reads
			if(k){
				beg = len - (e->off + edge_overlap_strgraph(g, node_id, k, i));
				end = len - e->off;
			} else {
				beg = e->off;
				end = e->off + edge_overlap_strgraph(g, node_id, k, i);
			}
			dir = e->dir ^ k;
			push_markv(marks, (mark_t){e->node_id, dir, 0, beg});
			push_markv(marks, (mark_t){e->node_id, dir, 1, end});
		}
	}
	sort_array(marks->buffer, marks->size, mark_t, a.pos > b.pos);
	clear_uuhash(hash);
	for(i=0;i<marks->size;i++){
		m = ref_markv(marks, i);
		U.key = m->node_id; U.val = m->dir;
		if(m->end){
			remove_uuhash(hash, U);
			reset_iter_uuhash(hash);
			while((u = ref_iter_uuhash(hash))){
				lval = ovl_uniq_long_id(m->node_id, u->key, m->dir ^ u->val);
				lptr = prepare_u64hash(closed_alns, lval, &exists);
				if(exists) continue;
				*lptr = lval;
				push_u8v(candidates, lval);
			}
		} else {
			put_uuhash(hash, U);
		}
	}
}

uint64_t abo_core_wtgbo(StringGraph *g, HZMAux *aux, float min_sm, int refine, int ncpu, u64hash *closed_alns, sgbiedgev *biedges, FILE *out){
	u8v *candidates;
	markv *marks;
	uuhash *hash;
	uint64_t ret, ncand;
	uint32_t node_id, i;
	def_counter(cnt);
	thread_preprocess(mgbo);
	candidates = init_u8v(1024);
	marks = init_markv(1024);
	hash = init_uuhash(1023);
	ret = 0;
	thread_beg_init(mgbo, ncpu);
	mgbo->g = g;
	mgbo->aux = init_hzmaux();
	copy_param_hzmaux(mgbo->aux, aux);
	mgbo->min_sm = min_sm;
	mgbo->refine = refine;
	mgbo->obj_id = 0xFFFFFFFFU;
	mgbo->ret = 0;
	thread_end_init(mgbo);
	beg_counter(cnt);
	ncand = 0;
	for(node_id=0;node_id<g->n_rd;node_id++){
		run_counter(cnt, 1000, stderr, 0);
		clear_u8v(candidates);
		clear_uuhash(hash);
		collect_anchor_candidates_wtgbo(g, node_id, closed_alns, candidates, marks, hash);
		if(candidates->size == 0) continue;
		ncand += candidates->size;
		for(i=0;i<candidates->size;i++){
			thread_wait_one(mgbo);
			if(mgbo->ret) ret += output_wtgbo(g, mgbo->obj_id, mgbo->qry_id, mgbo->dir, mgbo->aux, biedges, out);
			reset_hzmaux(mgbo->aux);
			mgbo->obj_id = candidates->buffer[i] >> 33;
			mgbo->qry_id = (candidates->buffer[i] >> 1) & 0xFFFFFFFFU;
			mgbo->dir    = candidates->buffer[i] & 0x01;
			thread_wake(mgbo);
		}
	}
	thread_wait_all(mgbo);
	end_counter(cnt, stderr);
	fprintf(stderr, "[%s] %llu candidates\n", date(), (unsigned long long)ncand);
	thread_beg_close(mgbo);
	if(mgbo->ret) ret += output_wtgbo(g, mgbo->obj_id, mgbo->qry_id, mgbo->dir, mgbo->aux, biedges, out);
	free_hzmaux(mgbo->aux);
	thread_end_close(mgbo);
	free_u8v(candidates);
	free_markv(marks);
	free_uuhash(hash);
	return ret;
}

int usage(){
	printf(
	"WTGBO: Overlapper based on overlap graph\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtgbo [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -i <string> Long reads sequences file(s), + *\n"
	" -b <string> Long reads retained region, often from wtobt, +\n"
	"             Format: read_name\\toffset\\tlength\\toriginal_len\n"
	" -j <string> Overlap file(s), + *\n"
	"             Format: reads1\\t+/-\\tlen1\\tbeg1\\tend1\\treads2\\t+/-\\tlen2\\tbeg2\\tend2\\tscore\n"
	" -L <string> Load pairs of read name from file, will avoid to calculate overlap them again, + [NULL]\n"
	" -s <int>    Minimum alignment score, [200]\n"
	" -m <float>  Minimum alignment identity, [0.6]\n"
	" -u <int>    Maximum margin of alignment, [100]\n"
	" -o <string> Output file of new overlaps, *\n"
	" -9 <string> Record pairs of sequences have beed aligned regardless of successful, including pairs from '-L'\n"
	"             Format: read1\\tread2\n"
	" -f          Force overwrite output file\n"
	" -c <int>    Minimum estimated coverage of edge to be trusted, [1]\n"
	"             edge coverage is calculated by counting overlaps that can replace this edge\n"
	" -Q          Use number of matches as alignment score\n"
	" -q <float>  Best score cutoff, say best overlap MUST have alignment score >= <-r> * read's best score [0.95]\n"
	" -H          Turn off homopolymer compression\n"
	" -z <int>    Smaller kmer size (z-mer), 5 <= <-z> <= %d, [10]\n"
	" -Z <int>    Filter high frequency z-mers, maybe repetitive, [100]\n"
	" -y <int>    Zmer window, [800]\n"
	" -R <int>    Minimum size of seeding region within zmer window, [200]\n"
	" -r <int>    Minimum size of total seeding region for zmer windows, [300]\n"
	" -l <int>    Maximum variant of uncompressed sizes between two matched hz-kmer, [2]\n"
	" -M <int>    Alignment penalty: match, [2]\n"
	" -X <int>    Alignment penalty: mismatch, [-5]\n"
	" -O <int>    Alignment penalty: insertion or deletion, [-3]\n"
	" -E <int>    Alignment penalty: gap extension, [-1]\n"
	" -T <int>    Alignment penalty: read end clipping, 0: distable HSP extension, otherwise set to -50 or other [-50]\n"
	" -w <int>    Minimum bandwidth, iteratively doubled to maximum [50]\n"
	" -W <int>    Maximum bandwidth, [3200]\n"
	" -n          Refine the alignment\n"
	" -N <int>    Max turns of iteration, [5]\n"
	"\n"
	, HZM_MAX_SEED_ZMER);
	return 1;
}

int main(int argc, char **argv){
	obj_desc_t wsg = markv_obj_desc;
	obj_desc_t ttt = wsg;
	wsg = ttt;
	StringGraph *g;
	sgbiedgev *biedges;
	u64hash *closed_alns;
	HZMAux *aux;
	FileReader *fr;
	Sequence *seq;
	sg_biedge_t *b;
	cplist *pbs, *ovls, *obss, *obts;
	char *output, *pairoutf;
	FILE *out, *pairout;
	unsigned long long n, nn, i, val;
	int c, edgecov_cutoff, overwrite, min_score, margin, mat_score, max_ext, max_iter, iter;
	int ncpu, w, W, ew, M, X, O, E, T, hz, zsize, kwin, kstep, ztot, zovl, zcut, kvar, refine;
	ztot = 0;
	int wushigang = ztot;
	int tmp = wushigang;
	wushigang = tmp;
	float best_score_cutoff, min_id;
	pbs = init_cplist(4);
	ovls = init_cplist(4);
	obss = init_cplist(4);
	obts = init_cplist(4);
	output = NULL;
	pairoutf = NULL;
	min_score = 200;
	min_id = 0.6;
	margin = 100;
	edgecov_cutoff = 1;
	mat_score = 0;
	best_score_cutoff = 0.95;
	overwrite = 0;
	max_ext = 0;
	max_iter = 5;
	ncpu = 1;
	w = 50;
	ew = 800;
	W = 3200;
	M = 2;
	X = -5;
	O = -3;
	E = -1;
	T = -50;
	hz = 1;
	zsize = 10;
	kwin = 800;
	kstep = 0;
	ztot = 300;
	zovl = 200;
	zcut = 100;
	kvar = 2;
	refine = 0;
	while((c = getopt(argc, argv, "hi:b:j:L:s:m:u:o:9:fQq:c:t:Hz:Z:y:l:r:R:w:e:W:M:X:O:E:T:nN:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'i': push_cplist(pbs, optarg); break;
			case 'b': push_cplist(obts, optarg); break;
			case 'j': push_cplist(ovls, optarg); break;
			case 'L': push_cplist(obss, optarg); break;
			case 's': min_score = atoi(optarg); break;
			case 'm': min_id = atof(optarg); break;
			case 'u': margin = atoi(optarg); break;
			case 'o': output = optarg; break;
			case '9': pairoutf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'Q': mat_score = 1; break;
			case 'q': best_score_cutoff = atof(optarg); break;
			case 'c': edgecov_cutoff = atoi(optarg); break;
			case 't': ncpu = atoi(optarg); break;
			case 'H': hz = 0; break;
			case 'z': zsize = atoi(optarg); break;
			case 'Z': zcut = atoi(optarg); break;
			case 'y': kwin = atoi(optarg); break;
			case 'l': kvar = atoi(optarg); break;
			case 'r': ztot = atof(optarg); break;
			case 'R': zovl = atof(optarg); break;
			case 'w': w = atoi(optarg); break;
			case 'e': ew = atoi(optarg); break;
			case 'W': W = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'O': O = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'T': T = atoi(optarg); break;
			case 'n': refine = 1; break;
			case 'N': max_iter = atoi(optarg); break;
			default: return usage();
		}
	}
	if(output == NULL) return usage();
	if(pbs->size == 0) return usage();
	if(ovls->size == 0) return usage();
	if(!overwrite && file_exists(output)){
		fprintf(stderr, "File exists! '%s'\n\n", output);
		return usage();
	}
	g = init_strgraph();
	g->min_score = min_score;
	g->min_id = min_id;
	g->max_ovl_margin = margin;
	g->mat_score = mat_score;
	fr = fopen_m_filereader(pbs->size, pbs->buffer);
	fprintf(stderr, "[%s] loading reads\n", date());
	seq = NULL;
	while(fread_seq(&seq, fr)){
		if(seq->qual.size == seq->seq.size * 7){ // f5q format
			push_read5q_strgraph(g, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size, seq->qual.string);
		} else push_read_strgraph(g, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
		if((g->n_rd % 10000) == 0){
			fprintf(stderr, "\r%u", g->n_rd); fflush(stderr);
		}
	}
	fclose_filereader(fr);
	fprintf(stderr, "\r[%s] Done, %u reads\n", date(), (unsigned)g->n_rd); fflush(stderr);
	if(obts->size){
		fprintf(stderr, "[%s] loading reads obt information\n", date());
		if((fr = fopen_m_filereader(obts->size, obts->buffer)) == NULL) exit(1);
		while((c = fread_table(fr)) != -1){
			if(fr->line->string[0] == '#') continue;
			if(c < 3) continue;
			set_read_clip_strgraph(g, get_col_str(fr, 0), atoi(get_col_str(fr, 1)), atoi(get_col_str(fr, 2)));
		}
		fclose_filereader(fr);
		fprintf(stderr, "[%s] Done\n", date()); fflush(stderr);
	} else {
		fprintf(stderr, "[%s] No obt information\n", date()); fflush(stderr);
	}
	generate_nodes_strgraph(g);
	aux = init_hzmaux();
	aux->zsize = zsize;
	aux->hz = hz;
	aux->zwin = kwin;
	aux->zstep = kstep;
	aux->zovl = zovl;
	aux->zmax = zcut;
	aux->zvar = kvar;
	aux->w = w;
	aux->W = W;
	aux->ew = ew;
	aux->rw = w;
	aux->M = M;
	aux->X = X;
	aux->I = O;
	aux->D = O;
	aux->E = E;
	aux->T = T;
	out = strcmp(output, "-")? fopen(output, "w") : stdout;
	closed_alns = init_u64hash(1023);
	if(obss->size){
		uint64_t naln;
		uint32_t pb1, pb2;
		if((fr = fopen_m_filereader(obss->size, obss->buffer)) == NULL) exit(1);
		fprintf(stderr, "[%s] loading pairs of read name that alreadly tested\n", date());
		naln = 0;
		while(fread_table(fr) != -1){
			if(fr->line->string[0] == '#') continue;
			if((naln % 10000) == 0){ fprintf(stderr, "\r%llu", (unsigned long long)naln); fflush(stderr); }
			naln ++;
			if((pb1 = kv_get_cuhash(g->rdname2id, get_col_str(fr, 0))) == 0xFFFFFFFFU) continue;
			if((pb2 = kv_get_cuhash(g->rdname2id, get_col_str(fr, 1))) == 0xFFFFFFFFU) continue;
			val = ovl_uniq_long_id(pb1, pb2, 0);
			put_u64hash(closed_alns, val);
		}
		fclose_filereader(fr);
		fprintf(stderr, "\r[%s] there were %llu existing tested pairs\n", date(), (unsigned long long)closed_alns->count);
	}
	biedges = NULL;
	iter = 0;
	while(iter < max_iter){
		iter ++;
		fprintf(stderr, "---------------------------\n");
		fprintf(stderr, "[%s] iteration %d\n", date(), iter); fflush(stderr);
		clear_sgedgev(g->edges);
		zeros_bitvec(g->node_status);
		zeros_bitvec(g->node_flags);
		zeros_bitvec(g->node_atts);
		if(iter == 1){
			fr = fopen_m_filereader(ovls->size, ovls->buffer);
			fprintf(stderr, "[%s] loading alignments\n", date()); fflush(stderr);
			biedges = load_overlaps_strgraph(g, fr, closed_alns);
			fclose_filereader(fr);
			fprintf(stderr, "[%s] Done\n", date()); fflush(stderr);
		} else {
			fprintf(stderr, "[%s] bulding edges\n", date()); fflush(stderr);
			for(i=0;i<g->nodes->size;i++){
				ref_sgnodev(g->nodes, i)->edge_cnts[0] = 0;
				ref_sgnodev(g->nodes, i)->edge_cnts[1] = 0;
			}
			for(i=0;i<biedges->size;i++){
				b = ref_sgbiedgev(biedges, i);
				if(ref_sgnodev(g->nodes, b->node_id[0])->edge_cnts[b->dir[0]] >= SG_MAX_EDGE) continue;
				if(ref_sgnodev(g->nodes, b->node_id[1])->edge_cnts[!b->dir[1]] >= SG_MAX_EDGE) continue;
				ref_sgnodev(g->nodes, b->node_id[0])->edge_cnts[b->dir[0]] ++;
				ref_sgnodev(g->nodes, b->node_id[1])->edge_cnts[!b->dir[1]] ++;
			}
			load_overlaps_core_strgraph(g, biedges);
			fprintf(stderr, "[%s] Done\n", date()); fflush(stderr);
		}
		fprintf(stderr, "[%s] calculating edge coverage ...\n", date());
		cal_edge_coverage_strgraph(g);
		n = remove_duplicate_edges_strgraph(g);
		fprintf(stderr, "[%s] removed %llu duplicate edges\n", date(), n);
		fprintf(stderr, "[%s] Done\n", date());
		n = mask_contained_reads_strgraph(g, NULL);
		fprintf(stderr, "[%s] masked %llu contained reads\n", date(), n);
		n = mask_low_cov_edge_strgraph(g, edgecov_cutoff);
		fprintf(stderr, "[%s] masked %llu low coverage (<%u) edges\n", date(), n, edgecov_cutoff);
		n = best_overlap_strgraph(g, best_score_cutoff);
		fprintf(stderr, "[%s] 'best_overlap' cut %llu non-best edges\n", date(), n);
		fprintf(stderr, "[%s] graph based overlapping\n", date());
		n = gbo_core_wtgbo(g, aux, max_ext, min_id, refine, ncpu, closed_alns, biedges, out);
		fprintf(stderr, "[%s] Done, %llu new overlaps\n", date(), n);
		nn = n;
		fprintf(stderr, "[%s] anchoring based overlapping\n", date());
		n = abo_core_wtgbo(g, aux, min_id, refine, ncpu, closed_alns, biedges, out);
		fprintf(stderr, "[%s] Done, %llu new overlaps\n", date(), n);
		nn += n;
		fflush(out);
		if(nn == 0) break;
	}
	if(strcmp(output, "-")) fclose(out);
	if(pairoutf){
		pairout = fopen(pairoutf, "w");
		uint64_t *ptr;
		uint32_t id1, id2;
		reset_iter_u64hash(closed_alns);
		while((ptr = ref_iter_u64hash(closed_alns))){
			id1 = (*ptr) >> 33;
			id2 = ((*ptr) & 0xFFFFFFFFU) >> 1;
			fprintf(pairout, "%s\t%s\n", g->rdnames->buffer[id1], g->rdnames->buffer[id2]);
		}
		fclose(pairout);
	}
	free_u64hash(closed_alns);
	free_sgbiedgev(biedges);
	free_hzmaux(aux);
	free_cplist(pbs);
	free_cplist(ovls);
	free_cplist(obss);
	free_cplist(obts);
	free_strgraph(g);
	return 0;
}

