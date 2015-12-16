*********
* wtzmo *
*********

The main overlaper of SMARTdenovo. Seed index -> Seed map -> Alignment.

wtzmo uses homopolymer compressed kmer as seeds. To trun off homopolymer
compression, see wtzmo -H 0. A seed contains [read_id, strand] in 4 bytes.
Please note that, wtzmo only keep the smaller between one kmer and its reverse
complementary, if equals, skips. Seeds are sorted by multiple value composed of
kmer, read_id and strand. A hashtable records {kmer_value =>
region_in_seed_array}. In implement, I scan the reads twice, first count kmer
in hastable and create a empty but reffered array, then fill the data. In the
concern of memory, wtzmo only index kmers tailed with C/A, halve the memory.
wtzmo -k 16, means the compressed kmer size is set to 16. wtzmo -K 500, the
upper K, means discarding high frequency (>500) kmer.

When query one read, wtzmo split the query into kmers. Then query each kmer
against the hashtable and seed array to find which [read_id,strand] cover
enough region on query sequence. Because seeds in a kmer are sorted by
read_id+strand, a process like merge-sort is used to calculate the covered size
for each existed read_id+strand. Here, wtzmo -A 500, means a heap records the
top 500 best covered reads. wtzmo -d 300 means a candidate must cover at least
300 bp uncompressed region on query. Suppose seed1 cover 1-21, seed2 cover
3-22, the sum is 1-22.

As we have candidates for one query, we can just use banded-SW to align them to
get the final alignments. However, it will spend too much time.

To speed it up, wtzmo build a similar seed index on the query sequence,
excepting the kmer_size is smaller, let call it zmer, wtzmo -z 10 is about
this. In this zmer-index, each zmer has its offset and strand. Filter high
frequency zmers by wtzmo -Z 100. Then match one candidate's zmers on
zmer-index. For now, we got pairs of matched seeds between candidate and query,
the next job is search the co-linear matched zmers.

wtzmo search the synteny within one window size specified by wtzmo -y 800. The
reason for why not search the synteny for whole overlapped region, is high rate
of INDELs make the distance between two nearby zmer high variable. wtzmo -R 200
has the similar meaning with wtzmo -d, but is it used in a window of 800bp.

When wtzmo matched lots of pairs of zmer-windows, it search coliner pairs in DP
manner. The final synteny of zmer-windows is called seed map. If a seed map
covers too small part of the query sequence, it will be discarded. It is
specified by wtzmo -r 300. Here, 300 means the solid zmer regions inside a
window, not the whole window. A candidate success in this step, will go to the
alignment procedure.

The alignment is splited into four parts.

a) matched zmer. just recover the compression, and add gaps.

b) within zmer window. Global alignment between two adject matched zmer.

c) Between two adject windows. Use banded global SW algorithm, the band size
   iteratively increase from wtzmo -w 50 to wtzmo -W 3200, until the alignment
   score is positive.

d) Extending of two ends. Local alignment. Its band width is specified by wtzmo
   -e 800.

wtzmo -T -50 means if the best score is bigger than the score at ending by less
than 50, will take the global alignment as result.

The order of alignment for candidates is sorted by their zmers covered size on
query. When wtzmo success to perform wtzmo -B 100 alignments, it will give up
to try more candidates. Use wtzmo -m 0.6 and wtmzo -s 200 to tell which
alignment is successful.

The order of alignment for candidates is also affected by estimation of
repetitive region. wtzmo calculates the zmer-window coverage on query, and
gives every position a weight, which is negatively correlated with coverage.
This weighting strategy avoids lots of useless computation. For example, a
query has a repeat in the middle, will recruit floods of false candidates that
contain the similar/same repeat. In this case, wtzmo will decrease the weight
of those flase candidates, and preffer to perform alignment on other
candidates. wtzmo -q 100 will set the weight of a position to zero when its
coverage reach to 100.

To get accurate CIGAR of alignment, add wtzmo -n. Refinement will be performed
based on the original alignment with flexible band width for each position.
That is wtzmo use global alignment on the alignment but with danymatic band
width. The basic band width is the same with wtzmo -w 50. When come to a gap of
size N, the bandwidth around it will increase N - d, d is the distance from
gap, the bandwidth increase by N inside the gap. wtzmo will trim the large
bandwidth into reasonable value based on DP alignment process.

To examine what happens during the alignment, use wtzmo -v, or -vv, -vvv,
-vvvv. The more -v will give more detailed information, catch them on the
STDERR.

If you have the read mask file, use wtzmo -b wt.cyc.clp to trim the reads
inside wtzmo. It will save the disk storage. wtzmo -J 1000 will filter all
reads have length less than 1000 bp. If you don't want to alignment a pair of
reads, add them in a file like 'already_aligned.pairs' and use wtzmo -L
already_aligned.pairs.

wtzmo is designed as a overlapper for de novo assembly. In the context of
String Graph, if one read is contained by another, it won't contribute in
graph. By default, when wtzmo find a candidate contained by query read, it will
mark the candidate as contained, and refuse to use it as query. Please note
that, the order of query sequence is sorted by read length DSC. But this
contained read can be used as candiate in other query. That is, we will know
which reads contain this shorter read, which will be informative in consensus
calling.

For huge genome, wtzmo can be run in low memory by wtzmo -G <N>. N is the
number of parts of kmer-index. Most of memory is spent on kmer-index. Thus,
split the kmer-index into N parts will reduce the memory into about 1/N. wtzmo
will record the candidates in additional memory and return to normal flow after
query all kmer-index parts.

To invoke wtzmo parallelly on multiple nodes, use -P <total_nodes> -p
<index_of_node>, e.g. -P 60 -p 0 on the first node of 60 nodes. It will take
the same memory on each node as only on one node.

Example: wtzmo -t 32 -i wt.fa -o wt.zmo.ovl

Ouput format: tab-delimited
 1, qry_name
 2, qry_strand: +/-
 3, qry_length
 4, qry_beg: 0-based
 5, qry_end: exclusive
 6, sbj_name
 7, sbj_strand
 8, sbj_length
 9, sbj_beg
10, sbj_end
11, score
12, idenity: 0.0 - 1.0
13, n_mat: number of matches
14, n_mis: mismatch
15, n_ins: insertion
16, n_del: deletion
17, cigar: CIGAR in SAM

*********
* wtobt *
*********

Trim reads base on overlaps. wtobt takes overlaps as input, trims high error
ending and chimeria. It try to retain max part of one read. First, it find a
max continous region that covered by other aligned reads. Then, it detect spurs
which one read get partial alignment on it, and counting how many reads cross
the spurs (as m) and how many reads get partial alignments (as n). If m no less
than half of the average depth at the spur, reject a chimeria. If m is bigger
than half of m, reject. otherwise, a chimra is detected. wtobt will retain the
larger part.

wtobt -i wt.fa -j wt.zmo.ovl -o wt.zmo.obt -c 2

wtobt output a read mask file with lines like: read_name offset length.

*********
* wtgbo *
*********

SMARTdenovo cannot find all of pairwise alignments as nearly all of other
aligners. It may miss key overlaps on graph. To rescue overlaps that are
important to assembly graph, wtgbo scans pairs of reads that might have
overlaps in: a) two reads have overlaps with the same other read, b) two reads
are connected by at most N steps, N is defined as 2. wtgbo try to align those
pairs in the same manner as wtzmo's candidates and query.

After generate new valid overlaps, wtgbo build a newer best overlap graph, and
infer new potential overlaps, until none new valid overlaps or max iterations.

wtgbo -t 32 -i wt.fa -j wt.zmo.ovl -o wt.zmo.gbo

wtgbo combine the algorithm of wtzmo and wtlay, its parameters like that of
wtzmo and wtlay. Its output is the same with wtzmo.

*********
* wtclp *
*********

The goal of wtclp is to maximize the total length of valid overlaps by trimming
reads or discarding reads (wtclp -F).

wtclp puts one read as reference, and tiles all reads having overlaps regards
of spurs. A function call_legal_overlaps_wtclp is used to calculate the length
of valid overlaps, all operations aim to maxmize the result of this function.

If -F is not specified, wtclp first clip high error ends by calculating
coverage. The threshold is set by wtclp -c 2. In discarding mode (wtclp -F), a
read is discarded or kept as whole.

wtclp detects all structure errors as chimera. Two algorithms are used: a)
depth depended, like wtobt; b) graph based. For (b, if one read connecting two
subgraphs by itself, it is thought to be structure error. wtclp checks whether
there is an alternative path formed by valid overlaps of tiled reads.

In NGS de novo assembly, people often plot kmer-distribution to estimate the
genome size and other features. TGS has very high error rate, it is not
possible to get genome size as usual. wtclp and wtobt use a new strategy to
estimate the genome size. As we get the total length of reads, if we know the
average depth, we can compute the genome size. The average depth is figured out
in the similar way as finding peak in NGS kmer-plot.

wtclp has the same output as wtobt. The reason of having both wtobt and wtclp
is to save the disk storage. Will discuss it in following section.

*********
* wtext *
*********

wtext is used in one kind of SMARTdenovo pipelines. It takes overlaps from
wtzmo and reads mask file from wtobt as input, and curates the overlaps. After
trimming reads, some alignments may be wrong, due to coordinate. wtext will try
to extend local alignment to a global alignment by wtext -T -100, like -T in
wtzmo.

wtext requires cigar of alignment, so that we need to keep cigars in the
alignment files. Cigars take most of the disk usage of alignment files, and
become hard in very big dataset. To save the disk usage, I introduced wtclp -F
mode. It either discards whole read or keeps it as orginal, and needn't to
invoke wtext to curate alignments.

*********
* wtlay *
*********

Roughly, wtlay implements BOG to generate layout of reads. It may take another
size of this document to discribe it. Just list some parameters.

 - wtlay -w 100. If an overlap is not end-to-end, but leaving N bp unaligned,
   and the N is no greater than 100, wtlay trust it as true overlap.

 - wtlay -c 1. Given a edge/overlap E between two node/reads A and B, the
   coverage of E is computed as how many edges draw out the same path as E. If
   the coverage of E is less than 1, will mask E as unreliable.

 - wtlay -r 0.95. which overlap is the best overlap? As high INDEL rate, wtlay
   doesn't trust the longest one on faith.  wtlay says the best overlap should
   have alignment score no less than 0.95 * <max score of its same strand>.

 - wtlay -q 0.4. I feel uneasy about it. wtlay doesn't merge bubbles, but cut
   up one path instead, which will leave islands.  To avoid to output them,
   wtlay filtered an unitig haing more than 40% of its length aligned on
   another untig.  It may bring the assembly size down.

 - wtlay -Q gCwgBgRURg. Having funs!

The prefix of output files is specified by wtlay -o <wt.lay>. 

wt.lay.utg is a fasta file, contains uncorrected sequences of unitigs.

wt.lay is a layout file, contains all the information needed by consensus caller.

>utg1
Y/N	read_name	strand	offset	length	sequence
...
Y/N: whether this read is used to build backbone or not. Reads in backbone should be not contained by others.
strand: for human read
offset: the offset to previous Y-starting read.
length: for huamn read
sequence: trimmed if need, reversed if need, the direct-in-use sequence used in calling consensus without any addtional operation.

wt.lay.<N>.dot is the graphviz source file.

*********
* wtcns *
*********

wtcns implemented the DAGCon algorithm described in HGAP paper. Alignment
algorithm is integrated into wtcns, thus doesn't need other alignment tool.
wtcns takes the layout file as input, and outputs consensus sequences of fasta
format.

The consensus sequences from wtcns are much accurater than PacBio reads, may
reach to 99.7%, but still cannot fit the need of genome assembly. If you have
other tools (e.g. Quiver) to improve the consensus sequences, please reduce the
number of iterations to save time, by wtcns -n 1 or -n 2.

**********
* wtcorr *
**********

Correcting long noisy reads on DBG from short accurate reads using k-mer
moving. The DBG contains a smaller kmer (e.g. k=25) for k-mer moving, and a
bigger kmer (e.g. k=41) to verify the path of kmer moving. The bigger kmers are
stored in counting-bloom-filter to save the memory.

Building the DBG
 wtcorr -t 32 -i read_1.fa -i read_2.fa -k 25 -K 41 -w sr.k25K41.dbg

Correcting
 wtcorr -t 32 -r sr.k25K41.dbg -c 5 -0 3 -1 2 -L -o corrected.fa raw.fa

wtcorr is very slow in large genome. It doesn't trust any matched kmer, and try
to start graph alignment from any possible matched kmer. Given a matched kmer
between PacBio read and DBG, it performs forward and backword kmer moving, and
merges the alignments. By default, it is local alignment mode. It will give
more than one corrected fragments. There may be more than hundreds of matched
kmers, there would be many kmers are random-matched (depends on the genome size
and PacBio errors). So that, when come to a kmer contained in prevois local
alignment, it may loss the true. I have worked out a algorithm which can start
alignments from all matched kmers, and use dynamic programming to avoid repeat
calculation. However, it is not implemented by now.

Whether error correction helps or not? Correcting pacbio reads using short
reads or itselves, MUST answer this question: Can we distinguish repetitive
sequences simplely? Here, repetitive is defined by calculation in the
background of high sequencing error, espically INDELs. For less repetitive
genomes, error correction will works well and simplify the process of assembly.
However, TGS is expected to solve complicated genomes, which always have flood
of repeats, will error correction still works well? In SMARTnodeovo, I choose
to assemble un-corrected long reads.

*********
* wtcyc *
*********

Align read against its reverse complementary sequence, to detect reads which
may miss its adpter. wtcyc generate a read mask file to discribe which part of
reads can be used in further analysis. User needn't process the raw reads file,
downstream programs can recongize this mask file.

 wtcyc -t 32 -i wt.fa -o wt.cyc.clp -a wt.cyc.info

wt.cyc.clp is the resulting mask file. wt.cyc.info is the alignments.

wtcyc is obselete, because it may clip the real DNA palindrome (exists in large
genome), and cause breakpionts on genome. I classfy this kind of structure
error as chimeric, and solve it in wtclp.
