#!/bin/bash

echo "Please modify the input file and delete this line"; exit

echo "==========================================="
cat $0
echo "==========================================="

NCPU=64
IDENTITY="0.1"
MARGIN=300
DOT_MATRIX="-U 256 -U 64 -U 512 -U 0.1 -U 0.05"

date

# comment the below line after generated wt.fa
wtpre -J 10000 input_raw_pacbio_seqs.fa >wt.fa

date

# overlap
wtzmo -t $NCPU -i wt.fa  -fo wt.dmo -k 16 $DOT_MATRIX -m $IDENTITY -A 500

date

# clip read
wtobt -i wt.fa -j wt.dmo -fo wt.dmo.obt -c 2 -m $IDENTITY -w $MARGIN

date

# overlap again
wtzmo -t $NCPU -i wt.fa -b wt.dmo.obt  -fo wt.dmo.ovl -k 16 $DOT_MATRIX -m $IDENTITY -A 500

date

# layout
wtlay -i wt.fa -b wt.dmo.obt -j wt.dmo.ovl -fo wt.dmo.lay -m $IDENTITY -w $MARGIN -r 0.95 -c 1

date

seq_n50.pl wt.dmo.lay.utg

#date

#wtcns -t $NCPU -i wt.dmo.lay -fo wt.dmo.cns.fa

#date

#seq_n50.pl wt.dmo.cns.fa
