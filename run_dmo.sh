#!/bin/bash

echo "Please modify the input file and delete this line"; exit

echo "==========================================="
cat $0
echo "==========================================="

NCPU=64
IDENTITY=0.1
MARGIN=300
DOT_MATRIX="-U 128 -U 64 -U 160 -U 1.0 -U 0.05"

date

# comment the below line after generated wt.fa
wtpre -J 10000 input_raw_pacbio_seqs.fa >wt.fa

date

# overlap
wtzmo -t $NCPU -i wt.fa  -fo wt.dmo.ovl -k 16 -z 10 -Z 16 $DOT_MATRIX -m $IDENTITY -A 1000

date

# clip read
wtclp -i wt.dmo.ovl -fo wt.dmo.obt -d 3 -FT -m $IDENTITY -k $MARGIN

date

# layout
wtlay -i wt.fa -b wt.dmo.obt -j wt.dmo.ovl -fo wt.dmo.lay -m $IDENTITY -w $MARGIN -r 0.95 -c 1

date

seq_n50.pl wt.dmo.lay.utg

#date

#wtcns -t $NCPU -i wt.dmo.lay -fo wt.dmo.cns.fa

#date

#seq_n50.pl wt.dmo.cns.fa
