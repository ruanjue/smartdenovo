#!/bin/bash

echo "Please modify the input file and delete this line"; exit

echo "==========================================="
cat $0
echo "==========================================="

NCPU=64
IDENTITY=0.6

date

# comment the below line after generated wt.fa
wtpre -J 10000 input_raw_pacbio_seqs.fa >wt.fa

date

# overlap
wtzmo -t $NCPU -i wt.fa  -fo - -k 16 -s 200 -m $IDENTITY | cut -f 1-16 >wt.zmo.ovl.short

date

# rescue overlaps on graph
wtgbo -t $NCPU -i wt.fa -j wt.zmo.ovl.short -fo - | cut -f 1-16 >wt.zmo.gbo.short

date

# clip read
wtclp -i wt.zmo.ovl.short -i wt.zmo.gbo.short -fo wt.zmo.obt -F -d 2

date

# layout
wtlay -i wt.fa -b wt.zmo.obt -j wt.zmo.ovl.short -j wt.zmo.gbo.short -fo wt.zmo.lay -s 200 -m $IDENTITY -R -r 1 -c 1

date

seq_n50.pl wt.zmo.lay.utg

#date

#wtcns -t $NCPU -i wt.zmo.lay -fo wt.zmo.cns.fa

#date

#seq_n50.pl wt.zmo.cns.fa
