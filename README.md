## Getting Started

```sh
# Download sample PacBio from the PBcR website
wget -O- http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz | tar zxf -
awk 'NR%4==1||NR%4==2' selfSampleData/pacbio_filtered.fastq | sed 's/^@/>/g' > reads.fa
# Install SMARTdenovo
git clone https://github.com/ruanjue/smartdenovo.git && (cd smartdenovo; make)
# Assemble (raw unitigs in wtasm.lay.utg; consensus unitigs: wtasm.cns)
smartdenovo/smartdenovo.pl reads.fa > wtasm.mak
make -f wtasm.mak
```

## Introduction

SMARTdenovo is a *de novo* assembler for PacBio and Oxford Nanopore (ONT)
data. It produces an assembly from all-vs-all raw read alignments without
an error correction stage. It also provides tools to generate accurate
consensus sequences, though a platform dependent consensus polish tools (e.g.
Quiver for PacBio or Nanopolish for ONT) are still required for higher
accuracy.

SMARTdenovo consists of several separate command line tools: **wtzmo** for read
overlapping, **wtgbo** to rescue missing overlaps, **wtclp** for identifying
low-quality regions and chimaera, and **wtcns** or **wtmsa** to produce better
unitig consensus. The `smartdenovo.pl` script provides a convenient interface
to call these programs in one go. If you do not care about the internal of
SMARTdenovo, you may simply run with:
```sh
/path/to/smartdenovo/smartdenovo.pl -p prefix reads.fa > prefix.mak
make -f prefix.mak
```
It calls other SMARTdenovo executables in the same directory containing
`smartdenovo.pl`. After assembly, the raw unitigs are reported in file
`prefix.lay.utg` and consensus unitigs in `prefix.cns`. If you want to know
more about how SMARTdenovo works in detail, please see [README-tools.md][rt].

[rt]: README-tools.md
