#!/bin/usr/sh
file1=$1
file2=$2
out1=$3
out2=$4
metric=$5
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6\
  -m 18 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT\
  -A GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\
  -o $out1 -p $out2 \
 $file1 $file2 > $metric
