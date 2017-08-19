#!/bin/usr/sh
file1=$1
out1=$2
metric=$3
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6  -m 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -o $out1  --too-short-output $file1.tooshort.fastq $file1  > $metric
