#!/bin/bash

ChIA-PET2 -g bwa_index/mm10.fa -b mm10.chrom.sizes -f reads_1.fastq -r reads_2.fastq -t 16 -A linkerA -B linkerB -m 0 -o output_chiapet2 -n data
