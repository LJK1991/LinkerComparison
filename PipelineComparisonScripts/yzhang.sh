#!/bin/bash

timefile="/media/draco/lucask/AlgorithmTest/yzhang/timelog.txt"

echo -e "START: $(date)\n" >> ${timefile}

conda activate yzhang
echo $CONDA_DEFAULT_ENV

split-seq all --fq1 /media/draco/lucask/AlgorithmTest/SRR6750041_1.fastq.gz \
              --fq2 /media/draco/lucask/AlgorithmTest/SRR6750041_2.fastq.gz \
              --output_dir /media/draco/lucask/AlgorithmTest/yzhang/ \
              --chemistry v1 \
              --nthreads 8 \
              --genome_dir /media/draco/lucask/AlgorithmTest/yzhang/mouse_genome/

conda deactivate
echo $CONDA_DEFAULT_ENV

echo -e "STOP: $(date)\n" >> ${timefile}

