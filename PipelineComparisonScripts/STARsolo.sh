#!/bin/bash

#getting the important options, namelue input files and the output directory, other options will be set in the script and should not change

while getopts f:r:o: flag
do
	case "$flag" in
		f) forward=${OPTARG};;
		r) reverse=${OPTARG};;
		o) output=${OPTARG};;
	esac
done

echo ${forward}
echo ${reverse}
echo ${output}


#STARsolo
echo "STEP1: STARsolo $(date)\n" > ${timefile}

STAR --runThreadN 8 \
	 --genomeDir /media/draco/lucask/genomes/mouse/M28_ALL/M28_ALL_STAR2-7-10 \
	 --readFilesIn ${forward} ${reverse} \
	 --readFilesCommand zcat \
	 --outSAMtype BAM SortedByCoordinate \
	 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
	 --sjdbGTFfile /media/draco/lucask/genomes/mouse/M28_ALL/genocde.vM28.basic.annotation.gtf \
	 --soloType CB_UMI_Complex \
	 --soloCBposition 0_10_0_17 0_48_0_55 0_86_0_93 \
	 --soloUMIposition 0_0_0_9 \
	 --soloCBwhitelist /media/draco/lucask/SPLTsq_Round1BC.txt /media/draco/lucask/SPLTsq_Round2BC.txt /media/draco/lucask/SPLTsq_Round3BC.txt \
	 --soloCBmatchWLtype EditDist_2 \
	 --soloBarcodeReadLength 0 \
	 --soloFeatures Gene GeneFull Velocyto \
	 --soloCellReadStats Standard \
	 --soloOutFileNames /media/draco/lucask/AlgorithmTest/STARsolo/ features.tsv barcodes.tsv matrix.mtx

#END
echo "END: $(date)" >> ${timefile}
