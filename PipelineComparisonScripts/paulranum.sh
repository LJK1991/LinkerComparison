#!/bin/bash

while getopts f:r:o: flag
do
	case "$flag" in
		f) forward=${OPTARG};;
		r) reverse=${OPTARG};;
		o) output=${OPTARG};;
	esac
done

timefile="${output}/timelog.txt"

echo -e "START: $(date)\n" > ${timefile}

conda activate scvelo
echo $CONDA_DEFAULT_ENV

/media/draco/lucask/SPLiT-Seq_demultiplexing/splitseqdemultiplex_0.2.1.sh -f ${forward} -r ${reverse} -o ${output} -m 1 -t 64000

conda deactivate

conda activate STAR
echo $CONDA_DEFAULT_ENV

pushd ${output}
echo $(pwd)

echo -e "ALIGN: $(date)\n" >> ${timefile}
STAR --runThreadN 8 \
     --readFilesIn MergedCells_1.fastq \
     --genomeDir /media/draco/lucask/genomes/mouse/M28_ALL/M28_ALL_STAR2-7-10/ \
     --outSAMtype BAM SortedByCoordinate

echo -e "featureCounts $(date)\n" >> ${timefile}
/home/lucask/subread-2.0.3-source/bin/featureCounts -F GTF \
						    -a /media/draco/lucask/genomes/mouse/M28_ALL/gencode.vM28.basic.annotation.gtf \
						    -t exon,gene,transcript,CDS,UTR \
						    -o gene_assigned \
						    -R BAM Aligned.sortedByCoord.out.bam \
						    -T 8
conda deactivate

conda activate scvelo
echo $CONDA_DEFAULT_ENV

echo -e "Sort and Index: $(date)\n" >> ${timefile}
/home/lucask/samtools/samtools sort -@ 8 Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
/home/lucask/samtools/samtools index -@ 8 assigned_sorted.bam

echo -e "Create count.mtx: $(date)\n" >> ${timefile}
/home/lucask/miniconda3/bin/umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --temp-dir /media/draco/lucask/SPLITCROP_Whole/ -I assigned_sorted.bam -S counts.tsv.gz

conda deactivate
popd

now=$(date)
echo -e "STOP: ${now}\n" >> ${timefile}
