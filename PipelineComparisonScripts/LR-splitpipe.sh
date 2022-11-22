#!/bin/bash

while getopts f:r:o: flag
do
        case "$flag" in
                f) forward=${OPTARG};;
                r) reverse=${OPTARG};;
                o) output=${OPTARG};;
        esac
done

timefile="${output}/timelog_two.txt"

#echo -e "START: $(date)\n" >> ${timefile}

conda activate scvelo
echo $CONDA_DEFAULT_ENV

#python /media/draco/lucask/splitpipe/splitpipe/demultiplex.py find_bcs -f ${reverse} -o ${output} -t 8 --chunksize 10**5 --verbosity 2 --delete_input

echo -e "Add bam tag: $(date)\n" >> ${timefile}

python /media/draco/lucask/splitpipe/splitpipe/add_bam_tag_fastq.py -f ${forward} -r ${output}/_seq_corrected_bcs.tsv -o ${output} --cutoff 100 --verbosity 2

conda deactivate

conda activate STAR
echo $CONDA_DEFAULT_ENV

pushd ${output}
echo $(pwd)

echo -e "ALIGN: $(date)\n" >> ${timefile}
STAR --runThreadN 8 \
     --readFilesIn ${output}/_barcoded_cells.fastq \
     --genomeDir /media/draco/lucask/genomes/mouse/M28_ALL/M28_ALL_STAR2-7-10/ \
     --outSAMtype BAM SortedByCoordinate

echo -e "featureCounts $(date)\n" >> ${timefile}
/home/lucask/subread-2.0.3-source/bin/featureCounts -F GTF \
                                                    -a /media/draco/lucask/genomes/mouse/M28_ALL_STAR2-7-10 \
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
/home/lucask/miniconda3/bin/umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv.gz

conda deactivate
popd

now=$(date)
echo -e "STOP: ${now}\n" >> ${timefile}

