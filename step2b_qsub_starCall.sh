#!/bin/bash

SAMTOOLS_call=/u/home/g/gandalm/project-geschwind/bin/samtools-1.3/samtools
STAR_call=/u/home/g/gandalm/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
jav=/u/local/apps/java/jre1.8.0_77/bin/java
pic=/u/home/g/gandalm/project-geschwind/bin/picard-2.6.2.jar
genomeFA=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
genomeDir=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Sequence/STAR_index
gtfFile=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
rootdir=/u/home/g/gandalm/project-geschwind/TSC_MIA_Silva
refFlat=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.refFlat.txt.gz
refBed=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.bed


name=$1

if [ ! -s $rootdir/data/STAR_bam/$name ]; then mkdir $rootdir/data/STAR_bam/$name; fi

#STAR Alignment to Genome
if [ ! -s ${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam ]; then
echo "starting STAR alignment"
$STAR_call \
--runThreadN 2 \
--genomeDir ${genomeDir} \
--outFileNamePrefix ${rootdir}/data/STAR_bam/$name/${name}. \
--readFilesCommand gunzip -c \
--readFilesIn ${rootdir}/data/fastq_merged/${name}_R1_001.fastq.gz ${rootdir}/data/fastq_merged/${name}_R2_001.fastq.gz \
--outSAMtype BAM Unsorted SortedByCoordinate
fi


#Samtools sorting and indexing of BAM file
$SAMTOOLS_call index ${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam

#PicardTools RNA alignment & quality metrics
$jav -Xmx4g -jar $pic ReorderSam \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.Aligned.sortedByCoord.out.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  REFERENCE=$genomeFA ## Reorder the .bam file according to the reference at hand

$jav -Xmx4g -jar $pic CollectAlignmentSummaryMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/alignment_stats.txt \
  ASSUME_SORTED=false \
  ADAPTER_SEQUENCE=null ## Collect alignment metrics if the file is not present

$jav -Xmx4g -jar $pic CollectRnaSeqMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/rnaseq_stats.txt \
  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
  REF_FLAT=$refFlat \
  ASSUME_SORTED=false ## Collect sequencing metrics if the file is not present

$jav -Xmx4g -jar $pic CollectGcBiasMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/gcbias_stats.txt \
  ASSUME_SORTED=false \
  CHART_OUTPUT=${rootdir}/data/STAR_bam/$name/gcbias_chart.pdf \
  SUMMARY_OUTPUT=${rootdir}/data/STAR_bam/$name/gcbias_summary.txt

$jav -Xmx4g -jar $pic CollectInsertSizeMetrics \
  REFERENCE_SEQUENCE=${genomeFA} \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/insert_size_metrics.txt \
  ASSUME_SORTED=false \
  HISTOGRAM_FILE=${rootdir}/data/STAR_bam/$name/insert_size_histogram.pdf

$jav -Xmx4g -jar $pic MarkDuplicates \
  INPUT=${rootdir}/data/STAR_bam/$name/${name}.reordered_reads.bam \
  METRICS_FILE=${rootdir}/data/STAR_bam/$name/duplication_stats.txt \
  ASSUME_SORTED=false \
  OUTPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads.bam \
  REMOVE_DUPLICATES=TRUE

$jav -Xmx4g -jar $pic SortSam \
  INPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads.bam \
  OUTPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads_sorted.bam \
  SORT_ORDER=coordinate ## Sort the de-duplicated file

$jav -Xmx4g -jar $pic BuildBamIndex \
  INPUT=${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads_sorted.bam

rm ${rootdir}/data/STAR_bam/$name/reordered_duplication_marked_reads.bam
rm ${rootdir}/data/STAR_bam/$name/reordered_reads.bam

