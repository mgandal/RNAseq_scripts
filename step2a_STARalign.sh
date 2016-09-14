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

cd ${rootdir}/data/fastq_merged

for file in *_R1_001.fastq.gz; do

name=`basename $file _R1_001.fastq.gz`
echo $name
qsub \
  -cwd -V -o ${rootdir}/code/log \
  -e ${rootdir}/code/log \
  -l h_rt=12:00:00,h_data=12G,highp -pe shared 8  \
  ${rootdir}/code/step2b_qsub_starCall.sh ${name}

done
