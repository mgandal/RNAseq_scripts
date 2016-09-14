#!/bin/bash

SAMTOOLS_call=/u/home/g/gandalm/project-geschwind/bin/samtools-1.3/samtools
STAR_call=/u/home/g/gandalm/bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
jav=/usr/bin/java
pic=/u/home/g/gandalm/project-geschwind/bin/picard-2.6.2.jar
genomeFA=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
genomeDir=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Sequence/STAR_index
gtfFile=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
rootdir=/u/home/g/gandalm/project-geschwind/TSC_MIA_Silva
refFlat=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.refFlat.txt.gz
refBed=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.bed


#qrsh -l h_rt=24:00:00,h_data=80G -pe shared 8
$STAR_call \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/STAR_index \
  --genomeFastaFiles /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa \
  --sjdbGTFfile /hp_shares/mgandal/datasets/refGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf \
  --sjdbOverhang 49
