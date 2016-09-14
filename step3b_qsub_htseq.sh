#!/bin/bash
export PYTHONPATH=$PYTHONPATH:~/bin/anaconda/bin/
htseq=/u/home/g/gandalm/bin/HTSeq-0.6.1p1/scripts/htseq-count
ref_gtf=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Annotation/gencode.vM10.annotation.gtf
ref_fa=/u/home/g/gandalm/project-geschwind/RefGenome/mmul10/GencodeM10/Sequence/GRCm38.p4.genome.fa
outdir=/u/home/g/gandalm/project-geschwind/TSC_MIA_Silva/data/HTSeq

file=$1
filenm=`basename $file`

${htseq} --format=bam --order=name --stranded=reverse --mode=union ${file} ${ref_gtf} > ${outdir}/$filenm
