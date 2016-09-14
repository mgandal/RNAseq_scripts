#!/bin/bash

bamdir=/u/home/g/gandalm/project-geschwind/TSC_MIA_Silva/data/STAR_bam/
codedir=/u/home/g/gandalm/project-geschwind/TSC_MIA_Silva/code/

cd $bamdir

for bam in */*Aligned.out.bam; do
  echo $bam
  qsub -V -cwd -e $codedir/log -o $codedir/log -l h_rt=04:00:00,h_data=8G,highp $codedir/step3b_qsub_htseq.sh $bamdir/$bam
done
