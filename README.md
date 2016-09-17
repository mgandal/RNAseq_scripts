# RNAseq_Pipeline


### Software Used
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) / [MultiQC](http://multiqc.info/) to compile
* [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
* [Samtools](http://www.htslib.org/doc/samtools.html)
* [HTSeq-Count](http://www-huber.embl.de/HTSeq/)
 * Install [anaconda distribution of python](https://www.continuum.io/downloads)
 * Install HTSeq-Count build from [here](https://pypi.python.org/pypi/HTSeq)
 * Install [pysam module](https://pypi.python.org/pypi/pysam/0.8.3) for python
* [PicardTools](https://broadinstitute.github.io/picard/)

### Workflow
1. Download .FASTQ files from Basespace / Synapse / SRA
2. QC on FastQ files using FastQC / MultiQC
3. Merge fastq files from same sample across lanes
4. Align to reference genome with STAR
5. Index, sort, QC on BAM using Picard / samtools
6. Quantify reads with HTSeq-Count
