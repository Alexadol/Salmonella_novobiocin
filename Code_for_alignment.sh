#!/bin/bash

#MAKE SURE THAT HISAT2 AND SAMTOOLS IS INSTALLED 

#You should create directory with all files with reads which should be aligned at the genome
PATH_TO_READS='./Reads/*'

mkdir SAM_GenBank
mkdir BAM_GenBank
mkdir SORTED_BAM_GenBank

#Further command can be used for creating genome index. If you need to create this - uncomment further row and make sure that genome file is in the directory where script is launched
#hisat2-build Salmonella_genome.fasta Salmonella_genome_gb 


for READ in $PATH_TO_READS; do
  hisat2  --dta -x Salmonella_genome_gb -U $READ -S ./SAM_GenBank/${READ##*/}.sam
  echo hisat OK
  echo ${READ##*/}
  samtools view -S -b ./SAM_GenBank/${READ##*/}.sam > ./BAM_GenBank/${READ##*/}.bam
  samtools sort ./BAM_GenBank/${READ##*/}.bam -o ./SORTED_BAM_GenBank/${READ##*/}.sorted.bam
  echo samtools OK
  

