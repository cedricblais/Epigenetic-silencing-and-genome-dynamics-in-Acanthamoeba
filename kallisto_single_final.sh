#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd



#The following if for RNA-seq data that are not paired-end reads

REFERENCE=<path to exome fasta>
READS=<path to fastq>
INDEX=<name of index file to be generated>
DIR=<name of output directory to be created>
SD=<standard deviation of the fragments>
LEN=<fragment length>

source activate kallisto

kallisto index -i $INDEX $REFERENCE

kallisto quant -i $INDEX -o $DIR --single -l $LEN -s $SD $READS

conda deactivate

