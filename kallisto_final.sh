#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd


REVERSE_READS=<reverse reads>
FORWARD_READS=<forward reads>
REFERENCE=<path to exome fasta>
INDEX=<name of index file to be generated>
DIR=<name of output directory to be created>


source activate kallisto

kallisto index -i $INDEX $REFERENCE

kallisto quant -i $INDEX $FORWARD_READS $REVERSE_READS -o $DIR

conda deactivate

