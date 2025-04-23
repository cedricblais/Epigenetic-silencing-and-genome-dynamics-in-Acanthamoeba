#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 15




#run the program HERE

source activate funannotate
source activate vardict-java
OUTNAME=<output name of your choice>
PROTEINS=<path to protein fasta>

/scratch2/software/interproscan-5.52-86.0/interproscan.sh -version

/scratch2/software/interproscan-5.52-86.0/interproscan.sh -i $PROTEINS --output-file-base $OUTNAME -cpu 15 -f xml,tsv --pathways -iprlookup --goterms

conda deactivate

