#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

#nanopolish version 0.14.0
source activate nanopolish

FASTQ=<fastq file>
ALIGNMENT=<fastq reads assembled to assembly (should be a sorted bam file)>
ASSEMBLY=<genome assembly>
OUTPUT=<output file prefix of your choice>
FAST5_DIR=<directory of fast5 files>

#we index all fast5 in the fast5 directory to the fastq. We don't need to feed those files to the script directly but they are needed, which is why we run this from the directory.
nanopolish index -d $FAST5_DIR $FASTQ

#We call methylation on our reads.
nanopolish call-methylation -t 10 -r $FASQ -b $ALIGNMENT -g $ASSEMBLY > $OUTPUT

#We calculate the methylation % at each genomic site
python calculate_methylation_frequency.py -i $OUTPUT -s > ${OUTPUT}_frequency.tsv

conda deactivate