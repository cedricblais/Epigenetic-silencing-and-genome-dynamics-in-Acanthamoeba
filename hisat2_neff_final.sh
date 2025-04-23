#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10




FASTA=<path to genome>
PREFIX=<prefix for output>
FORWARD_READS=<forward read fastq>
REV_READS=<reverse read fastq>
THREADS=10
THREADS_MINUS_ONE=9

source activate hisat2
hisat2-build -f $FASTA $PREFIX
hisat2 -q -p $THREADS --phred33 -x $PREFIX -1 $FORWARD_READS -2 $REV_READS -S $PREFIX.sam
source deactivate
 
/usr/bin/samtools view -@ $THREADS_MINUS_ONE -b -o $PREFIX.bam $PREFIX.sam
/usr/bin/samtools sort -@ $THREADS_MINUS_ONE -o $PREFIX.sorted.bam $PREFIX.bam
/usr/bin/samtools index -b $PREFX.sorted.bam $PREFIX.sorted.bam.bai
