#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 12


REFERENCE=<reference genome>
READS=<fastq reads>
PREFIX=<output prefix>
THREADS=12
THREADS_MINUS_ONE=11


source activate minimap2
minimap2 -a -Y -t $THREADS -x map-ont $REFERENCE $READS > $PREFIX.sam 
conda deactivate


/usr/bin/samtools view -@ $THREADS_MINUS_ONE -b -o $PREFIX.bam $PREFIX.sam
/usr/bin/samtools sort -@ $THREADS_MINUS_ONE -o $PREFIX.sorted.bam $PREFIX.bam
/usr/bin/samtools index -b $PREFIX.sorted.bam $PREFIX.sorted.bam.bai
