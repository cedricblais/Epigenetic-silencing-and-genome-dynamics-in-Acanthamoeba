#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o viral_recall.log
#$ -cwd
#$ -pe threaded 2

GENOME=<Path to assembly>
OUTPUT=<Output directory>


source activate python37-generic
python viralrecall.py -i $GENOME -p $OUTPUT -db marker -t 2
conda deactivate