#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 4



REF_NAME=<reference strain name>
QUE_NAME=<query strain name>
REFERENCE=<path to reference genome>
QUERY=<path to query genome>
BREAKLEN='1000'
OUTPUT_PREFIX=<prefix of your choice>

source activate mummer4
nucmer --version
nucmer $REFERENCE $QUERY --prefix=$OUTPUT_PREFIX --threads 4 --breaklen $BREAKLEN

DELTA_FILE=$OUTPUT_PREFIX'.delta'
dnadiff --d $DELTA_FILE
show-diff $DELTA_FILE > $DELTA_FILE'.show-diff'

conda deactivate
