#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 15


#diamond-custom version: diamond version 2.0.14

DATABASE=<blast formatted database>
QUERY=<protein queries>
OUTPUT=<output file>

diamond-custom blastp --threads 15 -d $DATABASE -q $QUERY -o $OUTPUT  --salltitles --max-target-seqs 500 --evalue 0.001 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp --more-sensitive --masking seg  




