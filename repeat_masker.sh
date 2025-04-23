
#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

cd $PWD


original_genome=<genome fasta>
root=<root for input file>
nthreads=10
RMdir=<output directory>

source activate repeatmasker+modeller-no-update

export PATH=/misc/scratch2/software/anaconda/envs/repeatmasker+modeller-no-update/custom/RepeatModeler-open-1.0.11:$PATH

echo 'Modeling and Masking genome'
BuildDatabase -name $root.index $original_genome
sleep 30s

RepeatModeler -pa $nthreads -database $root.index > $root.run.out 
RepeatMasker \
    -pa $nthreads \
    -s \
	-xsmall \
    -dir $RMdir \
    -a \
    -inv \
    -lib $root.index-families.fa \
     $original_genome

source deactivate
