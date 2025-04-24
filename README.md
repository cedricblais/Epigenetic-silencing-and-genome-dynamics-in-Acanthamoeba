# Bioinformatic analysis of viral DNA in Acanthamoeba str. C3 and Neff

This repository documents the bioinformatic tools used for the study "Epigenetic silencing and genome dynamics determine the fate of giant virus endogenizations in Acanthamoeba". Files associated with the code deposited here can be found at [10.5281/zenodo.15269025](https://doi.org/10.5281/zenodo.15269025).

## Viral sequence detection

Viral proteins were detected using a diamond BLASTp search of all predicted proteins and intergenic ORFs in Neff and C3 against a local copy of the NCBI database. `qsub diamond_blast.sh`

Taxonomy was assigned to hits using [acc2tax](https://github.com/richardmleggett/acc2tax). Hits to other Acanthamoeba proteins, as well as entries with taxonomy listed as 'NO TAXA FOUND', 'other entries' and 'Unknown' were filtered out using `grep -v ',Acanthamoeba,' $ | grep -v 'NO TAXA FOUND' | grep -v 'other entries' | grep -v 'Unknown'`. Proteins were retained as viral if either the top hit after filtering was viral, or if at least 5/10 top blast hits were viral.

Viral hallmark genes were detected using ViralRecall:

`qsub viral_recall_hallmark_genes.sh`

## Conservation and delineation of viral regions

The Neff and C3 genomes were aligned with mummer4, and the show-diff function was used to identify differences in the assemblies:
`qsub [mummer4_genome_alignment_final.sh](mummer4_genome_alignment_final.sh)`

Overlapping differences in the alignments were concatenated together to form 'divergent regions', using a custom script which accounts for translocations. Viral regions were defined as divergent regions which had at least one fully overlapping gene. Unconserved genes were defined as genes which fully overlapped with show-diff features corresponding to gaps in the alignment (GAP, BRK, JMP, INV, SEQ).
`python SCRIPT`

## Expression

To assess global expression levels, HISAT2 was used to align reads to the reference:

`qsub [hisat2_C3_final.sh](hisat2_C3_final.sh)`

The specific transcripts per million (TPM) values for all predicted genes and viral intergenic ORFs were calculated using kallisto:

`qsub kallisto_final.sh` for Neff.

`qsub kallisto_single_final.sh` for c3.

A t-test comparing the TPM values of conserved and unconserved viral and non-viral genes was performed:

`Rscript expression_violin.R`

## Methylation analysis

Methylation calls using nanopolish require aligning nanopore reads to a genome reference. This was done with minimap2:

`qsub minimap2_final.sh`

The resulting .bam file is then used as an input file for the nanopolish script:

`qsub call_methylation_nanopolish_final.sh`

A t-test comparing the methylation levels of CpG sites across different genomic contexts was then performed:

`Rscript genomic_context_methylation_violin.R`

## Mobile elements

Mobile elements were predicted using repeatmasker:

`qsub repeat_masker.sh`

Low_complexity and Simple_repeat predictions were removed from the repeatmasker .out file using `grep -v $file.out`. The file was then made R friendly by replacing sequences of multiple spaces with tabs and removing leading characters: `sed 's/  */ /g' file.out | cut -c2- | sed 's/ /       /g'`

A chi-square test was performed to compare the statistical significance of differences in mobile element composition of different genomic contexts:

`Rscript chi_square_repeat_masker.R`

## Functional annotation

The function of Neff proteins was predicted using interposcan:

`qsub interproscan.sh`

## Phylogenetic trees

Phylogenetic trees were generated using IQ-tree:

`qsub iqtree.sh`
