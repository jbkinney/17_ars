#!/usr/bin/env bash

# Example usage
#$ ./prep_genome.sh genome/s288c.fa 

# Get arguments
genome_file=$1 # Genome fasta file

# Index genome
echo "Indexing genome..." >&2
bwa index $genome_file

# Create .genome file (listing chromosome lengths)
echo "Computing chromosome lengths..." >&2
samtools faidx $genome_file
awk -v OFS='\t' {'print $1,$2'} $genome_file.fai > $genome_file.genome

echo "Done." >&2
