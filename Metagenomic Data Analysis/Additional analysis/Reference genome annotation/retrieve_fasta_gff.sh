#!/bin/bash

# Set variables for paths, select root genome ID according to genome of interest
input_fasta="/storage/homefs/ib21q318/UHGGv1/UHGG_reps.fasta"
output_dir="/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/prokka_annotation"
root_genome_id="GUT_GENOME000221"

# Count the number of sequences that start with the root genome ID
n=$(grep -c "^>${root_genome_id}_" $input_fasta)

# Loop through each contig and process it
for i in $(seq 1 $n)
do
    # Create the contig identifier (e.g., _1, _2, etc.)
    contig_id="${root_genome_id}_${i}"

    # Define the output FASTA file for the current contig
    output_fasta="${output_dir}/${contig_id}_reps.genes.fasta"

    # Extract sequences for the current contig and save to a new FASTA file
    awk -v contig_id="^>${contig_id}" '/^>/{p=0} $0 ~ contig_id{p=1} p' $input_fasta > $output_fasta

    # Run Prokka annotation on the extracted FASTA file
    prokka --outdir "${output_dir}/${contig_id}_prokka" --prefix "${contig_id}" $output_fasta
done

# Combine all GFF files into one
output_gff="${output_dir}/combined_annotations.gff"
cat ${output_dir}/GUT_GENOME000221_*_prokka/*.gff > $output_gff
