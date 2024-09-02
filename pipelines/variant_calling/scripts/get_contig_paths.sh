#!/bin/bash

# Directory where the contigs are stored
contigs_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/serotyping/contigs"

# Output TSV file
output_file="contigs.tsv"

# Write the header line
echo -e "sample_id\tcontig_path" > $output_file

# Loop over the subdirectories in contigs_dir
for sample_dir in "$contigs_dir"/*; do
    # Check if sample_dir is a directory
    if [[ -d $sample_dir ]]; then
        # Extract the sample ID from the directory name
        sample_id=$(basename "$sample_dir")

        # Loop over the .fa files in the subdirectory
        for contig_file in "$sample_dir"/*.fa; do
            # Check if contig_file is a file
            if [[ -f $contig_file ]]; then
                # Write a line to the TSV file
                echo -e "$sample_id\t$contig_file" >> $output_file
            fi
        done
    fi
done