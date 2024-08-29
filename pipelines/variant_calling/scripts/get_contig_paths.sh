#!/bin/bash

# Directory where the contigs are stored
contigs_dir="/Users/flg9/Desktop/Developer/grubaugh_lab/snakemake_workflows/GitHub/PGCOE_Bacseq/pipelines/variant_calling/contigs"

# Output TSV file
output_file="contigs.tsv"

# Write the header line
echo -e "sample_id\tcontig_path" > $output_file

# Loop over the files in contigs_dir
for contig_file in "$contigs_dir"/*.fa; do
    # Check if contig_file is a file
    if [[ -f $contig_file ]]; then
        # Extract the base sample name
        base_sample_name=$(basename "$contig_file" .contigs_velvet.fa)

        # Write a line to the TSV file
        echo -e "$base_sample_name\t$contig_file" >> $output_file
    fi
done