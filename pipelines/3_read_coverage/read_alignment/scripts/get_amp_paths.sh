#!/bin/bash

# Directory where the raw reads are stored
raw_reads_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/read_alignment/Github/seqs/amplicon/*"

# Output TSV file
output_file="amplicon_samples.tsv"

# Write the header line
echo -e "sample_id\tr1\tr2" > $output_file

# Loop over the directories in raw_reads_dir
for sample_dir in $(ls -d $raw_reads_dir); do
    # Extract the base sample name
    base_sample_name=$(basename $sample_dir)

    # Construct the paths to the R1 and R2 files
    r1_path=$(ls $sample_dir/*R1_001.fastq.gz)
    r2_path=$(ls $sample_dir/*R2_001.fastq.gz)

    # Check if both R1 and R2 files exist
    if [[ -e $r1_path && -e $r2_path ]]; then
        # If both files exist, write a line to the TSV file
        echo -e "$base_sample_name\t$r1_path\t$r2_path" >> $output_file
    fi
done