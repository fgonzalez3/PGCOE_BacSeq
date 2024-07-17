# This will create a CSV file with the paths to the SP GPSC fasta files

#!/bin/bash

# Directory where GPSCs are stored
GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/SP/"

# Output CSV file
output_file="SP_seqs.tsv"

# Write the header line
echo -e "sample_id\tseq_path" > $output_file

# Loop over the fasta files in GPSCs dir
for fasta_file in $GPSC_dir*.fasta; do
    # Extract the sample_id from the file name
    sample_id=$(basename $fasta_file .fasta)
    # Write a line to the CSV file with the path to the fasta file
    echo -e "$sample_id\t$fasta_file" >> $output_file
done

#!/bin/bash

# Directory where GPSCs are stored
GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/SP/"

# Output text file
output_file="SP_seqs_paths.txt"

# Loop over the fasta files in GPSCs dir and write their paths to the output file
for fasta_file in $GPSC_dir*.fasta; do
    echo "$fasta_file" >> $output_file
done