#!/bin/bash

# point to dir where seqs are stored
TB_GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/TB/"
SP_GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/SP/"

# specify output file
output_file="SP_TB_seqs.tsv"

# create header for output file
echo -e "sample_id\tseq_path\tspecies" > $output_file

# loop over files
for fasta_file in $TB_GPSC_dir*.fa; do
    sample_id=$(basename $fasta_file .fa)
    echo -e "$sample_id\t$fasta_file\tTB" >> $output_file
done

# Loop over the fasta files in SP GPSCs dir and add "SP" as species
for fasta_file in $SP_GPSC_dir*.fasta; do
    sample_id=$(basename $fasta_file .fasta)
    echo -e "$sample_id\t$fasta_file\tSP" >> $output_file
done