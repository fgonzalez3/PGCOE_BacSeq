import pandas as pd

configfile: "config/genome_annotation.yaml"

samples_df = pd.read_csv("tsv/SP_TB_seqs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
SEQ = {row.sample_id: {"sample_id": row.sample_id, "seq": row.seq_path, "bac": row.species} for row in samples_df.itertuples()}
sample_to_genus = {sample_id: details['bac'] for sample_id, details in SEQ.items()}

# Generate paths manually using a list comprehension
genus_sample_paths = []
for sample, genus in sample_to_genus.items():
    path = f"results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff"
    genus_sample_paths.append(path)

# Function to get genus for a sample from the TSV sequences
def get_genus_from_tsv(sample):
    if sample in SEQ:
        return SEQ[sample]['bac']
    else:
        raise ValueError(f"Sample {sample} not found in TSV file.")

# Function to get genus for a sample from the config sequences
def get_genus_from_config(sample):
    for genus in config['samples']:
        if sample in config['samples'][genus]:
            return genus
    raise ValueError(f"Sample {sample} not found in config file.")

# Function to get path for a sample from the TSV sequences
def get_path_from_tsv(sample):
    if sample in SEQ:
        return SEQ[sample]['seq']
    else:
        raise ValueError(f"Sample {sample} not found in TSV file.")

# Function to get path for a sample from the config sequences
def get_path_from_config(sample):
    for genus in config['samples']:
        if sample in config['samples'][genus]:
            return config['samples'][genus][sample]
    raise ValueError(f"Sample {sample} not found in config file.")

rule all:
    input:
        expand("results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff", genus=["SP", "TB"], sample=list(config["samples"]["SP"].keys()) + list(config["samples"]["TB"].keys())),
        genus_sample_paths,
        expand("results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff", sample=SAMPLES, genus=[get_genus_from_tsv(sample) for sample in SAMPLES]),
        expand("results/Fig2/{genus}/genome_annotation/roary/accessory.header.embl", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/accessory_binary_genes.fa", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/accessory_binary_genes.fa.newick", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/accessory_graph.dot", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/accessory.tab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/blast_identity_frequency.Rtab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/clustered_proteins", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/core_accessory.header.embl", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/core_accessory.tab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/core_accessory_graph.dot", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/core_alignment_header.embl", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/gene_presence_absence.Rtab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/number_of_conserved_genes.Rtab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/number_of_genes_in_pan_genome.Rtab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/number_of_new_genes.Rtab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/number_of_unique_genes.Rtab", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/pan_genome_reference.fa", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/summary_statistics.txt", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/core_gene_alignment.aln", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/gene_presence_absence.csv", genus=["SP", "TB"], sample=SAMPLES),
        expand("results/Fig2/{genus}/genome_annotation/roary/tree.newick", genus=["SP", "TB"], sample=SAMPLES)

rule prokka:
    """
    Add information regarding genes, their location, and other features using Prokka for both SP and TB
    """
    input:
        seq_path = lambda wildcards: get_path_from_tsv(wildcards.sample)
    params:
        species = lambda wildcards: SEQ[wildcards.sample]['bac']
    output:
        "results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff"
    conda: 
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/Fig2/{genus}/genome_annotation/logs/prokka/{sample}_prokka.log"
    shell:
        """
        genus=$(if [ "{wildcards.genus}" = "SP" ]; then echo "Streptococcus"; else echo "Mycobacterium"; fi)
        prokka --force --cpus {threads} --kingdom Bacteria --genus $genus \
        --outdir "results/Fig2/{wildcards.genus}/genome_annotation/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.seq_path} &> {log}
        """

rule prokka_ref_outgroup:
    """
    Run Prokka individually on our refseq and the outgroup omitted from the above sequence list
    """
    input:
        ref_outgroup = "seqs/refs_outgroups/{sample}.fasta"
    output:
        gff = "results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff"
    params:
        genus = lambda wildcards: get_genus_from_config(wildcards.sample)  
    conda:  
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/Fig2/{genus}/genome_annotation/logs/{sample}_prokka.log"
    shell:
        """
        genus=$(if [ "{params.genus}" = "SP" ]; then echo "Streptococcus"; else echo "Mycobacterium"; fi)
        prokka --force --cpus {threads} --kingdom Bacteria --genus $genus \
        --outdir "results/Fig2/{genus}/genome_annotation/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.ref_outgroup} &> {log}
        """

rule roary:
    """
    Run pangenome analysis using Roary
    """
    input:
        gff_files=expand("results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff", genus=["SP", "TB"], sample=SAMPLES)
    output:
        "results/Fig2/{genus}/genome_annotation/roary/accessory.header.embl",
        "results/Fig2/{genus}/genome_annotation/roary/accessory_binary_genes.fa",
        "results/Fig2/{genus}/genome_annotation/roary/accessory_binary_genes.fa.newick",
        "results/Fig2/{genus}/genome_annotation/roary/accessory_graph.dot",
        "results/Fig2/{genus}/genome_annotation/roary/accessory.tab",
        "results/Fig2/{genus}/genome_annotation/roary/blast_identity_frequency.Rtab",
        "results/Fig2/{genus}/genome_annotation/roary/clustered_proteins",
        "results/Fig2/{genus}/genome_annotation/roary/core_accessory.header.embl",
        "results/Fig2/{genus}/genome_annotation/roary/core_accessory.tab",
        "results/Fig2/{genus}/genome_annotation/roary/core_accessory_graph.dot",
        "results/Fig2/{genus}/genome_annotation/roary/core_alignment_header.embl",
        "results/Fig2/{genus}/genome_annotation/roary/gene_presence_absence.Rtab",
        "results/Fig2/{genus}/genome_annotation/roary/number_of_conserved_genes.Rtab",
        "results/Fig2/{genus}/genome_annotation/roary/number_of_genes_in_pan_genome.Rtab",
        "results/Fig2/{genus}/genome_annotation/roary/number_of_new_genes.Rtab",
        "results/Fig2/{genus}/genome_annotation/roary/number_of_unique_genes.Rtab",
        "results/Fig2/{genus}/genome_annotation/roary/pan_genome_reference.fa",
        "results/Fig2/{genus}/genome_annotation/roary/summary_statistics.txt",
        "results/Fig2/{genus}/genome_annotation/roary/core_gene_alignment.aln",
        "results/Fig2/{genus}/genome_annotation/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="results/Fig2/{genus}/genome_annotation/roary"
    threads: 4
    log:
        "results/Fig2/{genus}/genome_annotation/logs/roary/roary.log"
    shell:
        """
        rm -r results/Fig2/{genus}/genome_annotation/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas
    """
    input:
        aln="results/Fig2/{genus}/genome_annotation/roary/core_gene_alignment.aln"
    output:
        nwk="results/Fig2/{genus}/genome_annotation/roary/tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "results/Fig2/{genus}/genome_annotation/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """