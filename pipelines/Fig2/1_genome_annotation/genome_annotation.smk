import pandas as pd

configfile: "config/genome_annotation.yaml"

# load tsv
samples_df = pd.read_csv("tsv/SP_TB_seqs.tsv", sep='\t')

# create dictionary for tsv
sample_info = {row['sample_id']: row for index, row in samples_df.iterrows()}

# define out paths
output_paths = []
roary_output = []
fasttree_output = []
log_paths = []

for sample in samples_df['sample_id']:
    genus = sample_info[sample]['genus']
    if genus == "Streptococcus":
        output_paths.append(f"results/Fig2/SP/genome_annotation/prokka/{sample}/{sample}.gff")
        log_paths.append(f"results/Fig2/SP/genome_annotation/logs/prokka/{sample}_prokka.log")
        roary_output.append(f"results/Fig2/SP/genome_annotation/roary/core_gene_alignment.aln")
        fasttree_output.append(f"results/Fig2/SP/genome_annotation/roary/tree.newick")
    elif genus == "Mycobacterium":
        output_paths.append(f"results/Fig2/TB/genome_annotation/prokka/{sample}/{sample}.gff")
        log_paths.append(f"results/Fig2/TB/genome_annotation/logs/prokka/{sample}_prokka.log")
        roary_output.append(f"results/Fig2/TB/genome_annotation/roary/core_gene_alignment.aln")
        fasttree_output.append(f"results/Fig2/TB/genome_annotation/roary/tree.newick")

rule all:
    input:
        output_paths,
        expand("results/SP/prokka_roary/prokka/{sample}/{sample}.gff", sample=config["samples"]["SP"]),
        expand("results/TB/prokka_roary/prokka/{sample}/{sample}.gff", sample=config["samples"]["TB"]), 
        roary_output, 
        fasttree_output
       

rule prokka:
    input:
        seq_path=lambda wildcards: sample_info[wildcards.sample]['seq_path']
    params:
        genus=lambda wildcards: sample_info[wildcards.sample]['genus']
    output:
        gff="results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff"
    conda:
        "envs/prokka.yaml"
    threads: 4
    log:
        log="results/Fig2/{genus}/genome_annotation/logs/prokka/{sample}_prokka.log"
    shell:
        """
        echo "Debug: Processing sample {wildcards.sample} as genus {params.genus}"
        prokka --force --cpus {threads} --kingdom Bacteria --genus {params.genus} \
        --outdir "$(dirname {output.gff})" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.seq_path} &> {log.log}
        """

rule prokka_ref_outgroup:
    """
    Run Prokka individually on our refseq and the outgroup omitted from the above sequence list
    """
    input:
        ref_outgroup=lambda wildcards: config['samples'][wildcards.genus][wildcards.sample]
    output:
        "results/{genus}/prokka_roary/prokka/{sample}/{sample}.gff"
    conda:  
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/{genus}/prokka_roary/logs/prokka/{sample}_prokka.log"
    shell:
        """
        prokka --force --cpus {threads} --kingdom Bacteria --genus {genus} \
        --outdir "results/{genus}/prokka_roary/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.ref_outgroup} &> {log}
        """

rule roary:
    input:
        gff_files=expand("results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff", genus=["SP", "TB"], sample=samples_df['sample_id'])
    output:
        "results/Fig2/{genus}/genome_annotation/roary/core_gene_alignment.aln"
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