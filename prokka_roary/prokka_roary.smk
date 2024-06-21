import pandas as pd

configfile: "config/SP_prokka_roary.yaml"

samples_df = pd.read_csv("SP_GPSCs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONSENSUS = {row.sample_id: {"sample_id": row.sample_id, "cons_seq": row.seq_path} for row in samples_df.itertuples()}


rule all:
    input:
        expand("results/SP/prokka_roary/prokka/{sample}/{sample}.gff", sample=SAMPLES),
        "results/SP/prokka_roary/roary/accessory.header.embl",
        "results/SP/prokka_roary/roary/accessory_binary_genes.fa",
        "results/SP/prokka_roary/roary/accessory_binary_genes.fa.newick",
        "results/SP/prokka_roary/roary/accessory_graph.dot",
        "results/SP/prokka_roary/roary/accessory.tab",
        "results/SP/prokka_roary/roary/blast_identity_frequency.Rtab",
        "results/SP/prokka_roary/roary/clustered_proteins",
        "results/SP/prokka_roary/roary/core_accessory.header.embl",
        "results/SP/prokka_roary/roary/core_accessory.tab",
        "results/SP/prokka_roary/roary/core_accessory_graph.dot",
        "results/SP/prokka_roary/roary/core_alignment_header.embl",
        "results/SP/prokka_roary/roary/gene_presence_absence.Rtab",
        "results/SP/prokka_roary/roary/number_of_conserved_genes.Rtab",
        "results/SP/prokka_roary/roary/number_of_genes_in_pan_genome.Rtab",
        "results/SP/prokka_roary/roary/number_of_new_genes.Rtab",
        "results/SP/prokka_roary/roary/number_of_unique_genes.Rtab",
        "results/SP/prokka_roary/roary/pan_genome_reference.fa",
        "results/SP/prokka_roary/roary/summary_statistics.txt",
        "results/SP/prokka_roary/roary/core_gene_alignment.aln",
        "results/SP/prokka_roary/roary/gene_presence_absence.csv",
        "results/SP/prokka_roary/roary/GPSC_pneumo_tree.newick"

rule prokka:
    """
    Add information regarding genes, their location, and other features using Prokka
    """
    input:
        prokka_input=lambda wildcards: CONSENSUS[wildcards.sample]["cons_seq"]
    output:
        "results/SP/prokka_roary/prokka/{sample}/{sample}.gff"
    conda: 
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/SP/prokka_roary/logs/prokka/{sample}_prokka.log"
    shell:
        """
        prokka --force --cpus {threads} --kingdom Bacteria --genus Streptococcus \
        --outdir "results/SP/prokka_roary/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.prokka_input} &> {log}
        """

rule roary:
    """
    Run pangenome analysis using Roary
    """
    input:
        gff_files=expand("results/SP/prokka_roary/prokka/{sample}/{sample}.gff", sample=SAMPLES)
    output:
        "results/SP/prokka_roary/roary/accessory.header.embl",
        "results/SP/prokka_roary/roary/accessory_binary_genes.fa",
        "results/SP/prokka_roary/roary/accessory_binary_genes.fa.newick",
        "results/SP/prokka_roary/roary/accessory_graph.dot",
        "results/SP/prokka_roary/roary/accessory.tab",
        "results/SP/prokka_roary/roary/blast_identity_frequency.Rtab",
        "results/SP/prokka_roary/roary/clustered_proteins",
        "results/SP/prokka_roary/roary/core_accessory.header.embl",
        "results/SP/prokka_roary/roary/core_accessory.tab",
        "results/SP/prokka_roary/roary/core_accessory_graph.dot",
        "results/SP/prokka_roary/roary/core_alignment_header.embl",
        "results/SP/prokka_roary/roary/gene_presence_absence.Rtab",
        "results/SP/prokka_roary/roary/number_of_conserved_genes.Rtab",
        "results/SP/prokka_roary/roary/number_of_genes_in_pan_genome.Rtab",
        "results/SP/prokka_roary/roary/number_of_new_genes.Rtab",
        "results/SP/prokka_roary/roary/number_of_unique_genes.Rtab",
        "results/SP/prokka_roary/roary/pan_genome_reference.fa",
        "results/SP/prokka_roary/roary/summary_statistics.txt",
        "results/SP/prokka_roary/roary/core_gene_alignment.aln",
        "results/SP/prokka_roary/roary/gene_presence_absence.csv"
    conda:
        "envs/roary_2.yaml"
    params:
       outdir="results/SP/prokka_roary/roary"
    threads: 4
    log:
        "results/SP/prokka_roary/logs/roary/roary.log"
    shell:
        """
        rm -r results/SP/prokka_roary/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas
    """
    input:
        aln="results/SP/prokka_roary/roary/core_gene_alignment.aln"
    output:
        nwk="results/SP/prokka_roary/roary/GPSC_pneumo_tree.newick"
    conda:
        "envs/roary_2.yaml"
    log:
        "results/SP/prokka_roary/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """

#rule generate_visualizations:
    #"""
    #Generate some visualizations to check that Roary output is fine
    #"""
    #input:
        #tree = "results/SP/roary/GPSC_pneumo_tree.newick", 
        #csv = "results/SP/roary/gene_presence_absence.csv"
    #output:
        #"pangenome_frequency.png",
        #"pangenome_matrix.png",
        #"pangenome_pie.png"
    #conda:
        #"pipelines/prokka_roary/envs/roary.yaml"
    #shell:
        #"""
        #wget -O roary_plots.py https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py
        #cp {input.tree} .
        #cp {input.csv} .
        #python roary_plots.py GPSC_pneumo_tree.newick gene_presence_absence.csv
        #"""