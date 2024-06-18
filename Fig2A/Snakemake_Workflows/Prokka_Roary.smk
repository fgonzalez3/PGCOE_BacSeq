import pandas as pd

configfile: "config/SP_prokka_roary.yaml"

samples_df = pd.read_csv("SP_GPSCs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONSENSUS = {row.sample_id: {"sample_id": row.sample_id, "cons_seq": row.cons_seq} for row in samples_df.itertuples()}

rule all:
    input:
        expand("data/results/prokka_annotated_GPSCs/{sample}/{sample}.gff", sample=config["samples"]),
        "data/results/roary/accessory.header.embl",
        "data/results/roary/accessory_binary_genes.fa",
        "data/results/roary/accessory_binary_genes.fa.newick",
        "data/results/roary/accessory_graph.dot",
        "data/results/roary/accessory.tab",
        "data/results/roary/blast_identity_frequency.Rtab",
        "data/results/roary/clustered_proteins",
        "data/results/roary/core_accessory.header.embl",
        "data/results/roary/core_accessory.tab",
        "data/results/roary/core_accessory_graph.dot",
        "data/results/roary/core_alignment_header.embl",
        "data/results/roary/gene_presence_absence.Rtab",
        "data/results/roary/number_of_conserved_genes.Rtab",
        "data/results/roary/number_of_genes_in_pan_genome.Rtab",
        "data/results/roary/number_of_new_genes.Rtab",
        "data/results/roary/number_of_unique_genes.Rtab",
        "data/results/roary/pan_genome_reference.fa",
        "data/results/roary/summary_statistics.txt",
        "data/results/roary/core_gene_alignment.aln",
        "data/results/roary/gene_presence_absence.csv",
        "data/results/roary/GPSC_pneumo_tree.newick", 
        "pangenome_frequency.png", 
        "pangenome_matrix.png",
        "pangenome_pie.png"
        

rule prokka:
    """
    Add information regarding genes, their location, and other features using Prokka
    """
    input:
        prokka_input = lambda wildcards: config["samples"][wildcards.sample]
    output:
        "data/results/prokka_annotated_GPSCs/{sample}/{sample}.gff"
    conda: 
        "envs/prokka.yaml"
    threads: 4
    log:
        "data/results/logs/prokka_annotated_GPSCs/{sample}_prokka.log"
    shell:
        """
        prokka --force --cpus {threads} --kingdom Bacteria --genus Streptococcus \
        --outdir "data/results/prokka_annotated_GPSCs/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.prokka_input} &> {log}
        """

rule roary:
    """
    Run pangenome analysis using Roary
    """
    input:
        gff_files=expand("data/results/prokka_annotated_GPSCs/{sample}/{sample}.gff", sample=config["samples"])
    output:
        "data/results/roary/accessory.header.embl",
        "data/results/roary/accessory_binary_genes.fa",
        "data/results/roary/accessory_binary_genes.fa.newick",
        "data/results/roary/accessory_graph.dot",
        "data/results/roary/accessory.tab",
        "data/results/roary/blast_identity_frequency.Rtab",
        "data/results/roary/clustered_proteins",
        "data/results/roary/core_accessory.header.embl",
        "data/results/roary/core_accessory.tab",
        "data/results/roary/core_accessory_graph.dot",
        "data/results/roary/core_alignment_header.embl",
        "data/results/roary/gene_presence_absence.Rtab",
        "data/results/roary/number_of_conserved_genes.Rtab",
        "data/results/roary/number_of_genes_in_pan_genome.Rtab",
        "data/results/roary/number_of_new_genes.Rtab",
        "data/results/roary/number_of_unique_genes.Rtab",
        "data/results/roary/pan_genome_reference.fa",
        "data/results/roary/summary_statistics.txt",
        "data/results/roary/core_gene_alignment.aln",
        "data/results/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="data/results/roary"
    threads: 4
    log:
        "data/results/logs/roary/roary.log"
    shell:
        """
        rm -r data/results/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file 
    """
    input:
        aln = "data/results/roary/core_gene_alignment.aln"
    output:
        nwk = "data/results/roary/GPSC_pneumo_tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "data/results/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """

rule generate_visualizations:
    """
    Generate some visualizations using our newick
    """
    input:
        tree = "data/results/roary/GPSC_pneumo_tree.newick", 
        csv = "data/results/roary/gene_presence_absence.csv"
    output:
        "pangenome_frequency.png",
        "pangenome_matrix.png",
        "pangenome_pie.png"
    conda:
        "envs/roary.yaml"
    shell:
        """
        wget -O roary_plots.py https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py
        cp {input.tree} .
        cp {input.csv} .
        python roary_plots.py GPSC_pneumo_tree.newick gene_presence_absence.csv
        """
