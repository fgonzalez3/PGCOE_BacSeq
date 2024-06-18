import pandas as pd

configfile: "config/SP_prokka_roary.yaml"

samples_df = pd.read_csv("SP_reps.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONSENSUS = {row.sample_id: {"sample_id": row.sample_id, "cons_seq": row.cons_seq} for row in samples_df.itertuples()}

rule all:
    input:
        "parnas/reps/SP/roary/accessory.header.embl",
        "parnas/reps/SP/roary/accessory_binary_genes.fa",
        "parnas/reps/SP/roary/accessory_binary_genes.fa.newick",
        "parnas/reps/SP/roary/accessory_graph.dot",
        "parnas/reps/SP/roary/accessory.tab",
        "parnas/reps/SP/roary/blast_identity_frequency.Rtab",
        "parnas/reps/SP/roary/clustered_proteins",
        "parnas/reps/SP/roary/core_accessory.header.embl",
        "parnas/reps/SP/roary/core_accessory.tab",
        "parnas/reps/SP/roary/core_accessory_graph.dot",
        "parnas/reps/SP/roary/core_alignment_header.embl",
        "parnas/reps/SP/roary/gene_presence_absence.Rtab",
        "parnas/reps/SP/roary/number_of_conserved_genes.Rtab",
        "parnas/reps/SP/roary/number_of_genes_in_pan_genome.Rtab",
        "parnas/reps/SP/roary/number_of_new_genes.Rtab",
        "parnas/reps/SP/roary/number_of_unique_genes.Rtab",
        "parnas/reps/SP/roary/pan_genome_reference.fa",
        "parnas/reps/SP/roary/summary_statistics.txt",
        "parnas/reps/SP/roary/core_gene_alignment.aln",
        "parnas/reps/SP/roary/gene_presence_absence.csv",
        "parnas/reps/SP/roary/GPSC_pneumo_tree.newick"

rule roary:
    """
    Run pangenome analysis using Roary
    """
    input:
        gff_files=expand("parnas/reps/SP/gffs/{sample}.gff", sample=SAMPLES)
    output:
        "parnas/reps/SP/roary/accessory.header.embl",
        "parnas/reps/SP/roary/accessory_binary_genes.fa",
        "parnas/reps/SP/roary/accessory_binary_genes.fa.newick",
        "parnas/reps/SP/roary/accessory_graph.dot",
        "parnas/reps/SP/roary/accessory.tab",
        "parnas/reps/SP/roary/blast_identity_frequency.Rtab",
        "parnas/reps/SP/roary/clustered_proteins",
        "parnas/reps/SP/roary/core_accessory.header.embl",
        "parnas/reps/SP/roary/core_accessory.tab",
        "parnas/reps/SP/roary/core_accessory_graph.dot",
        "parnas/reps/SP/roary/core_alignment_header.embl",
        "parnas/reps/SP/roary/gene_presence_absence.Rtab",
        "parnas/reps/SP/roary/number_of_conserved_genes.Rtab",
        "parnas/reps/SP/roary/number_of_genes_in_pan_genome.Rtab",
        "parnas/reps/SP/roary/number_of_new_genes.Rtab",
        "parnas/reps/SP/roary/number_of_unique_genes.Rtab",
        "parnas/reps/SP/roary/pan_genome_reference.fa",
        "parnas/reps/SP/roary/summary_statistics.txt",
        "parnas/reps/SP/roary/core_gene_alignment.aln",
        "parnas/reps/SP/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="parnas/reps/SP/roary"
    threads: 4
    log:
        "parnas/reps/SP/logs/roary/roary.log"
    shell:
        """
        rm -r parnas/reps/SP/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file 
    """
    input:
        aln = "parnas/reps/SP/roary/core_gene_alignment.aln"
    output:
        nwk = "parnas/reps/SP/roary/GPSC_pneumo_tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "parnas/reps/SP/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """
