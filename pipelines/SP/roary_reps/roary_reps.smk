import pandas as pd

configfile: "config/SP_prokka_roary.yaml"

samples_df = pd.read_csv("SP_reps.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONSENSUS = {row.sample_id: {"sample_id": row.sample_id, "cons_seq": row.cons_seq} for row in samples_df.itertuples()}

rule all:
    input:
        config["file_paths"]

rule roary:
    """
    Run pangenome analysis using Roary on Parnas reps
    """
    input:
        gff_files=expand("parnas/reps/SP/gffs/{sample}.gff", sample=SAMPLES)
    output:
        "results/SP/roary_reps/accessory.header.embl",
        "results/SP/roary_reps/accessory_binary_genes.fa",
        "results/SP/roary_reps/accessory_binary_genes.fa.newick",
        "results/SP/roary_reps/accessory_graph.dot",
        "results/SP/roary_reps/accessory.tab",
        "results/SP/roary_reps/blast_identity_frequency.Rtab",
        "results/SP/roary_reps/clustered_proteins",
        "results/SP/roary_reps/core_accessory.header.embl",
        "results/SP/roary_reps/core_accessory.tab",
        "results/SP/roary_reps/core_accessory_graph.dot",
        "results/SP/roary_reps/core_alignment_header.embl",
        "results/SP/roary_reps/gene_presence_absence.Rtab",
        "results/SP/roary_reps/number_of_conserved_genes.Rtab",
        "results/SP/roary_reps/number_of_genes_in_pan_genome.Rtab",
        "results/SP/roary_reps/number_of_new_genes.Rtab",
        "results/SP/roary_reps/number_of_unique_genes.Rtab",
        "results/SP/roary_reps/pan_genome_reference.fa",
        "results/SP/roary_reps/summary_statistics.txt",
        "results/SP/roary_reps/core_gene_alignment.aln",
        "results/SP/roary_reps/gene_presence_absence.csv"
    conda:
        "pipelines/roary_reps/envs/roary.yaml"
    params:
       outdir="results/SP/roary_reps"
    threads: 4
    log:
        "results/SP/logs/roary_reps/roary.log"
    shell:
        """
        rm -r results/SP/roary_reps
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas reps
    """
    input:
        aln = "results/SP/roary_reps/core_gene_alignment.aln"
    output:
        nwk = "results/SP/roary_reps/GPSC_pneumo_tree.newick"
    conda:
        "pipelines/roary_reps/envs/roary.yaml"
    log:
        "results/SP/logs/FastTree_reps/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """

rule fastani_refseq_reps:
    """
    Run FastANI between our RefSeq and our Parnas reps 
    """
    input:
        queryseqs = "pipelines/SP/roary_reps/parnas_reps.txt",
        ref = "pipelines/SP/prokka_roary/GPSCs/NC_017592.fasta"
    output:
        "results/SP/fastani/fastani.out"
    conda:
        "pipelines/roary_reps/envs/fastani.yaml"
    log:
        "results/SP/logs/fastani/fastani.log"
    shell:
        """
        fastANI --ql {input.queryseqs} -r {input.ref} -o {output} >> {log} 2>&1
        """