import pandas as pd

configfile: "config/TB_roary_reps.yaml"

samples_df = pd.read_csv("tsv/TB_reps.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
SEQ = {row.sample_id: {"sample_id": row.sample_id, "seq": row.seq_path} for row in samples_df.itertuples()}

rule all:
    input:
        config["file_paths"]

rule roary:
    """
    Run pangenome analysis using Roary on Parnas reps
    """
    input:
        gff_files=expand("/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/prokka_roary/results/TB/prokka_roary/prokka/{sample}/{sample}.gff", sample=SAMPLES)
    output:
        "results/TB/roary_reps/roary/accessory.header.embl",
        "results/TB/roary_reps/roary/accessory_binary_genes.fa",
        "results/TB/roary_reps/roary/accessory_binary_genes.fa.newick",
        "results/TB/roary_reps/roary/accessory_graph.dot",
        "results/TB/roary_reps/roary/accessory.tab",
        "results/TB/roary_reps/roary/blast_identity_frequency.Rtab",
        "results/TB/roary_reps/roary/clustered_proteins",
        "results/TB/roary_reps/roary/core_accessory.header.embl",
        "results/TB/roary_reps/roary/core_accessory.tab",
        "results/TB/roary_reps/roary/core_accessory_graph.dot",
        "results/TB/roary_reps/roary/core_alignment_header.embl",
        "results/TB/roary_reps/roary/gene_presence_absence.Rtab",
        "results/TB/roary_reps/roary/number_of_conserved_genes.Rtab",
        "results/TB/roary_reps/roary/number_of_genes_in_pan_genome.Rtab",
        "results/TB/roary_reps/roary/number_of_new_genes.Rtab",
        "results/TB/roary_reps/roary/number_of_unique_genes.Rtab",
        "results/TB/roary_reps/roary/pan_genome_reference.fa",
        "results/TB/roary_reps/roary/summary_statistics.txt",
        "results/TB/roary_reps/roary/core_gene_alignment.aln",
        "results/TB/roary_reps/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="results/TB/roary_reps/roary"
    threads: 4
    log:
        "results/TB/roary_reps/logs/roary/roary.log"
    shell:
        """
        rm -r results/TB/roary_reps/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas reps
    """
    input:
        aln = "results/TB/roary_reps/roary/core_gene_alignment.aln"
    output:
        nwk = "results/TB/roary_reps/roary/GPSC_pneumo_tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "results/TB/roary_reps/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """

rule fastani_refseq_reps:
    """
    Run FastANI between our RefSeq and our Parnas reps 
    """
    input:
        queryseqs = "txt/TB_parnas_reps.txt",
        ref = "seqs/TB_reps/H37Rv.fasta"
    output:
        "results/TB/roary_reps/fastani/fastani.out"
    conda:
        "envs/fastani.yaml"
    log:
        "results/TB/roary_reps/logs/fastani/fastani.log"
    shell:
        """
        fastANI --ql {input.queryseqs} -r {input.ref} -o {output} >> {log} 2>&1
        """