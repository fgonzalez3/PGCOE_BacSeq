import pandas as pd

configfile: "config/SP_pangenome.yaml"

samples_df = pd.read_csv("tsv/SP_allseqs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
SEQ = {row.sample_id: {"sample_id": row.sample_id, "seq": row.seq_path} for row in samples_df.itertuples()} 

rule all:
    input:
        expand("results/{target}/prokka/{sample}/{sample}.gff", target=config["target"],sample=SAMPLES),
        expand("results/{target}/roary/core_gene_alignment.aln", target=config["target"]),
        expand("results/{target}/roary/gene_presence_absence.csv", target=config["target"]),
        expand("results/{target}/roary/tree.newick", target=config["target"])

rule prokka:
    """
    Obtain gff files needed for Roary
    """
    input:
        inputseqs=lambda wildcards: SEQ[wildcards.sample]["seq"]
    output:
        "results/{target}/prokka/{sample}/{sample}.gff"
    params:
        target=config["target"]
    conda: 
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/{target}/logs/prokka/{sample}_prokka.log"
    shell:
        """
        prokka --force --cpus {threads} --kingdom Bacteria --genus {params.target} \
        --outdir "results/{params.target}/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.inputseqs} &> {log}
        """

rule roary:
    """
    Obtain core gene alignment and gene presence absence file
    """
    input:
        gff_files=expand("results/{target}/prokka/{sample}/{sample}.gff", sample=SAMPLES, target=config["target"])
    output:
        "results/{target}/roary/core_gene_alignment.aln",
        "results/{target}/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="results/{target}/roary", 
       target=config["target"]
    threads: 4
    log:
        "results/{target}/logs/roary/roary.log"
    shell:
        """
        rm -r results/{params.target}/prokka_roary/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas
    """
    input:
        aln="results/{target}/roary/core_gene_alignment.aln"
    output:
        nwk="results/{target}/roary/tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "results/{target}/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """

rule fastani:
    """
    Run FastANI between our RefSeq and our queryseqs 
    """
    input:
        queryseqs = "txt/SP_allseqs_paths.txt",
        ref = "seqs/SP/NC_017592.fasta"
    output:
        "results/{target}/fastani/fastani.out"
    conda:
        "envs/fastani.yaml"
    log:
        "results/{target}/logs/fastani/fastani.log"
    shell:
        """
        fastANI --ql {input.queryseqs} -r {input.ref} -o {output} >> {log} 2>&1
        """