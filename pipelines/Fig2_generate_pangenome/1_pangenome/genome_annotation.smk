import pandas as pd

configfile: "config/SP_pangenome.yaml"

samples_df = pd.read_csv("tsv/SP_allseqs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
SEQ = {row.sample_id: {"sample_id": row.sample_id, "seq": row.seq_path} for row in samples_df.itertuples()} 

rule all:
    input:
        expand("results/{genus}/prokka/{sample}/{sample}.gff", genus=config["genus"],sample=SAMPLES),
        expand("results/{genus}/roary/core_gene_alignment.aln", genus=config["genus"]),
        expand("results/{genus}/roary/gene_presence_absence.csv", genus=config["genus"]),
        expand("results/{genus}/roary/tree.newick", genus=config["genus"])

rule prokka:
    """
    Obtain gff files needed for Roary
    """
    input:
        inputseqs=lambda wildcards: SEQ[wildcards.sample]["seq"]
    output:
        "results/{genus}/prokka/{sample}/{sample}.gff"
    params:
        genus=config["genus"]
    conda: 
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/{genus}/logs/prokka/{sample}_prokka.log"
    shell:
        """
        prokka --force --cpus {threads} --kingdom Bacteria --genus {params.genus} \
        --outdir "results/{params.genus}/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.inputseqs} &> {log}
        """

rule roary:
    """
    Obtain core gene alignment and gene presence absence file
    """
    input:
        gff_files=expand("results/{genus}/prokka/{sample}/{sample}.gff", sample=SAMPLES, genus=config["genus"])
    output:
        "results/{genus}/roary/core_gene_alignment.aln",
        "results/{genus}/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="results/{genus}/roary", 
       genus=config["genus"]
    threads: 4
    log:
        "results/{genus}/logs/roary/roary.log"
    shell:
        """
        rm -r results/{params.genus}/prokka_roary/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas
    """
    input:
        aln="results/{genus}/roary/core_gene_alignment.aln"
    output:
        nwk="results/{genus}/roary/tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "results/{genus}/logs/FastTree/FastTree.log"
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
        "results/{genus}/fastani/fastani.out"
    conda:
        "envs/fastani.yaml"
    log:
        "results/{genus}/logs/fastani/fastani.log"
    shell:
        """
        fastANI --ql {input.queryseqs} -r {input.ref} -o {output} >> {log} 2>&1
        """