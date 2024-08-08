import pandas as pd

configfile: "config/SP_AMR.yaml"

samples_df = pd.read_csv("tsv/SP_amplicon_isolate_samples.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/shovill/{sample}/contigs.fa", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/contigs.gfa", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/shovill.corrections", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/shovill.log", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/spades.fasta", sample=SAMPLES, genera=config["genera"]), 
        "results/{genera}/abricate/abricate.tsv", genera=config["genera"],   
        expand("results/{genera}/abricate/{sample}.txt", sample=SAMPLES, genera=config["genera"])

rule shovill:
    """
    Assemble reads using Shovill 
    """ 
    input:
        r1=lambda wildcards: READS[wildcards.sample]["r1"],
        r2=lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        "results/{genera}/shovill/{sample}/contigs.fa", 
        "results/{genera}/shovill/{sample}/contigs.gfa", 
        "results/{genera}/shovill/{sample}/shovill.corrections", 
        "results/{genera}/shovill/{sample}/shovill.log", 
        "results/{genera}/shovill/{sample}/spades.fasta"
    params:
        genera=config["genera"]
    conda:
        "envs/shovill.yaml"
    shell:
        """
        shovill --outdir results/{params.genera}/shovill/{wildcards.sample} --R1 {input.r1} --R2 {input.r2} --force
        """

rule abricate:
    """
    Run in-silico AMR prediction
    """
    input:
        "results/{genera}/shovill/{sample}/contigs.fasta"
    output:
        "results/{genera}/abricate/abricate.tsv",
        "results/{genera}/abricate/{sample}.txt"
    conda:
        "envs/abricate.yaml"
    shell:
        """
        abricate {input} --minid 75 --mincov 75
        """