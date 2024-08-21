import yaml

configfile: "config/SP_AMR.yaml"

SAMPLES = list(config['samples'].keys())
READS = config['samples']

rule all:
    input:
        expand("results/{genera}/shovill/{sample}/contigs.fa", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/contigs.gfa", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/shovill.corrections", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/shovill.log", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/spades.fasta", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/abricate/{sample}/abricate.tsv", sample=SAMPLES, genera=config["genera"])

rule shovill:
    """
    Assemble reads using Shovill 
    """ 
    input:
        fwd=lambda wildcards: READS[wildcards.sample]["r1"],
        rev=lambda wildcards: READS[wildcards.sample]["r2"]
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
        shovill --outdir results/{params.genera}/shovill/{wildcards.sample} --R1 {input.fwd} --R2 {input.rev} --ram 100 --force
        """

rule abricate:
    """
    Run in-silico AMR prediction
    """
    input:
        "results/{genera}/shovill/{sample}/contigs.fa"
    output:
        "results/{genera}/abricate/{sample}/abricate.tsv"
    params:
        genera=config["genera"]
    conda:
        "envs/abricate.yaml"
    shell:
        """
        abricate {input} --minid 75 --mincov 75 > {output}
        """