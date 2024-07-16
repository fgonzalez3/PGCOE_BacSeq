import yaml

with open('config/SP_genome_annotation.yaml', 'r') as file:
    config = yaml.safe_load(file)

full_paths = {filename: config['base_path'] + filename for filename in config['samples']}

rule all:
    input:
        expand("results/{genus}/fastani/fastani.out", genus=config["genus"])
       
rule prokka:
    """
    Obtain gff files needed for Roary
    """
    input:
        inputseqs = lambda wildcards: full_paths[wildcards.sample]
    output:
        gff="results/{genus}/prokka/{wildcards.sample}/{wildcards.sample}.gff"
    params:
        genus=config["genus"]
    conda:
        "envs/prokka.yaml"
    threads: 4
    log:
        log="results/{genus}/logs/prokka/{wildcards.sample}_prokka.log"
    shell:
        """
        echo "Processing sample {wildcards.sample} as genus {params.genus}"
        prokka --force --cpus {threads} --kingdom Bacteria --genus {params.genus} \
        --outdir "$(dirname {output.gff})" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.inputseqs} &> {log}
        """

rule roary:
    """
    Obtain core gene alignment and gene presence absence file
    """
    input:
        expand("results/Fig2/{genus}/genome_annotation/prokka/{sample}/{sample}.gff", genus=config["genus"], sample=config["samples"])
    output:
        "results/{genus}/roary/presence_absence.csv",
        "results/{genus}/roary/core_gene_alignment.aln"
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
        rm -rf results/{params.genus}/roary 
        roary -e -n -v -p {threads} -f {params.outdir} {input}
        """
    
rule fastree:
    """
    Generate newick file using core gene alignment from Roary
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
    Run FastANI between refseq and queryseqs
    """
    input:
        queryseqs = [config["base_path"] + sample for sample in config["samples"]],
        ref = config["base_path"] + config["ref"]
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