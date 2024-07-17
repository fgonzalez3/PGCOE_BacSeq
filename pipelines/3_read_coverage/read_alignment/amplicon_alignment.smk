import pandas as pd

configfile: "config/SP_read_aln.yaml"

samples_df = pd.read_csv("tsv/SP_amplicons.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

SUBSETS = set(config['SUBSETS'])

rule all:
    input:
        expand("results/{genera}/ref_index/indexed_ref.0123", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.amb", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.ann", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.bwt.2bit.64", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.pac", genera=config["genera"]),
        expand("results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/coverage_subsets/{sample}/{sample}_{subset}_aligned_sorted.bam", sample=SAMPLES, subset=SUBSETS, genera=config["genera"]),
        expand("results/{genera}/coverage_subsets/{sample}/{sample}_{subset}_aligned_sorted.csv", sample=SAMPLES, subset=SUBSETS, genera=config["genera"])  

rule bwa_build:
    """
    Create index of our refseq
    """
    input:
        ref=config["coordinate_sequence"]
    output:
        "results/{genera}/ref_index/indexed_ref.0123",
        "results/{genera}/ref_index/indexed_ref.amb",
        "results/{genera}/ref_index/indexed_ref.ann",
        "results/{genera}/ref_index/indexed_ref.bwt.2bit.64",
        "results/{genera}/ref_index/indexed_ref.pac"
    params:
        genera=config["genera"]
    log:
        "results/{genera}/logs/bwa_build/build.log"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        mkdir -p results/{genera}/ref_index/
        bwa-mem2 index -p results/{genera}/ref_index/indexed_ref {input.ref} > {log}
        """

rule bwa:
    """
    Align raw reads to our indexed refseq
    """
    input:
        fwd=lambda wildcards: READS[wildcards.sample]["r1"],
        rev=lambda wildcards: READS[wildcards.sample]["r2"],
        idx=[
            "results/{genera}/ref_index/indexed_ref.0123",
            "results/{genera}/ref_index/indexed_ref.amb",
            "results/{genera}/ref_index/indexed_ref.ann",
            "results/{genera}/ref_index/indexed_ref.bwt.2bit.64",
            "results/{genera}/ref_index/indexed_ref.pac"
        ]
    output:
        bam="results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam"
    conda:
        "envs/read_aln.yaml"
    shell: 
        """
        bwa-mem2 mem -t 32 -p results/{genera}/ref_index/indexed_ref {input.fwd} {input.rev} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.bam} 
        """

rule coverage_subsets:
    """
    Obtain coverage for untrimmed amplicon sequences in subsets
    """
    input:
        "results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam"
    params:
        subset_vals = lambda wildcards: float(wildcards.subset)
    output:
        subset_mapping = "results/{genera}/coverage_subsets/{sample}/{sample}_{subset}_aligned_sorted.bam", 
        subset_csv = "results/{genera}/coverage_subsets/{sample}/{sample}_{subset}_aligned_sorted.csv"
    log:
        "results/{genera}/logs/coverage_subsets/{sample}_{subset}_aligned_sorted.log"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools view -@ 4 -b -F 4 -F 2048 -s {params.subset_vals} {input} > {output.subset_mapping} && samtools coverage {output.subset_mapping} > {output.subset_csv} 2>> {log}
        """