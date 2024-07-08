import pandas as pd

configfile: "config/read_aln.yaml"

samples_df = pd.read_csv("tsv/mNGS_samples.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2, "r1_output": row.r1_output, "r2_output": row.r2_output} for row in samples_df.itertuples()}
SUBSETS = config["SUBSETS"]  

rule all:
    input:
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.0123",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.amb",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.ann",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.bwt.2bit.64",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.pac",
        expand("results/SP/mNGS_aln/trim_adapters/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES),
        expand("results/SP/mNGS_aln/trim_adapters/{sample}/{sample}_val_2.fq.gz", sample=SAMPLES),
        expand("results/SP/mNGS_aln/bwa/{sample}/{sample}_aligned_sorted.bam", sample=SAMPLES), 
        expand("results/SP/mNGS_aln/coverage_subsets/{sample}/{sample}_{subset}_trimmed_aligned_sorted.bam", sample=SAMPLES, subset=SUBSETS), 
        expand("results/SP/mNGS_aln/coverage_subsets/{sample}/{sample}_{subset}_trimmed_aligned_sorted.csv", sample=SAMPLES, subset=SUBSETS)

rule trim_adapters:
    """
    Trim adapter sequences and QC
    """
    input:
        r1=lambda wildcards: READS[wildcards.sample]["r1"],
        r2=lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        r1 = "results/SP/mNGS_aln/trim_adapters/{sample}/{sample}_val_1.fq.gz",
        r2 = "results/SP/mNGS_aln/trim_adapters/{sample}/{sample}_val_2.fq.gz"
    conda:
        "envs/trim_galore.yaml"
    shell:
        """
        trim_galore --paired --fastqc --length 100 --cores 4 {input.r1} {input.r2} --output_dir results/SP/mNGS_aln/trim_adapters/{wildcards.sample} \
        --basename {wildcards.sample}
        """

rule bwa_build:
    """
    Create index of our consensus sequence that was used as PrimalScheme reference
    """
    input:
        ref=config["coordinate_sequence"]
    output:
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.0123",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.amb",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.ann",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.bwt.2bit.64",
        "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.pac"
    log:
        "results/SP/mNGS_aln/logs/bwa_build/build.log"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        mkdir -p results/SP/mNGS_aln/bwa_build/ref_index/
        bwa-mem2 index -p results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref {input.ref} > {log}
        """

rule bwa:
    """
    Align mNGS reads to our indexed reference sequence
    """
    input:
        fwd=lambda wildcards: READS[wildcards.sample]["r1"],
        rev=lambda wildcards: READS[wildcards.sample]["r2"],
        idx=[
            "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.0123",
            "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.amb",
            "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.ann",
            "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.bwt.2bit.64",
            "results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref.pac"
        ]
    output:
        aligned_sorted_bam = "results/SP/mNGS_aln/bwa/{sample}/{sample}_aligned_sorted.bam"
    conda:
        "envs/read_aln.yaml"
    shell: 
        """
        bwa-mem2 mem -t 32 -p results/SP/mNGS_aln/bwa_build/ref_index/indexed_ref {input.fwd} {input.rev} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aligned_sorted_bam} 
        """

rule coverage_subsets:
    """
    Obtain coverage stats for mNGS read alignments
    """
    input:
        "results/SP/mNGS_aln/bwa/{sample}/{sample}_aligned_sorted.bam"
    params:
        mapping_subset = lambda wildcards: float(wildcards.subset)
    output:
        subset_mapping = "results/SP/mNGS_aln/coverage_subsets/{sample}/{sample}_{subset}_trimmed_aligned_sorted.bam", 
        subset_csv = "results/SP/mNGS_aln/coverage_subsets/{sample}/{sample}_{subset}_trimmed_aligned_sorted.csv"
    log:
        "results/SP/mNGS_aln/logs/coverage_subsets/{sample}_{subset}_trimmed_aligned_sorted.log"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools view -@ 4 -b -F 4 -F 2048 -s {params.mapping_subset} {input} > {output.subset_mapping} && samtools coverage {output.subset_mapping} > {output.subset_csv} 2>> {log}
        """