import pandas as pd

configfile: "config/SP_variant_calling.yaml"

samples_df = pd.read_csv("tsv/SP_amplicon_samples.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/ref_index/indexed_ref.0123", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.amb", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.ann", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.bwt.2bit.64", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.pac", genera=config["genera"]),
        expand("results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/variant_calls/{sample}/{sample}_variants.vcf.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/variant_calls/{sample}/{sample}_filtered_variants.vcf.gz", sample=SAMPLES, genera=config["genera"]),
        expand("reslts/{genera}/snippy/{sample}_snps.vcf", sample=SAMPLES, genera=config["genera"])

rule bwa_build:
    """
    Create index of our refseq
    """
    input:
        ref=config["ref"]
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
        "envs/variant_calling.yaml"
    shell:
        """
        mkdir -p results/{params.genera}/ref_index/
        bwa-mem2 index {input.ref} -p results/{params.genera}/ref_index/indexed_ref > {log} 2>&1
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
        "envs/variant_calling.yaml"
    params:
        genera=config["genera"]
    shell: 
        """
        bwa-mem2 mem -t 32 results/{params.genera}/ref_index/indexed_ref {input.fwd} {input.rev} > results/{params.genera}/alignments/{wildcards.sample}/{wildcards.sample}_aligned.sam
        samtools view -h -b -F 4 -F 2048 results/{params.genera}/alignments/{wildcards.sample}/{wildcards.sample}_aligned.sam | samtools sort -o {output.bam}
        """

rule variant_calling:
    input:
        bam="results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam"
    output:
        vcf="results/{genera}/variant_calls/{sample}/{sample}_variants.vcf.gz"
    conda:
        "envs/variant_calling.yaml"
    params:
        genera=config["genera"], 
        ref=config["ref"]
    shell:
        """
        bcftools mpileup -Ou -f {params.ref} {input.bam} -o temp.bcf
        bcftools call --ploidy 1 -vc -Oz -o {output.vcf} temp.bcf
        tabix -p vcf {output.vcf}
        rm temp.bcf
        """

rule variant_filter:
    """
    Filter variants based on quality and depth
    """
    input:
        vcf="results/{genera}/variant_calls/{sample}/{sample}_variants.vcf.gz"
    output:
        filtered_vcf="results/{genera}/variant_calls/{sample}/{sample}_filtered_variants.vcf.gz"
    conda:
        "envs/variant_calling.yaml"
    params:
        genera=config["genera"]
    shell:
        """
        bcftools filter -O z -o {output.filtered_vcf} -i 'QUAL>10 & DP>10' {input.vcf}
        """