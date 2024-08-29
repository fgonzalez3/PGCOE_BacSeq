import pandas as pd

configfile: "config/SP_variant_calling.yaml"

samples_df = pd.read_csv("tsv/contigs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONTIGS = {row.sample_id: row.contig_path for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/variant_calls/{sample}/{sample}_variants.vcf.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/variant_calls/{sample}/{sample}_filtered_variants.vcf.gz", sample=SAMPLES, genera=config["genera"])

rule minimap2_align:
    """
    Align contigs to refseq using Minimap2
    """
    input:
        contigs=lambda wildcards: CONTIGS[wildcards.sample],
        ref=config["ref"]
    output:
        bam="results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam"
    conda:
        "envs/variant_calling.yaml"
    params:
        genera=config["genera"]
    shell: 
        """
        minimap2 -ax asm5 {input.ref} {input.contigs} > results/{params.genera}/alignments/{wildcards.sample}/{wildcards.sample}_aligned.sam
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