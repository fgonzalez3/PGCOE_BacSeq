configfile: "config/ref_prokka_roary.yaml"

rule all:
    input:
        "s_pneumo_data/data/results/primer_fa/pos_primers.bed", 
        "s_pneumo_data/data/results/primer_fa/neg_primers.bed",
        "s_pneumo_data/data/results/primer_fa/pos.scheme.primer.fasta",
        "s_pneumo_data/data/results/primer_fa/neg.scheme.primer.fasta",
        expand("s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.1.bt2", sample=config["samples"]),
        expand("s_pneumo_data/data/results/bowtie/{sample}_primers_pos.bam", sample=config["samples"]),
        expand("s_pneumo_data/data/results/bowtie/{sample}_primers_neg.bam", sample=config["samples"]),
        expand("s_pneumo_data/data/results/samtools_merge/{sample}_merged.bam", sample=config["samples"]),
        expand("s_pneumo_data/data/results/samtools_depth/{sample}.depth", sample=config["samples"]), 
        expand("s_pneumo_data/data/results/depth_vis/{sample}_depth.png", sample=config["samples"]), 
        expand("s_pneumo_data/data/results/coverage_pc/{sample}_coverage.csv", sample=config["samples"])

rule split_bed:
    """
    Use awk to extract positive and negative strand primers 
    """
    input:
        "s_pneumo_data/input/scheme.primer.bed"
    output:
        pos_pr="s_pneumo_data/data/results/primer_fa/pos_primers.bed", 
        neg_pr="s_pneumo_data/data/results/primer_fa/neg_primers.bed"
    shell:
        """
        awk '$6 == "+" {{print}}' {input} > {output.pos_pr} 
        awk '$6 == "-" {{print}}' {input} > {output.neg_pr} 
        """

rule get_fasta:
    """
    Convert BED to FASTA for downstream alignment 
    """
    input:
        pos_bed = "s_pneumo_data/data/results/primer_fa/pos_primers.bed", 
        neg_bed = "s_pneumo_data/data/results/primer_fa/neg_primers.bed"
    output:
        bedfasta_pos = "s_pneumo_data/data/results/primer_fa/pos.scheme.primer.fasta", 
        bedfasta_neg = "s_pneumo_data/data/results/primer_fa/neg.scheme.primer.fasta"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/scripts/get_indiv_fasta.py {input.pos_bed} {input.neg_bed} {output.bedfasta_pos} {output.bedfasta_neg}
        """

rule bowtie2_build:
    """
    Create index of GPSCs
    """
    input:
        queryseqs = lambda wildcards: config["samples"][wildcards.sample]
    output:
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.1.bt2",
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.2.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.3.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.4.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.rev.1.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.rev.2.bt2"
    conda:
        "envs/bowtie_samtools.yaml"
    log:
        "s_pneumo_data/data/results/logs/{sample}_bowtie_build.log"
    shell:
        """
        bowtie2-build -f {input.queryseqs} s_pneumo_data/data/results/ref_index/{wildcards.sample}_indexed_ref >> {log} 2>&1
        """

rule bowtie_primers:
    """
    Align primers to GPSCs
    """
    input:
        pos = rules.get_fasta.output.bedfasta_pos, 
        neg = rules.get_fasta.output.bedfasta_neg,
        idx = [
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.1.bt2",
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.2.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.3.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.4.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.rev.1.bt2", 
        "s_pneumo_data/data/results/ref_index/{sample}_indexed_ref.rev.2.bt2"  
        ]
    output:
        aln_pos = "s_pneumo_data/data/results/bowtie/{sample}_primers_pos.bam",
        aln_neg = "s_pneumo_data/data/results/bowtie/{sample}_primers_neg.bam"
    conda:
        "envs/bowtie_samtools.yaml"
    log:
        "s_pneumo_data/data/results/logs/bowtie_{sample}.log"
    shell:
        """
        bowtie2 -x s_pneumo_data/data/results/ref_index/{wildcards.sample}_indexed_ref -f {input.pos} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aln_pos} >> {log} 2>&1
        bowtie2 -x s_pneumo_data/data/results/ref_index/{wildcards.sample}_indexed_ref -f {input.neg} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aln_neg} >> {log} 2>&1
        """

rule samtools_merge:
    """
    Merge BAM files for visualization 
    """
    input:
        fwd=rules.bowtie_primers.output.aln_pos, 
        rev=rules.bowtie_primers.output.aln_neg
    output:
        merged_bam="s_pneumo_data/data/results/samtools_merge/{sample}_merged.bam"
    conda:
        "envs/bowtie_samtools.yaml"
    shell:
        """
        samtools merge -o {output.merged_bam} {input.fwd} {input.rev}
        """

rule depth:
    """
    Calculate depth for visualization
    """
    input:
        primer_alns = rules.samtools_merge.output.merged_bam
    output:
        primer_depth = "s_pneumo_data/data/results/samtools_depth/{sample}.depth"
    conda:
        "envs/bowtie_samtools.yaml"
    shell:
        """
        samtools depth -a {input.primer_alns} > {output.primer_depth}
        """

rule visualize_depth:
    """
    Visualize depth for each primer alignment
    """
    input:
        py_input = rules.depth.output.primer_depth
    output:
        "s_pneumo_data/data/results/depth_vis/{sample}_depth.png"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/scripts/plot_cov_primers.py {input.py_input} {output} 
        """

rule coverage:
    """
    Get % coverage for each primer/GPSC pair
    """
    input:
        "s_pneumo_data/data/results/samtools_merge/{sample}_merged.bam"
    output:
        "s_pneumo_data/data/results/coverage_pc/{sample}_coverage.csv"
    conda:
        "envs/bowtie_samtools.yaml"
    shell:
        """
        samtools coverage {input} > {output}
        """