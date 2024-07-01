configfile: "config/aln_stats.yaml"

rule all:
    input:
        "results/SP/primer_fa/fwd_primers.bed", 
        "results/SP/primer_fa/rev_primers.bed"
        "results/SP/primer_fa/fwd.scheme.primer.fasta",
        "results/SP/primer_fa/rev.scheme.primer.fasta",
        expand("results/SP/bwa_index/{sample}_indexed_ref.0123", sample=SAMPLES),
        expand("results/SP/bwa_index/{sample}_indexed_ref.amb", sample=SAMPLES),
        expand("results/SP/bwa_index/{sample}_indexed_ref.ann", sample=SAMPLES),
        expand("results/SP/bwa_index/{sample}_indexed_ref.bwt.2bit.64", sample=SAMPLES),
        expand("results/SP/bwa_index/{sample}_indexed_ref.pac", sample=SAMPLES),
        expand("results/SP/bwa_align/{sample}_primers_fwd.bam", sample=config["samples"]),
        expand("results/SP/bwa_align/{sample}_primers_rev.bam", sample=config["samples"]),
        expand("results/SP/samtools_merge/{sample}_merged.bam", sample=config["samples"]),
        expand("results/SP/samtools_depth/{sample}.depth", sample=config["samples"]), 
        expand("results/SP/depth_vis/{sample}_depth.png", sample=config["samples"]), 
        expand("results/SP/coverage_pc/{sample}_coverage.csv", sample=config["samples"])

rule split_bed:
    """
    Use awk to extract fwd and rev primers 
    """
    input:
        "pipelines/alignment_stats/scheme.primer.bed"
    output:
        fwd_pr="results/SP/primer_fa/fwd_primers.bed", 
        rev_pr="results/SP/primer_fa/rev_primers.bed"
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
        fwd_bed = "results/SP/primer_fa/fwd_primers.bed", 
        rev_bed = "results/SP/primer_fa/rev_primers.bed"
    output:
        bedfasta_fwd = "results/SP/primer_fa/fwd.scheme.primer.fasta", 
        bedfasta_rev = "results/SP/primer_fa/rev.scheme.primer.fasta"
    conda:
        "pipelines/alignment_stats/envs/biopython.yaml"
    shell:
        """
        python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/pipelines/alignment_stats/scripts/get_indiv_fasta.py {input.pos_bed} {input.neg_bed} {output.bedfasta_pos} {output.bedfasta_neg}
        """

rule bwa_build:
    """
    Create index of GPSCs
    """
    input:
        queryseqs = lambda wildcards: config["samples"][wildcards.sample]
    output:
        "results/SP/bwa_index/{sample}_indexed_ref.0123",
        "results/SP/bwa_index/{sample}_indexed_ref.amb",
        "results/SP/bwa_index/{sample}_indexed_ref.ann",
        "results/SP/bwa_index/{sample}_indexed_ref.bwt.2bit.64",
        "results/SP/bwa_index/{sample}_indexed_ref.pac"
    conda:
        "pipelines/alignment_stats/envs/read_aln.yaml"
    log:
        "results/SP/logs/{sample}_bwa_build.log"
    shell:
        """
        mkdir -p results/SP/ref_index/
        bwa-mem2 index -p results/SP/ref_index {input.queryseqs} > {log}
        """

rule bwa_align:
    """
    Align primers to GPSCs
    """
    input:
        fwd = rules.get_fasta.output.bedfasta_fwd, 
        rev = rules.get_fasta.output.bedfasta_rev,
        idx = [
        "results/SP/bwa_index/{sample}_indexed_ref.0123",
        "results/SP/bwa_index/{sample}_indexed_ref.amb",
        "results/SP/bwa_index/{sample}_indexed_ref.ann",
        "results/SP/bwa_index/{sample}_indexed_ref.bwt.2bit.64",
        "results/SP/bwa_index/{sample}_indexed_ref.pac"  
        ]
    output:
        aln_fwd = "results/SP/bwa_align/{sample}_primers_fwd.bam",
        aln_rev = "results/SP/bwa_align/{sample}_primers_rev.bam"
    conda:
        "pipelines/alignment_stats/envs/read_aln.yaml"
    log:
        "results/SP/logs/bwa_align/bwa_{sample}.log"
    shell:
        """
        bwa-mem2 mem -t 32 -p results/SP/ref_index/{wildcards.sample}_indexed_ref {input.fwd} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aln_fwd} >> {log} 2>&1
        bwa-mem2 mem -t 32 -p results/SP/ref_index/{wildcards.sample}_indexed_ref {input.rev} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aln_rev} >> {log} 2>&1
        """

rule samtools_merge:
    """
    Merge BAM files for visualization 
    """
    input:
        fwd=rules.bwa_align.output.aln_fwd, 
        rev=rules.bwa_align.output.aln_rev
    output:
        merged_bam="results/SP/samtools_merge/{sample}_merged.bam"
    conda:
        "pipelines/alignment_stats/envs/read_aln.yaml"
    shell:
        """
        samtools merge -o {output.merged_bam} {input.fwd} {input.rev}
        """

rule depth:
    """
    Calculate depth for visualization
    """
    input:
        primer_alns=rules.samtools_merge.output.merged_bam
    output:
        primer_depth="results/SP/samtools_depth/{sample}.depth"
    conda:
        "pipelines/alignment_stats/envs/read_aln.yaml"
    shell:
        """
        samtools depth -a {input.primer_alns} > {output.primer_depth}
        """

rule visualize_depth:
    """
    Visualize depth for each primer alignment
    """
    input:
        py_input=rules.depth.output.primer_depth
    output:
        "results/SP/depth_vis/{sample}_depth.png"
    conda:
        "pipelines/alignment_stats/envs/biopython.yaml"
    shell:
        """
        python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/pipelines/alignment_stats/scripts/plot_cov_primers.py {input.py_input} {output} 
        """

rule coverage:
    """
    Get % coverage for each primer/GPSC pair
    """
    input:
        "results/SP/samtools_merge/{sample}_merged.bam"
    output:
        "results/SP/coverage_pc/{sample}_coverage.csv"
    conda:
        "pipelines/alignment_stats/envs/bowtie_samtools.yaml"
    shell:
        """
        samtools coverage {input} > {output}
        """