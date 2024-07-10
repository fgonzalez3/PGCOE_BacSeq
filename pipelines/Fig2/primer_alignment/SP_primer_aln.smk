configfile: "config/SP_primer_aln.yaml"

rule all:
    input:
        "results/SP/primer_aln/primer_fa/fwd_primers.bed", 
        "results/SP/primer_aln/primer_fa/rev_primers.bed", 
        "results/SP/primer_aln/primer_fa/fwd.scheme.primer.fasta",
        "results/SP/primer_aln/primer_fa/rev.scheme.primer.fasta",
        expand("results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.1.bt2", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.2.bt2", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.3.bt2", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.4.bt2", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.1.bt2", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.2.bt2", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_align/{sample}_primers_fwd.bam", sample=config["samples"]),
        expand("results/SP/primer_aln/bowtie_align/{sample}_primers_rev.bam", sample=config["samples"]),
        expand("results/SP/primer_aln/samtools_depth_indiv_primers/{sample}_fwd.depth", sample=config["samples"]), 
        expand("results/SP/primer_aln/samtools_depth_indiv_primers/{sample}_rev.depth", sample=config["samples"]), 
        expand("results/SP/primer_aln/samtools_merge/{sample}_merged.bam", sample=config["samples"]),
        expand("results/SP/primer_aln/samtools_depth/{sample}.depth", sample=config["samples"]), 
        #expand("results/SP/primer_aln/depth_vis/{sample}_depth.png", sample=config["samples"]), 
        expand("results/SP/primer_aln/coverage_pc/{sample}_coverage.csv", sample=config["samples"])

rule split_bed:
    """
    Use awk to extract fwd and rev primers 
    """
    input:
        "SP.scheme.primer.bed"
    output:
        fwd_pr="results/SP/primer_aln/primer_fa/fwd_primers.bed", 
        rev_pr="results/SP/primer_aln/primer_fa/rev_primers.bed"
    shell:
        """
        awk '$6 == "+" {{print}}' {input} > {output.fwd_pr} 
        awk '$6 == "-" {{print}}' {input} > {output.rev_pr} 
        """

rule get_fasta:
    """
    Convert BED to FASTA for downstream alignment 
    """
    input:
        fwd_bed = "results/SP/primer_aln/primer_fa/fwd_primers.bed", 
        rev_bed = "results/SP/primer_aln/primer_fa/rev_primers.bed"
    output:
        bedfasta_fwd = "results/SP/primer_aln/primer_fa/fwd.scheme.primer.fasta", 
        bedfasta_rev = "results/SP/primer_aln/primer_fa/rev.scheme.primer.fasta"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/primer_alignment/scripts/get_indiv_fasta.py {input.fwd_bed} {input.rev_bed} {output.bedfasta_fwd} {output.bedfasta_rev}
        """

rule bowtie_build:
    """
    Create index of GPSCs
    """
    input:
        queryseqs = lambda wildcards: config['samples'][wildcards.sample]
    output:
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.1.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.2.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.3.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.4.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.1.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.2.bt2"
    conda:
        "envs/read_aln.yaml"
    log:
        "results/SP/primer_aln/logs/bowtie_build/{sample}_bowtie_build.log"
    shell:
        """
        bowtie2-build -f {input.queryseqs} results/SP/primer_aln/bowtie_index/{wildcards.sample}/{wildcards.sample}_indexed_ref > {log}
        """

rule bowtie_align:
    """
    Align primers to GPSCs
    """
    input:
        fwd = rules.get_fasta.output.bedfasta_fwd, 
        rev = rules.get_fasta.output.bedfasta_rev,
        idx = [
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.1.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.2.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.3.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.4.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.1.bt2",
        "results/SP/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.2.bt2" 
        ]
    output:
        aln_fwd = "results/SP/primer_aln/bowtie_align/{sample}_primers_fwd.bam",
        aln_rev = "results/SP/primer_aln/bowtie_align/{sample}_primers_rev.bam"
    conda:
        "envs/read_aln.yaml"
    log:
        "results/SP/primer_aln/logs/bowtie_align/bowtie_{sample}.log"
    shell:
        """
        bowtie2 -x results/SP/primer_aln/bowtie_index/{wildcards.sample}/{wildcards.sample}_indexed_ref -f {input.fwd} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aln_fwd} 2> {log}
        bowtie2 -x results/SP/primer_aln/bowtie_index/{wildcards.sample}/{wildcards.sample}_indexed_ref -f {input.rev} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aln_rev} 2> {log}
        """

rule fwd_rev_depth:
    """
    Calculate depth for fwd and rev primers separately
    """
    input:
        fwd_primer_alns = rules.bowtie_align.output.aln_fwd,
        rev_primer_alns = rules.bowtie_align.output.aln_rev 
    output:
        fwd_primer_depth = "results/SP/primer_aln/samtools_depth_indiv_primers/{sample}_fwd.depth",
        rev_primer_depth = "results/SP/primer_aln/samtools_depth_indiv_primers/{sample}_rev.depth"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools depth -a {input.fwd_primer_alns} > {output.fwd_primer_depth}
        samtools depth -a {input.rev_primer_alns} > {output.rev_primer_depth}
        """

rule samtools_merge:
    """
    Merge BAM files for visualization 
    """
    input:
        fwd=rules.bowtie_align.output.aln_fwd, 
        rev=rules.bowtie_align.output.aln_rev
    output:
        merged_bam="results/SP/primer_aln/samtools_merge/{sample}_merged.bam"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools merge -o {output.merged_bam} {input.fwd} {input.rev}
        """

rule depth:
    """
    Calculate depth for fwd and rev primers in single BAM input
    """
    input:
        primer_alns=rules.samtools_merge.output.merged_bam
    output:
        primer_depth="results/SP/primer_aln/samtools_depth/{sample}.depth"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools depth -a {input.primer_alns} > {output.primer_depth}
        """

#rule visualize_depth:
    #"""
    #Visualize depth for each primer alignment
    #"""
    #input:
        #py_input=rules.depth.output.primer_depth
    #output:
        #"results/SP/primer_aln/depth_vis/{sample}_depth.png"
    #conda:
        #"envs/biopython.yaml"
    #shell:
        #"""
        #python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/primer_alignment/scripts/plot_cov_primers.py {input.py_input} {output} 
        #"""

rule coverage:
    """
    Get % coverage for each primer/GPSC pair
    """
    input:
        "results/SP/primer_aln/samtools_merge/{sample}_merged.bam"
    output:
        "results/SP/primer_aln/coverage_pc/{sample}_coverage.csv"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools coverage {input} > {output}
        """
