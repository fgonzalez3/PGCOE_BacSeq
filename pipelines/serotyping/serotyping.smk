import pandas as pd 

configfile: "config/SP_sero.yaml"

samples_df = pd.read_csv("tsv/SP_isolates.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/shovill/{sample}/contigs.fa", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/contigs.gfa", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/shovill.corrections", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/shovill.log", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/shovill/{sample}/spades.fasta", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/mlst/{sample}_mlst.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/poppunk/GPSC_assignments_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/poppunk/GPSC_assignments.dists.pkl", genera=config["genera"]),
        expand("results/{genera}/poppunk/GPSC_assignments.h5", genera=config["genera"]),
        expand("results/{genera}/poppunk/GPSC_assignments.dists.npy", genera=config["genera"]),
        expand("results/{genera}/poppunk/GPSC_assignments_external_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/poppunk/GPSC_assignments_unword_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/microreact/microreact_microreact_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/microreact/microreact_core_NJ.nwk", genera=config["genera"]),
        expand("results/{genera}/microreact/microreact_perplexity50.0_accessory_mandrake.dot", genera=config["genera"])

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
    shell:
        """
        source activate /home/flg9/.conda/envs/shovill
        shovill --outdir results/{params.genera}/shovill/{wildcards.sample} --R1 {input.fwd} --R2 {input.rev} --force
        """

rule mlst:
    """
    In-silico MLST assignment
    """
    input:
        "results/{genera}/shovill/{sample}/contigs.fa"
    output:
        csv = "results/{genera}/mlst/{sample}_mlst.csv"
    conda:
        "envs/mlst.yaml"
    shell:
        """
        mlst --csv {input} > {output.csv}
        """

rule poppunk:
    """
    GPSC assignment
    """
    input:
        queryseqs = "txt/qfile_contigs.txt"
    output:
        "results/{genera}/poppunk/GPSC_assignments_clusters.csv",
        "results/{genera}/poppunk/GPSC_assignments.dists.pkl",
        "results/{genera}/poppunk/GPSC_assignments.h5",
        "results/{genera}/poppunk/GPSC_assignments.dists.npy",
        "results/{genera}/poppunk/GPSC_assignments_external_clusters.csv",
        "results/{genera}/poppunk/GPSC_assignments_unword_clusters.csv"
    params:
        genera=config["genera"]
    conda:
        "envs/poppunk.yaml"
    shell:
        """
        poppunk_assign --db GPS_v9 --external-clustering GPS_v9_external_clusters.csv \
        --output results/{params.genera}/poppunk --query {input.queryseqs} --threads 8 --update-db
        """

rule microreact:
    """
    Visualize GPSCs with Microreact
    """
    input: 
        "results/{genera}/poppunk/GPSC_assignments_clusters.csv"
    output:
        "results/{genera}/microreact/microreact_microreact_clusters.csv",
        "results/{genera}/microreact/microreact_core_NJ.nwk",
        "results/{genera}/microreact/microreact_perplexity50.0_accessory_mandrake.dot"
    params:
        genera=config["genera"]
    conda:
        "envs/poppunk.yaml"
    shell:
        """
        poppunk_visualise --ref-db GPS_v9 --query-db results/{params.genera}/poppunk --output results/{params.genera}/microreact --previous-clustering {input} --perplexity 50 --microreact --threads 8
        """