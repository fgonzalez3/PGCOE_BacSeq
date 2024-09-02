import pandas as pd 

configfile: "config/SP_sero.yaml"

samples_df = pd.read_csv("tsv/contigs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONTIGS = {row.sample_id: row.contig_path for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/mlst/{sample}_mlst.csv", genera=config["genera"], sample=SAMPLES),
        expand("results/{genera}/poppunk/poppunk_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/poppunk/poppunk.dists.pkl", genera=config["genera"]),
        expand("results/{genera}/poppunk/poppunk.h5", genera=config["genera"]),
        expand("results/{genera}/poppunk/poppunk.dists.npy", genera=config["genera"]),
        expand("results/{genera}/poppunk/poppunk_external_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/poppunk/poppunk_unword_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/microreact/microreact_microreact_clusters.csv", genera=config["genera"]),
        expand("results/{genera}/microreact/microreact_core_NJ.nwk", genera=config["genera"]),
        expand("results/{genera}/microreact/microreact_perplexity50.0_accessory_mandrake.dot", genera=config["genera"])

rule mlst:
    """
    In-silico MLST assignment
    """
    input:
        "contigs/{sample}/contigs.fa"
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
        queryseqs = "qfile_contigs.txt"
    output:
        "results/{genera}/poppunk/poppunk_clusters.csv",
        "results/{genera}/poppunk/poppunk.dists.pkl",
        "results/{genera}/poppunk/poppunk.h5",
        "results/{genera}/poppunk/poppunk.dists.npy",
        "results/{genera}/poppunk/poppunk_external_clusters.csv",
        "results/{genera}/poppunk/poppunk_unword_clusters.csv"
    params:
        genera=config["genera"]
    shell:
        """
        source activate /home/flg9/.conda/envs/poppunk
        poppunk_assign --db GPS_v8_ref --external-clustering GPS_v8_external_clusters.csv --output results/{params.genera}/poppunk --query {input.queryseqs} --threads 8 --update-db
        """

rule microreact:
    """
    Visualize GPSCs with Microreact
    """
    input: 
        "results/{genera}/poppunk/poppunk_clusters.csv"
    output:
        "results/{genera}/microreact/microreact_microreact_clusters.csv",
        "results/{genera}/microreact/microreact_core_NJ.nwk",
        "results/{genera}/microreact/microreact_perplexity50.0_accessory_mandrake.dot"
    params:
        genera=config["genera"]
    shell:
        """
        source activate /home/flg9/.conda/envs/poppunk
        poppunk_visualise --ref-db GPS_v8_ref --query-db results/{params.genera}/poppunk --output results/{params.genera}/microreact --previous-clustering {input} --perplexity 50 --microreact --threads 8
        """