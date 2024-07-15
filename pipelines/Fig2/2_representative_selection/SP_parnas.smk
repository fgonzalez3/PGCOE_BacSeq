configfile: "config/SP_parnas.yaml"
RADIUSES = set(config['radius'])

rule all:
    input:
        "results/SP/parnas/estimated_number_representatives.csv",
        "results/SP/parnas/pneumo_representative_strains.tre", 
        "results/SP/parnas/Estimated_SPneumo_Reps_Needed.jpg", 
        "results/SP/parnas/Pneumo_subtree.tre",
        expand("results/SP/parnas/{rad}/colors_tree.csv", rad = RADIUSES),
        expand("results/SP/parnas/{rad}/{rad}_pneumo_cluster_reps.tab", rad = RADIUSES), 
        expand("results/SP/parnas/{rad}/colors_formatted.csv", rad = RADIUSES),
        expand("results/SP/parnas/{rad}/{rad}_Pneumo_subtree.tre", rad = RADIUSES)

rule estimate_representative_number:
    """
    This estimates how many representatives is sufficient to maximally cover diversity
    """
    input:
        tree = config["tree"]
    output:
        scores = "results/SP/parnas/estimated_number_representatives.csv",
        representatives = "results/SP/parnas/representatives.txt",
        colortree = "results/SP/parnas/pneumo_representative_strains.tre", 
        subtree1 = "results/SP/parnas/Pneumo_subtree.tre" 
    params: 
        max_num = config["max_num"]
    shell:
        """
        parnas -t {input.tree} -n {params.max_num} \
        --diversity {output.scores} > {output.representatives} \
        --color {output.colortree} \
        --subtree {output.subtree1} 
        """

rule create_diversity_figure:
    """
    This creates a figure using the diversity covered data in the prior step
    """
    input:
        scores = rules.estimate_representative_number.output.scores
    output:
        figure = "results/SP/parnas/Estimated_SPneumo_Reps_Needed.jpg"
    shell:
        r"""
        Rscript "scripts/estimate_reps_needed.R" --input {input.scores} --output {output.figure}
        """

rule colortrees_with_radius:
    """
    Find a minimum number of representatives that cover most of the diversity across our tree
    """
    input:
        tree = config["tree"] 
    params:
        radius = "{rad}"
    output:
        colors = "results/SP/parnas/{rad}/colors_tree.csv",
        reps = "results/SP/parnas/{rad}/{rad}_pneumo_cluster_reps.tab", 
        subtree2 = "results/SP/parnas/{rad}/{rad}_Pneumo_subtree.tre"
    shell:
        """
        parnas -t {input.tree} --cover --radius {params.radius} \
        --subtree {output.subtree2} \
        --clusters {output.colors} > {output.reps} \
        """

rule combine_colortrees:
    """
    Combine reps with numerous radiuses from prior rule 
    """
    input:
        colors = rules.colortrees_with_radius.output.colors, 
        reps = rules.colortrees_with_radius.output.reps
    output:
        formatted = "results/SP/parnas/{rad}/colors_formatted.csv"
    shell:
        r"""
        Rscript scripts/colorcombiner.R \
        --representatives-file {input.reps} \
        --colors-file {input.colors} \
        --output {output}
        """