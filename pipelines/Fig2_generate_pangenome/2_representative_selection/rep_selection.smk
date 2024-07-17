configfile: "config/SP_parnas.yaml"
RADIUSES = set(config['radius'])

rule all:
    input:
        expand("results/{target}/parnas/representatives.txt", target=config["target"]), 
        expand("results/{target}/parnas/{rad}/colors_tree.csv", rad = RADIUSES, target=config["target"]),
        expand("results/{target}/parnas/{rad}/{rad}_cluster_reps.tab", rad = RADIUSES, target=config["target"]), 
        expand("results/{target}/parnas/{rad}/colors_formatted.csv", rad = RADIUSES, target=config["target"]),
        expand("results/{target}/parnas/{rad}/{rad}_subtree.tre", rad = RADIUSES, target=config["target"])

rule estimate_representative_number:
    """
    This estimates how many representatives is sufficient to maximally cover diversity
    """
    input:
        tree = config["tree"]
    output:
        scores = "results/{target}/parnas/estimated_number_representatives.csv",
        representatives = "results/{target}/parnas/representatives.txt",
        colortree = "results/{target}/parnas/representative_strains.tre", 
        subtree1 = "results/{target}/parnas/subtree.tre" 
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
        figure = "results/{target}/parnas/Estimated_Reps_Needed.jpg"
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
        colors = "results/{target}/parnas/{rad}/colors_tree.csv",
        reps = "results/{target}/parnas/{rad}/{rad}_cluster_reps.tab", 
        subtree2 = "results/{target}/parnas/{rad}/{rad}_subtree.tre"
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
        formatted = "results/{target}/parnas/{rad}/colors_formatted.csv"
    shell:
        r"""
        Rscript scripts/colorcombiner.R \
        --representatives-file {input.reps} \
        --colors-file {input.colors} \
        --output {output}
        """