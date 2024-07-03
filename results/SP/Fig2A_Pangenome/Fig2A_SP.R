# remove any remaining environments
rm(list = ls())

library(ggplot2)
library(ggtree)
library(tidytree)
library(ape)
library(phylotools)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(dplyr)

# define path to pneumo tree
homewd= '/Users/flg9/Desktop/Developer/grubaugh_lab/snakemake_workflows/pangenomics_test/'
setwd(paste0(homewd))

# load in tree
sp_tree <-  read.tree(file = paste0(homewd, "GPSC_pneumo_tree.newick"))

# root with outgroup
rooted_sp_tree <- midpoint(sp_tree)

# rename tree tips 
rooted_sp_tree$tip.label
rooted_sp_tree$tip.label <- gsub("_reference", "", rooted_sp_tree$tip.label)
rooted_sp_tree$tip.label

# Now, rooted_sp_tree$tip.label will have "_reference" removed from the names
print(rooted_sp_tree$tip.label)

# take a quick look at the tree
rooted_sp_ggtree <- ggtree(rooted_sp_tree) +
  geom_nodelab(aes(label=""), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='black')

rooted_sp_ggtree

# load in matrix data
roary <- read.csv("gene_presence_absence_copy.csv", header = T)

# transform gene presence absence data into a matrix 
roary[roary != ""] <- 1
roary[is.na(roary)] <- 0
roary <- as.data.frame(sapply(roary, as.numeric))

# sort the matrix by the sum of strains presence 
roary_sorted <- roary[order(rowSums(roary), decreasing = TRUE), ]

# sort roary matrix data by phylogenetic tips 
#tip_labels <- rooted_sp_tree$tip.label
#roary_sorted <- roary_sorted[, tip_labels] # extracting this directly from ggtree plot is tricky, so hardcoding this

# hard code order of matrix 
order <- c("GPSC3", "GPSC37", "GPSC22", "GPSC17", "GPSC4", "GPSC32", 
           "GPSC8", "GPSC25", "GPSC34", "GPSC26", "NC_017592", "GPSC15", "GPSC5", 
           "JYGP01")

# rename colnames from roary df and sort
colnames(roary_sorted) <- gsub("_reference", "", colnames(roary_sorted))
roary_sorted <- roary_sorted[, order]

# assign color palette
color_palette <- colorRampPalette(c("white", "#93a3b1"))(2)  # adjust for binary data

# define breaks 
breaks <- c(-0.5, 0.5, 1.5)

# generate heatmap with the above breaks
matrix_plot <- pheatmap(t(roary_sorted),
         color = color_palette,
         breaks = breaks,  # Use the custom breaks
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         legend = F, 
         show_rownames = F)

# combine the tree and matrix to see how it looks
combined_plot <- plot_grid(rooted_sp_ggtree, matrix_plot$gtable, ncol = 2, align = 'h')

combined_plot

# now read in the fastani output
fastani <- read.delim('fastani.out', header = F)

# change up the structure of this file 
fastani$V1 <- sub(".*(GPSC\\d+)_reference\\.fasta", "\\1", fastani$V1)
fastani$V1 <- sub(".*(JYGP01)\\.fasta", "\\1", fastani$V1)
fastani$V1 <- sub(".*(NC_017592)\\.fasta", "\\1", fastani$V1)

fastani$V2 <- sub(".*(NC_017592)\\.fasta", "\\1", fastani$V2)

colnames(fastani) <- c("GPSC", "Reference", "ANI", "Seq_Fragments", "Ortho_Matches")

# now plot ANI across gene matrix 
head(fastani)

# normalise
fastani$Normalized_ANI <- (fastani$ANI - min(fastani$ANI)) / (max(fastani$ANI) - min(fastani$ANI))

# create color pallete for ANI vals
ani_colors <- colorRampPalette(c("blue", "red"))(100)

# map normalized values to colors
fastani$ANI_Color <- ani_colors[as.integer(fastani$Normalized_ANI * 99) + 1]

# create a data frame for row annotations
row_annotations <- data.frame(ANI_Color = fastani$ANI_Color)
rownames(row_annotations) <- fastani$GPSC 

