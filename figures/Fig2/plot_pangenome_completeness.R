rm(list = ls())

library(ggtree)
library(ape)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(phytools)
library(patchwork)

SP_treefile       <- "./SP_tree.newick"
TB_treefile       <- "./TB_tree.newick"
SP_amps           <- "./SP_amp_positions.csv"
TB_amps           <- "./TB_amp_positions.csv"
SP_assemblyfile   <- "./SP_assembly_lengths.csv"
TB_assemblyfile   <- "./TB_assembly_lengths.csv"
SP_fastani        <- "./SP_fastani.out"
TB_fastani        <- "./TB_fastani.out"
SP_pangenomefile  <- "./SP_gene_presence_absence.csv"
TB_pangenomefile  <- "./TB_gene_presence_absence.csv"
SP_outfile  <- "./SP_combined_plot.png"
TB_outfile  <- "./TB_combined_plot.png"


#SP_outgroup="JYGP01"
#TB_outgroup= "NC_019950"


# tree plot ---------------------------------------------------------------

############# SP##################
# load in tree
SP_tree <-  read.tree(SP_treefile)
SP_tree$tip.label <- gsub("_reference", "", SP_tree$tip.label)

# root with midpoint
#tree <- ladderize(root(tree,outgroup))
SP_tree <- midpoint_root(SP_tree)

treesize = max(SP_tree$edge.length)

treeplot <- ggtree(SP_tree,ladderize=F) +
  #geom_nodelab(aes(label=""), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align=T, offset = treesize*0.01, linetype="dotted", 
              linesize = 0.4,size=3) + 
  scale_x_continuous(expand=c(0,treesize*0.35))
treeplot

#get tree order 
istip <- SP_tree$edge[,2]<=length(SP_tree$tip.label)
treeorder <- SP_tree$tip.label[SP_tree$edge[istip,2]]

# pangenome representation ------------------------------------------------

# load in matrix data
pangenome <- read.csv(pangenomefile, header = T)

pangenome <- pangenome %>% 
  select(!c(2,3,6,7:14)) %>% 
  rename_with(function(x) {gsub("_reference","",x)}) %>% 
  rename_with(function(x) {gsub("\\.\\.","_",x)}) %>% 
  rename_with(tolower,c(1:3)) %>% 
  mutate(across(!c(1:3), function(x) {x != ""})) %>%
  mutate(gindex = seq.int(nrow(pangenome))) %>% 
  relocate(gindex) %>% 
  pivot_longer(!c(1:4),names_to = "strain",values_to="present") %>% 
  mutate(strain = factor(strain,ordered=T,levels=treeorder))


# read in divergence / fastani output
mash <- read.table(mash, sep="\t",header = F,
                       col.names = c("ref","strain","mashdist","pval","matching-hashes")) %>%
  mutate(across(c(ref,strain),function(x) {gsub(".fasta","",gsub(".*\\/","",x,perl=T))})) %>%
  mutate(across(c(ref,strain),function(x) {gsub("_reference","",x)}))

fastani <- read.table(fastani, sep="\t",header = F,
                                  col.names = c("strain","ref","identity","fragments","matches")) %>%
  mutate(across(c(strain,ref),function(x) {gsub(".fasta","",gsub(".*\\/","",x,perl=T))})) %>%
  mutate(across(c(strain,strain),function(x) {gsub("_reference","",x)}))


# plot pangenome / divergence heatmap 
pangenome_fastani <- merge(pangenome, fastani)

panplot <- ggplot(subset(pangenome_mash, present), aes(x = gindex, y = strain, fill = mashdist)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "darkred", mid = "lightgrey", midpoint = 0.05, limits = c(0.00, 0.1), name = "seq identity") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("pangenome representation")
panplot

panplot2 <- ggplot(subset(pangenome_fastani,present),aes(x=gindex,y=strain,fill=identity)) + 
  geom_tile() + 
  scale_fill_gradient(high="darkred", low="white",limits=c(80,100),name="seq identity") + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  theme(legend.position="bottom",
        panel.background=element_blank(),
        panel.grid = element_blank(),
        axis.title=element_blank(),axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("pangenome representation")
panplot2


# read amplicon mappings --------------------------------------------------
assembly_lengths <- read.csv("SP_assmebly_lengths.csv", header=FALSE, col.names=c("strain", "start", "end")) %>%
  mutate(strain = factor(strain, ordered=TRUE, levels=treeorder)) %>%
  mutate(y = as.numeric(strain) - 0.5)

ampmap <- read.csv(fwd_or_rev_amplicons, col.names = c("strain", "start", "end", "panel")) %>%
  mutate(strain = factor(strain, ordered=TRUE, levels=treeorder),
         panel_numeric = as.numeric(factor(panel, levels = c("fwd", "rev")))) %>%
  mutate(y = as.numeric(strain) - 1 + (panel_numeric - 1) * 0.5)

ampplot <- ggplot(ampmap, aes(xmin=start, xmax=end, ymin=y+0.05, ymax=y+0.45)) + 
  geom_rect() + 
  geom_segment(data=assembly_lengths, aes(x=start, xend=end, y=y, yend=y)) + 
  scale_y_continuous(breaks=seq(0.5, length(treeorder), 1), labels=treeorder,
                     limits=c(0, length(treeorder)), expand=c(0, 0)) + 
  xlab("mb") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("amplicon coverage")

ampplot


# assemble combined plots -------------------------------------------------
treeplot + panplot + ampplot + plot_layout(widths=c(2,3,3))
ggsave(outfile, width=350,height=150,dpi=400,units="mm")


treeplot + panplot + panplot2 + plot_layout(widths=c(2,3,3))
ggsave(outfile, width=350,height=150,dpi=400,units="mm")



