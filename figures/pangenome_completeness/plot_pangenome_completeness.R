library(ggtree)
library(ape)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
#library(tidytree)
library(phytools)
#library(pheatmap)
#library(cowplot)

treefile       <- "./SP_ml_phylogeny.newick"
ampliconfile   <- "./SP_amplicon_positions.csv"
divergencefile <- "./SP_fastani.out"
pangenomefile  <- "./SP_gene_presence_absence_copy.csv"


#outgroup="JYGP01"


# tree plot ---------------------------------------------------------------

# load in tree
tree <-  read.tree(treefile)
tree$tip.label <- gsub("_reference", "", tree$tip.label)

# root with midpoint
#tree <- ladderize(root(tree,outgroup))
tree <- midpoint_root(tree)

treeplot <- ggtree(tree,ladderize=F) +
  #geom_nodelab(aes(label=""), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align=T, linetype="dotted", linesize = 0.1, size = 2.5)
treeplot

#get tree order 
istip <- tree$edge[,2]<=length(tree$tip.label)
treeorder <- tree$tip.label[tree$edge[istip,2]]

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
divtable <- read.table(divergencefile, sep="\t",header = F,
                       col.names = c("strain","ref","identity","fragments","matches")) %>% 
                  mutate(across(c(strain,ref),function(x) {gsub(".fasta","",gsub(".*\\/","",x,perl=T))})) %>%
                  mutate(across(c(strain,ref),function(x) {gsub("_reference","",x)}))


# plot pangenome / divergence heatmap 
pangenome <- merge(pangenome,divtable)
panplot <- ggplot(subset(pangenome,present),aes(x=gindex,y=strain,fill=identity)) + 
  geom_tile() + 
  scale_fill_gradient(high="darkred",low="white",limits=c(80,100),name="seq identity") + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(legend.position="bottom",
        panel.background=element_blank(),panel.grid = element_blank(),
        axis.title=element_blank(),axis.text.x=element_blank())
panplot



# assemble combined plots -------------------------------------------------


treeplot | panplot


ampmap <- read.table("SP_amplicon_positions.csv",sep="\t")

