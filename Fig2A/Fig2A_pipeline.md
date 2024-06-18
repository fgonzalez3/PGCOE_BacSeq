# S. pneumoniae pangenome pipeline

Input data Needed:
1. Serotype 3 reference sequence used to create amplicon scheme
   a. [NC_017592](https://www.ncbi.nlm.nih.gov/nuccore/NC_017592.1) 
3. S. mitis outgroup for phylogenetic tree
   a. [JYGP01](https://www.ncbi.nlm.nih.gov/nuccore/NZ_JYGP00000000.1)
5. Query sequences (n = 35) representing majority diversity (62%) on GPS database
   a. These were pulled from [Gladestone et al. 2019](https://www.thelancet.com/article/S2352-3964(19)30259-2/fulltext#%20). 

Pangenome pipeline:
1. All genomes collected for this analysis were complete, so the first step was to run [Prokka](https://github.com/tseemann/prokka) whose .gff files would serve as input for [Roary](https://github.com/sanger-pathogens/Roary).
    a. Scripts and needed files can be found here. 
