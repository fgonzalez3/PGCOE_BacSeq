

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(argparse)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(viridis)
library(sysfonts)
library(forcats)

coverageinfo_tb <- read_delim('/Users/chaneykalinich/Documents/PGCoE/github/combinedcoverage.tsv',delim='\t')
datasheet_tb <- read_delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pgcoe_pipeline/data/PGCOE_Datasheet - M.tuberculosis_Seq.csv',delim=',')%>%
  select("Seq_ID","Original_ID","Sample_Type","Sample_source","collection_date","CT","starting_quant","Template_dilution","Primer_Conc","NGS_Prep_Method",
         "amp_date","NGS_primers","NGS_Run_ID","NGS_Date")
coverageinfo_spn_amp <- read_delim('/Users/chaneykalinich/Documents/PGCoE/SP_alignment_13JUN2024/cov_stats/amplicon_seq_13JUN2024.csv',delim=',') %>%
  select("sample","subsample","X.rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq",
         "Original_ID","Sample_Type","Sample_Dilution","PCR_Primer_Scheme","Primer_Conc","NGS_Prep_Method","NGS_Run_Date",
         "Ct1","Ct2")%>%
  rename("rname"=X.rname)
coverageinfo_spn_mngs <- read_delim('/Users/chaneykalinich/Documents/PGCoE/SP_alignment_13JUN2024/cov_stats/mNGS_seq_13JUN2024.csv',delim=',')%>%
  select("sample","subsample","X.rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq",
         "Original_ID","Sample_Type","Sample_Dilution","PCR_Primer_Scheme","Primer_Conc","NGS_Prep_Method","NGS_Run_Date",
         "Ct1","Ct2")%>%
  rename("rname"=X.rname)
datasheet_spn <- read_delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pgcoe_pipeline/data/PGCOE_Datasheet - S.pneumo_Seq.csv',delim=',')%>%
  select("Seq_ID","Original_ID","Sample_Type","Template_Dilution","PCR_Primer-Scheme","Primer_Conc","NGS_Prep_Method","NGS_Run_ID",
         "PCR_Primer-Scheme","NGS_Run_Date")
qpcr_data <- read_delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pgcoe_pipeline/data/PGCOE_Datasheet - qPCR.csv',delim=',')%>%
  mutate(FAM_Ct=ifelse(is.na(FAM_Ct),40,FAM_Ct)) %>%
  select(!Template_Dilution)

data_combined_tb <- coverageinfo_tb %>%
  left_join(datasheet_tb, by=c('sample'='Seq_ID')) %>%
  mutate(pathogen="TB")
data_combined_spn <- coverageinfo_spn_amp %>%
  full_join(coverageinfo_spn_mngs) %>%
  mutate(sample=str_c("Yale-",sample)) %>%
  mutate(NGS_Run_Date=as.Date(NGS_Run_Date,"%Y-%m-%d"))%>%
  mutate(NGS_Prep_Method=as.factor(NGS_Prep_Method))%>%
  left_join(datasheet_spn,by=c('sample'='Seq_ID')) %>%
  select("sample","subsample","rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq",
         "Original_ID.y","Sample_Type.y","Template_Dilution","PCR_Primer_Scheme","Primer_Conc.y","NGS_Prep_Method.y","NGS_Run_Date.y",
         "Ct1","Ct2","NGS_Run_ID")%>%
  #left_join(datasheet_spn,by=c('sample'='Seq_ID',"Sample_Type","Primer_Conc","NGS_Prep_Method","NGS_Run_Date")) #%>%
  #rename(c("NGS_primers"='PCR_Primer_Scheme',"NGS_Date"='NGS_Run_Date',"Template_dilution"="Sample_Dilution")) %>%
  left_join(qpcr_data,by=c("Original_ID.y"="Sample_ID")) %>% #keep value from qPCR spreadsheet if there is one, otherwise check for previous one
  mutate(CT_touse = ifelse(!is.na(FAM_Ct),FAM_Ct,Ct1)) %>%
  rename(c("Original_ID"="Original_ID.y","Sample_Type"="Sample_Type.y","Template_dilution"="Template_Dilution","NGS_primers"="PCR_Primer_Scheme",
           "CT"="CT_touse","Primer_Conc"="Primer_Conc.y","NGS_Prep_Method"="NGS_Prep_Method.y","NGS_Date"="NGS_Run_Date.y"))%>%
  mutate(pathogen="SPn")

runs_touse <- c("CS006","CS009","CS015","CS018","CS019","TB001","TB003")

data_combined <- data_combined_tb %>%
  full_join(data_combined_spn,by=c("sample","subsample","rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq",
                                   "Original_ID","Sample_Type","CT","Template_dilution","Primer_Conc","NGS_Prep_Method","NGS_primers","NGS_Run_ID","NGS_Date",
                                   "pathogen")) %>%
  select("sample","subsample","rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq",
         "Original_ID","Sample_Type","CT","Template_dilution","Primer_Conc","NGS_Prep_Method","NGS_primers","NGS_Run_ID","NGS_Date",
         "pathogen","Sample_source","collection_date","starting_quant","amp_date") %>%
  mutate(NGS_primers=as.factor(NGS_primers)) %>%
  mutate(Primer_Conc=as.factor(Primer_Conc)) %>%
  filter(Primer_Conc %in% c("0","100uM","200uM")) %>%
  mutate(subsample=as.factor(subsample)) %>%
  filter(NGS_Run_ID %in% runs_touse)

#statsheaders <- c('sample','subsample','rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')

coverageplot <- data_combined %>%
  mutate(Template_dilution=replace_na(Template_dilution,1)) %>%
  #mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  filter(Sample_Type!='Water'&subsample==1) %>%
  ggplot() +
  geom_point(aes(x=Template_dilution,y=coverage,color=Original_ID, shape=Primer_Conc)) +
  scale_x_log10() +
  facet_grid(.~pathogen) +
  #theme_minimal_hgrid()+
  #theme_stata() +
  theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0)) +
  labs(x="Fold dilution from template DNA",
       y="Coverage of Reference (%)") +
  scale_color_discrete(name="Sample")
plot(coverageplot)

#fix facet labels
amp.labs<-c("With amplification","Without amplification")
names(amp.labs) <- c("TBv1.1","None")
colors <- c("#532989","#7F9DD5")
coverageplot_ct <-data_combined %>%
  #mutate(Template_dilution=replace_na(Template_dilution,1)) %>%
  filter(Sample_Type!='Water') %>%
  filter(subsample==1) %>%
  #filter(subsample %in% c(1,0.5,0.25,0.05)) %>%
  #mutate(NGS_primers=replace_na(NGS_primers,'None')) %>%
  #mutate(NGS_primers=factor(NGS_primers,levels=c("TBv1.1","None"))) %>%
  mutate(primersused=fct_collapse(Primer_Conc,
               "Without Amplification"="0",
               other_level = "With Amplification"))%>%
  #mutate(ct2=replace_na(CT,15))
  #mutate(logdil=log10(Template_dilution)+15) %>%
  #mutate(ct2=log10(Template_dilution)+15) %>%
  ggplot() +
  geom_hline(aes(yintercept=80)) +
  geom_smooth(aes(x=CT,y=coverage,group=primersused,color=primersused),method='loess') +
  geom_point(aes(x=CT,y=coverage,color=primersused),size=5) +
  #scale_x_log10() +
  facet_grid(rows=vars(pathogen))  +
  #facet_grid(rows=vars(NGS_primers),labeller=labeller(NGS_primers=amp.labs)) +
  labs(x='Cycle threshold',y="Genome \ncoverage (%)",color="") +
  theme_cowplot(font_size=20,font_family = 'sans',rel_small=15/20) +
  theme(legend.position=c(.1,.25)) +
  #scale_color_viridis(discrete=TRUE,option="C",labels=c("With Amplification","Without Amplification"))
  scale_color_manual(values=colors)
  #geom_hline(aes(yintercept=80,linetype="dashed")) +
  #geom_smooth(aes(x=CT,y=coverage),method='loess')
plot(coverageplot_ct)
ggsave2('/Users/chaneykalinich/Documents/PGCoE/ctvcov_tb003.tiff',coverageplot_ct,width=9.32,height=3.07,units="in")

coverageplot_sq <- data_combined %>%
  #mutate(Template_dilution=replace_na(Template_dilution,1)) %>%
  filter(Original_ID!='NTC') %>%
  filter(subsample==1.0) %>%
  #filter(subsample %in% c(1,0.5,0.25,0.05)) %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None')) %>%
  mutate(NGS_primers=factor(NGS_primers,levels=c("TBv1.1","None"))) %>%
  ggplot() +
  geom_hline(aes(yintercept=80)) +
  geom_smooth(aes(x=starting_quant,y=coverage,group=NGS_primers,color=NGS_primers),method='loess') +
  geom_point(aes(x=starting_quant,y=coverage,color=NGS_primers),size=5) +
  scale_x_log10(n.breaks=8) +
  scale_y_continuous(limits=c(0,100))+
  #facet_grid(cols=vars(subsample), rows=vars(NGS_primers))  +
  #facet_grid(rows=vars(NGS_primers),labeller=labeller(NGS_primers=amp.labs)) +
  labs(x='Starting quantity',y="Genome \ncoverage (%)",color="") +
  theme_cowplot(font_size=20,font_family = 'sans',rel_small=15/20) +
  background_grid()+
  theme(legend.position=c(.7,.25)) +
  #scale_color_viridis(discrete=TRUE,option="C",labels=c("With Amplification","Without Amplification"))
  scale_color_manual(values=colors,labels=c("With Amplification","Without Amplification"))
plot(coverageplot_sq)


foldchange <- data_combined %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  pivot_wider(names_from="NGS_primers",values_from='coverage',id_cols=c(Original_ID,CT,subsample,Template_dilution)) %>%
  filter(Original_ID %in% c('111-10712-18','111-5264-18')) %>%
  ggplot() +
  geom_linerange(aes(ymin=None,ymax=TBv1.1,x=CT)) +
  geom_point(aes(x=CT,y=None,color=Original_ID)) +
  geom_point(aes(x=CT,y=TBv1.1,color=Original_ID)) +
  facet_wrap('subsample') +
  labs(x='Ct',y="Coverage of H37Rv Reference (%)") 
plot(foldchange)


foldchange_reads<- data_combined %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  pivot_wider(names_from="NGS_primers",values_from='numreads',id_cols=c(Original_ID,CT,subsample,Template_dilution)) %>%
  filter(Original_ID %in% c('111-10712-18','111-5264-18')) %>%
  mutate(changereads=TBv1.1/None) %>%
  ggplot() +
  geom_linerange(aes(ymin=None,ymax=TBv1.1,x=CT)) +
  geom_point(aes(x=CT,y=None,color=Original_ID)) +
  geom_point(aes(x=CT,y=TBv1.1,color=Original_ID)) +
  facet_wrap('subsample') +
  labs(x='Ct',y="Number of reads aligned to H37Rv reference") +
  scale_y_log10()
plot(foldchange_reads)

datacheck <- data_combined %>%
  filter(subsample==1.0) %>%
  select(sample, numreads, covbases,coverage, meandepth,meanbaseq,Original_ID,starting_quant,Primer_Conc)
