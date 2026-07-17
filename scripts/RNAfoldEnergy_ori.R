rm(list=ls())
gc()
library(ggplot2)
options(bedtools.path = "/opt/homebrew/bin/")
library(bedtoolsr)

setwd('/Users/pollet/home1/PolyploidiX/Li_Zhen/Retrozyme_project/Retrozyme_bioinfo_Zhen/git_repository/retrozyme/retrozyme_data')

Retrozyme_rnafold=function(species){
  retrozyme_tblfile=paste('Retrozymes_detection/', species, '/retrozyme.tbl', sep = '')
  if(!file.exists(retrozyme_tblfile)){
    return(c())
  }
  retrozyme_tbl=read.csv2(retrozyme_tblfile,  sep = '\t', header = F, stringsAsFactors = F)
  colnames(retrozyme_tbl)=c('Chrm', 'start', 'stop','rpt', 'consensus_size', 'strand', 'type', 'rtzname', 'motifs', 'active_count', 'active_rpt', 'cons', 'hhr', 'family', 'appearence')
  #retrozyme_tbl=retrozyme_tbl[which(retrozyme_tbl$appearence>=2), ]
  rtzname_list=retrozyme_tbl$rtzname
  
  rnafold_df=as.data.frame(do.call('rbind' , lapply(rtzname_list, function(x){
    rnafold_file=paste('Retrozymes_detection/', species, '/HMM_cluster/', x, '.hmm.clust', sep = '')
    rnafold_df=read.csv2(rnafold_file, sep = '\t', header = F, stringsAsFactors = F)
    rnafold_values=as.numeric(rnafold_df$V4)
    return(c(x, mean(rnafold_values)))
  })), stringsAsFactors = F)
  
  print(species)
  
  rnafold_df$V2=as.numeric(rnafold_df$V2)
  colnames(rnafold_df)=c('rtzname', 'rnafold')
  family_count_df=as.data.frame(table(retrozyme_tbl$family), stringsAsFactors = F)
  colnames(family_count_df)=c('family', 'family_count')
  retrozyme_tbl=merge(retrozyme_tbl, rnafold_df)
  retrozyme_tbl=merge(retrozyme_tbl, family_count_df)
  retrozyme_tbl$species=species
  repeat_summary_file=paste('Retrozymes_detection/', species, '/repeat_summary.tbl', sep = '')
  
  write.table(file = repeat_summary_file, x = retrozyme_tbl, col.names = T, 
              sep = '\t', quote = F, row.names = F)
  return(retrozyme_tbl[, c(1, 2, 3, 16, 17, 18)])
}

All_rnafold_df=as.data.frame(do.call('rbind', lapply(list.files('Retrozymes_detection/'), Retrozyme_rnafold)), 
                             stringsAsFactors = F)

All_rnafold_df$species=sapply(All_rnafold_df$species, function(x){gsub('_', ' ', x)})
species_groupinfo_df=read.csv2('species_classname.txt', stringsAsFactors = F, header = F, sep = '\t')
colnames(species_groupinfo_df)=c('species', 'clades')
head(species_groupinfo_df)

All_rnafold_df=merge(All_rnafold_df, species_groupinfo_df)

head(All_rnafold_df)
unique(All_rnafold_df$clades)

All_rnafold_df$clades=factor(All_rnafold_df$clades, levels = c(" Aves", " Reptilia"," Mammalia", " Amphibia", 
                                                               " Sarcopterygii", " Actinopterygii", " Chondrichthyes", 
                                                               " Insecta", " Anthozoa"))

ggplot(All_rnafold_df, aes(y=rnafold, fill=clades))+
  geom_histogram(binwidth = 0.5,color='black')+
  facet_grid(.~clades, scales = 'free_x')+
  theme(legend.position = 'none')+
  scale_x_continuous(n.breaks = 3)+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  scale_fill_manual(values=rev(c("#2ca02c", "#f7b6d2","#9467bd" ,"#1f77b4", "#8c564b", 
                                 "#d62728", "#e377c2", "#7f7f7f" ,"#ff7f0e")))
ggsave('rnafold_distribution.png', width = 12, height = 6)

ggplot(All_rnafold_df, aes(x=rnafold))+
  geom_histogram(binwidth = 0.5,color='black', fill='white')

ggsave('rnafold_distribution_sum.png', width = 6, height = 6)

head(All_rnafold_df)

cor.test(All_rnafold_df$rnafold, All_rnafold_df$family_count)

?cor.test
ggplot(All_rnafold_df, aes(x=rnafold, y=family_count))+
  geom_jitter()+
  stat_smooth(method = 'lm')+
  scale_y_log10()

ggsave('rnafold_appearence.png', width = 6, height = 6)


