rm(list=ls())
gc()

library(ggplot2)
library(bedtoolsr)
options(bedtools.path = "/home/zhenli/bioinfo/bedtools2/bin/")


setwd("~/remote2/Toutatis_backup/Retrozyme/Retrozymes_detection/Xenopus_tropicalis/")

RNA_depth_level=function(bamfile){
  #bamfile='~/remote2/EST/aln_minimap.bam'
  #bamfile = "~/remote2/RNAseq/Xenopus_tropicalis/SRR19548800.bam"
  active_retrozyme_df=read.csv2('active_retrozyme.bed', stringsAsFactors = F, header = F, sep = '\t')
  colnames(active_retrozyme_df)=c('chrmid', 'start', 'stop', 'name', 'rnafold', 'strand', 
                            'family', 'rpt', 'monomer_length','unitlength')

  active_retrozyme_w50_df=bedtoolsr::bt.makewindows(b = active_retrozyme_df, w = 50)
  colnames(active_retrozyme_w50_df)=c('chrmid', 'start', 'stop')
  active_retrozyme_w50_df = bt.intersect(a=active_retrozyme_w50_df, b=active_retrozyme_df, wo=T)[,c(1,2,3,7,8,9,10,11,12,13,14)]
  write.table( active_retrozyme_w50_df, 'active_retrozyme_w50.tbl',quote = F, sep = "\t", col.names = F, row.names = F)
  
  colnames(active_retrozyme_w50_df)=c('chrmid', 'start', 'stop', 'name', 'rnafold', 'strand', 
                                      'family', 'rpt', 'monomer_length','unitlength', 'intersection')

  coverage_df = bt.coverage(a = active_retrozyme_w50_df, b = bamfile, counts = T, split = T)
  colnames(coverage_df)=c('chrmid', 'start', 'stop', 'name', 'rnafold', 'strand', 
                          'family', 'rpt', 'monomer_length','unitlength', 'intersection', 'depth')
  return(coverage_df)
  #summary_depth=coverage_df[,c('intersection', 'depth')]
  #summary_depth=as.data.frame(table(summary_depth$depth), stringsAsFactors = F)
  #summary_depth$Var1=factor(summary_depth$Var1, levels = as.character(as.numeric(summary_depth$Var1)))
  #return(summary_depth)
}
est_depth=RNA_depth_level(bamfile = '~/remote2/EST/aln_minimap.bam')
head(est_depth)

RNA_seq_depth=read.csv2('active_retrozyme_w50.rnaseq', stringsAsFactors = F, header=F, sep="\t")
colnames(RNA_seq_depth)=c('chrmid', 'start', 'stop', 'name', 'rnafold', 'strand', 
                          'family', 'rpt', 'monomer_length','unitlength', 'intersection', 'depth')
head(RNA_seq_depth)
RNA_depth=data.frame(name=RNA_seq_depth$name, EST=est_depth$depth, RNAseq=RNA_seq_depth$depth)
RNA_depth$id=row.names(RNA_depth)
RNA_depth=tidyr::gather(RNA_depth, key='type', value = 'depth', c(2,3))
head(RNA_depth)

ggplot(RNA_depth, aes(x=name, y=depth, color=type))+
  geom_jitter(alpha=0.8)+
  coord_flip()+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(colour = "black",size=1.2))+
  theme(axis.text.y = element_text(size=12,face = "bold",colour = "black",family="arial"),
        axis.text.x = element_text(size=12,face = "bold",colour = "black",family="arial"),
        axis.title =  element_text(size=15,face = "bold",colour = "black",family="arial"),
        legend.position = 'right', 
        legend.text = element_text(size=15,face = "bold",colour = "black",family="arial"),
        legend.title = element_text(size=15,face = "bold",colour = "black",family="arial"))+
  ylab('Depth')+
  xlab('Name')+
  scale_color_manual(values=c("#1f77b4", "#9467bd"))

ggsave('RNAseq_est.jpeg', width = 8, height = 10)

pi_egg=RNA_depth_level(bamfile = '~/remote2/Small_RNA/PiRNA_xtro/SRR033660.egg.PiRNA.bam')
pi_egg$stage='egg'
pi_embryo_s8=RNA_depth_level(bamfile = '~/remote2/Small_RNA/PiRNA_xtro/SRR505561.embryo_s8.bam')
pi_embryo_s8$stage='embryo s8'
pi_embryo_s10=RNA_depth_level(bamfile = '~/remote2/Small_RNA/PiRNA_xtro/SRR505562.embryo_s10.bam')
pi_embryo_s10$stage='embryo s10'
pi_embryo_s18=RNA_depth_level(bamfile = '~/remote2/Small_RNA/PiRNA_xtro/SRR505563.embryo_s18.bam')
pi_embryo_s18$stage='embryo s18'

piRNA_total_df = rbind(pi_egg, pi_embryo_s8, pi_embryo_s10, pi_embryo_s18)

piRNA_total_df$stage=factor(piRNA_total_df$stage, levels = c('egg', 'embryo s8', 'embryo s10', 'embryo s18'))

head(piRNA_total_df)
ggplot(piRNA_total_df, aes(x=name, y=depth, color=stage))+
  geom_jitter(alpha=0.8)+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(colour = "black",size=1.2))+
  theme(axis.text.y = element_text(size=12,face = "bold",colour = "black",family="arial"),
        axis.text.x = element_text(size=12,face = "bold",colour = "black",family="arial"),
        axis.title =  element_text(size=15,face = "bold",colour = "black",family="arial"),
        legend.position = 'right', 
        legend.text = element_text(size=15,face = "bold",colour = "black",family="arial"),
        legend.title = element_text(size=15,face = "bold",colour = "black",family="arial"))+
  ylab('Depth')+
  xlab('Name')+
  coord_flip()+
  scale_color_manual(values=c("#1f77b4", "#9467bd", "#2ca02c", "#FFD966"))

ggsave('piRNA_depth.jpeg', width = 8, height = 10)


