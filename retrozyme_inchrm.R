library(here)
library(ggplot2)
options(bedtools.path = "/opt/homebrew/bin/")
library(bedtoolsr)




coverage_df_func=function(wkdir){
  retrozyme_tbl_file <- paste(wkdir,"repeat_summary.tbl",sep="/")
  retrozyme_tbl=read.csv2(retrozyme_tbl_file,  sep = '\t', header = T, stringsAsFactors = F)
  retrozyme_tbl$rnafold=as.numeric(retrozyme_tbl$rnafold)
  active_family=unique(retrozyme_tbl[which(retrozyme_tbl$rnafold<=-5), 'family'])

  chrm_info_file=paste(wkdir,"chrminfo.txt",sep="/")
  chrm_info=read.csv(chrm_info_file, sep = '\t', header = F, stringsAsFactors = F)
  colnames(chrm_info)=c('chrmid', 'chrmname')

  genome_window_df_file <- paste(wkdir,"genomewindow_w100k_s50k.bed",sep="/")
  genome_window_df=read.csv(genome_window_df_file, stringsAsFactors = F, header = F, sep = '\t')
  colnames(genome_window_df)=c('chrmid', 'start', 'end')

  genome_window_df=genome_window_df[which(genome_window_df$chrmid %in% chrm_info$chrmid),]

  ## To classify monomers and multimers
  monomer_multimer_split=function(cluster_file){
  cluster_bed = read.csv2(cluster_file, header = F, sep = '\t', stringsAsFactors = F)
  colnames(cluster_bed)=c('chrmid', 'start', 'end', 'identity', 'coverage', 'strand', 'cluster_id')
  cluster_idlist=cluster_bed$cluster_id
  cluster_iddf=as.data.frame(table(cluster_idlist), stringsAsFactors = F)
  monomer_idlist=cluster_iddf[which(cluster_iddf$Freq==1), "cluster_idlist"]
  dimer_idlist=cluster_iddf[which(cluster_iddf$Freq==2), "cluster_idlist"]
  multimer_idlist=cluster_iddf[which(cluster_iddf$Freq>2), "cluster_idlist"]
  monomer_bed=cluster_bed[which(cluster_bed$cluster_id %in% monomer_idlist), ]
  dimer_bed=cluster_bed[which(cluster_bed$cluster_id %in% dimer_idlist), ]
  multimer_bed= cluster_bed[which(cluster_bed$cluster_id %in% multimer_idlist), ]
  return(list(monomer_bed, dimer_bed, multimer_bed))
}

  monomer_total_bed=data.frame()
  dimer_total_bed=data.frame()
  multimer_total_bed=data.frame()

  ## To merge all monomer dataframe and multimer frames, respectively.
  for(clustername in active_family){
  clusterfile=paste(wkdir,"/Clusters/", clustername, ".clust.bed6", sep = "")
  classify_list = monomer_multimer_split(clusterfile)
  monomer_total_bed=rbind(monomer_total_bed, classify_list[[1]])
  dimer_total_bed=rbind(dimer_total_bed, classify_list[[2]])
  multimer_total_bed=rbind(multimer_total_bed, classify_list[[3]])
  }

  monomer_total_bed=monomer_total_bed[order(monomer_total_bed$chrmid, monomer_total_bed$start), ]
  dimer_total_bed=dimer_total_bed[order(dimer_total_bed$chrmid, dimer_total_bed$start), ]
  multimer_total_bed=multimer_total_bed[order(multimer_total_bed$chrmid, multimer_total_bed$start), ]

  multimer_coverage_df=bt.coverage(genome_window_df, multimer_total_bed)
  colnames(multimer_coverage_df)=c('chrmid', 'start', 'stop', 'counts', 'length', 'windowsize', 'multimer')

  monomer_coverage_df=bt.coverage(genome_window_df, monomer_total_bed)
  colnames(monomer_coverage_df)=c('chrmid', 'start', 'stop', 'counts', 'length', 'windowsize', 'monomer')

  dimer_coverage_df=bt.coverage(genome_window_df, dimer_total_bed)
  colnames(dimer_coverage_df)=c('chrmid', 'start', 'stop', 'counts', 'length', 'windowsize', 'dimer')

  coverage_total_df=cbind(multimer_coverage_df[,c(1,2,3,6,7)], monomer_coverage_df[,7], dimer_coverage_df[,7])
  colnames(coverage_total_df)[c(6, 7)]=c('monomer', 'dimer')

  coverage_total_df=tidyr::gather(coverage_total_df, key='type', value='Coverage', 
                                  c('monomer', 'dimer','multimer'))
  coverage_total_df$Coverage=100*coverage_total_df$Coverage

  coverage_total_df[which(coverage_total_df$type=='multimer'), 'Coverage']=-1* coverage_total_df[which(coverage_total_df$type=='multimer'), 'Coverage']

  coverage_total_df$type=factor(coverage_total_df$type, levels = c('monomer', 'dimer', 'multimer'))
  coverage_total_df=merge(coverage_total_df, chrm_info)
  coverage_total_df$chrmname=factor(coverage_total_df$chrmname, levels = chrm_info$chrmname)
  coverage_total_df$start=coverage_total_df$start/1000000
  return(coverage_total_df)
}

Xl_dir<-here("retrozyme_data","Retrozymes_detection","Xenopus_laevis")
Xt_dir<-here("retrozyme_data","Retrozymes_detection","Xenopus_tropicalis")
Xb_dir<-here("retrozyme_data","Retrozymes_detection","Xenopus_borealis")

Xl_coverage_df=coverage_df_func(Xl_dir)
Xt_coverage_df=coverage_df_func(Xt_dir)
Xb_coverage_df=coverage_df_func(Xb_dir)

Xl_coverage_df$species='XL'
Xb_coverage_df$species='XB'
Xt_coverage_df$species='XT'
Xenopus=rbind(Xl_coverage_df, Xb_coverage_df, Xt_coverage_df)
head(Xenopus)
head(Xt_coverage_df)

Xenopus$label=paste(Xenopus$chrmname, Xenopus$species, sep = '-')
levels(factor(Xenopus$label))

chrmbasenm = c("Chr1L", 'Chr2L', 'Chr3L', 'Chr4L', 'Chr5L', 'Chr6L', 'Chr7L', 
               'Chr8L', 'Chr9_10L', "Chr1S", 'Chr2S', 'Chr3S', 'Chr4S', 'Chr5S', 
               'Chr6S', 'Chr7S', 'Chr8S', 'Chr9_10S')
Xt_chrmname=c("Chr1", 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10')

c(paste(chrmbasenm, c('XL'), sep = '-'),
  paste(chrmbasenm, c('XB'), sep = '-'),
  paste(Xt_chrmname, c('XT')))

Xenopus$label=factor(Xenopus$label, levels = c(paste(chrmbasenm, c('XL'), sep = '-'),
                                               paste(chrmbasenm, c('XB'), sep = '-'),
                                               paste(Xt_chrmname, c('XT'), sep = '-')))

plot_xenopus <- ggplot(Xenopus, aes(x=start, y=Coverage, color=type))+
  geom_line(linewidth=0.6)+
  theme_bw()+
  scale_x_continuous(name = "Position (Mb)",n.breaks = 3) +
  scale_y_continuous(name = "Coverage %") +
  facet_wrap(.~label, ncol = 9, scales = 'free_x')+
  labs(title = "Retrozyme coverage across Xenopus chromosomes")+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  coord_cartesian(ylim=c(-2,2))+
  theme(legend.position = "right",legend.title = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, family = 'arial'))+
  theme(panel.spacing = unit(0.8, "lines"),
        axis.text.x = element_text(hjust = 0.8))+
  scale_color_manual(values = c('#B99B6B', '#698269', '#913175'))


plot_xenopus_laevis <- Xenopus %>% filter(species=="XL") %>% ggplot(aes(x=start, y=Coverage, color=type))+
  geom_line(linewidth=0.6)+
  theme_bw()+
  scale_x_continuous(name = "Position (Mb)",n.breaks = 3) +
  scale_y_continuous(name = "Coverage %") +
  facet_wrap(.~label, ncol = 9, scales = 'free_x')+
  labs(title = "Retrozyme coverage across Xenopus laevis chromosomes")+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  coord_cartesian(ylim=c(-2,2))+
  theme(legend.position = "right",legend.title = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, family = 'arial'))+
  theme(panel.spacing = unit(0.8, "lines"),
        axis.text.x = element_text(hjust = 0.8))+
  scale_color_manual(values = c('#B99B6B', '#698269', '#913175'))
plot_xenopus_laevis
Xl_chrm_distribution<-here("img","Xl_chrm_distribution.png")
ggsave(Xl_chrm_distribution, width = 14, height = 10)



plot_xenopus_borealis <- Xenopus %>% filter(species=="XB") %>% ggplot(aes(x=start, y=Coverage, color=type))+
  geom_line(linewidth=0.6)+
  theme_bw()+
  scale_x_continuous(name = "Position (Mb)",n.breaks = 3) +
  scale_y_continuous(name = "Coverage %") +
  facet_wrap(.~label, ncol = 9, scales = 'free_x')+
  labs(title = "Retrozyme coverage across Xenopus borealis chromosomes")+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  coord_cartesian(ylim=c(-2,2))+
  theme(legend.position = "right",legend.title = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, family = 'arial'))+
  theme(panel.spacing = unit(0.8, "lines"),
        axis.text.x = element_text(hjust = 0.8))+
  scale_color_manual(values = c('#B99B6B', '#698269', '#913175'))
plot_xenopus_borealis
Xb_chrm_distribution<-here("img","Xb_chrm_distribution.png")
ggsave(Xb_chrm_distribution, width = 14, height = 10)



plot_xenopus_tropicalis <- Xenopus %>% filter(species=="XT") %>% ggplot(aes(x=start, y=Coverage, color=type))+
  geom_line(linewidth=0.6)+
  theme_bw()+
  scale_x_continuous(name = "Position (Mb)",n.breaks = 3) +
  scale_y_continuous(name = "Coverage %") +
  facet_wrap(.~label, ncol = 9, scales = 'free_x')+
  labs(title = "Retrozyme coverage across Xenopus tropicalis chromosomes")+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  coord_cartesian(ylim=c(-2,2))+
  theme(legend.position = "right",legend.title = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, family = 'arial'))+
  theme(panel.spacing = unit(0.8, "lines"),
        axis.text.x = element_text(hjust = 0.8))+
  scale_color_manual(values = c('#B99B6B', '#698269', '#913175'))

plot_xenopus_tropicalis
Xt_chrm_distribution<-here("img","Xt_chrm_distribution.png")
ggsave(Xt_chrm_distribution, width = 14, height = 10)


