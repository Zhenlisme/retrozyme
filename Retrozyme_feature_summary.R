rm(list=ls())
gc()
library(ggplot2)
library(bedtoolsr)
options(bedtools.path = "/home/zhenli/bioinfo/bedtools2/bin/")

setwd('/home/zhenli/remote2/Toutatis_backup/Retrozyme/')


## To obtain summary information
Retrozyme_summary=function(species){
  print(species)
  retrozyme_tblfile=paste('Retrozymes_detection/', species, '/repeat_summary.tbl', sep = '')
  if(!file.exists(retrozyme_tblfile)){
    return(c(species, 0, 0, 0, 0))
  }
  retrozyme_tbl=read.csv2(retrozyme_tblfile,  sep = '\t', header = T, stringsAsFactors = F)
  #colnames(retrozyme_tbl)=c('Chrm', 'start', 'stop','rpt', 'consensus_size', 'strand', 'type', 'rtzname', 'motifs', 'active_count', 'active_rpt', 'cons', 'hhr', 'family', 'appearence')
  retrozyme_tbl$rnafold=as.numeric(retrozyme_tbl$rnafold)
  retrozyme_tbl=retrozyme_tbl[which(retrozyme_tbl$rnafold<=-5), ]
  active_family=unique(retrozyme_tbl$family)
  FAMILY_NUM=length(active_family)
  ACTIVE_COUNT=sum(as.numeric(retrozyme_tbl$active_count))
  if(ACTIVE_COUNT==0){
    return(c(species, ACTIVE_COUNT, 0, 0, 0))
  }
  head(retrozyme_tbl)
  RPT_LIMIT=max(as.numeric(retrozyme_tbl$rpt))
  monomer_proportion_list=c()
  for(clust in active_family){
    clusterfile=paste('Retrozymes_detection/', species,'/Clusters/', clust, '.clust.bed6', sep = '')
    cluster_bed = read.csv2(clusterfile, header = F, sep = '\t', stringsAsFactors = F)
    colnames(cluster_bed)=c('chrmid', 'start', 'end', 'identity', 'coverage', 'strand', 'cluster_id')
    cluster_idlist=cluster_bed$cluster_id
    cluster_iddf=as.data.frame(table(cluster_idlist), stringsAsFactors = F)
    monomer_idlist=cluster_iddf[which(cluster_iddf$Freq==1), "cluster_idlist"]
    multimer_idlist=cluster_iddf[which(cluster_iddf$Freq>=2), "cluster_idlist"]
    
    monomer_bed=cluster_bed[which(cluster_bed$cluster_id %in% monomer_idlist), ]
    multimer_bed= cluster_bed[which(cluster_bed$cluster_id %in% multimer_idlist),]
    monomer_count = length(row.names(monomer_bed))
    multimer_count=length(row.names(multimer_bed))
    monomer_poportion=round(monomer_count/(monomer_count+multimer_count), 2)
    monomer_proportion_list=c(monomer_proportion_list, monomer_poportion)
  }
  MONOMER_PROPORTION=median(monomer_proportion_list)   ## Use median instad of mean to represent the monomerproportion.
  return(c(species, ACTIVE_COUNT, RPT_LIMIT, FAMILY_NUM, MONOMER_PROPORTION))
}

Retrozyme_summary_df=as.data.frame(do.call('rbind', lapply(list.files('Retrozymes_detection/'), Retrozyme_summary)), stringsAsFactors = F)
Retrozyme_summary_df
colnames(Retrozyme_summary_df)=c('species', 'active_count', 'rpt_limit', 'family_num', 'monomer_proportion')

write.table(file = 'repeat_summary.tbl', x = Retrozyme_summary_df, col.names = T, sep = '\t', quote = F, row.names = F)

### To obtain consensus length varies, repeat count varies, cluster number varies
Retrozyme_VARY=function(species){
  VARY_df=data.frame()
  retrozyme_tblfile=paste('Retrozymes_detection/', species, '/repeat_summary.tbl', sep = '')
  if(!file.exists(retrozyme_tblfile)){
    return(c('0', 0,0,0,0,species))
  }
  retrozyme_tbl=read.csv2(retrozyme_tblfile,  sep = '\t', header = T, stringsAsFactors = F)
  retrozyme_tbl$rnafold=as.numeric(retrozyme_tbl$rnafold)
  retrozyme_tbl=retrozyme_tbl[which(retrozyme_tbl$rnafold<=-5), ]
  active_family=unique(retrozyme_tbl$family)
  ACTIVE_COUNT=sum(as.numeric(retrozyme_tbl$active_count))
  if(ACTIVE_COUNT==0){
    return(c('0', 0,0,0,0,species))
  }
  CONS_VARY_df=aggregate(as.numeric(retrozyme_tbl$consensus_size), list(retrozyme_tbl$family), median)
  colnames(CONS_VARY_df)=c('family', 'cons_length')
  CLUSTER_VARY_df=aggregate(as.numeric(retrozyme_tbl$appearence), list(retrozyme_tbl$family), median)
  colnames(CLUSTER_VARY_df)=c('family', 'appearence')
  
  RPTLIMIT_VARY_df=aggregate(as.numeric(retrozyme_tbl$rpt), list(retrozyme_tbl$family), max)
  colnames(RPTLIMIT_VARY_df)=c('family', 'rpt_limit')
  
  monomer_proportion_list=c()
  for(family in active_family){
    clusterfile=paste('Retrozymes_detection/', species,'/Clusters/', family, '.clust.bed6', sep = '')
    cluster_bed = read.csv2(clusterfile, header = F, sep = '\t', stringsAsFactors = F)
    colnames(cluster_bed)=c('chrmid', 'start', 'end', 'identity', 'coverage', 'strand', 'cluster_id')
    cluster_idlist=cluster_bed$cluster_id
    cluster_iddf=as.data.frame(table(cluster_idlist), stringsAsFactors = F)
    monomer_idlist=cluster_iddf[which(cluster_iddf$Freq==1), "cluster_idlist"]
    multimer_idlist=cluster_iddf[which(cluster_iddf$Freq>=2), "cluster_idlist"]
    
    monomer_bed=cluster_bed[which(cluster_bed$cluster_id %in% monomer_idlist), ]
    multimer_bed= cluster_bed[which(cluster_bed$cluster_id %in% multimer_idlist),]
    monomer_count = length(row.names(monomer_bed))
    multimer_count=length(row.names(multimer_bed))
    monomer_poportion=round(monomer_count/(monomer_count+multimer_count), 2)
    monomer_proportion_list=c(monomer_proportion_list, c(family, monomer_poportion))
  }
  monomer_proportion_df = as.data.frame(matrix(monomer_proportion_list, ncol = 2, byrow = T), stringsAsFactors = F)
  colnames(monomer_proportion_df)=c('family', 'monomer_proportion')
  VARY_df=merge(merge(CLUSTER_VARY_df, CONS_VARY_df), monomer_proportion_df)
  VARY_df=merge(VARY_df, RPTLIMIT_VARY_df)
  VARY_df$species=species  
  return(VARY_df)
}

Retrozyme_VARY_df=as.data.frame(do.call('rbind', lapply(list.files('Retrozymes_detection/'), Retrozyme_VARY)), stringsAsFactors = F)
write.table(file = 'repeat_vary.tbl', x = Retrozyme_VARY_df, col.names = T, sep = '\t', quote = F, row.names = F)

head(Retrozyme_VARY_df)

ggplot(Retrozyme_VARY_df, aes(x=as.numeric(monomer_proportion), y=as.numeric(rpt_limit)))+
  geom_point()+
  facet_wrap(.~species, ncol = 5, scales = 'free')
