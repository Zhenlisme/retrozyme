rm(list=ls())
gc()


library(bedtoolsr)
library(ggplot2)
options(bedtools.path = "~/bioinfo/bedtools2/bin/")

setwd('~/remote2/Toutatis_backup/Retrozyme/Retrozymes_detection/Xenopus_tropicalis/')
genome_size_df=read.csv2('genomesize.tbl', header = F, stringsAsFactors = F, sep = '\t')
head(genome_size_df)

genome_size_df=rbind(
  data.frame(V1=paste(genome_size_df$V1, '+',sep = ':'),
             V2=genome_size_df$V2, stringsAsFactors = F),
  data.frame(V1=paste(genome_size_df$V1, '-',sep = ':'),
             V2=genome_size_df$V2, stringsAsFactors = F))

genome_size_df=genome_size_df[order(genome_size_df$V1, genome_size_df$V2),]

gtf <- rtracklayer::import('~/remote2/genomedb/Xenopus_laevis/genomic.gtf')

gtf_df=as.data.frame(gtf)
head(gtf_df)

gene_df=gtf_df[gtf_df$type=='transcript', c(1,2,3,10,4,5)]
gene_df$seqnames=paste(gene_df$seqnames, gene_df$strand, sep = ':')
gene_df=gene_df[order(gene_df$seqnames, gene_df$start),]

intergenes_df=bt.complement(gene_df,g = genome_size_df, L = T)
intergenes_df$gene_id='intergenic'
colnames(intergenes_df)=c('seqnames','start','end','gene_id')
intergenes_df$seqnames=as.character(intergenes_df$seqnames)

exon_df=gtf_df[gtf_df$type=='exon', c(1,2,3,10,4,5)]
exon_df$seqnames=paste(exon_df$seqnames, exon_df$strand, sep = ':')
exon_df$seqnames=paste(exon_df$seqnames, exon_df$gene_id, sep = '__')
exon_df=exon_df[order(exon_df$seqnames, exon_df$start),]
exon_df=bt.merge(exon_df)
head(exon_df)

ano_df=as.data.frame(do.call('rbind', sapply(as.character(exon_df$V1), 
                    function(x){strsplit(x,split='__')})), stringsAsFactors = F)
head(ano_df)
exon_df=cbind(ano_df, exon_df)[,c(1,4,5,2)]
colnames(exon_df)=c('seqnames','start','end','gene_id')
exon_df=exon_df[order(exon_df$seqnames, exon_df$start),]

exonic_intergenic_df=rbind(exon_df, intergenes_df, stringsAsFactors = F)
exonic_intergenic_df=exonic_intergenic_df[order(exonic_intergenic_df$seqnames, exonic_intergenic_df$start), ]

head(exonic_intergenic_df)

intronic_df=bt.complement(exonic_intergenic_df, g=genome_size_df, L = T)
intronic_df=unique(bt.intersect(intronic_df, gene_df, wo = T, f = 1)[,c(1,2,3,7)])
colnames(intronic_df)=c('seqnames', 'start', 'end', 'gene_id')
intronic_df$seqnames=as.character(intronic_df$seqnames)
#intronic_df=aggregate(intronic_df$gene_id, list(intronic_df$seqnames, intronic_df$start, intronic_df$end), 
#          function(x){stringi::stri_paste(x, collapse = ', ')})

#colnames(intronic_df)=c('seqnames', 'start', 'end', 'gene_id')
intronic_df=intronic_df[order(intronic_df$seqnames, intronic_df$start),]

nc_info_df=as.data.frame(do.call('rbind',sapply(intronic_df$seqnames, function(x){strsplit(x, split = ':')})),
              stringsAsFactors = F)
intronic_df$type='intronic'
head(intronic_df)
intronic_df=cbind(nc_info_df, intronic_df, stringsAsFactors = F)[,c(1,4,5,6,7,2)]
colnames(intronic_df)=c('chrid', 'start', 'end', 'gene_id', 'type', 'strand')

exon_df$type='exon'
nc_info_df=as.data.frame(do.call('rbind',sapply(exon_df$seqnames, function(x){strsplit(x, split = ':')})),
                         stringsAsFactors = F)
head(cbind(nc_info_df, exon_df, stringsAsFactors = F))

exon_df=cbind(nc_info_df, exon_df, stringsAsFactors = F)[,c(1,4,5,6,7,2)]
colnames(exon_df)=c('chrid', 'start', 'end', 'gene_id', 'type', 'strand')
exon_df=exon_df[order(exon_df$chrid, exon_df$start),]

nc_info_df=as.data.frame(do.call('rbind',sapply(intergenes_df$seqnames, function(x){strsplit(x, split = ':')})),
                         stringsAsFactors = F)
head(cbind(nc_info_df, intergenes_df, stringsAsFactors = F))

intergenes_df=cbind(nc_info_df, intergenes_df, stringsAsFactors = F)[,c(1,4,5,6,6,2)]
colnames(intergenes_df)=c('chrid', 'start', 'end', 'gene_id', 'type', 'strand')
head(intergenes_df)

Anno_df=rbind(intergenes_df, exon_df, intronic_df, stringsAsFactors=F)
Anno_df=Anno_df[order(Anno_df$chrid, Anno_df$start), ]


write.table(Anno_df, file = '/home/zhenli/remote2/genomedb/Xenopus_laevis/genic_intergenic_region.tbl', 
            sep = '\t', col.names = T, row.names = F, quote = F)

#########################################################
setwd('~/remote2/Toutatis_backup/Retrozyme/Retrozymes_detection/Xenopus_tropicalis/')

Anno_df=read.csv2('~/remote2/genomedb/Xenopus_tropicalis/genic_intergenic_region.tbl', 
                  sep = '\t', header = T, stringsAsFactors = F)
head(Anno_df)
Anno_df=Anno_df[order(Anno_df$chrid, as.numeric(Anno_df$start)), ]
genome_size_df=read.csv2('genomesize.tbl', header = F, stringsAsFactors = F, sep = '\t')
head(genome_size_df)
retrozyme_df=read.csv2('repeat_summary.tbl', header = T, stringsAsFactors = F, sep = '\t')
retrozyme_df[which(as.numeric(retrozyme_df$rnafold)>=-5), ]

active_family=unique(retrozyme_df[which(as.numeric(retrozyme_df$rnafold)<=-5), 'family'])
active_family

retrozyme_df=data.frame()

for(clust in active_family){
  clusterfile=paste('./Clusters/', clust, '.clust.bed6', sep = '')
  cluster_bed = read.csv2(clusterfile, header = F, sep = '\t', stringsAsFactors = F)
  colnames(cluster_bed)=c('chrmid', 'start', 'end', 'identity', 'coverage', 'strand', 'cluster_id')
  cluster_bed$retro_family=clust
  retrozyme_df=rbind(retrozyme_df, cluster_bed)
}

retrozyme_df=retrozyme_df[, c(1,2,3,8,4,6)]
retrozyme_df=retrozyme_df[order(retrozyme_df$chrmid, retrozyme_df$start), ]
###bedtools fisher test 

bt.fisher(retrozyme_df, Anno_df[which(Anno_df$type=='exon'),], genome_size_df)
bt.fisher(retrozyme_df, Anno_df[which(Anno_df$type=='intergenic'),], genome_size_df)
bt.fisher(retrozyme_df, Anno_df[which(Anno_df$type=='intronic'),], genome_size_df)
###

retrozyme_coverage_df=bt.coverage(Anno_df, retrozyme_df, s = T)

head(retrozyme_coverage_df)
retrozyme_exon_df=retrozyme_coverage_df[which(retrozyme_coverage_df$V8!=0 & retrozyme_coverage_df$V5=='exon'), ]
write.table(retrozyme_exon_df, file = 'retrozyme_in_exon.bed', sep = '\t', col.names = F, row.names = F, quote = F)

region_total_length_df=aggregate(retrozyme_coverage_df$V9, list(retrozyme_coverage_df$V5), sum)
region_retro_length_df=aggregate(retrozyme_coverage_df$V8, list(retrozyme_coverage_df$V5), sum)
colnames(region_total_length_df)=c('region', 'region_length')
colnames(region_retro_length_df)=c('region', 'retro_length')

coverage_summary=merge(region_retro_length_df, region_total_length_df)
coverage_summary$coverage=100*coverage_summary$retro_length/coverage_summary$region_length
row.names(coverage_summary)=coverage_summary$region
coverage_summary

table_function=function(summary_tbl, region='exon', alternative='less'){
  retro_total=sum(summary_tbl$retro_length)
  retro_region=summary_tbl[region, 'retro_length']
  
  region_total=sum(summary_tbl$region_length)
  region_length=summary_tbl[region, 'region_length']
  null_pro=region_length/region_total
  prop_test=prop.test(retro_region, retro_total, p=null_pro, alternative = alternative)
  return(prop_test)
  #nonretro_region=region_length-retro_region
  #retro_nonregion=sum(summary_tbl$retro_length)-retro_region
  
  #nonretro_non_region=sum(summary_tbl[setdiff(summary_tbl$region,c(region)), 'region_length'])-
  # sum(summary_tbl[setdiff(summary_tbl$region,c(region)), 'retro_length'])
  
  #retro_matrix=matrix(c(retro_region, nonretro_region, retro_nonregion, nonretro_non_region), byrow = T, nrow = 2)
  #colnames(retro_matrix)=c('retro', 'nonretro')
  #row.names(retro_matrix)=c(region, paste('non', region, sep = '-'))
  #print(retro_matrix)
  #chisq_test=chisq.test(retro_matrix, correct = F)
  #return(chisq_test)
}

table_function(coverage_summary, region = 'exon', alternative = 'less')
table_function(coverage_summary, region = 'intergenic', alternative = 'greater')
table_function(coverage_summary, region = 'intronic', alternative = 'less')

## retrozymes are significantly excluded from exon and intronic region , but enriched in intergenic region.
coverage_summary


ggplot(coverage_summary, aes(x=region, y=coverage, fill=region))+
  geom_bar(stat = 'identity', width = 0.7)+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=20,family="arial"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  theme(axis.text = element_text(size=20,face = "bold",colour = "black",family="arial"),
        text = element_text(size=20,face = "bold",colour = "black",family="arial"))+
  ylab('Fold')+
  scale_fill_manual(values=c("#1f77b4", "#9467bd", "#2ca02c"))

ggsave("Genomic_region_proportion.png",width = 8,height = 8)




coverage_summary=tidyr::gather(coverage_summary, key=program, value = length, c(2,3))

ggplot(coverage_summary, aes(x=program, y=length, fill=region))+
  geom_bar(stat = 'identity', width = 0.7, position = position_dodge2())+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=20,family="arial"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom")+
  theme(axis.text = element_text(size=20,face = "bold",colour = "black",family="arial"),
        text = element_text(size=20,face = "bold",colour = "black",family="arial"),
        legend.text = element_text(face = "bold",colour = "black",family="arial"))+
  ylab('Length (nt)')+
  scale_x_discrete(labels=c(c('Genomic', 'Retrozyme')))+
  scale_fill_manual(values=c("#1f77b4", "#9467bd", "#2ca02c"))+
  scale_y_log10()

ggsave("Genomic_region.png",width = 8,height = 8)

