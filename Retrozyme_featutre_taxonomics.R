rm(list=ls())
gc()

library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(ggstar)
BiocManager::install("ggtreeExtra")

setwd('../Retrozyme/')
genome_size_df=read.csv2('Genome.size.tbl', stringsAsFactors = F, header = F, sep = '\t')
colnames(genome_size_df)=c('species', 'genomesize')

retrozyme_summary = read.csv2('repeat_summary.tbl', stringsAsFactors = F, header = T, sep = '\t')

retrozyme_summary=merge(retrozyme_summary, genome_size_df)

retrozyme_summary$genomesize=as.numeric(retrozyme_summary$genomesize)/1000000000
retrozyme_summary$species=sapply(retrozyme_summary$species, function(x){gsub('_', ' ', x)})


species_groupinfo_df=read.csv2('species_classname.txt', stringsAsFactors = F, header = F, sep = '\t')
colnames(species_groupinfo_df)=c('species', 'clades')
retrozyme_summary=merge(retrozyme_summary, species_groupinfo_df)

tree=read.tree("tree.nwk")
tree$tip.label=sapply(tree$tip.label, function(x){gsub('_', ' ', x)})
groupInfo <- split(species_groupinfo_df$species, species_groupinfo_df$clades)
tree_plot=ggtree(groupOTU(tree, groupInfo), size=1, aes(color=group, face='bold'))
species_name = rev(get_taxa_name(tree_plot))
color_range = rep(c("#2ca02c", "#f7b6d2","#9467bd" ,"#1f77b4", "#8c564b", 
                    "#d62728", "#e377c2", "#7f7f7f" ,"#ff7f0e"), 
                  c(1,3,1, 7, 1, 7, 6, 5, 4))
species_name
color_range

#########################################Retrozyme varies ####################
retrozyme_vary=read.csv2('repeat_vary.tbl', stringsAsFactors = F, header = T, sep = '\t')
retrozyme_vary$species=sapply(retrozyme_vary$species, function(x){gsub('_', ' ', x)})
retrozyme_vary=merge(retrozyme_vary, species_groupinfo_df)
head(retrozyme_vary)
retrozyme_vary=retrozyme_vary[order(retrozyme_vary$clades), ]
head(retrozyme_vary)


retrozyme_vary$species=factor(retrozyme_vary$species, levels = species_name)

ggplot(retrozyme_vary, aes(y=species, x=as.numeric(cons_length), color=clades))+
  geom_jitter(size=2)+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=20,family="arial"),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  theme(axis.text.y = element_text(size=20,face = "bold.italic",colour = color_range,family="arial"),
        axis.text.x = element_text(size=20,face = "bold",colour = 'black',family="arial"),
        text = element_text(size=20,face = "bold",colour = "black",family="arial"))+
  xlab('Monomer size')+
  scale_x_log10()+
  scale_color_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))

ggsave("MonomerSize_range.png",width = 12,height =12)


ggplot(retrozyme_vary, aes(y=species, x=as.numeric(rpt_limit), color=clades))+
  geom_jitter(size=4)+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=20,family="arial"),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  theme(axis.text.y = element_text(size=20,face = "bold.italic",colour = color_range,family="arial"),
        axis.text.x = element_text(size=20,face = "bold",colour = 'black',family="arial"),
        text = element_text(size=20,face = "bold",colour = "black",family="arial"))+
  xlab('Repeat limits')+
  scale_x_log10()+
  scale_color_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                              "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))

ggsave("RepeatLimits_range.png",width = 12,height = 12)

#######################################################

Genome_size_plot=ggtree(groupOTU(tree, groupInfo), size=1, aes(color=group))+geom_tiplab(size=4, family='arial')+
  scale_color_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                              "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b", "black"))+
  guides(color=F)+
  geom_fruit(data=retrozyme_summary, geom=geom_bar, 
             mapping=aes(y=species, x=sqrt(as.numeric(genomesize)), fill=clades), 
             pwidth=0.4, stat="identity", offset = 0.5, orientation="y", 
             axis.params=list(axis="x", text.size=4, nbreak=3, title='Genome size',family = 'arial',
                              title.size=5, title.height =0.06, hjust=0.5, vjust=1),
             grid.params=list())+
  scale_fill_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", 
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b", "black"),
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))



Active_Count=Genome_size_plot+geom_fruit(data=retrozyme_summary, geom=geom_bar, 
                          mapping=aes(y=species, x=log(as.numeric(active_count)+1, 10), fill=clades), 
                          pwidth=0.4, stat="identity", 
                          orientation="y", 
                          axis.params=list(
                            axis="x", text.size=4, nbreak=3, title='Active count',title.size=5, family = 'arial',
                            title.height =0.06, hjust=0.5, vjust=1),  
                          grid.params=list())

head(retrozyme_summary)

Repeatlimit_vary=Active_Count+geom_fruit(data=retrozyme_summary, geom=geom_bar,
                                         mapping=aes(y=species, x=log(as.numeric(rpt_limit)+1,10), fill=clades), 
                                         pwidth=0.4, stat='identity',
                                         orientation="y",
                                         axis.params=list(axis="x", text.size=4, title='Repeat limit',nbreak=3,family = 'arial',
                                                          title.size=5, title.height =0.06, hjust=0.5, vjust=1),
                                         grid.params=list())


#Monomer_size=Repeatlimit_vary+geom_fruit(data=retrozyme_summary, geom=geom_bar, 
#                            mapping=aes(y=species, x=as.numeric(cons_length), fill=clades), 
#                            pwidth=0.4,
#                            orientation="y",
#                            axis.params=list(axis="x", text.size=4, nbreak=4, title='Monomer size',
#                                             title.size=5, title.height =0.06, hjust=0.5, vjust=1),
#                            grid.params=list(),
#                            lwd=.1,
#                            outlier.size=0.5,
#                            outlier.stroke=0.08,
#                            outlier.shape=21)

family_num = Repeatlimit_vary+geom_fruit(data=retrozyme_summary, geom=geom_bar, 
             mapping=aes(y=species, x=log(as.numeric(family_num)+1, 10), fill=clades), 
             pwidth=0.4, stat="identity",
             orientation="y", 
             axis.params=list(
               axis="x", text.size=4, nbreak=3, title='Families',title.size=5, title.height =0.06, family = 'arial',
               hjust=0.5, vjust=1),  
             grid.params=list())



head(retrozyme_summary)
monomer_proportion=family_num+
  geom_fruit(data=retrozyme_summary, geom=geom_bar, 
             mapping=aes(y=species, x=as.numeric(monomer_proportion), fill=clades), 
             pwidth=0.4, stat="identity",
             orientation="y", 
             axis.params=list(axis="x", text.size=4, nbreak=3, title='Monomer proportion',family = 'arial',
                              title.size=5, title.height =0.06, hjust=0.5, vjust=1), 
             grid.params=list())+
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size = 15), 
        legend.text=element_text(size=10),
        legend.spacing.y = unit(0.5, "cm"),
        legend.position = 'bottom')
monomer_proportion

ggsave("Feature_accross_trees.png",width = 20,height =10)

