rm(list=ls())
gc()

library(here)
library(cowplot)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(ggstar)


PLE_prodf=read.csv2(here("retrozyme_data","PLE.summary.tbl"), stringsAsFactors = F, header = F, sep = '\t')
colnames(PLE_prodf)=c('species', 'PLE', 'count')
PLE_prodf$PLE=sapply(PLE_prodf$PLE, function(x){gsub('RT_', '', x)})
PLE_prodf$PLE=sapply(PLE_prodf$PLE, function(x){gsub('GIY_', '', x)})
PLE_prodf$PLE_type=sapply(PLE_prodf$PLE, function(x){
  split_str = strsplit(x, split = "/")[[1]]
  RT = split_str[1]
  GIY = split_str[2]
  if(RT==GIY){
    return(RT)
  }else{
    return('NA')
  }
})
PLE_prodf=PLE_prodf[which(PLE_prodf$PLE_typ != 'NA'), ]
PLE_prodf$species=sapply(PLE_prodf$species, function(x){gsub('_', ' ', x)})

species_groupinfo_df=read.csv2(here("retrozyme_data","species_classname.txt"), stringsAsFactors = F, header = F, sep = '\t')
colnames(species_groupinfo_df)=c('species', 'clades')
head(species_groupinfo_df)

PLE_prodf=merge(PLE_prodf, species_groupinfo_df, all.y = T)

groupInfo <- split(species_groupinfo_df$species, species_groupinfo_df$clades)


tree=read.tree(here("retrozyme_data","tree.nwk"))
tree$tip.label=sapply(tree$tip.label, function(x){gsub('_', ' ', x)})


PLE_prodf_sum = aggregate(PLE_prodf$count, list(PLE_prodf$species), sum)
colnames(PLE_prodf_sum)=c('species', 'sum')
head(PLE_prodf_sum)
PLE_prodf=merge(PLE_prodf, PLE_prodf_sum)
PLE_prodf$proportion=PLE_prodf$count/PLE_prodf$sum

p=ggtree(groupOTU(tree, groupInfo), size=1, aes(color=group, face='bold'))+
  #geom_tiplab(size=4, family='arial')+
  scale_color_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                              "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b", "black"))+
  theme_tree2(legend.position = 'none')+
  theme(axis.text.x = element_text(face = 'bold', family='arial', colour = 'black'),
        axis.line.x = element_line(linewidth=0.5))

species_name = rev(get_taxa_name(p))
color_range = rep(c("#2ca02c", "#f7b6d2","#9467bd" ,"#1f77b4", "#8c564b", 
                    "#d62728", "#e377c2", "#7f7f7f" ,"#ff7f0e"), 
                  c(1,3,1, 7, 1, 7, 6, 5, 4))

PLE_prodf$species=factor(PLE_prodf$species, levels = species_name)

pbar=ggplot(PLE_prodf, aes(x=count, y=species, fill=PLE_type))+
  geom_bar(stat = 'identity')+
  theme(panel.border = element_rect(colour = "black",size=0.5))+
  theme(axis.title = element_blank(), 
        axis.ticks.y =element_blank(),
        axis.text.y = element_text(family = 'arial', size = 10, face = 'bold.italic', hjust=0, colour= color_range),
        axis.text.x = element_text(family = 'arial', size = 10, face = 'bold', hjust=0.5, colour = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(face = 'bold'))+
  scale_fill_manual(values=c("#0E8388", "#F7C04A","#617143",  "#DF7857"))+
  scale_x_log10()

cowplot::plot_grid(p,pbar, rel_widths = c(1, 1.8), vjust = c(1, 10))

PLE_proportion_image = here("img","PLE_proportion.png")
ggsave(PLE_proportion_image, width = 8, height = 8)

genome_size_df=read.csv2(here("retrozyme_data","Genome.size.tbl"), stringsAsFactors = F, header = F, sep = '\t')
colnames(genome_size_df)=c('species', 'genomesize')
retrozyme_summary = read.csv2(here("retrozyme_data","repeat_summary.tbl"), stringsAsFactors = F, header = T, sep = '\t')

retrozyme_summary=merge(retrozyme_summary, genome_size_df)

retrozyme_summary$genomesize=as.numeric(retrozyme_summary$genomesize)/1000000000
retrozyme_summary$species=sapply(retrozyme_summary$species, function(x){gsub('_', ' ', x)})
retrozyme_summary$species=factor(retrozyme_summary$species, levels = species_name)
retrozyme_summary=merge(retrozyme_summary, species_groupinfo_df)
head(retrozyme_summary)

p_genome_size = ggplot(retrozyme_summary, aes(x=genomesize, y=species, fill = clades))+
  geom_bar(stat = 'identity')+
  scale_x_sqrt()+
  theme(panel.border = element_rect(colour = "black",size=0.5))+
  theme(axis.title = element_blank(), 
        axis.ticks.y =element_blank(),
        axis.text.y = element_text(family = 'arial', size = 10, face = 'bold.italic', hjust=0, colour= color_range),
        axis.text.x = element_text(family = 'arial', size = 10, face = 'bold', hjust=0.5, colour = 'black'),
        legend.position = 'none')+
  scale_fill_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))

head(retrozyme_summary)

p_active_count = ggplot(retrozyme_summary, aes(x=active_count, y=species, fill = clades))+
  geom_bar(stat = 'identity')+
  scale_x_log10()+
  theme(panel.border = element_rect(colour = "black",size=0.5))+
  theme(axis.title = element_blank(), 
        axis.ticks.y =element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = 'arial', size = 10, face = 'bold', hjust=0.5, colour = 'black'),
        legend.position = 'none')+
  scale_fill_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))

p_rpt_limit = ggplot(retrozyme_summary, aes(x=as.numeric(rpt_limit), y=species, fill = clades))+
  geom_bar(stat = 'identity')+
  scale_x_log10()+
  theme(panel.border = element_rect(colour = "black",size=0.5))+
  theme(axis.title = element_blank(), 
        axis.ticks.y =element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = 'arial', size = 10, face = 'bold', hjust=0.5, colour = 'black'),
        legend.position = 'none')+
  scale_fill_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))

p_family_num = ggplot(retrozyme_summary, aes(x=as.numeric(family_num), y=species, fill = clades))+
  geom_bar(stat = 'identity')+
  scale_x_log10()+
  theme(panel.border = element_rect(colour = "black",size=0.5))+
  theme(axis.title = element_blank(), 
        axis.ticks.y =element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = 'arial', size = 10, face = 'bold', hjust=0.6, colour = 'black'),
        legend.position = 'none')+
  scale_fill_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))

p_monomer_proportion = ggplot(retrozyme_summary, aes(x=as.numeric(monomer_proportion)*100, y=species, fill = clades))+
  geom_bar(stat = 'identity')+
  theme(panel.border = element_rect(colour = "black",size=0.5))+
  theme(axis.title = element_blank(), 
        axis.ticks.y =element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = 'arial', size = 10, face = 'bold', hjust=0.5, colour = 'black'),
        legend.position = 'none')+
  scale_fill_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                             "#f7b6d2", "#e377c2", "#7f7f7f", "#8c564b"))+
  scale_x_continuous(n.breaks = 3)
p_monomer_proportion

legend_share <- get_legend(p_monomer_proportion + 
    guides(fill = guide_legend(ncol = 1)) +
    theme(legend.position = "right", legend.title = element_blank(),
          legend.text = element_text(face = 'bold', family = 'arial'),
          axis.text = element_text(face = 'bold', family = 'arial')))

cowplot::plot_grid(p, p_genome_size, p_active_count, p_rpt_limit,
                   p_family_num, p_monomer_proportion, legend_share, 
                   rel_widths = c(1, 2, 1, 1, 1, 1, 0.8), 
                   ncol = 7)

ggsave('retrozyme_feature_character.jpeg', width = 15, height = 8)

