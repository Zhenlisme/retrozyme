library(here)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(ggstar)


#groupInfo <- split(species_groupinfo_df$species, species_groupinfo_df$clades)


tree=read.tree(here("retrozyme_data","pLTR","XT","pltr_plus_retrozyme.phylip_phyml_tree.txt"))
tree$node.label=sapply(tree$node.label, function(x){
  if(stringi::stri_length(x)==0){
    return(x)
  }
  else{
    bootstrap_value=round(as.numeric(x), 1)
    if(bootstrap_value>=0.7){
      return(as.character(bootstrap_value))
    }else{
      return(c(''))
    }
  }
})

groupInfo=as.data.frame(do.call('rbind', lapply(as.vector(tree$tip.label), function(x){
  species_label=length(grep('^X', x))
  if(species_label!=0){
    family=strsplit(x, '-')[[1]][1]
    family=paste('Rtz', family, sep = '_')
  }else{
    ple_label=length(grep('^Neptune|^Penelope', x))
    if(ple_label!=0){
      family='PLE-like'
    }else{
      family='Unknown TE'
      #rep_label=length(grep('^REP4', x))
      #if(rep_label!=0){
      #  family='REP4_XT'
      #}else{
      #  family='Other'
      #}
    }}
  return(c(x, family))
  })), stringsAsFactors = F)
colnames(groupInfo)=c('tename', 'group')

groupInfo <- split(groupInfo$tename, groupInfo$group)

ggtree(groupOTU(tree, groupInfo), size=.7, aes(color=group), layout = 'circular')+
  geom_nodelab(geom = 'text', color='black', size=2)+
  scale_color_manual(values=c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                              "#8c564b", "black"))+
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size = 15, family = 'arial'), 
        legend.text=element_text(size=10, family = 'arial'),
        legend.spacing.y = unit(0.1, "cm"),
        legend.position = 'right')

Rtz_PLE_tree_image = here("img","Rtz_PLE.tree.png")
ggsave(Rtz_PLE_tree_image,width = 8,height =8)

