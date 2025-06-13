
library(ggplot2)
rnafold_value_Amx=read.csv2('~/Toutatis_backup/Retrozyme/Retrozymes_detection/Ambystoma_mexicanum/Rtz353_331.rnafold',
                            stringsAsFactors = F, header = F)

p=ggplot(rnafold_value_Amx, aes(x=as.numeric(V1)))+
  geom_histogram(binwidth = 0.5,color='black')+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  xlab('rnafold')

ggsave('Rtz331_353_Ax_rnafold.png',p, width = 6, height = 6)
