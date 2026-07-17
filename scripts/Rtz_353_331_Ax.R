library(here)
library(ggplot2)
rnafold_value_Amx=read.csv2(here("retrozyme_data","Retrozymes_detection","Ambystoma_mexicanum","Rtz353_331.rnafold"),
                            stringsAsFactors = F, header = F)

p=ggplot(rnafold_value_Amx, aes(x=as.numeric(V1)))+
  geom_histogram(binwidth = 0.5,color='black')+
  theme(axis.text = element_text(size=12,face = "bold",colour = "black", family = 'arial'),
        text = element_text(size=12,face = "bold",colour = "black", family = 'arial'))+
  xlab('rnafold') +
  ggtitle("RNAfold MFE values for the axolotl retrozyme Rtz331_353_Ax")

Rtz331_353_Ax_rnafold_image = here("img","Rtz331_353_Ax_rnafold.png")
ggsave(Rtz331_353_Ax_rnafold_image, p, width = 7, height = 6)
