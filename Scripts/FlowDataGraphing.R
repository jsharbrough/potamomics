#### Flow Cytometry graphing #
library(tidyverse)
#setwd("Documents/Biology/Sharbrough_Lab/Projects/Potamomics/Flow Cytometry/")
d<-read.table('Pest_Pant_CRBC_20230127-20230213.txt',header=T,sep='\t')
pdf("Pest-31_Y2-01_CRBC_20230127-20230213.DAPI_Fluorescence.pdf", width = 3.25,height=2.44)
ggplot(data=d,mapping=aes(x=DAPI_A,fill=Species)) + geom_density(alpha=0.25) + scale_x_continuous(name="DAPI Fluorescence", limits=c(50,325), breaks=c(50,100,150,200,250,300)) + theme(panel.grid=element_line(colour='black'),panel.background=element_rect(fill='transparent',colour='black'),axis.ticks=element_line(colour='black',size=0.5),axis.text=element_text(colour='black'),legend.position=c(0.15,0.75),legend.title=element_blank()) + scale_y_continuous(name="Density") + scale_fill_manual(values=c("darkgrey", "cornflowerblue", "forestgreen"))
dev.off()
