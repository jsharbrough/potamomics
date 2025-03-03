#load libraries
library(tidyverse)

#load data
setwd('~/Documents/Biology/Sharbrough_Lab/Projects/Potamomics/Genome Duplication/Pant-uncollapsed_collapsed')
d <- read.table('Pant_collapsed-uncollapsed.1-to-2_Orthogroups.contigPositions.txt',header=T,sep='\t')
geneD <- read.table('Pant_collapsed-uncollapsed.1-to-2_Orthogroups.scaffoldGenePositions.txt',header=T,sep='\t')


#plot synteny across the genome
scaffold_palettes <- geneD %>% group_by(hapScaffold) %>% summarise(OGIDs = list(unique(OGID))) %>% rowwise() %>% mutate(palette = list(scales::hue_pal()(length(OGIDs)))) %>% ungroup()
geneD <- geneD %>% left_join(scaffold_palettes, by = "hapScaffold") %>% rowwise() %>% mutate(fill_color = palette[which(OGID == OGIDs)]) %>% ungroup()
pdf("synteny_plot.pdf",width=50,height=210.6445313)
ggplot(data = d) +  facet_wrap(~hapScaffold, ncol = 32, scales = "free_x") +  geom_rect(data = geneD,           mapping = aes(xmin = contigScaffoldGeneStart,                          xmax = contigScaffoldGeneStop,                          ymin = Y - 0.125,                          ymax = Y + 0.125,                          fill = fill_color),            color = 'transparent') +  geom_segment(data = d, mapping = aes(x = x, y = newY, xend = xend, yend = newY)) +  geom_segment(data = geneD,               mapping = aes(x = contigScaffoldGeneMiddle,                             y = Y - 0.13,                             xend = scaffoldGene_scaffoldGeneMiddle,                             yend = Y - 0.375),               lwd = 0.05) +  theme(panel.background = element_rect(fill = 'transparent', colour = 'black'),        panel.grid = element_blank(),        legend.position = 'NA',        axis.text.y = element_blank(),        axis.ticks.y = element_blank(),        axis.title.y = element_blank()) +  scale_y_continuous(limits = c(-.7, 0.7)) +  geom_text(data = d, mapping = aes(x = -300000, y = -0.5, label = "uncollapsed\ncontigs"), angle = 90) +  geom_text(data = d, mapping = aes(x = -300000, y = 0, label = "collapsed\nscaffold"), angle = 90) +  geom_text(data = d, mapping = aes(x = -300000, y = 0.5, label = "uncollapsed\ncontigs"), angle = 90) +  scale_x_continuous(name = "Nucleotide Position (bp)", limits = c(-400000, NA)) + scale_fill_identity()
dev.off()

#subset for scaffold 14
gene_scf14<-geneD[which(geneD$hapScaffold=='scaffold14'),]
scf14<-d[which(d$hapScaffold=='scaffold14'),]
scf14_1to1<-read.table('Pant_collapsed-uncollapsed.1-to-1_Orthogroups.scaffoldGenePositions.txt',header=T,sep='\t')
potSpeciesSingleCopyGenes<-read.table('Pant_poly_Pkait_Pest.1-to-1_Orthogroups.scaffoldGenePositions.txt',header=T,sep='\t')
potSpeciesContigs<-read.table('Pant_poly_Pkait_Pest.Orthogroups.contigPositions.txt',header=T,sep='\t')
potSpeciesDoubleCopyGenes<-read.table('Pant_poly_Pkait_Pest.2-to-1_Orthogroups.scaffoldGenePositions.txt',header=T,sep='\t')
pdf("synteny_plot.scaffold14.pdf",width=6.25,height=9)
ggplot() + geom_segment(data=potSpeciesSingleCopyGenes,mapping=aes(x=contigScaffoldGeneMiddle,xend=partnerGene_scaffoldGeneMiddle,yend=partnerGene_line_Y,y=line_Yend),col='darkgrey',lwd=0.05) + geom_rect(data=potSpeciesSingleCopyGenes,mapping=aes(xmin=contigScaffoldGeneStart,xmax=contigScaffoldGeneStop,ymin=Y-0.0625,ymax=Y+0.0625),fill="black",col='transparent') + geom_rect(data=scf14_1to1,mapping=aes(xmin=contigScaffoldGeneStart,xmax=contigScaffoldGeneStop,ymin=Y-0.0625,ymax=Y+0.0625),fill="black",col='transparent') + geom_segment(data=scf14_1to1,mapping=aes(x=contigScaffoldGeneMiddle,xend=partnerGene_geneContigScaffoldMiddle,y=line_Y,yend=line_Yend),col='darkgrey',lwd=0.05) + geom_rect(data=gene_scf14,mapping=aes(xmin=contigScaffoldGeneStart,xmax=contigScaffoldGeneStop, ymin=Y-0.125,ymax=Y+0.125,fill=OGID),col='transparent') + geom_segment(data=scf14,mapping=aes(x=x,y=newY,xend=xend,yend=newY)) + geom_segment(data=potSpeciesContigs,mapping=aes(x=x,y=newY,xend=xend,yend=newY)) + geom_segment(data=gene_scf14,mapping=aes(x=contigScaffoldGeneMiddle,y=-0.13,xend=partnerGene_contigScaffoldGeneMiddle,yend=-0.375, col=OGID),lwd=0.2) + geom_segment(data=gene_scf14,mapping=aes(x=contigScaffoldGeneMiddle,y=0.13,xend=scaffoldGene_scaffoldGeneMiddle,yend=0.375, col=OGID),lwd=0.2) + theme(panel.background=element_rect(fill='transparent',colour='black'),panel.grid=element_blank(),legend.position='NA',axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(colour='black')) + scale_y_continuous(limits=c(NA,NA)) + annotate("text",x=-300000,y=c(-0.5,0,0.5),label=c("uncollapsed\ncontigs","uncollapsed\ncontigs","collapsed\nscaffold"),angle=90) + scale_x_continuous(name="Collapsed Assembly Scaffold 14 Nucleotide Position (bp)", limits=c(-400000,NA)) + geom_rect(data=potSpeciesDoubleCopyGenes,mapping=aes(xmin=contigScaffoldGeneStart,xmax=contigScaffoldGeneStop, ymin=Y-0.125,ymax=Y+0.125,fill=OGID),col='transparent') + geom_segment(data=potSpeciesDoubleCopyGenes,mapping=aes(x=contigScaffoldGeneMiddle,y=line_Yend,xend=partnerGene_scaffoldGeneMiddle,yend=partnerGene_line_Y, col=OGID),lwd=0.2)
dev.off()

