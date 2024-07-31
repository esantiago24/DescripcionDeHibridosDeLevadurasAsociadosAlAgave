library(ggplot2)
library(GGally)
library(readxl)
library(reshape)
library(dplyr)
library(RColorBrewer)

setwd("~/Desktop/Genomic_Sciences/CuartoAño/Trees/")
distances<-read.csv("./metadata/DistsnGroupsOnly.csv",header = TRUE,sep = ",")
colnames(distances)[1] = "Familia"
colnames(distances)[2] = "Distancia a S. cerevisiae"
colnames(distances)[3] = "Distancia a S. paradoxus"
distancesNoNatHybs4Thesis<-distances[c(-64,-65),]

distancesNoNatHybs4Thesis$Familia<-replace(distancesNoNatHybs4Thesis$Familia, distancesNoNatHybs4Thesis$Familia!=14,1)
summary(distancesNoNatHybs)
#varSACE<-var(distancesNoNatHybs4Thesis$`Distancia a S. cerevisiae`)
#varSAPA<-var(distancesNoNatHybs4Thesis$`Distancia a S. paradoxus`)
#sdSACE<-sd(distancesNoNatHybs4Thesis$`Distancia a S. cerevisiae`)
#sdSAPA<-sd(distancesNoNatHybs4Thesis$`Distancia a S. paradoxus`)

parallel_plot<-ggparcoord(distancesNoNatHybs4Thesis,columns=c(2,3),groupColumn = 1,showPoints = TRUE, alphaLines = 0.5,scale = "globalminmax") + xlab("") + ylab("Distancia filogenética (sustituciones por sitio)") + scale_color_brewer(palette="Dark2")
#ggparcoord(distances,columns=c(2,3),groupColumn = 1,showPoints = TRUE, alphaLines = 0.3,scale = "globalminmax")
parallel_plot

distancesNoNatHybs<-melt(distancesNoNatHybs)

ggplot(distancesNoNatHybs, aes(x=value, fill=variable)) + geom_density(alpha=.5) + xlab("Distancia filogenética (sustituciones por sitio)") + ylab("Densidad") + theme(legend.position = "none") + 
  annotate("text", x = 0.17, y = 5, size = 7, label = "S. cerevisiae", fontface = "italic", color="#F8766D")+
  geom_segment(aes(x = 0.03, y = 0, xend = 0.03, yend = 10.5)) +
  annotate("text", x = 0.08, y = 25, size = 7, label = "S. paradoxus", fontface = "italic", color="#00BFC4")+
  geom_segment(aes(x = 0.04, y =0, xend = 0.04, yend = 27))

density_plot<-ggplot(distancesNoNatHybs, aes(x=x)) + xlab("Distancia filogenética (sustituciones por sitio)") + ylab("Densidad") +
                                       xlim(0,0.25)+ylim(-30,30)+ 
                                       geom_density(aes(x=`Distancia a S. cerevisiae`, y=..density..),fill="#F8766D")+ 
                                       annotate("text", x = 0.19, y = 25, size = 7, label = "S. cerevisiae", fontface = "italic", color="#F8766D")+
                                       geom_segment(aes(x = 0.03, y = 0, xend = 0.03, yend = 10.5))+
                                       geom_density(aes(x=`Distancia a S. paradoxus`, y=-..density..),fill="#00BFC4")+
                                       annotate("text", x = 0.19, y = -10, size = 6, label = "S. paradoxus", fontface = "italic", color="#00BFC4")+
                                       geom_segment(aes(x = 0.04, y =0, xend = 0.04, yend = -27))

parallel_plot
density_plot
