library(ggplot2)
library(dplyr)
library(ggrepel)
library(presize)

setwd("~/Desktop/Genomic_Sciences/CuartoAño/Mitochondria/NuclearMitCorrelation")
file<-read.csv("./HybsAndNewGroups.csv",header=TRUE)
file<-file[,c(1,6,7,8,9)]
file<-na.omit(file)

A<-function(x) x*100
file[2:5]<-apply(file[,2:5],2,A)

#Scatterplot of mitochondrial vs nuclear proportions
ggplot(file) + ggtitle("Relación entre genoma nuclear vs genoma mitocondrial de S. cerevisiae") + geom_text_repel(aes(x=X.SACE, y=X.mtSACE, color=New.Hybrid.Group.ESG, label=New.Hybrid.Group.ESG), max.overlaps = 50 ,point.size = 5, segment.size=0.25) + ylim(0,100) + xlab("Porcentaje de lecturas mapeadas a genoma nuclear") + ylab("Porcentaje de lecturas mapeadas a genoma mitocondrial") + theme(plot.title = element_text(hjust = 0.8,size=14),legend.position = "None", axis.text = element_text(size=15)) + scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) + scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) + geom_hline(yintercept=c(33,50,66), linetype="dashed", color = "gray", linewidth=0.6) + geom_vline(xintercept=c(33,50,66), linetype="dashed", color = "gray", linewidth=0.6)

#Visual inspection of normality via QQplot and histogram for nuclear genome
qqnorm(file$X.SACE, pch = 1, frame = FALSE)
qqline(file$X.SACE, col = "steelblue", lwd = 2)

x <- file$X.SACE
h<-hist(x, breaks=10, col="red", xlab="SACE",
        main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

#Shapiro test to assess normality
shapiro.test(file$X.SACE)
#Data DO NOT seem to follow normality

########################################################

#Visual inspection of normality via QQplot and histogram for mitochondrial genome
qqnorm(file$X.mtSACE, pch = 1, frame = FALSE)
qqline(file$X.mtSACE, col = "steelblue", lwd = 2)

x <- file$X.mtSACE
h<-hist(x, breaks=20, col="red", xlab="SACE",
        main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

#Shapiro test to asses normality
shapiro.test(file$X.mtSACE)
#Data DO NOT seem to follow normality

##############################################################################################3
#Due to data not following normality, a non-parametric test will be used

file4Corr<-read.csv("./HybsAndNewGroups.csv",header=TRUE)
file4Corr<-file4Corr[,c(1,6,7,8,9)]
file4Corr<-na.omit(file4Corr)

A<-function(x) x*100
file4Corr[2:5]<-apply(file4Corr[,2:5],2,A)

file4Corr<-file4Corr %>% group_by(New.Hybrid.Group.ESG) %>% slice_sample(n=1)


#Spearman Correlation Test
res<-cor.test(file4Corr$X.SACE,file4Corr$X.mtSACE, method="spearman")
cat("Correlation coefficient:",res$estimate, "with a p-value of",  res$p.value )

#Calculate confidence interval
prec_cor(r=0.68, n = 17, conf.level = 0.95, method = 'spearman')

ggplot(file4Corr) + ggtitle("Porcentaje de reads mapeadas a S. cerevisiae") + annotate (geom = "text", x= 55 , y = 5, label="rho = 0.68, int. confianza = (0.24, 0.88), p-valor = 0.00229",size=5) + geom_text_repel(aes(x=X.SACE, y=X.mtSACE, color=New.Hybrid.Group.ESG,label=New.Hybrid.Group.ESG),max.overlaps = 30 ,point.size = 3) + ylim(0,100) + xlab("Genoma Nuclear") + ylab("Genoma Mitocondrial") + theme(plot.title = element_text(hjust = 0.5),legend.position = "None") + scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) + scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))


