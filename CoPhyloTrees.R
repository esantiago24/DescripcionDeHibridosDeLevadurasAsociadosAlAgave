library(ape)
library(readxl)
library(phytools)
library(devEMF)
setwd("~/Desktop/Genomic_Sciences/CuartoAño/Trees/AllHybTrees/")

#Read Parental trees
SACE_tree<-read.tree("./RAxML_bipartitionsBranchLabels.Matrix_SNPs_SACE_from_CONC_gt_onlySNPs_filtered_missing_10_plus_SRR800827_recode_min4_AllHybSites_noN.tree")
SAPA_tree<-read.tree("./RAxML_bipartitionsBranchLabels.Matrix_SNPs_SAPA_from_CONC_gt_onlySNPs_filtered_missing_10_plus_YMX506H02_recode_min4_AllHybSites_noN.tree")
AllParentals<-read_excel("../20230516_HybridParentals_SACEnSAPA_noN.xlsx", sheet="AllParentals")
SACE_Parentals<-read.csv("../metadata/SampleSheet_SACE199_SaceForHybTrees.csv")
SACE_Parentals<-SACE_Parentals[,c("Nombre_Fastq","Project","Country","Geographic_location")]
SAPA_parentals<-read.csv("../metadata/SampleSheet_SAPA88_SpBandSpD.csv")
SAPA_parentals<-SAPA_parentals[,c("Nombre_Fastq","Project")]

association<-cbind(AllParentals$Closest_SACE_ID,AllParentals$Closest_SAPA_ID) #Create a dataframe to indicate the associations between the closest parentals of a given hybrid
Cophylo<-cophylo(SACE_tree,SAPA_tree,assoc = association,rotate=TRUE)
#str(Cophylo$trees)

########################This chunk was used to label the tips of the trees according to the labels in the samplesheet. It can be omitted##############################
#SACEtips2plot_YMX<-which(Cophylo$trees[[1]]$tip.label %in% SACE_Parentals[SACE_Parentals$Project == "YMX",]$Nombre_Fastq)
#SACEtips2plot_FG<-which(Cophylo$trees[[1]]$tip.label %in% SACE_Parentals[SACE_Parentals$Project == "JS-FG",]$Nombre_Fastq)
#SACEtips2plot_TEQ<-which(Cophylo$trees[[1]]$tip.label %in% SACE_Parentals[SACE_Parentals$Project == "Lachance1995",]$Nombre_Fastq)
#SACEtips2plot_AMIX<-which(Cophylo$trees[[1]]$tip.label %in% SACE_Parentals[SACE_Parentals$Project == "AmMix2",]$Nombre_Fastq)
#SAPAtips2plot_YMX<-which(Cophylo$trees[[2]]$tip.label %in% SAPA_parentals[SAPA_parentals$Project == "SACEandSAPA_Trees_NatEnvMezcal",]$Nombre_Fastq)
#SAPAtips2plot_SpB<-which(Cophylo$trees[[2]]$tip.label %in% SAPA_parentals[SAPA_parentals$Project == "SpB",]$Nombre_Fastq)
#SAPAtips2plot_SpD<-which(Cophylo$trees[[2]]$tip.label %in% SAPA_parentals[SAPA_parentals$Project == "SpD1",]$Nombre_Fastq)

#pdf(file="CoPhylo_Rotated_wTips_nLabels.pdf")
#plot(Cophylo,fsize=c(.4,.5),type="phylogram",show.tip.label=FALSE)
#tiplabels.cophylo(tip=SACEtips2plot_YMX,frame="none",bg="grey",cex=.7,col="red")
#tiplabels.cophylo(tip=SACEtips2plot_FG,frame="none",bg="grey",cex=.7,col="blue")
#tiplabels.cophylo(tip=SACEtips2plot_TEQ,frame="none",bg="grey",cex=.7,col="yellow")
#tiplabels.cophylo(tip=SACEtips2plot_TEQ,frame="none",bg="grey",cex=.7,col="green")


#tiplabels.cophylo(tip=SAPAtips2plot_YMX,frame="none",bg="grey",cex=.7,col="purple",which = "right")
#tiplabels.cophylo(tip=SAPAtips2plot_SpB,frame="none",bg="grey",cex=.7,col="orange",which = "right")

#dev.off()
#####################################################################################################################################################################
pdf(file="CoPhylo_Rotated_woTips_nLabels.pdf")
plot(Cophylo,fsize=c(.4,.5),type="phylogram",show.tip.label=FALSE)
dev.off()
