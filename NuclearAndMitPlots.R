library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(gridExtra)
setwd("~/Desktop/Genomic_Sciences/CuartoAño/Mitochondria/plots/HybGroups/")
file<-read.csv("./Hybs_AllPercentagesAndGroup.csv")
file<-file %>% drop_na()
path_out <- "/home/esantiago/Desktop/Genomic_Sciences/CuartoAño/Mitochondria/plots/NuclearAndMit"
groups=split(file,file$ID)
groups<- lapply(seq_along(groups), function(x) as.data.frame(groups[[x]])[,1:5]) 
#grouplist<-c("1","2a","2b","3","4a","4b","5","6a","6b","7","8","9a","9b","Undefined")

#write.csv(file, paste0(path_out,"/Hybs_AllPercentagesAndGroup.csv"), row.names=FALSE)

for(i in 1:length(groups)){  #for cycle to create the n number of plots
  currentgroup=as.data.frame(groups[i])
  currentgroup=currentgroup %>% rename(mtSACE = X.mtSACE , mtSAPA = X.mtSAPA, SACE = X.SACE, SAPA = X.SAPA)
  currentgroup_Nuclear=currentgroup[,1:3]
  currentgroup_Mit=currentgroup[,c(1,4,5)]
  frmtdcurrentgroup_Nuclear<-melt(currentgroup_Nuclear,var.id="ID")
  frmtdcurrentgroup_Mit<-melt(currentgroup_Mit,var.id="ID")
  colnames(frmtdcurrentgroup_Nuclear)<-c("Sample","Reference","FractionOfMappedReads")
  colnames(frmtdcurrentgroup_Mit)<-c("Sample","Reference","FractionOfMappedReads")
  
  frmtdcurrentgroup_Nuclear<-frmtdcurrentgroup_Nuclear %>% left_join(select(file,ID,Group),by=c("Sample"="ID"))
  frmtdcurrentgroup_Mit<-frmtdcurrentgroup_Mit %>% left_join(select(file,ID,Group),by=c("Sample"="ID"))
  hsize=3
  frmtdcurrentgroup_Mit <- frmtdcurrentgroup_Mit%>%
    mutate(gap=hsize)
  
  p <- ggplot(frmtdcurrentgroup_Nuclear,aes(x = Sample,y = FractionOfMappedReads,fill=Reference)) +
    geom_bar(stat = "identity", width = 0.3)+
    #geom_text(aes(size=9,label=sprintf(FractionOfMappedReads,fmt='%#.2f')),position=position_stack(vjust = 0.5))+
    theme(axis.text.x = element_blank(), axis.text.y = element_blank() , axis.title = element_blank()) + theme(legend.position = "none") + theme_void()
  #svgimg(i,p)   
  
  mit<-ggplot(frmtdcurrentgroup_Mit,aes(x = gap,y = FractionOfMappedReads,fill=Reference)) +
    geom_col()+
    coord_polar(theta="y") + xlim(c(0.2,hsize+0.5)) +
    #geom_text(aes(size=10,label=sprintf(FractionOfMappedReads,fmt='%#.2f')),position=position_stack(vjust = 0.5))+
    labs(title=paste("            Group: ",frmtdcurrentgroup_Mit$Group,"_",frmtdcurrentgroup_Mit$Sample,sep=""))  + theme_void() + theme(legend.position = "none") 
  
  graph<-arrangeGrob(mit,p,ncol=2)
  name <- paste(path_out,"/NuclearAndMitDNA_",currentgroup_Mit$ID,".png",sep="")
  ggsave(name,plot = graph ,width = 14, height = 12, units = "cm",bg = "white", dpi = 300)
  print(name)
}
