library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Missing csv file.", call.=FALSE)
} else {
  path =  args[1]
  Q = args[2]
}

path_out <- "/mnt/Timina/lmorales/Public/ymez/data/figures/04_coverage/mtReferencePercentage"
file = read.csv(path,header=T)

cols=ncol(file)
file=file[,-(cols)] #Remove the Total column since it will not be used for plotting
N = 8
chunks <- split(1:length(file$Sample), ceiling(seq_along(file$Sample)/N)) #Divide the total number of samples into N groups. This will be used to plot N samples per graph.
n_graphs <- ceiling(length(file$Sample)/N)

for(i in 2:(dim(file)[2]-1)){ #Compute the operation to obtain the percentage of reads that correspond to SACE and SAPA, respectively.
  for(j in 1:dim(file)[1]){
    file[j,i] <- file[j,i]/file[j,dim(file)[2]]
  }
}

svgimg <- function(n,graph){ #Function to save the plot as .png
  name <- paste(path_out,"/percentage_reference_mitDNA_",Q,"_",n,".png",sep="")
  if(Q == "C_p"){name <- paste(path_out,"/percentage_reference_mitDNA_",n,".png",sep="")}
  print(name)
  png(name, width = 17, height = 12, units = "cm",bg = "white", res = 300)
  print(graph)
  dev.off()
}

for(i in 1:n_graphs){  #for cycle to create the n number of plots
  file=file[,-(cols-1)]
  mltp<-melt(file[chunks[[i]],],id.var="Sample")
  colnames(mltp) <- c("Sample","Reference","FractionOfMappedReads")
  p <- ggplot(mltp,aes(x = Sample,y = FractionOfMappedReads,fill = Reference)) +
    geom_bar(stat = "identity")+
    theme(text = element_text(size=10), axis.text.x = element_text(angle =90, hjust = 1))
  svgimg(i,p)
}
