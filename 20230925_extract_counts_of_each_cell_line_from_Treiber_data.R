# this script is to extract spectrum counts of each cell line from Treiber's data
# 2023/09/25 made

# set function for converting NA to 0
na.zero <-function(x){
  for (i in 1:ncol(x)) {
    x[,i] <-ifelse(is.na(x[,i]),0,x[,i])
  }
  return(x)
  }

# set funtion for calculating pecentage
cal.percent <-function(x,data){
  
  if(x==1){
    total <-apply(data, 1, sum)
    
    for (i in 1:nrow(data)) {
      data[i,] <-data[i,]/total[i]*100
    }}

  if(x==2){
    total <-apply(data, 2, sum)
    
    for (i in 1:ncol(data)) {
      data[,i] <-data[,i]/total[i]*100
    }
  }
  return(data)
  }

# import Treiber's mass spectrum data
# this data is located at ""
setwd("C:/Rdata/Treiber_data")
treiber.counts <-read.csv("Treiber_spectrum_counts_hit.csv",header = T,sep=",",check.names = F,stringsAsFactors = F)

# extract name of cell lines 
cell.lines <-unique(treiber.counts$`cell line`)

# convert NA to 0
treiber.counts <-na.zero(treiber.counts)

data <-treiber.counts[treiber.counts[,2]=="HEPG2",c(8:79)]
miRNA <-cal.percent(1,data)
RBP <-cal.percent(2,data)

# arrange treiber's data for each cell line
library(reshape2)
for (i in 1:length(cell.lines)) {
  # generate file name
  file.name <-paste0("pecentage_of_Treiber_count_in_",cell.lines[i],".txt")
  
  # extract only rows of a certain cell line
  cell.counts <-treiber.counts[treiber.counts$`cell line`==cell.lines[i],]
  
  # calculate percentage of count
  RBP.counts <-cal.percent(1,cell.counts[,c(8:79)])
  miRNA.counts <-cal.percent(2,cell.counts[,c(8:79)])
  
  # bind column of RBP
  RBP.counts[,73] <-cell.counts[,4]
  colnames(RBP.counts)[73] <-"gene symbol"
  miRNA.counts[,73] <-cell.counts[,4]
  colnames(miRNA.counts)[73] <-"gene symbol"
  
  # arrange extracted data
  RBP.counts.melt <-melt(RBP.counts,id.vars="gene symbol")
  miRNA.counts.melt <-melt(miRNA.counts,id.vars="gene symbol")
  
  # merge two percentages
  colnames(RBP.counts.melt) <-c("RBP","miRNA","sum_RBP_per_total")
  colnames(miRNA.counts.melt) <-c("RBP","miRNA","sum_miRNA_per_total")
  merged.melt <-merge(RBP.counts.melt,miRNA.counts.melt,by=c("RBP","miRNA"))
  merged.melt[,5] <-(merged.melt[,3]+merged.melt[,4])/2
  colnames(merged.melt)[5] <-"average"
  
  # output data
  write.table(merged.melt,file.name,sep = "\t",row.names = F,quote = F)
}
