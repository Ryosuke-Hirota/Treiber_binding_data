#this script is to generate table about treiber's spectrum counts by mass spectrometry. 

# set function for converting NA to 0
na.zero <-function(x){
  for (i in 1:ncol(x)) {
    x[,i] <-ifelse(is.na(x[,i]),0,x[,i])
  }
  return(x)
}

# set function for calculating each percentage
cal.percent1 <-function(data){
  
    total <-apply(data, 1, sum)
    
    for (i in 1:nrow(data)) {
      data[i,] <-data[i,]/total[i]*100
    }
    
    mean.percent <-apply(data,2,mean)
    
  return(mean.percent)
}

cal.percent2 <-function(data){
  
  total <-apply(data, 1, sum)
  sum.total <-sum(total)
  
  miRNA.total <-apply(data, 2,sum)
  
  percent <-1:72
  
  for (i in 1:length(miRNA.total)) {
    
   percent[i] <-miRNA.total[i]/sum.total*100
}
  
  return(percent)
}

# import mass spectrum data from Treber's paper
# this data is located at "https://github.com/Ryosuke-Hirota/Treiber_binding_data"
setwd("C:/Rdata/Treiber_data")
treiber.counts <-read.csv("Treiber_spectrum_counts_hit.csv",header = T,stringsAsFactors = F,check.names = F)

# convert NA to 0
treiber.counts <-na.zero(treiber.counts)

# extract names of RBP and miRNA
RBP <-unique(treiber.counts[,4])
miRNA <-colnames(treiber.counts)[8:79]

# make empty data frame and name column and row
sm1 <-as.data.frame(matrix(nrow = 180,ncol = 72))
rownames(sm1) <-RBP
colnames(sm1) <-miRNA
sm2 <-as.data.frame(matrix(nrow = 180,ncol = 72))
rownames(sm2) <-RBP
colnames(sm2) <-miRNA

#sum percentages by RBP
for (i in 1:length(RBP)) {
  RBP.counts <-subset(treiber.counts,treiber.counts[,4]==RBP[i])
  mean.percent <-cal.percent1(RBP.counts[,c(8:79)])
  total.percent <-cal.percent2(RBP.counts[,c(8:79)])
  
  sm1[i,] <-mean.percent
  sm2[i,] <-total.percent
  }

# make empty data frame
final.sm <-as.data.frame(matrix(nrow = 180*72,ncol=5))
colnames(final.sm) <-c("RBP","miRNA","mean_percent","total_percent","average")

# calculate value for heatmap
for (i in 1:nrow(sm1)) {
  for (k in 1:ncol(sm1)) {
   
    final.sm[(i-1)*72+k,1] <-rownames(sm1)[i] 
    final.sm[(i-1)*72+k,2] <-colnames(sm1)[k]
    final.sm[(i-1)*72+k,3] <-sm1[i,k]
    final.sm[(i-1)*72+k,4] <-sm2[i,k]
    final.sm[(i-1)*72+k,5] <-(final.sm[(i-1)*72+k,3]+final.sm[(i-1)*72+k,4])/2
     
  }}


# output value used in heatmap
write.table(final.sm,"treiber_heatmap_score.txt",sep="\t",row.names = F,quote = F)

# arrange data frame for drawing heatmap
final.sm1 <-final.sm
final.sm1[,1] <-factor(final.sm1[,1],levels =rev(RBP) )
final.sm1[,2] <-factor(final.sm1[,2],levels =miRNA )
final.sm1[final.sm1[,5]<3,5] <-0
final.sm1[final.sm1[,5]>=33,5] <-33

# draw heatmap
heatmap1 <- ggplot(final.sm1, aes(x = miRNA, y = RBP, fill = average))+ 
            geom_tile()+
            scale_fill_gradientn("Value", colours = c("white","lightskyblue","blue4"))+
            theme(axis.text.x = element_text(size = 1, colour = "black",angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size = 1, colour = "black"))
heatmap1
# ggsave("Treiber_heatmap.pdf",heatmap1,width = ,height = ,dpi=)