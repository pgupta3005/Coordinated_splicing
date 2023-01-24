#!/apps/codes/R/4.1.1/bin/Rscript

#loading libraries (if any)
library(dplyr)

#########################
### to format samples ###
#########################

#sample_names is defined in step2_exonchunks.R

b<-read.table("sample_names", sep="\t", header=F)
var_names<-b[,1]
file_names<-b[,2]
#species_names<-b[,3]

#number of replicates per sampletype
replicates <- c(5,5,2,2,12)

#specify column with structural category info
strc_cat <- "structural_category"
#specify txid column name
col_txid<-"associated_transcript"

#reading .csv or .txt files 
for(s in 1:length(var_names))
{
  x <- read.table(paste(var_names[s], ".csv", sep=""), sep=",", header=TRUE)
  x <- x[x[,strc_cat]== "full-splice_match",]
  assign(var_names[s], x)
}

#####################
### Main function ###
#####################

for(s in 1:length(var_names))
{
  rep<-replicates[s]
  rep_pre<-0

  x<-get(var_names[s])
  newdf<-x
  
  #in newdf replace NAs with 0
  newdf[is.na(newdf)]<-0
  
  #and make a total or counts ka column
  if(rep>1)
    newdf$counts<-rowSums(newdf[,4:(3+rep)])
  if(rep==1)
    newdf$counts<-newdf[,4]
  
  #assign(var_names[s], newdf)
  write.table(newdf, file_names[s], sep = "\t", quote=F)
  
}

#saveRDS(var_names, "var_names.RDS")
#saveRDS(var_names2, "var_names2.RDS")

print("done")
