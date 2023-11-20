#!$rscript

#loading libraries
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tibble)
library(parallel)

#######################################
### read files and take user inputs ###
#######################################

## There should be a tab separated file named 'sample_names' with
## the first col is a short name for the sample
## the second col is the filename of the counts 
## the third col is the name of the reference gtf (genome annotation file) with extension

a<-read.table("sample_names", sep="\t", header=F)
species_names<-a[,1]     ##define species name
reference_names<-a[,3]   #define filenames for reference

#Accept number of threads from user as the first argument after the scriptname

num_threads<-as.integer(detectCores()/4) ##default value

tmp<-commandArgs(trailingOnly = TRUE)
if(length(tmp)==1)
  num_threads<-as.numeric(tmp[1])

#print(length(tmp))
#print(tmp)

#####################################
### PART A: create txdb object(s) ###
#####################################

#save txdb objects
for(i in 1:length(species_names))    #length same as that of reference name
{
  x<-makeTxDbFromGFF(reference_names[i])
  assign(paste0("txdb_", species_names[i]), x)
  saveDb(x, paste0("txdb_", species_names[i]))
}

#to load txdb objects
#for(i in species_names)
#  assign(paste0("txdb_", i), loadDb(paste0("txdb_", i)))


####################################################################
### PART B: create exon-wise and trasncript-wise granges objects ###
####################################################################

for(i in species_names)
{
  x<-get(paste0("txdb_", i))
  
  y<-exons(x, columns=c("EXONNAME","TXNAME","GENEID"))
  tmp<-as.data.frame(y)
  y$new_EXONNAME<-paste(tmp$GENEID, tmp$start, tmp$end, sep="-")
  assign(paste0("exon_", i), y)
  saveRDS(y, paste0("exon_", i, ".RDS"))
  
  y<-transcripts(x, columns=c("TXNAME","EXONNAME","GENEID"))
  y$temp_EXONNAME<-paste(y$GENEID, y$EXONNAME, sep="-")
  assign(paste0("tx_", i), y)
}

##############################################################
### PART C: change to new_EXONNAME in both granges objects ###
##############################################################

func_ID_converter=function(x){
  
  l<-unlist(x)
  n<-strsplit(l,"-")
  g<-n[[1]][1]             #[[anything]][1] is geneid
  
  for(j in 1:length(n))
  {
    e<-df_exon_name[df_exon_name$EXONNAME==n[[j]][2],"new_EXONNAME"]
    
    if(length(e)==1)
      l[j]<-e
    else if(length(e)>1)
      l[j]<-df_exon_name[(df_exon_name$EXONNAME==n[[j]][2] & unlist(df_exon_name$GENEID)==g),"new_EXONNAME"]
  }
  
  return(l)
  
}

func_add_newEXONNAME=function(tx_name, exon_name, core_num){
  
  df_exon_name<<-as.data.frame(exon_name)     #important assign function
  df_tx_name<<-as.data.frame(tx_name)         #important assign function
  newID<-mclapply(df_tx_name[,"temp_EXONNAME"], func_ID_converter, mc.cores=core_num)
  tx_name$new_EXONNAME<-CharacterList(newID)
  return(tx_name)
  
}

for(i in species_names)
{
  tx<-get(paste0("tx_", i))
  exon<-get(paste0("exon_", i))
  out<-func_add_newEXONNAME(tx, exon, num_threads)
  print(i)
  assign(paste0("tx_", i), out)
  saveRDS(out, paste0("tx_", i, ".RDS"))
}
