#!$rscript

#loading libraries
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tibble)
library(parallel)

###########################################
### read files and initialize variables ###
###########################################

b<-read.table("sample_names", sep="\t", header=F)
var_names<-b[,1]
#file_names<-b[,2]
#species_names<-b[,3]

num_threads<-as.integer(detectCores()/8) ##default

tmp<-commandArgs(trailingOnly = TRUE)
if(length(tmp)==1)
  num_threads<-as.numeric(tmp[1])

for(s in 1:length(var_names))
{
  assign(paste(var_names[s], "_used", sep=""), readRDS(paste(var_names[s], "_used.RDS", sep="")))
  assign(paste("df_exon_", var_names[s], sep=""), readRDS(paste("df_exon_", var_names[s], "_final.RDS", sep="")))
  assign(paste("constitutive_EXONCHUNK_", var_names[s], sep=""), readRDS(paste("constitutive_EXONCHUNK_", var_names[s], ".RDS",sep="")))
  assign(paste("non_constitutive_EXONCHUNK_", var_names[s], sep=""), readRDS(paste("non_constitutive_EXONCHUNK_", var_names[s], ".RDS",sep="")))
}

################################################################################
### define function to get relative and absolute position of exonchunk pairs ###
################################################################################

RelPos=function(x){
  temp=strsplit(x,"x")
  gene=temp[[1]][1]
  exon1<-temp[[1]][2]
  exon2<-temp[[1]][3]
  
  #should we use exon_full so that any alternative exon that is entirely absent from our dataset is also considered while calculating position? 
  #very tough to accommodate here :(
  
  tmp<-subset(df_exon, GENEID==gene)
  list_of_n<-as.integer(gsub(".*_","",tmp$EXONCHUNK))

  n1<-as.integer(gsub(".*_","",exon1))
  n2<-as.integer(gsub(".*_","",exon2))
  
  u=n2-n1
  
  ### category var for relative position ###

  if(u==1) {
    category<-"adjacent"
  } else if(u>1) {
      c1<-which(df_exon$EXONCHUNK==exon1)
      c2<-which(df_exon$EXONCHUNK==exon2)
      sub_exon<-df_exon[c1:c2,]
    
      c<-constitutive_EXONCHUNK[constitutive_EXONCHUNK %in% sub_exon$EXONCHUNK]
      c<-length(unique(unlist(c)))
      if(c==0)
        category<-"distant1"    #no constitutive exon 
      if(c!=0)
        category<-"distant2"    #one or more constitutive exon
  }
 
  ### exon_pos var for absolute position ###

  if(n1==1){
    exon_pos1<-"start"
  } else {
    exon_pos1<-"internal"
  }

  if(n2==length(unique(list_of_n))){
    exon_pos2<-"end"
  } else {
    exon_pos2<-"internal"
  }

  return(paste(category, exon_pos1, exon_pos2, sep="_"))
}

#####################
### Main function ###
#####################

for(s in 1:length(var_names))
{
  print(paste("loading", var_names[s]))
  
  df_exon<-get(paste("df_exon_", var_names[s], sep=""))
  constitutive_EXONCHUNK<-get(paste("constitutive_EXONCHUNK_", var_names[s], sep=""))
  non_constitutive_EXONCHUNK<-get(paste("non_constitutive_EXONCHUNK_", var_names[s], sep=""))
  df_exon_non_constitutive<-subset(df_exon, EXONCHUNK %in% non_constitutive_EXONCHUNK)
  gene_non<-unique(unlist(df_exon_non_constitutive$GENEID))
  
  ############################################
  ### Part A: Creating pairs of exonchunks ###
  ############################################
  
  table1<-matrix(nrow=0, ncol=3)
  colnames(table1)<-c("GENEID","EXONCHUNK1","EXONCHUNK2")
  
  for(i in 1:length(gene_non))
  {
    #the gene must have more than one non_constitutive exon (in df_exon) 
    
    a<-subset(df_exon, unlist(GENEID) %in% gene_non[i] & EXONCHUNK %in% non_constitutive_EXONCHUNK)

    if(length(unique(unlist(a$EXONCHUNK)))>1)
    {
      x<-combn(unique(unlist(a$EXONCHUNK)), 2)
      y<-replicate(ncol(x),gene_non[i])
      temp<-cbind(y,t(x))
      table1<-rbind(table1, temp)
    }
  }
  
  table1<-as.data.frame(table1)
  table1$combine<-paste(table1$GENEID, table1$EXONCHUNK1, table1$EXONCHUNK2, sep="x")
  
  print(paste("made", nrow(table1), "exon pairs from", length(unique(table1$GENEID)), "genes"))
  
  ##############################################################################
  ### Part B: Determining relative and absolute positions of exonchunk pairs ###
  ##############################################################################

  table1<-table1 %>% add_column(relpos=NA, abspos1=NA, abspos2=NA)
  
  output<-mclapply(table1$combine, RelPos, mc.cores=num_threads)
  output<-strsplit(unlist(output), "_")
  
  table1$relpos<-sapply(output, "[[", 1)
  table1$abspos1<-sapply(output, "[[", 2)
  table1$abspos2<-sapply(output, "[[", 3)
  
  print("assigned relative and absolute positions for all exon pairs")
  
  saveRDS(table1, paste("table_", var_names[s], ".RDS", sep=""))
  
}
