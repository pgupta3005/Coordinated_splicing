#!$rscript

#loading libraries
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tibble)

###########################################
### read files and initialise variables ###
###########################################

b<-read.table("sample_names", sep="\t", header=F)
var_names<-b[,1]
#file_names<-b[,2]
#species_names<-b[,3]

col_txid<-"associated_transcript"
col_counts<-"counts"

for(s in 1:length(var_names))
{
  assign(paste(var_names[s], "_used", sep=""), readRDS(paste(var_names[s], "_used.RDS", sep="")))
  assign(paste("df_tx_", var_names[s], sep=""), readRDS(paste("df_tx_", var_names[s], "_final.RDS",sep="")))
  assign(paste("df_exon_", var_names[s], sep=""), readRDS(paste("df_exon_", var_names[s], "_final.RDS", sep="")))
}

#####################
### Main function ###
#####################

for(s in 1:length(var_names))
{
  x<-get(paste(var_names[s], "_used", sep=""))
  
  df_tx<-get(paste("df_tx_", var_names[s], sep=""))
  df_exon<-get(paste("df_exon_", var_names[s], sep=""))
  
  print(paste("loading", var_names[s]))
  
  list_of_gene<-unique(unlist(df_exon$GENEID))        #same as  list_of_gene2<-unique(unlist(df_tx$GENEID))
  #list_of_new_EXONNAME<-unique(unlist(df_exon$new_EXONNAME))
  list_of_EXONCHUNK<-unique(unlist(df_exon$EXONCHUNK))  
  
  list_of_constitutive_EXONCHUNK<-list()
  
  ### for every gene ###
  
  ##################################################################################
  ### using full-length vs truncated tx info to identify alt vs const exonchunks ###
  ##################################################################################

  for(i in 1:length(list_of_gene))
  {
    df_exon_gene1<-subset(df_exon, unlist(GENEID) %in% list_of_gene[i])
    df_tx_gene1<-subset(df_tx, unlist(GENEID) %in% list_of_gene[i])
    tx_all<-unique(unlist(df_tx_gene1$TXNAME))   
    #there will be transcripts that were not found in our sample, but that is okay
    #however, this contains full and fragment both - so change the txs_all set for each exon within the loop
    
    exonchunk_all<-unique(df_exon_gene1$EXONCHUNK)            
    
    for(j in 1:length(exonchunk_all))
    {
      df_exon_gene1_exon1<-subset(df_exon_gene1, EXONCHUNK==exonchunk_all[j])
      tx_of_exon1<-unique(unlist(df_exon_gene1_exon1$TXNAME))
      
      #define a new total tx set for every exon
      inc_tx<-list()
      
      for(k in 1:length(tx_all))
      {
        
        if(df_tx_gene1[df_tx_gene1$TXNAME==tx_all[k],"TSS"]=="full"){
          inc_tx<-append(inc_tx, tx_all[k])
        } else {
          n_of_start<-which(exonchunk_all==df_tx_gene1[df_tx_gene1$TXNAME==tx_all[k],"TSS"])
          n_of_end<-which(exonchunk_all==df_tx_gene1[df_tx_gene1$TXNAME==tx_all[k],"TTS"])
        
          if((j>=n_of_start & j<=n_of_end) | (j<=n_of_start & j>=n_of_end)){
            inc_tx<-append(inc_tx, tx_all[k])
          }
        }
      }
      
      inc_tx<-unlist(inc_tx)
      
      total_1<-sum(x[x[, col_txid] %in% inc_tx, col_counts])
      total_2<-sum(x[x[, col_txid] %in% tx_of_exon1, col_counts])
      
      if(total_2>(0.95*total_1))
        list_of_constitutive_EXONCHUNK<-append(list_of_constitutive_EXONCHUNK, exonchunk_all[j])
    }
    #print(i)
  }
  
  list_of_constitutive_EXONCHUNK<-unique(unlist(list_of_constitutive_EXONCHUNK))
  saveRDS(list_of_constitutive_EXONCHUNK, paste("constitutive_EXONCHUNK_",var_names[s], ".RDS", sep=""))
  
  g<-unique(gsub("_.*","",list_of_constitutive_EXONCHUNK))
  
  l1<-length(list_of_constitutive_EXONCHUNK)
  l2<-length(g)
  
  ####################################################################
  print(paste("found", l1, "constitutive exonchunks in", l2, "genes"))
  
  list_of_non_constitutive_EXONCHUNK<-list_of_EXONCHUNK[!(list_of_EXONCHUNK %in% list_of_constitutive_EXONCHUNK)]
  saveRDS(list_of_non_constitutive_EXONCHUNK, paste("non_constitutive_EXONCHUNK_",var_names[s],".RDS",sep=""))
  
  ng<-unique(gsub("_.*","",list_of_non_constitutive_EXONCHUNK))
  
  l3<-length(list_of_non_constitutive_EXONCHUNK)
  l4<-length(ng)
  
  ########################################################################
  print(paste("found", l3, "non-constitutive exonchunks in", l4, "genes"))
  
  #roughly 1:3 or 1:2 ratio of non-constitutive to constitutive exonchunks
}
