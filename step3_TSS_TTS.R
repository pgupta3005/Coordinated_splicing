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

for(s in 1:length(var_names))
{
  assign(paste("df_tx_", var_names[s], sep=""), readRDS(paste("df_tx_", var_names[s], ".RDS",sep="")))
  assign(paste("df_exon_", var_names[s], sep=""), readRDS(paste("df_exon_", var_names[s], ".RDS", sep="")))
}

#####################
### Main function ###
#####################

for(s in 1:length(var_names))
{
  df_tx<-get(paste("df_tx_", var_names[s], sep=""))
  df_exon<-get(paste("df_exon_", var_names[s], sep=""))
  
  print(paste("working with", var_names[s]))
  
  all_genes<-unique(unlist(df_tx$GENEID))
  #all_genes<-unique(unlist(df_exon$GENEID))      #should be same
  
  ############################################################################################
  ### Part A: Mark which FSM transcripts are possibly a fragment of another FSM transcript ###
  ############################################################################################

  print("identifying FSM transcripts that may be possible fragments of another FSM")
  
  short_tx<-list()
  
  for(i in 1:length(all_genes))
  {
    tx_all<-unique(df_tx[unlist(df_tx$GENEID)==all_genes[i],"TXNAME"])
    #print(length(tx_all))
    #there is a chance that some FSM got removed due to all NA exons
    #after this fragment calculation, we'll remove genes with only one full FSM transcript
    
    if(length(tx_all)>1)
    {
      l1<-as.data.frame(t(combn(unique(tx_all), 2)))
    
      for(j in 1:nrow(l1))
      {
        #compare longer %in% shorter tx wrt EXONCHUNK
        e1<-unlist(df_tx[df_tx$TXNAME==l1[j,1],"EXONCHUNK"])
        e2<-unlist(df_tx[df_tx$TXNAME==l1[j,2],"EXONCHUNK"])
        
        #print(e1)
        #print(e2)
        
        if(length(e1)>=length(e2)) {
          if(all(diff(which(e1 %in% e2))==1) && length(intersect(e1,e2))==length(e2))
            short_tx<-append(short_tx,l1[j,2])
        } else {
          if(all(diff(which(e2 %in% e1))==1) && length(intersect(e1,e2))==length(e1))
            short_tx<-append(short_tx,l1[j,1])
        }
      }
    }
  }
  
  short_tx<-unique(unlist(short_tx))
  t1<-length(unique(df_tx$TXNAME))
  
  print(paste("found", length(short_tx), "transcripts to be possible fragments out of", t1, "transcripts"))
  
  #############################################################
  ### Part B: demarcate by adding two columns : TSS and TTS ###
  #############################################################

  #for full transcripts, add "full" for both columns
  #for fragmented transcripts, write the first and last exonchunk names, resp.
  
  df_tx<-df_tx %>% add_column(TSS=NA, TTS=NA)
  
  for(i in 1:nrow(df_tx))
  {
    tx_name<-df_tx[i,"TXNAME"]
    if(tx_name %in% short_tx)
    {
      e<-unlist(df_tx[i,"EXONCHUNK"])
      df_tx[i,"TSS"]<-e[1]
      df_tx[i,"TTS"]<-e[length(e)]
      
      if(e[1]=="NA" | e[length(e)]=="NA")
      {
        df_tx[i,"TSS"]<-"full"
        df_tx[i,"TTS"]<-"full"
      }
      
    } else {
      df_tx[i,"TSS"]<-"full"
      df_tx[i,"TTS"]<-"full"
    }
  }
  
  ####################################################################
  ### Part C: Subset to remove genes with only one full transcript ###
  ####################################################################

  #this removal must be done for both df_tx and df_exon --> then overwrite old files
  
  df_tx_tmp<-subset(df_tx, TSS=="full")
  tb<-as.data.frame(table(unlist(df_tx_tmp$GENEID)))
  g1<-nrow(tb)
  
  tb<-subset(tb, Freq!=1)
  g2<-nrow(tb)
  new_genes<-tb$Var1
  df_tx<-subset(df_tx, unlist(GENEID) %in% new_genes)
  t2<-length(unique(df_tx$TXNAME))
  df_exon<-subset(df_exon, unlist(GENEID) %in% new_genes)
  
  print(paste("retaining", g2, "genes out of", g1, "and", t2, "transcripts out of", t1))
  
  #losing roughly 15% genes as they contain only one FSM (plus any fragment(s) or no fragments)
  
  saveRDS(df_tx, paste("df_tx_", var_names[s], "_final.RDS", sep=""))
  saveRDS(df_exon, paste("df_exon_", var_names[s], "_final.RDS", sep=""))
  
}


