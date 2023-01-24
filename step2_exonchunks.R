#!shebang

#loading libraries
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tibble)
library(parallel)

###########################################
### read files and initialize variables ###
###########################################

## There should be a tab separated file named 'sample_names' with 
#the first col is a short name for the sample 
#the second col is the name of the tab separated counts file
#the third col is the short name of the corresponding reference

b<-read.table("sample_names", sep="\t", header=F)
var_names<-b[,1]
file_names<-b[,2]
species_names<-b[,3]   #these until now are repeated, but will be made unique for every sample

num_threads<-as.integer(detectCores()/8) ##default

tmp<-commandArgs(trailingOnly = TRUE)
if(length(tmp)==1)
  num_threads<-as.numeric(tmp[1])

for(s in 1:length(var_names))
{
  assign(var_names[s], read.table(file_names[s], header=T, sep="\t"))
  tmp<-readRDS(paste0("tx_", species_names[s], ".RDS"))
  assign(paste0("tx_",var_names[s]), tmp)
  assign(paste0("df_tx_",var_names[s]), as.data.frame(tmp))
  tmp<- readRDS(paste0("exon_", species_names[s], ".RDS"))
  assign(paste0("exon_",var_names[s]),tmp)
  assign(paste0("df_exon_",var_names[s]), as.data.frame(tmp))
}

rm(tmp)

col_txid<-"associated_transcript"
col_counts<-"counts"

#####################################################
### Define function to add EXONCHUNK info to df_tx###
#####################################################

#this will run once for every sample - but it's quick as the number of transcripts per sample is limited

IDconverterCHUNK=function(y){
  l<-unlist(y)
  n<-strsplit(l,"-")
  g<-n[[1]][1]    #[[anything]][1] is geneid
  for(j in 1:length(n))
  {
    e<-df_exon[(df_exon$EXONNAME==n[[j]][2] & unlist(df_exon$GENEID)==g),"EXONCHUNK"] 
    #next if else is to accommodate certain exons that were removed via 5% cutoff
    if(length(e)!=0){
      l[j]<-unique(e)
    } else {
      l[j]<-"NA"
    }
  }
  return(l)
}


#####################
### Main function ###
#####################

for(s in 1:length(var_names))
{
  x<-get(var_names[s])
  
  print(paste(var_names[s],"loaded"))
  
  df_tx<-get(paste0("df_tx_",var_names[s]))
  df_exon<-get(paste0("df_exon_",var_names[s]))
  
  list_of_tx<-unique(x[,col_txid])
  print(paste("total transcripts are", length(list_of_tx)))
  
  df_tx<-subset(df_tx, TXNAME %in% list_of_tx)
  print(paste("number of rows in df_tx are", nrow(df_tx)))
  
  list_of_gene<-unique(unlist(df_tx$GENEID))
  print(paste("total genes are", length(list_of_gene)))
  
  ##################################################################################################
  ### Part A: remove genes with single transcripts and where less than 10 molecules are obtained ###
  ##################################################################################################
  
  print("removing genes with single transcripts and where less than 10 molecules are obtained")
  
  tb<-as.data.frame(table(unlist(df_tx$GENEID)))
  tb<-tb %>% add_column(gene_tot=NA)

  for(i in 1:nrow(tb))
  {
    g<-tb[i,1]
    tx<-df_tx[unlist(df_tx$GENEID)==g, "TXNAME"]
    c<-sum(x[x[,col_txid] %in% tx, col_counts])
    tb[i,3]<-c
  }

  tb<-subset(tb, Freq>1 & gene_tot>=10)
  list_of_gene<-as.character(tb$Var1)  
  
  print(paste("new total genes are", length(list_of_gene), sep=" "))
  
  df_tx<-subset(df_tx, unlist(GENEID) %in% list_of_gene)

  list_of_tx<-unique(df_tx$TXNAME)
  print(paste("new transcripts are", length(list_of_tx), sep=" "))
  
  #starting with new_EXONNAME directly
  list_of_exons<-unique(unlist(df_tx$new_EXONNAME))
  print(paste("total newexonname are", length(list_of_exons), sep=" "))
  
  df_exon<-subset(df_exon, new_EXONNAME %in% list_of_exons)
  print(paste("nrow in df_exon are", nrow(df_exon), sep=" "))
  
  df_exon<-df_exon %>% add_column(tot_exon=NA, tot_gene=NA)
  
  #####################################################################################
  ### Part B: remove exons which are present in less upto 5% molecules of that gene ###
  #####################################################################################
  
  print(paste("removing exons which are present in only upto", "5%", "molecules - very slow"))
  
  for(i in 1:length(list_of_exons))
  {
    e<-list_of_exons[i]
    ind<-which(df_exon$new_EXONNAME==e)
    tx<-unique(unlist(df_exon[ind, "TXNAME"]))
    gene<-unique(unlist(df_exon[ind, "GENEID"]))
    c_exon<-sum(x[x[,col_txid] %in% tx, col_counts])
    tx_all<-unique(unlist(df_exon[unlist(df_exon[,"GENEID"])==gene, "TXNAME"]))
    c_gene<-sum(x[x[,col_txid] %in% tx_all, col_counts])
    df_exon[ind,"tot_exon"]<-c_exon
    df_exon[ind,"tot_gene"]<-c_gene
  }
  
  #maybe can save here to help troubleshoot later
  
  df_exon<-subset(df_exon, tot_exon>=(0.05*tot_gene))    #this cutoff helps removes the longest weirdest exons
  
  print(paste("new nrow in df_exon are", nrow(df_exon), sep=" "))
  
  list_of_exons<-unique(unlist(df_exon$new_EXONNAME))
  print(paste("new newexonname are", length(list_of_exons), sep=" "))
  
  list_of_gene<-unique(unlist(df_exon$GENEID))
  print(paste("new genes are", length(list_of_gene), sep=" "))
  
  ################################################
  ### Part C: removing mono- and di-exon genes ###
  ################################################
  
  print("removing mono- and di-exon genes")
  
  tg<-as.data.frame(table(unlist(df_exon$GENEID)))
  tg<-subset(tg, Freq>2)

  list_of_gene<-as.character(tg$Var1)
  print(paste("new genes are", length(list_of_gene), sep=" "))
  
  df_exon<-subset(df_exon, unlist(GENEID) %in% list_of_gene)
  print(paste("new nrow in df_exon are", nrow(df_exon), sep=" "))
  
  list_of_exons<-unique(unlist(df_exon$new_EXONNAME))
  print(paste("new newexonname are", length(list_of_exons), sep=" "))
  
  ##############################################################
  ### Part D: defining EXONCHUNKS and adding info to df_exon ###
  ##############################################################

  print("defining exonchunks based on this set of new_EXONNAME - slow")

  #saveRDS(df_exon, paste("df_exon_", var_names[s], "_counts.RDS", sep=""))
  #this is just intermediate save if need be

  #define EXONCHUNK based on common 3' and/or 5' site(s)

  df_exon<- df_exon %>% add_column(EXONCHUNK=NA)

  for(i in 1:length(list_of_gene))
  {
    row_num<-which(unlist(df_exon$GENEID)==list_of_gene[i])
    n<-length(row_num)
    blacklist<-list()
    count<-1
    #print(i)
    start_list<-df_exon[row_num,"start"]
    end_list<-df_exon[row_num,"end"]
  
    for(j in 1:n)
    {
      #k<-row_num[j]
      if(!(j %in% blacklist))
      {
        start_j<-start_list[j]
        end_j<-end_list[j]
      
        df_exon[row_num[j], "EXONCHUNK"]<-paste(list_of_gene[i],count,sep="_")   #here
      
        count=count+1
      
        if(j<n)
        {
          for(l in (j+1):n)
          {
            #m<-row_num[l]
            start_l<-start_list[l]
            end_l<-end_list[l]
            for(p in 1:(l-1))
            {
              start_p<-start_list[p]
              end_p<-end_list[p]
            
              if(start_l==start_p | end_l==end_p)
              {
              df_exon[row_num[l],"EXONCHUNK"]<-df_exon[row_num[p],"EXONCHUNK"]
              blacklist<-append(blacklist,l)
              break
              }
            }
          }
        }
      }
    }
  }

  list_of_EXONCHUNK<-unique(unlist(df_exon$EXONCHUNK))  
  print(paste("length of EXONCHUNK is", length(list_of_EXONCHUNK)))
 
  #####################################################
  ### Part E: remove genes with upto two EXONCHUNKS ###
  #####################################################

  blk_list_genes<-list()
  
  for(i in df_exon$GENEID)
  {
    exonchunk_count<-length(unique(subset(df_exon, GENEID==i)$EXONCHUNK))
    if(exonchunk_count<=2)
      blk_list_genes<-append(blk_list_genes, i)
  }

  blk_list_genes<-unlist(blk_list_genes)
  list_of_gene<-list_of_gene[!(list_of_gene %in% blk_list_genes)]

  print(paste("new genes are", length(list_of_gene)))

  df_exon<-subset(df_exon, unlist(GENEID) %in% list_of_gene)
  print(paste("new nrow in df_exon are", nrow(df_exon)))

  list_of_exonchunk<-unique(unlist(df_exon$EXONCHUNK))
  print(paste("new exonchunks are", length(list_of_exonchunk)))

  saveRDS(df_exon, paste0("df_exon_", var_names[s], ".RDS"))
  print(paste0("df_exon_",var_names[s],".RDS saved"))
  
  ##############################################
  ### Part F: adding EXONCHUNK info to df_tx ###
  ##############################################

  df_tx<-subset(df_tx, unlist(GENEID) %in% list_of_gene)
 
  print("adding exonchunk info to df_tx")
  newID<-mclapply(df_tx[,"temp_EXONNAME"], IDconverterCHUNK, mc.cores=num_threads)   #quite fast
  df_tx$EXONCHUNK<-newID

  #################################################################
  ### Part G: removing transcripts with only "NA" exonchunks(s) ###
  #################################################################
  
  blk_ind<-list()
  
  for(i in 1:nrow(df_tx))
  {
    l<-unlist(df_tx[i,"EXONCHUNK"])
    if(sum(l=="NA")==length(l))
      blk_ind<-append(blk_ind,i)
  }

  if(length(unlist(blk_ind))!=0)
    df_tx<-df_tx[-unlist(blk_ind),]
  
  #removing "NA" from exonchunk columns
  
  #for(i in 1:nrow(df_tx))
  #{
    #l<-unlist(df_tx[i,"EXONCHUNK"])
    #n<-which(l=="NA")
    #if(length(n)!=0)
      #df_tx[i,"EXONCHUNK"]<-list(list(l[-n]))
  #}
  
  #######################################################
  ### Part H: removing genes with 1 non-NA transcript ###
  #######################################################

  blk_list_genes<-list()

  for(i in df_tx$GENEID)
  {
    tx_count<-length(unique(subset(df_tx, GENEID==i)$TXNAME))
    if(tx_count==1)
      blk_list_genes<-append(blk_list_genes, i)
  }

  blk_list_genes<-unlist(blk_list_genes)
  list_of_gene<-list_of_gene[!(list_of_gene %in% blk_list_genes)]

  print(paste("new genes are", length(list_of_gene)))

  
  ################################################################################
  ### Part I: removing genes same number of exonchunks in transcripts and gene ###
  ################################################################################
  
  blk_list_genes<-list()
  
  for(i in list_of_gene)
  {
    exon_df_sub<-subset(df_exon, GENEID == i)
    tx_df_sub<-subset(df_tx, GENEID == i)
    
    exonchunk_all<-unique(exon_df_sub$EXONCHUNK)
    exonchunk_txs_new<-unlist(tx_df_sub$EXONCHUNK)
    
    if(length(exonchunk_all)==length(exonchunk_txs_new)) 
      blk_list_genes<-append(blk_list_genes, i)
  }
  
  blk_list_genes<-unlist(blk_list_genes)
  list_of_gene<-list_of_gene[!(list_of_gene %in% blk_list_genes)]
  
  print(paste("new genes are", length(list_of_gene)))
  
  df_tx<-subset(df_tx, unlist(GENEID) %in% list_of_gene)
  
  print(paste("new nrow in df_tx are", nrow(df_tx)))

  saveRDS(df_tx, paste0("df_tx_", var_names[s], ".RDS"))
  print(paste0("df_tx_", var_names[s], ".RDS saved"))

  x<-x[(x[,col_txid] %in% df_tx$TXNAME),]
  saveRDS(x, paste0(var_names[s], "_used.RDS"))
  print(paste0(var_names[s], "_used.RDS saved"))

}

