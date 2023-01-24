#!shebang

#loading libraries
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tibble)
library(parallel)

###########################################
### read files and initialise variables ###
###########################################

a<-read.table("ref_names", sep="\t", header=F)
species_names<-a[,1]     ##define species name
reference_names<-a[,2] 

annot<-""  #default value
num_threads<-as.integer(detectCores()/8) ##default

tmp<-commandArgs(trailingOnly = TRUE)

if(length(tmp)>=1)
  num_threads<-as.numeric(tmp[1])

if(length(tmp)>1)
  annot<-tmp[2]  #re-assigned only if user has provided this

#use gtf files to get gene names for gene ids
if(annot=="annot")
{
  for(i in 1:length(species_names))    #length same as that of reference name
  {
    x<-rtracklayer::import(reference_names[i])
    x<-as.data.frame(x)
    x_sub<-dplyr::select(x, "gene_id", "gene_name")
    x_sub<-unique.data.frame(x_sub)
    assign(paste0("gtf_", species_names[i]), x_sub)
  }
}

b<-read.table("sample_names", sep="\t", header=F)
var_names<-b[,1]
#file_names<-b[,2]
species_names<-b[,3] ###reassigning species_name variable imp

for(s in 1:length(var_names))
{
  assign(paste(var_names[s], "_used", sep=""), readRDS(paste(var_names[s], "_used.RDS", sep="")))
  assign(paste("table1_", var_names[s], sep=""), readRDS(paste("table_", var_names[s], ".RDS", sep="")))
  assign(paste("df_exon_", var_names[s], sep=""), readRDS(paste("df_exon_", var_names[s], "_final.RDS", sep="")))
  assign(paste("df_tx_", var_names[s], sep=""), readRDS(paste("df_tx_", var_names[s], "_final.RDS", sep="")))
}

col_counts <- "counts"
col_txid <- "associated_transcript"

########################
### Define functions ###
########################

final_func=function(x){
  temp<-unlist(strsplit(x, "x"))
  gene<-temp[[1]]
  exon1<-temp[[2]]
  exon2<-temp[[3]]
  
  exonchunk_all<-unique(unlist(df_exon[df_exon$GENEID == gene, "EXONCHUNK"]))  
  e1_n<-which(exonchunk_all==exon1)
  e2_n<-which(exonchunk_all==exon2)
      
  tx_all<-unique(unlist(df_tx[df_tx$GENEID == gene, "TXNAME"]))  
  df_tx_sub<-subset(df_tx, TXNAME %in% tx_all)
  #head(df_tx_sub)
  
  #from all_tx, those tx containing the positions for both exons must only be used ahead
  
  #change the total_set, here all_tx
  #create a whitelist for exon1
  #create a whitelist for exon2
  #only txs common to these whitelists must be retained ultimately in all_tx
  
  inc_tx<-list()
     
  for(k in 1:length(tx_all))
  {
    if(df_tx_sub[df_tx_sub$TXNAME==tx_all[k],"TSS"]=="full"){
      inc_tx<-append(inc_tx, tx_all[k])
    } else {
      n_of_start<-which(exonchunk_all==df_tx_sub[df_tx_sub$TXNAME==tx_all[k],"TSS"])
      n_of_end<-which(exonchunk_all==df_tx_sub[df_tx_sub$TXNAME==tx_all[k],"TTS"])
    
      if(((e1_n>=n_of_start & e1_n<=n_of_end) | (e1_n<=n_of_start & e1_n>=n_of_end)) & ((e2_n>=n_of_start & e2_n<=n_of_end) | (e2_n<=n_of_start & e2_n>=n_of_end))){
        inc_tx<-append(inc_tx, tx_all[k])
      }
    }
  }
      
  inc_tx<-unlist(inc_tx)
  
  #then proceed normally
  
  tx_of_exon1<-intersect(unique(unlist(df_exon[df_exon$EXONCHUNK == exon1, "TXNAME"])), inc_tx)  
  non_tx_of_exon1<-inc_tx[!(inc_tx %in% tx_of_exon1)]
  
  tx_of_exon2<-intersect(unique(unlist(df_exon[df_exon$EXONCHUNK == exon2, "TXNAME"])), inc_tx)  
  non_tx_of_exon2<-inc_tx[!(inc_tx %in% tx_of_exon2)]
  
  a<-intersect(tx_of_exon1, tx_of_exon2)
  b<-intersect(tx_of_exon1, non_tx_of_exon2)
  c<-intersect(non_tx_of_exon1, tx_of_exon2)
  d<-intersect(non_tx_of_exon1, non_tx_of_exon2)
  
  empty_string<-""
  
  val_a<-sum(y[y[,col_txid] %in% a, col_counts])  
  val_b<-sum(y[y[,col_txid] %in% b, col_counts])
  val_c<-sum(y[y[,col_txid] %in% c, col_counts])
  val_d<-sum(y[y[,col_txid] %in% d, col_counts])
  
  total<-val_a + val_b + val_c + val_d
  
  co_inc<-val_a/(val_a + val_b + val_c + val_d)
  
  sig<-(val_a + val_d)/(val_a + val_b + val_c + val_d)
  
  mat<-matrix(c(val_a, val_b, val_c, val_d), nrow=2, ncol=2)
  f<-fisher.test(mat)
  p<-f$p.value
  
  if(val_a==0)
    val_a<-0.5
  if(val_b==0)
    val_b<-0.5
  if(val_c==0)
    val_c<-0.5
  if(val_d==0)
    val_d<-0.5
  
  lor<-log2((val_a*val_d)/(val_b*val_c))

  exonchunk_new<-unlist(df_tx_sub$EXONCHUNK)
  if(length(exonchunk_all)==length(exonchunk_new)) {
    take_ahead  <- "no"
  } else {
    take_ahead <- "yes"
  }

  empty_string<-paste(empty_string, total, co_inc, sig, p, lor, take_ahead, sep="_")
  
  return(empty_string)
}

########################################
### Main function: filling the table ###
########################################

for(s in 1:length(var_names))
{
  print(paste("loading", var_names[s]))
  
  y<-get(paste0(var_names[s], "_used"))
  table1<-get(paste0("table1_", var_names[s]))
  df_exon<-get(paste0("df_exon_", var_names[s]))
  df_tx<-get(paste0("df_tx_", var_names[s]))
  
  #############################################
  print("computing coordinated splicing stats")
  
  output<-mclapply(table1[,"combine"], final_func, mc.cores=num_threads)
  sp<-strsplit(unlist(output), "_")
  head(sp)
  
  new_mat<-matrix(nrow=nrow(table1), ncol=7)
  colnames(new_mat)<-c("total","co_inc","sigma","pval","LOR","cor_pval","take_ahead")
  table1<-cbind(table1,new_mat)
  
  table1$total<-sapply(sp,"[[",2)
  table1$co_inc<-sapply(sp,"[[",3)
  table1$sigma<-sapply(sp,"[[",4)
  table1$pval<-sapply(sp,"[[",5)
  table1$LOR<-sapply(sp,"[[",6)
  table1$take_ahead<-sapply(sp,"[[",7)

  #print(as.data.frame(table1))
  
  table1<-subset(table1, take_ahead=="yes")

  #remove the take_ahead col
  table1<-table1[,-14] 

  #changing the columns from character to numeric datatype (check!!)
  for(i in 8:12)
    table1[,i]<-as.numeric(table1[,i])

  ###################################################
  print("performing false discovery rate correction")

  table1_a<-subset(table1, relpos=="adjacent")
  table1_d1<-subset(table1, relpos=="distant1")
  table1_d2<-subset(table1, relpos=="distant2")

  #calculating corrected p-values
  table1_a[,"cor_pval"]<-p.adjust(table1_a[,"pval"], method = "BY")
  table1_d1[,"cor_pval"]<-p.adjust(table1_d1[,"pval"], method = "BY")
  table1_d2[,"cor_pval"]<-p.adjust(table1_d2[,"pval"], method = "BY")

  table2<-rbind(table1_a, table1_d1, table1_d2)  
  
  ###################################################
  
  if(annot=="annot")
  {
    print("adding meta info and saving two tables per sample")
    tmp<-get(paste0("gtf_", species_names[s]))
    table2<-left_join(table2, tmp, by=c("GENEID"="gene_id"))
  }
  
  write.csv(table2, paste0("filled_", var_names[s], ".csv"), quote=F)

  table3<-subset(table2, cor_pval<0.05)
  
  print(paste(nrow(table3), "exonchunk pairs (from", length(unique(table3$GENEID)), "genes) out of", nrow(table2), "exonchunk pairs (from", length(unique(table2$GENEID)), "genes) show coordinated splicing"))
  
  #saveRDS(table3, paste("significant_", var_names[s], ".RDS", sep=""))
  
}
