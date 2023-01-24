#!$rscript

list_of_CRAN_packages<-c("dplyr", "tibble", "parallel", "BiocManager")
list_of_bioc_packages<-c("rtracklayer", "GenomicFeatures")

new_CRAN_packages<-list_of_CRAN_packages[!(list_of_CRAN_packages %in% installed.packages()[,"Package"])]
new_bioc_packages<-list_of_bioc_packages[!(list_of_bioc_packages %in% installed.packages()[,"Package"])]

#if(length(new_CRAN_packages)>0)
#  install.packages(new_CRAN_packages, dependencies=T, repos='http://cran.us.r-project.org')

if(length(new_bioc_packages)>0)
  BiocManager::install(new_bioc_packages, dependencies=T)

lapply(c(list_of_CRAN_packages, list_of_bioc_packages), require, character.only = TRUE)
