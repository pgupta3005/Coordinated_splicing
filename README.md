# Coordinated_splicing

## Pre-requisites (should be present in the working directory in addition to the scripts)

### Files to be made

#### 1. **Tab-separated count matrices** for every sample 
Two columns are mandatory with *'isoform'* and *'counts'* column names. Transcripts under the 'isoform' column should be present in the GTF file.
#### 2. **Transcriptome GTF file** 
Such a file is typically generated as a part of trnascriptome reconstruction from long-read (using Stringtie or FLAIR). It should have entries for *'gene'*, *'exon'*, and *'transcript'* as the third column; the meta-data column (9th) must have gene_id, tx_id and exon_id wherever relevant. Additional tools/steps such as gffread, agat and sed may be required to convert a transcriptome GTF file into script-compatible format (transform_gtf.sh can be used as a reference). An example GTF file is uploaded. 
#### 3. **One tab-separated (with no header) index file** 
It is named *'sample_names'* with the first column as a short-hand notation of the samples (to be used as a prefix for outputs of all steps), the second column as the filename of the count matrix (including the extension, if any), and the third column as the filename of the transcriptome GTF file (including the extension). An example index file is provided.

### R packages to be installed
The following packages can be installed before use or via step0_install.R or while running wrapper_script.R directly
Installed either in the default R lib or in a new conda package [the setup must remain when running the wrapper_script.sh]

```{r}
install.packages("dplyr")
install.packages("tibble")
install.packages("parallel")
install.packages("BiocManager")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")
```

### Running the scripts

```{r}
bash wrapper_script.sh arg1 arg2 arg3
nohup bash wrapper_script.sh arg1 arg2 arg3 > log_out &
```

arg1 is the number of threads while genenrating 2 index files for each genome GTF used (step1) [default takes one-fourth of the available threads]

agr2 is the number of threads for other subsequent steps (steps 2/5/6) [default takes one-eighth of the available threads]

arg3 can be empty or *annot* if the GTF file contains gene_name can that info should be added to the final output .csv file

As the intermediate files can be big in size, it is recommended to run the pipeline on a server/HPC.

## Scripts and their significance

### Step 0: Installing some basic CRAN and BiocManager packages

### Step 1: Indexing the genome GTFs
Outputs --> txdb_\<sample> and tx_\<sample>.RDS and exon_\<sample>.RDS

### Step 2: Filtering
Outputs --> df_exon_\<sample>.RDS and df_tx_\<sample>.RDS 

### Step 3: Demarcating true full-length vs truncated trasncripts
Outputs --> df_exon_\<sample>_final.RDS and df_tx_\<sample>_final.RDS 

### Step 4: Identifying alternative and constitutive exons
Outputs --> constitutive_\<sample>.RDS and non_constitutive_\<sample>.RDS

### Step 5: Making exon pairs and assigning positions
Outputs --> table_\<sample>.RDS

### Step 6: Final calculations
Outputs --> filled_\<sample>.tsv, filled_\<sample>.RDS, significant_\<sample>.RDS 
