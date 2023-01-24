# Coordinated_splicing

## Pre-requisites (should be present in the working directory in addition to the six scripts 1-6)

### Files to be made

#### 1. **Tab-separated count matrices** for every sample 
Two columns are mandatory with *'associated_transcript'* and *'counts'* column names. Transcripts under the 'associated_transcript' column should be present in the GTF file.
#### 2. **Genome annotation GTF file** 
Should have entries for *'gene'*, *'exon'*, and *'transcript'* as the third column; the meta-data column (9th) must have gene_id, tx_id and exon_id wherever relevant. An example GTF file is uploaded.
#### 3. **Two tab-separated (with no header) index files** 
(A) One called *'ref_names'* with the first column as a short-hand notation of the reference(s) used and the second column as the name of the GTF file (in the working directory, with the extension), and (B) the other called *'sample_names'* with the first column as a short-hand notation for the samples used, the second column as the name of the count matrix (with entension, if any), and the third column as short-hand notation of their respective genomes. An example for each index file is provided.

### R pakcages to be installed

```{r}
install.packages("dplyr")
install.packages("tibble")
install.packages("parallel")  #if not installed by default
install.packages("BiocManager")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")
```

### Running the scripts

## shebang line ? - change to where the packages are installed

#### Run it as - possible options

```{r}
Rscript stepn_xyz.R
Rscript stepn_xyz.R arg
nohup Rscript stepn_xyz.R > log.out &
nohup Rscript stepn_xyz.R arg > log.out &
```
#### Steps 1, 2, 5, and 6 have components that can be parallelized. As defaults, these scripts take one-fourth, one-eighth, one-tenth and one-tenth of the available cores, repectively. If the user wishes to use lesser or more threads, they can specify the number of threads as the arg for these scripts.

#### Step 6 also allows adding **gene name** annotation from the gtf file (if present) in the final output. Default: no information is added. To add this info, add *'annot'* as the second argument after specifiying the number of threads for script 6.

## Scripts and their significance

### Step 1

### Step 2

### Step 3

### Step 4

### Step 5

### Step 6




