!#!/bin/bash

source ~/miniconda3/bin/activate env1

INPUT_DIR=$PWD
OUTPUT_DIR=$PWD

# ${i}.collapsed_classification.filtered_lite.gtf is the output from FLAIR
# ${i}_3.gtf can be used for this pipeline

while read i
do
gffread -E --keep-genes ${i}.collapsed_classification.filtered_lite.gtf -o- > ${i}_1.gtf
docker run --ipc=host --gpus all \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0 agat_convert_sp_gff2gtf.pl \
--in "${INPUT_DIR}"/${i}_1.gtf --out "${OUTPUT_DIR}"/${i}_2.gtf
sed 's/ID/exon_id/g' ${i}_2.gtf > ${i}_3.gtf
done < filenames
