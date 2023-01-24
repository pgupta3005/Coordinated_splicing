#!/bin/bash

threads1=$1
threads2=$2
annot=$3

if [ "$1" == "--help" ]; then
  echo "Usage: bash wrapper_script.sh 32 10 annot"
  echo "Arguments in order"
  echo "1: Number of threads for step 1 [default: one-fourth of the available threads]"
  echo "2: Number of threads for steps 2, 5 and 6 [default: one-tenth of the available threads]"
  echo "3: Whether to add gene name annot to final output [default: <empty> ] [or write <annot>]"
  exit 0
fi

rscript=$(which Rscript)

sed -i "s/\/usr\/bin\/Rscript/\$rscript/g" step0_install.R
sed -i "s/\/usr\/bin\/Rscript/\$rscript/g" step1_txdb_obj.R
sed -i 's/\/usr\/bin\/Rscript/$rscript/g' step2_exonchunks.R
sed -i 's/\/usr\/bin\/Rscript/$rscript/g' step3_TSS_TTS.R
sed -i 's/\/usr\/bin\/Rscript/$rscript/g' step4_const_alt_exons.R
sed -i 's/\/usr\/bin\/Rscript/$rscript/g' step5_pairs.R
sed -i 's/\/usr\/bin\/Rscript/$rscript/g' step6_calculations.R

Rscript step0_install.R
Rscript step1_txdb_obj.R $threads1
Rscript	step2_exonchunks.R $threads2
Rscript	step3_TSS_TTS.R
Rscript	step4_const_alt_exons.R
Rscript	step5_pairs.R $threads2
Rscript	step6_calculations.R $threads2 $annot
