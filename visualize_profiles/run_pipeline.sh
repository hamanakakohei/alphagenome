#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

export ALPHAGENOME_API_KEY="your_api_key"


VARIANT="chr1:1000:G>T"
CURIE_LIST=data/curies.txt
OUT=out.png


scripts/01.py \
  --variant $VARIANT \
  --curie_list $CURIE_LIST \
  --out $OUT \
  --tx_longest \
  --tx_protein_coding \
  --tx_high_support \
  --anno_variant \
  --anno_tss \
  --centered_on_variant \
  --plt_bp_margin 10000 
  #--gene_centered_on SERPINB7 \
  
