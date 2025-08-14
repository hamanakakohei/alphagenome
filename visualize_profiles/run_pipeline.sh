#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

export ALPHAGENOME_API_KEY=""


VARIANT="chr1:1000:G>T"
CURIE_LIST=data/curies.txt
REGION="chr1:1000-2000"
OUT=results/01/out.${REGION}.png
#MARGIN=5000
#GENE=aaa


mkdir -p results/01

scripts/01.py \
  --variant $VARIANT \
  --curie_list $CURIE_LIST \
  --out $OUT \
  --anno_variant \
  --plt_interval $REGION
  #--tx_high_support \
  #--tx_longest \
  #--tx_protein_coding \
  #--plt_bp_margin $MARGIN 
  #--gene_centered_on $GENE \
  #--centered_on_variant \
  #--anno_tss \
