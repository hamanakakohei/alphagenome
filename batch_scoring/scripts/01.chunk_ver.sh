#!/usr/bin/env bash
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate alphagenome
source ~/.bash_profile


ANALYSIS=$1   #chunk_000
VCF=data/vcf_chunks_1000/${ANALYSIS}.vcf
SCORE_DIR=results/01/${ANALYSIS}


mkdir -p logs
mkdir -p results/01
mkdir -p results/01/${ANALYSIS}


# 1：VCFを与えてスコアを計算する
scripts/01.chunk_ver.py \
  --vcf_file $VCF \
  --out_dir $SCORE_DIR \
  > logs/1.${ANALYSIS}.txt 2>&1
