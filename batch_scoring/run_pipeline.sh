#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
export ALPHAGENOME_API_KEY="your_api_key"


VCF=data/variants.vcf

AG_SCORE=results/01/scores.txt

GENCODE_GTF=data/gencode.v47.primary_assembly.annotation.gtf
GENCODE_RDS=data/gencode.v47.primary_assembly.annotation.rds
SCORE_TYPE=raw_score # quantile_score
PLOT_CHR=chr1
PLOT_STA=1000
PLOT_END=100000
PLOT_OUT=results/02/out.${PLOT_CHR}.${PLOT_STA}.${PLOT_END}.${SCORE_TYPE}.png


# 1
conda activate alphagenome
mkdir -p results/01

scripts/01.py \
  --vcf_file $VCF \
  --out $AG_SCORE \
  > logs/1.txt 2>&1

awk -F"\t" 'NR==1 || $15=="CL:0000312" || $15=="CL:1001606" || $15=="CL:2000092"' results/01/scores.txt > results/01/scores.some_cells.txt


# 2
conda activate gviz
mkdir -p results/02

scripts/02.R \
  --input results/01/scores.some_cells.txt \
  --gtf $GENCODE_GTF \
  --txdb_db $GENCODE_RDS \
  --chr $PLOT_CHR \
  --start $PLOT_STA \
  --end $PLOT_END \
  --score_type $SCORE_TYPE \
  --output $PLOT_OUT \
  > logs/2.txt 2>&1
  #--input $AG_SCORE \
