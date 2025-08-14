#!/usr/bin/env bash
set -euo pipefail

export ALPHAGENOME_API_KEY=""

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

VCF=data/aaa.vcf.gz

SCORES=results/01/scores.txt
CURIES=data/ontology_curies.txt
SCORES_OF_CURIES=results/02/scores.curies.txt

GENCODE_GTF=data/gencode.v47.primary_assembly.annotation.gtf
GENCODE_RDS=data/gencode.v47.primary_assembly.annotation.rds
SCORE_TYPE=raw_score
#SCORE_TYPE=quantile_score
PLOT_CHR=chr1
PLOT_STA=1000 
PLOT_END=2000
PLOT_OUT=results/03/out.${PLOT_CHR}.${PLOT_STA}.${PLOT_END}.${SCORE_TYPE}.png

mkdir -p logs


# 1：VCFを与えてスコアを計算する
mkdir -p results/01

scripts/01.py \
  --vcf_file VCF \
  --out $SCORES \
  > logs/1.txt 2>&1


# 2：スコアを興味のある細胞種に絞る
awk -F"\t" '
  NR==1 { print; next }       
  FNR==NR { curie[$1]; next } 
  $15 in curie                
  ' $CURIES $SCORES > $SCORES_OF_CURIES


# 3：図示する
conda activate gviz
mkdir -p results/03

scripts/03.R \
  --input $SCORES_OF_CURIES \
  --gtf $GENCODE_GTF \
  --txdb_db $GENCODE_RDS \
  --chr $PLOT_CHR \
  --start $PLOT_STA \
  --end $PLOT_END \
  --score_type $SCORE_TYPE \
  --output $PLOT_OUT \
  > logs/3.txt 2>&1
