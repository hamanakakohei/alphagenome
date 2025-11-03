#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
source ~/.bash_profile


ANALYSIS=aaa

VCF=data/aaa.vcf

SCORES=results/01/scores.${ANALYSIS}.txt
SCORES_FILT=results/02/scores.${ANALYSIS}.filtered.txt

GENCODE_GTF=data/gencode.v47.primary_assembly.annotation.gtf
GENCODE_RDS=data/gencode.v47.primary_assembly.annotation.rds
#SCORE_TYPE=raw_score
SCORE_TYPE=quantile_score
PLOT_CHR=chr1
PLOT_STA=1000 
PLOT_END=2000
PLOT_OUT=results/03/out.${ANALYSIS}.${PLOT_CHR}.${PLOT_STA}.${PLOT_END}.${SCORE_TYPE}.png

mkdir -p logs


# 1：VCFを与えてスコアを計算する
conda activate alphagenome
mkdir -p results/01

scripts/01.py \
  --vcf_file $VCF \
  --out $SCORES \
  > logs/1.${ANALYSIS}.txt 2>&1


# 2：スコアを興味のある細胞種に絞る
conda activate misc
scripts/02.py \
  -i $SCORES \
  -o $SCORES_FILT \
  --filter gene_name SERPINB7 \
  --filter gene_strand + \
  --filter track_strand + \
  --filter ontology_curie CL:0000312 NaN \
  --filter output_type RNA_SEQ SPLICE_SITES \
  > logs/2.${ANALYSIS}.txt 2>&1
  #--include-nan


# 3：図示する
conda activate gviz
mkdir -p results/03

scripts/03.R \
  --input $SCORES_FILT \
  --gtf $GENCODE_GTF \
  --txdb_db $GENCODE_RDS \
  --don_acc_max \
  --chr $PLOT_CHR \
  --start $PLOT_STA \
  --end $PLOT_END \
  --score_type $SCORE_TYPE \
  --output $PLOT_OUT \
  > logs/3.${ANALYSIS}.txt 2>&1
