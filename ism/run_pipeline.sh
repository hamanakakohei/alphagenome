#!/usr/bin/env bash
set -euo pipefail

source ~/.bash_profile
source ~/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome


# 1用
VARIANT="chr1:1000:G>T"
VARIANT_4_FILENAME=$(tr ':>' '_' <<< "$VARIANT")

PRED_CONTEXT=1MB
ISM_INTERVAL=512

#MODALITY=ATAC
#MODALITY=CONTACT_MAPS
#MODALITY=DNASE
#MODALITY=CHIP_TF
#MODALITY=CHIP_HISTONE
#MODALITY=CAGE
#MODALITY=PROCAP
MODALITY=RNA_SEQ
#MODALITY=SPLICE_SITES
#MODALITY=SPLICE_SITE_USAGE
#MODALITY=SPLICE_JUNCTIONS

OUT_ISM=results/01/out.${VARIANT_4_FILENAME}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.pkl


# 2用
# ISM結果をプロットする対象の実験を選ぶためのフィルター条件
# GENEはRNA_SEQの時のみ、他の時は外す
CURIE=CL:0000312
STRAND=+
GENE=aaa

OUT_PLOT=results/02/out.${VARIANT_4_FILENAME}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.${CURIE}.png
ISM_HEATMAP_MAT=results/02/ism_detail.${VARIANT_4_FILENAME}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.${CURIE}.pkl

LOG1=logs/01/${VARIANT_4_FILENAME}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.log
LOG2=logs/02/${VARIANT_4_FILENAME}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.${CURIE}.log


# 1
mkdir -p logs/01
scripts/01.py \
  --variant $VARIANT \
  --context_width $PRED_CONTEXT \
  --modality $MODALITY \
  --out $OUT_ISM \
  --ism_width_around_variant $ISM_INTERVAL \
  > $LOG1 2>&1


# 2
mkdir -p logs/02
scripts/02.py \
  --in_ism_res $OUT_ISM \
  --ontology_curie $CURIE \
  --strand $STRAND \
  --gene_name $GENE \
  --out_plot $OUT_PLOT \
  --out_ism_heatmap_mat $ISM_HEATMAP_MAT \
  > $LOG2 2>&1
  #--width $W \
  #--height_scale $H \
  #--name $NAME \
  #--data_source encode \
  #--endedness \
  #--genetically_modified \


# 3
OUT=results/03/ism_heatmap.png
scripts/03.py \
  --ism_heatmap_mat $ISM_HEATMAP_MAT \
  --out $OUT \
  #--start 1000 \
  #--end 1100 \
  #--width $W \
  #--height $H
