#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alphagenome

export ALPHAGENOME_API_KEY="your_api_key"

VARIANT="chr1:1000:G>T"
VARIANT2=$(tr ':>' '_' <<< "$VARIANT")

PRED_CONTEXT=1MB

#MODALITY=ATAC
#MODALITY=CONTACT_MAPS
#MODALITY=DNASE
#MODALITY=CHIP_TF
#MODALITY=CHIP_HISTONE
MODALITY=CAGE
#MODALITY=PROCAP
#MODALITY=RNA_SEQ
#MODALITY=SPLICE_SITES
#MODALITY=SPLICE_SITE_USAGE
#MODALITY=SPLICE_JUNCTIONS

ISM_INTERVAL=256
OUT_ISM=results/01/out.${VARIANT2}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.pkl

CURIE=CL:0000312
#CURIE=CL:1001606
#CURIE=CL:2000092
STRAND=+
#STRAND=plus
#STRAND=minus
GENE=SERPINB7
OUT_PLOT=results/02/out.${VARIANT2}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.${CURIE}.png

LOG1=logs/01.${VARIANT2}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.log
LOG2=logs/02.${VARIANT2}.${MODALITY}.${PRED_CONTEXT}.${ISM_INTERVAL}.${CURIE}.log

## 1
#mkdir -p logs
#
#scripts/01.py \
#  --variant $VARIANT \
#  --context_width $PRED_CONTEXT \
#  --modality $MODALITY \
#  --out $OUT_ISM \
#  --ism_width_around_variant $ISM_INTERVAL \
#  > $LOG1 2>&1
#  #--ism_region xxxx \


# 2
scripts/02.py \
  --ism_result $OUT_ISM \
  --out $OUT_PLOT \
  --ontology_curie $CURIE \
  --strand $STRAND \
  > $LOG2 2>&1
  #--gene_name $GENE \
  #--name $NAME \



#--data_source encode \
#--endedness \
#--genetically_modified \
