#!/usr/bin/env bash
#
# 1. batch/chunk ag annotation
# 2. APIなので毎回表のパターンが変わるかもなので、想定通りかチェック
# 3. 各バリアントに付く大量のアノテを適当なルールでサマライズする
set -euo pipefail


# 1
CHUNK_LIST=data/chunks.list

while read CHUNK; do
    bash scripts/01.chunk_ver.sh $CHUNK
done < $CHUNK_LIST


# 2
eval "$(conda shell.bash hook)"
conda activate misc

python3 scripts/02_validate_col_assumptions.py \
    --in random_ag_res.list


# 3
AG_RES_LIST=data/ag_res.list

while read AG_RES; do
    AG_SUMSTA=
    python3 scripts/03_summry_ag_table_of_each_variant.py \
        --in $AG_RES \
        --out $AG_SUMSTA
done < $AG_RES_LIST
