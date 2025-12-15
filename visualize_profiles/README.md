## 使い方

run_pipeline.sh
'''
VARIANT="chr1:1000:G>T"
CURIE_LIST=data/curies.txt
REGION="chr1:1000-2000"
OUT=results/01/out.${REGION}.png
'''
を編集して実行する。
curies.txtはtrack_metadata.txtからてきとうに興味のあるIDを使う。

## To do
- 引数で指定するようにする？
