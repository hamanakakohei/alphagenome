## 使い方

run_pipeline.sh内の
```bash
VARIANT="chr1:1000:G>T"
CURIE_LIST=data/curies.txt
REGION="chr1:1000-2000"
OUT=results/01/out.${REGION}.png
```
を編集して実行する。
curies.txtはtrack_metadata.txtからてきとうに興味のあるontology_curieのリストを入れる。

## To do
- 引数で指定するようにする？
