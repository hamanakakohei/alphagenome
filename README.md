# alphagenome
To do:
- 1 snp x 全データ（全細胞 x 全モダリティ）を表で保存する
- VCF x 全データ（全細胞 x 全モダリティ）を
-- 表で保存しつつ、
-- X軸はゲノムポジでスコアを各モダリティごとにプロットする
- ハプロタイプは
-- Y軸にRef/Altハプロタイプのスコアを水平ダッシュ線で参照として示しつつ、
-- 上の図と全く同じようにする

ontology_curieの値を1つ与えて、以下について全てのトラックを表示する
rna_seq
cage
procap
atac
dnase
chip_histone
splice_sites
splice_site_usage
splice_junctions
出来れば実験データをダウンロードして並べて示す

各モダリティのデータ（トラック）が、細胞種（ontology_curie列）で一意に決まるか？
- atac：各ontology_curieに一つだけ

| Assay title | data\_source | endedness | genetically\_modified | strand | count |
| ----------- | ------------ | --------- | --------------------- | ------ | ----- |
| ATAC-seq    | encode       | paired    | False                 | .      | 166   |
| ATAC-seq    | encode       | single    | False                 | .      | 1     |

- cage：ontology_curieにによってはLQhCAGEとhCAGEの両方ある

| Assay title | data\_source | strand | count |
| ----------- | ------------ | ------ | ----- |
| hCAGE       | fantom       | +      | 258   |
| hCAGE       | fantom       | -      | 258   |
| LQhCAGE     | fantom       | +      | 15    |
| LQhCAGE     | fantom       | -      | 15    |

  
- dnase：各ontology_curieに一つだけ

| Assay title | data\_source | endedness | genetically\_modified | strand | count |
| ----------- | ------------ | --------- | --------------------- | ------ | ----- |
| DNase-seq   | encode       | paired    | False                 | .      | 197   |
| DNase-seq   | encode       | single    | False                 | .      | 108   |


- rna_seq：各ontology_curieに多くのデータがある

| Assay title        | data\_source | endedness | genetically\_modified | strand | count |
| ------------------ | ------------ | --------- | --------------------- | ------ | ----- |
| total RNA-seq      | encode       | paired    | False                 | +      | 135   |
| total RNA-seq      | encode       | paired    | False                 | -      | 135   |
| polyA plus RNA-seq | encode       | paired    | False                 | +      | 75    |
| polyA plus RNA-seq | encode       | paired    | False                 | -      | 75    |
| total RNA-seq      | encode       | single    | False                 | +      | 61    |
| total RNA-seq      | encode       | single    | False                 | -      | 61    |
| polyA plus RNA-seq | gtex         | paired    | False                 | .      | 54    |
| polyA plus RNA-seq | encode       | single    | False                 | .      | 38    |
| polyA plus RNA-seq | encode       | paired    | False                 | .      | 31    |
| total RNA-seq      | encode       | paired    | False                 | .      | 2     |


- chip_histone：ontology_curie x histone_markペアに1つずつデータがある

| Assay title      | data\_source | endedness | genetically\_modified | strand | count |
| ---------------- | ------------ | --------- | --------------------- | ------ | ----- |
| Histone ChIP-seq | encode       | single    | False                 | .      | 1046  |
| Histone ChIP-seq | encode       | paired    | False                 | .      | 70    |

- chip_tf：genetically_modified==Falseならontology_curie x transcription_factorで一つに決まる（SEやPEがあるが、両方はない）。genetically_modifiedは色んな種類があり、name列で見分ける？

| Assay title | data\_source | endedness | genetically\_modified | strand | count |
| ----------- | ------------ | --------- | --------------------- | ------ | ----- |
| TF ChIP-seq | encode       | single    | False                 | .      | 595   |
| TF ChIP-seq | encode       | single    | True                  | .      | 559   |
| TF ChIP-seq | encode       | paired    | False                 | .      | 288   |
| TF ChIP-seq | encode       | paired    | True                  | .      | 175   |

- splice_sites：

| name     | strand | count |
| -------- | ------ | ----- |
| acceptor | +      | 1     |
| acceptor | -      | 1     |
| donor    | +      | 1     |
| donor    | -      | 1     |


- splice_site_usage：strand x Assay title (polyA+ / total) x ontology_curie x data_source (gtex / encode)なら一つに決まる。encode total: 195; encode polyA+: 118; gtex: 54

| Assay title        | data\_source | strand | count |
| ------------------ | ------------ | ------ | ----- |
| total RNA-seq      | encode       | -      | 195   |
| total RNA-seq      | encode       | +      | 195   |
| polyA plus RNA-seq | encode       | -      | 118   |
| polyA plus RNA-seq | encode       | +      | 118   |
| polyA plus RNA-seq | gtex         | .      | 54    |
| polyA plus RNA-seq | gtex         | -      | 54    |

- splice_junctions：なぜかAssay titleが設定されていない、代わりにnameを見ないといけない。367データあるが、恐らく上のencode total: 195; encode polyA+: 118; gtex: 54だろう。

- procap：6 ontology_curie x 2 strand = 12データのみ
