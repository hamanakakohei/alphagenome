# alphagenome
To do:
- 1 snp x 全データ（全細胞 x 全モダリティ）を表で保存する
- VCF x 全データ（全細胞 x 全モダリティ）を
-- 表で保存しつつ、
-- X軸はゲノムポジでスコアを各モダリティごとにプロットする
- ハプロタイプは
-- Y軸にRef/Altハプロタイプのスコアを水平ダッシュ線で参照として示しつつ、
-- 上の図と全く同じようにする

各モダリティのデータ（トラック）が、細胞種（ontology_curie列）で一意に決まるか？
- atac：各ontology_curieに一つだけ
| Assay title | data\_source | endedness | genetically\_modified | strand | count |
| ----------- | ------------ | --------- | --------------------- | ------ | ----- |
| ATAC-seq    | encode       | paired    | False                 | .      | 166   |
| ATAC-seq    | encode       | single    | False                 | .      | 1     |

- cage：ontology_curieにによってはLQhCAGEとhCAGEの両方ある
Assay title  data_source  strand
hCAGE        fantom       +       258
hCAGE        fantom       -       258
LQhCAGE      fantom       +        15
LQhCAGE      fantom       -        15
  
- dnase：各ontology_curieに一つだけ
Assay title  data_source  endedness  genetically_modified  strand
DNase-seq    encode       paired     False                 .       197
                          single     False                 l       108
  
- rna_seq：各ontology_curieに多くのデータがある
strand  Assay title         data_source  endedness  genetically_modified
+       total RNA-seq       encode       paired     False                   135
-       total RNA-seq       encode       paired     False                   135
+       polyA plus RNA-seq  encode       paired     False                    75
-       polyA plus RNA-seq  encode       paired     False                    75
+       total RNA-seq       encode       single     False                    61
-       total RNA-seq       encode       single     False                    61
.       polyA plus RNA-seq  gtex         paired     False                    54
                            encode       single     False                    38
                                         paired     False                    31
        total RNA-seq       encode       paired     False                     2
        
- chip_histone：ontology_curie x histone_markペアに1つずつデータがある
strand  Assay title       data_source  endedness  genetically_modified
.       Histone ChIP-seq  encode       single     False                   1046
                                       paired     False                     70

- chip_tf：genetically_modified==Falseならontology_curie x transcription_factorで一つに決まる（SEやPEがあるが、両方はない）。genetically_modifiedは色んな種類があり、name列で見分ける？
Assay title  data_source  endedness  genetically_modified  strand
TF ChIP-seq  encode       single     False                 .       595
                                     True                          559
                          paired     False                         288
                                     True                          175

- splice_sites：
name      strand
acceptor  +         1
          -         1
donor     +         1
          -         1
  
- splice_site_usage：strand x Assay title (polyA+ / total) x ontology_curie x data_source (gtex / encode)なら一つに決まる。encode total: 195; encode polyA+: 118; gtex: 54
strand  Assay title         data_source
-       total RNA-seq       encode         195
+       total RNA-seq       encode         195
-       polyA plus RNA-seq  encode         118
+       polyA plus RNA-seq  encode         118
                            gtex            54
-       polyA plus RNA-seq  gtex            54
  
- splice_junctions：なぜかAssay titleが設定されていない、代わりにnameを見ないといけない。367データあるが、恐らく上のencode total: 195; encode polyA+: 118; gtex: 54だろう。

- procap：6 ontology_curie x 2 strand = 12データのみ
