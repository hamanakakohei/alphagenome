# ism

ismのロゴ絵のY軸がどう計算されているか、わかりやすい説明がない気がするが、
https://www.alphagenomedocs.com/api/generated/alphagenome.interpretation.ism.ism_matrix.html#alphagenome.interpretation.ism.ism_matrix
にある通り、score[position, base] - mean(score[position, :])になる。
つまり、ある塩基について3種類のaltアレルのraw scoreを3で割ってマイナスにした値になる。
raw scoreは自分で選んだaggregation function次第。
また、ismパイプラインで計算したスコアは、予測に使うインターバルの中心の変異のみbatch_scoringパイプラインで計算したスコアと一致するが、
周辺の変異はスコアが微妙に異なる。
これは、ismパイプラインでは、予測に使うインターバルは固定されており、変異毎にそれを中心としたインターバルを使わないせいだと思う。


In-silico mutagenesis（01.py）の注意点
- variant_scorersがデフォルトではREFとALTの予測結果の比を出すので、差が欲しければ--aggregation_typeを変える
- 参照：https://www.alphagenomedocs.com/variant_scoring.html


図示（02.py）の注意点
- ISM結果を一つの実験に絞るために、track_metadata.txtを見てフィルター条件を考える
- もしくは複数の実験を平均したりしてもよいかも
