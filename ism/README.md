# ism

In-silico mutagenesis（01.py）の注意点
- variant_scorersがデフォルトではREFとALTの予測結果の比を出すので、差が欲しければ--aggregation_typeを変える
- 参照：https://www.alphagenomedocs.com/variant_scoring.html


図示（02.py）の注意点
- ISM結果を一つの実験に絞るために、track_metadata.txtを見てフィルター条件を考える
- もしくは複数の実験を平均したりしてもよいかも
