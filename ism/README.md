# ism

スコアリング（01.py）の注意点
- variant_scorersがデフォルトではREFとALTの予測結果の比を出すので、差が欲しければ変える
- 参照：https://www.alphagenomedocs.com/variant_scoring.html


図示（03.R）の注意点
- quantile_scoreはlog10(1 - quantile_score)したうえで、元の正負に合わせている
- 小さなスコアを過大に見せないため、最小でも-0.1~0.1の範囲をY軸として表示している
