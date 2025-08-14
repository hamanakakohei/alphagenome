# batch_scoring

事前にVCFをフィルターするコマンド例：
```bash
bcftools view \
  -i 'TYPE="snp" && QUAL>5 && FILTER="PASS" && INFO/DP>9 && FORMAT/GQ>5 && GT="1|0"' \
  -o $FILTERED_VCF \
  -Ov \
  $VCF \
  > logs/0.txt 2>&1
```


スコアリング（01.py）の注意点
- variant_scorersがREFとALTの予測結果の比を出す
- 参照：https://www.alphagenomedocs.com/variant_scoring.html


図示（03.R）の注意点
- quantile_scoreはlog10(1 - quantile_score)したうえで、元の正負に合わせている
- 小さなスコアを過大に見せないため、最小でも-0.1~0.1の範囲をY軸として表示している
