# batch_scoring

事前にVCFをフィルタリングするコマンド例：
```bash
bcftools view \
  -i 'TYPE="snp" && QUAL>5 && FILTER="PASS" && INFO/DP>9 && FORMAT/GQ>5 && GT="1|0"' \
  -o $FILTERED_VCF \
  -Ov \
  $VCF \
  > logs/0.txt 2>&1
```
