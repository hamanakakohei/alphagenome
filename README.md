# alphagenome

## 使い方
1. 環境変数を`~/.bash_profile`に加える：
```bash
export ALPHAGENOME_API_KEY="you_api_key"
```

2. alphagenome仮想環境を作る：
```bash
mamba create -n alphagenome python=3.11 matplotlib numpy pandas seaborn pillow
conda activate alphagenome
pip install -U alphagenome
```

3. gviz仮想環境はこう：
- https://github.com/hamanakakohei/gviz

4. utilsリポを`~/github/`以下に置く：
```bash
mkdir -p ~/github
git clone https://github.com/hamanakakohei/utils/ ~/github/utils
```

参考：track_metadata.txt
- alphagenomeで予測できる実験データの一覧
- これと同じもの：
  https://www.alphagenomedocs.com/colabs/tissue_ontology_mapping.html

## To do
- ハプロタイプで解析できるようにする
 （Y軸にRef/Altハプロタイプのスコアを水平ダッシュ線で参照として示しつつ、各バリアントの影響を示す）

## 注意点
- cage：ontology_curieにによってはLQhCAGEとhCAGEの両方ある
- splice_junctions：なぜかAssay titleが設定されていない、代わりにnameを見ないといけない。

## License and Disclaimer

This project uses code from [google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome), 
licensed under the Apache License, Version 2.0. You may obtain a copy of the license at: https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, this project is distributed 
on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
