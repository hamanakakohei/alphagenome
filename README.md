# alphagenome

track_metadata.txt
- alphagenomeで予測できる実験データの一覧
- これと同じもの：
  https://www.alphagenomedocs.com/colabs/tissue_ontology_mapping.html

To do
- ハプロタイプで解析できるようにする
 （Y軸にRef/Altハプロタイプのスコアを水平ダッシュ線で参照として示しつつ、各バリアントの影響を示す）

注意点
- cage：ontology_curieにによってはLQhCAGEとhCAGEの両方ある
- splice_junctions：なぜかAssay titleが設定されていない、代わりにnameを見ないといけない。

## License and Disclaimer

This project is based on code from 
[google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome), 
licensed under the Apache License, Version 2.0. You may obtain a copy of the license at: https://www.apache.org/licenses/LICENSE-2.0
Some examples and documentation are adapted from the original repository and are licensed 
under the Creative Commons Attribution 4.0 International License (CC-BY 4.0).

### Modifications

This repository includes modifications compared to the original alphagenome code:
- Added support for XYZ data format
- Changed the training loop
- Simplified the configuration system

Unless required by applicable law or agreed to in writing, this project is distributed 
on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.

