#!/usr/bin/env python3
"""
Plot ISM results for a specific ontology CURIE.
参考: https://colab.research.google.com/github/google-deepmind/alphagenome/blob/main/colabs/quick_start.ipynb
"""


# --- 描画パラメーター定義 ---
TITLE = ''
FIG_WIDTH = 35


# --- モジュール読み込み ---
import argparse
from pathlib import Path
import logging
import pickle
import numpy as np

# ロギング設定
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# alphagenome 関連
from alphagenome.interpretation import ism
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt


def extract_one_experiment(adata, filters_var=None, filters_obs=None):
    if filters_var is None: # 冗長に見えるが、ミュータブルなデフォルト引数はだめなので
        filters_var = {}
    if filters_obs is None:
        filters_obs = {}

    # 列フィルター
    mask_var = np.ones(len(adata.var), dtype=bool)
    for col, val in filters_var.items():
        if val is not None:
            mask_var &= (adata.var[col] == val)

    # 行フィルター
    mask_obs = np.ones(len(adata.obs), dtype=bool)
    for col, val in filters_obs.items():
        if val is not None:
            mask_obs &= (adata.obs[col] == val)

    # 両方のフィルターを適用
    values = adata.X[np.ix_(mask_obs, mask_var)]
    print( adata.var[mask_var].iloc[0] )
    print( adata.obs[mask_obs].iloc[0] )
    assert values.size == 1, f"Unexpected match count for filters: var={filters_var}, obs={filters_obs}"
    return values.flatten()[0]


def main():
    # --- 引数パーサー設定 ---
    parser = argparse.ArgumentParser(description="Plot ISM results")
    parser.add_argument("--ism_result", type=Path, default=Path("results/01/out.pkl"))
    parser.add_argument("--out", type=Path)

    # ISM結果をフィルタするための引数（一部のみ指定すれば良い）
    filter_cols = [
        "name", # "Assay title"は空白があるので使えないので、代わりに"name"を使う
        "ontology_curie",
        "strand",
        "data_source",
        "endedness",
        "genetically_modified",
        "assay_title",
    ]
    
    for col in filter_cols:
        parser.add_argument(f"--{col}", type=str)

    filter_rows = [
        "gene_id",
        "gene_name",
        "gene_type",
        # "strand", # これはfilter_cols（adata.var）にもあるので、どうすべきか、、
    ]
    
    for col in filter_rows:
        parser.add_argument(f"--{col}", type=str)

    args = parser.parse_args()

    # 出力ファイルパスの設定
    if args.out is None:
        args.out = Path(f"results/02/{args.ism_result.stem}.{args.ontology_curie}.png")

    # --- ISM結果の読み込み ---
    with open(args.ism_result, "rb") as f:
        d = pickle.load(f)
        variant_scores = d["variant_scores"]
        ism_interval = d["ism_interval"]
        modality = d["modality"]

    # --- ISM結果を一つのTrack結果に絞る ---
    filters_var = {col: getattr(args, col) for col in filter_cols}
    filters_obs = {col: getattr(args, col) for col in filter_rows}

    ism_result = ism.ism_matrix(
        [ extract_one_experiment(x[0], filters_var, filters_obs) for x in variant_scores ],
        variants = [ v[0].uns['variant'] for v in variant_scores ],
    )

    logging.info(f"ISM plot data shape: {ism_result.shape}")

    # --- プロット ---
    plot = plot_components.plot(
        [
            plot_components.SeqLogo(
                scores = ism_result,
                scores_interval = ism_interval,
                ylabel = f'ISM {str(ism_interval)} {args.ontology_curie} {modality}',
            )
        ],
        interval = ism_interval,
        title = TITLE,
        fig_width = FIG_WIDTH,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plot.figure.savefig(args.out, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()


#def extract_one_experiment(adata, **filters):
#    mask = np.ones(len(adata.var), dtype=bool)
#    for col, val in filters.items():
#        if val is not None:
#            mask &= (adata.var[col] == val)
#
#    values = adata.X[:, mask]
#    assert values.size == 1, f"Unexpected match count for filters: {filters}"
#    return values.flatten()[0]
