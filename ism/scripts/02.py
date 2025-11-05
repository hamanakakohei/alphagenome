#!/usr/bin/env python3
"""
Plot ISM results for a specific ontology CURIE.
参考: https://colab.research.google.com/github/google-deepmind/alphagenome/blob/main/colabs/quick_start.ipynb
"""

# To do:
# - Y軸のタイトルを改行する、visualize_profiles/scriptを参考にできる


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
    #print( adata.var[mask_var].iloc[0] )
    #print( adata.obs[mask_obs].iloc[0] )
    assert values.size == 1, f"Unexpected match count for filters: var={filters_var}, obs={filters_obs}"
    return values.flatten()[0]


def main():
    parser = argparse.ArgumentParser(description="01のISM結果を指定した条件でフィルターしつつ、base毎プロットを描きつつ、置換毎heatmap用の表を保存する")
    parser.add_argument("--in_ism_res", type=Path, default=Path("results/01/out.pkl"))
    parser.add_argument("--out_plot", type=Path)
    parser.add_argument("--out_ism_heatmap_mat", type=Path, help="バリアントごとのスコア表の出力先pkl")
    parser.add_argument("--title", type=str, default="")
    parser.add_argument("--width", type=int, default=20)
    parser.add_argument("--height_scale", type=float, default=1.0)

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


    # --- ISM結果の読み込み ---
    with open(args.in_ism_res, "rb") as f:
        d = pickle.load(f)
        variant_scores = d["variant_scores"]
        ism_interval = d["ism_interval"]
        modality = d["modality"]

    
    # --- ISM結果を一つのTrack結果に絞りプロット用にbase毎スコア表にする ---
    filters_var = {col: getattr(args, col) for col in filter_cols}
    filters_obs = {col: getattr(args, col) for col in filter_rows}

    ism_detail = [ extract_one_experiment(x[0], filters_var, filters_obs) for x in variant_scores ]
    variants = [ v[0].uns['variant'] for v in variant_scores ]

    ism_result = ism.ism_matrix(
        ism_detail,
        variants = variants,
    )
    
    # ついでにISM結果をヒートマップ用に置換毎スコア表にして保存する
    with open(args.out_ism_heatmap_mat, "wb") as f:
        pickle.dump({"ism_detail":ism_detail, "variants":variants}, f)    


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
        title = args.title,
        fig_width = args.width,
        fig_height_scale = args.height_scale,
    )

    args.out_plot.parent.mkdir(parents=True, exist_ok=True)
    plot.figure.savefig(args.out_plot, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()

