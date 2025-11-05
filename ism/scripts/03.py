#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pickle
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description="Plot ISM results")
parser.add_argument("--ism_heatmap_mat", type=Path, help="バリアントごとのスコア表pkl")
parser.add_argument("--out", type=Path)
parser.add_argument("--width", type=float, default=20)
parser.add_argument("--height", type=float, default=10)
parser.add_argument("--start", type=int, help="描画する範囲の開始座標", required=False)
parser.add_argument("--end", type=int, help="描画する範囲の終了座標", required=False)
args = parser.parse_args()


# === データをDataFrame化 ===
with open(args.ism_detail, "rb") as f:
    data = pickle.load(f)

variants = data["variants"]
scores = data["ism_detail"]

data = []
for v, s in zip(variants, scores):
    data.append({
        "chrom": v.chromosome,
        "pos": v.position,
        "ref": v.reference_bases,
        "alt": v.alternate_bases,
        "score": s
    })
df = pd.DataFrame(data)


# === 指定範囲でフィルタリング ===
if args.start is not None:
  df = df[df["pos"] >= args.start]
if args.end is not None:
  df = df[df["pos"] <= args.end]


# === 可視化のための整形 ===
bases = ["A", "C", "G", "T"]
positions = sorted(df["pos"].unique())

# ポジション×塩基の行列を作る
heatmap_data = pd.DataFrame(index=bases, columns=positions, dtype=float)

for _, row in df.iterrows():
    heatmap_data.loc[row["alt"], row["pos"]] = row["score"]


# === カラーマップと描画 ===
plt.figure(figsize=(args.width, args.height))

vmin = -max(abs(df["score"].min()), df["score"].max())  # 中心を0、赤-青で対称に
vmax =  max(abs(df["score"].min()), df["score"].max())

sns.heatmap(
    heatmap_data,
    cmap="RdBu_r",  # 赤＝プラス、青＝マイナス
    center=0,
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Score"},
    square=False,
    vmin=vmin,
    vmax=vmax,
)

# 各ポジションごとのリファレンス塩基を取得し、灰色で上書き
for pos in positions:
    ref_base = df.loc[df["pos"] == pos, "ref"].iloc[0]
    plt.gca().add_patch(plt.Rectangle(
        (positions.index(pos), bases.index(ref_base)), 1, 1,
        color="black", ec="white", lw=0.5
    ))

#plt.xlabel("Position")
#plt.ylabel("Base")

plt.tight_layout()
plt.savefig(args.out, dpi=300)
plt.close()
