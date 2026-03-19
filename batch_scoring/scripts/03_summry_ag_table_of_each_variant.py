#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np


def compute_stats(series: pd.Series) -> pd.Series:
    s = series.dropna()

    def nth_largest(s, n):
        if len(s) < n:
            return np.nan
        return s.nlargest(n).iloc[-1]

    def nth_smallest(s, n):
        if len(s) < n:
            return np.nan
        return s.nsmallest(n).iloc[-1]

    return pd.Series({
        "median": s.median(),
        "max": s.max(),
        "min": s.min(),
        "top5": nth_largest(s, 5),
        "bottom5": nth_smallest(s, 5),
        "top10": nth_largest(s, 10),
        "bottom10": nth_smallest(s, 10),
    })


def process(df, query, group_cols, label):
    sub = df.query(query).copy()
    if sub.empty:
        return pd.DataFrame()

    res = (
        sub.groupby(group_cols)["quantile_score"]
        .apply(compute_stats)
        .reset_index()
    )

    res["category"] = label
    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", required=True, help="AGアノテーション結果表")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')

    configs = [
        # simple
        ("output_type == 'ATAC'",       ["variant_id"], "ATAC"),
        ("output_type == 'CAGE'",       ["variant_id"], "CAGE"),
        ("output_type == 'CHIP_HISTON'", ["variant_id"], "CHIP_HISTON"),
        ("output_type == 'CHIP_TF'",    ["variant_id"], "CHIP_TF"),
        ("output_type == 'DNASE'",      ["variant_id"], "DNASE"),
        ("output_type == 'PROCAP'",     ["variant_id"], "PROCAP"),

        # histone細分化
        ("output_type == 'CHIP_HISTON' and histone_mark == 'H3K27ac'", ["variant_id"], "H3K27ac"),
        ("output_type == 'CHIP_HISTON' and histone_mark == 'H3K9ac'",  ["variant_id"], "H3K9ac"),
        ("output_type == 'CHIP_HISTON' and histone_mark == 'H3K4me3'", ["variant_id"], "H3K4me3"),

        # geneあり
        ("output_type == 'RNA_SEQ'",           ["variant_id", "gene_id", "gene_name"], "RNA_SEQ"),
        ("output_type == 'SPLICE_JUNCTIONS'",  ["variant_id", "gene_id", "gene_name"], "SPLICE_JUNCTIONS"),
        ("output_type == 'SPLICE_SITES'",      ["variant_id", "gene_id", "gene_name"], "SPLICE_SITES"),
        ("output_type == 'SPLICE_SITE_USAGE'", ["variant_id", "gene_id", "gene_name"], "SPLICE_SITE_USAGE"),
    ]

    results = []
    for query, group_cols, label in configs:
        res = process(df, query, group_cols, label)
        if not res.empty:
            results.append(res)

    final_df = pd.concat(results, ignore_index=True)
    final_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
