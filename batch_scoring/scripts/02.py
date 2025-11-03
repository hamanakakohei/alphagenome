#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter dataframe by column-value pairs (multiple allowed)."
    )
    parser.add_argument("-i", required=True, type=Path)
    parser.add_argument("-o", required=True, type=Path)
    parser.add_argument(
        "--filter",
        nargs="+",
        action="append",
        metavar=("COLUMN", "VALUES..."),
        help="Filtering condition: e.g. --filter colA a b c --filter colB x y",
    )
    parser.add_argument(
        "--include-nan",
        action="store_true",
        help="If set, keep rows where the filtered column is missing (NaN) in addition to matches.",
    )
    return parser.parse_args()

def main():
    args = parse_args()
    df = pd.read_csv(args.i, sep="\t", low_memory=False)

    # 複数列のフィルタ条件を順に適用
    if args.filter:
        for f in args.filter:
            col = f[0]
            values = f[1:]
            want_nan = any(v.lower() in {"nan", "na"} for v in values)

            # 欠損も含めて残すかどうか
            if args.include_nan or want_nan:
                missing_values = {".", ""} # "NA"
                mask = (
                    df[col].isin(values)
                    | df[col].isna()
                    | df[col].isin(missing_values)
                )
            else:
                mask = df[col].isin(values)

            df = df[mask]

    # 保存
    args.o.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.o, sep="\t", index=False)

if __name__ == "__main__":
    main()
