#!/usr/bin/env python3
"""
Run ISM scoring for a given variant using alphagenome.
参考: https://colab.research.google.com/github/google-deepmind/alphagenome/blob/main/colabs/quick_start.ipynb

-  CenterMaskScorer(        requested_output=ATAC,               width=501,   aggregation_type=DIFF_LOG2_SUM)
-  CenterMaskScorer(        requested_output=DNASE,              width=501,   aggregation_type=DIFF_LOG2_SUM)
-  CenterMaskScorer(        requested_output=CHIP_TF,            width=501,   aggregation_type=DIFF_LOG2_SUM)
-  CenterMaskScorer(        requested_output=CHIP_HISTONE,       width=2001,  aggregation_type=DIFF_LOG2_SUM)
-  CenterMaskScorer(        requested_output=CAGE,               width=501,   aggregation_type=DIFF_LOG2_SUM)
-  CenterMaskScorer(        requested_output=PROCAP,             width=501,   aggregation_type=DIFF_LOG2_SUM)
-  GeneMaskLFCScorer(       requested_output=RNA_SEQ)
-  PolyadenylationScorer()
-  GeneMaskSplicingScorer(  requested_output=SPLICE_SITES,       width=None)
-  GeneMaskSplicingScorer(  requested_output=SPLICE_SITE_USAGE,  width=None)
-  SpliceJunctionScorer()
-  ContactMapScorer()
-  GeneMaskActiveScorer(    requested_output=RNA_SEQ)
-  CenterMaskScorer(        requested_output=ATAC,               width=501,   aggregation_type=ACTIVE_SUM)
-  CenterMaskScorer(        requested_output=DNASE,              width=501,   aggregation_type=ACTIVE_SUM)
-  CenterMaskScorer(        requested_output=CHIP_TF,            width=501,   aggregation_type=ACTIVE_SUM)
-  CenterMaskScorer(        requested_output=CHIP_HISTONE,       width=2001,  aggregation_type=ACTIVE_SUM)
-  CenterMaskScorer(        requested_output=CAGE,               width=501,   aggregation_type=ACTIVE_SUM)
-  CenterMaskScorer(        requested_output=PROCAP,             width=501,   aggregation_type=ACTIVE_SUM)]
"""

import sys
import argparse
from pathlib import Path
import logging
import pickle

# ロギング設定
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# 自作モジュール読み込み
github_path = Path.home() / "github"
sys.path.append(str(github_path))
from utils.others import get_api_key

# alphagenome 関連
from alphagenome.data import genome
from alphagenome.interpretation import ism
from alphagenome.models import dna_client, variant_scorers


def get_context_width(key: str):
    return getattr(dna_client, f"SEQUENCE_LENGTH_{key}")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--api_key", type=str, help="ALPHAGENOME_API_KEY 環境変数でも可")
    parser.add_argument("--variant", type=str, required=True, help="例: chr22:36201698:A>C")
    parser.add_argument("--context_width", type=str, choices=["2KB", "16KB", "100KB", "500KB", "1MB"], required=True)
    parser.add_argument("--modality", type=str, choices=variant_scorers.RECOMMENDED_VARIANT_SCORERS.keys(), required=True)
    parser.add_argument("--scorer_class", type=str, help="例: CenterMaskScorer")
    parser.add_argument("--aggregation_type", type=str, help="例: DIFF_LOG2_SUM, DIFF_MEAN")
    parser.add_argument("--width", type=int, help="variant_scorerのスコアリング幅")
    parser.add_argument("--out", type=Path, default=Path("results/01/out.pkl"))

    # ISM領域指定
    width_group = parser.add_mutually_exclusive_group(required=True)
    width_group.add_argument("--ism_width_around_variant", type=int)
    width_group.add_argument("--ism_region", type=str)

    args = parser.parse_args()


    # モデル読み込み
    api_key = get_api_key(args.api_key)
    dna_model = dna_client.create(api_key)

    # 変異オブジェクト生成
    variant = genome.Variant.from_str(args.variant)

    # プレディクション時のコンテキスト幅設定
    seq_length = getattr(dna_client, f"SEQUENCE_LENGTH_{args.context_width}")
    sequence_interval = variant.reference_interval.resize(seq_length)

    # ISM領域設定
    if args.ism_width_around_variant:
        ism_interval = variant.reference_interval.resize(args.ism_width_around_variant)
    else:
        chrom, sta_end = args.ism_region.split(":")
        sta, end = map(int, sta_end.split("-"))
        ism_interval = genome.Interval(chrom, sta, end)

    # variant_scorer 設定
    default_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[args.modality]
    
    if hasattr(default_scorer, "aggregation_type") and hasattr(default_scorer, "width"):
         ScorerClass = getattr(
             variant_scorers,
             args.scorer_class or default_scorer.__class__.__name__
         )
         
         agg_type = getattr(
             variant_scorers.AggregationType,
             args.aggregation_type or default_scorer.aggregation_type.name
         )
         
         width = args.width or default_scorer.width
   
         variant_scorer = ScorerClass(
             requested_output=dna_client.OutputType[args.modality],
             width=width,
             aggregation_type=agg_type,
         )
    else:
        logging.warning(
            f"{default_scorer.__class__.__name__} は aggregation_type や width を持たないため、"
            f"デフォルトのインスタンスを使用します"
        )
        variant_scorer = default_scorer


    # ISMスコア計算
    variant_scores = dna_model.score_ism_variants(
        interval = sequence_interval,
        ism_interval = ism_interval,
        variant_scorers = [variant_scorer],
    )

    logging.info(f"variant_scores length (ism_interval x 3): {len(variant_scores)}")
    logging.info(f"one variant ISM data shape: {variant_scores[0][0].X.shape}")


    # 保存
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wb") as f:
        pickle.dump({
	    "variant_scores": variant_scores, 
	    "ism_interval": ism_interval,
	    "modality": args.modality
	}, f)


if __name__ == "__main__":
    main()
