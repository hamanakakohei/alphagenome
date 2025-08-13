#!/usr/bin/env python3
# 以下をコピペしただけ：
# https://www.alphagenomedocs.com/colabs/batch_variant_scoring.html
# 
# 元はgoogle colab用のスクリプトなので、ローカル用に適宜コメントアウト、追加している

import sys
from pathlib import Path
import argparse

# 自作モジュールを追加
github_path = Path.home() / "github"
sys.path.append(str(github_path))
from utils.vcf import read_vcf_as_df_4_alphagenome
from utils.others import get_api_key

def parse_args():
    parser = argparse.ArgumentParser(description="変異をVCFで与えて、AlphaGenomeの予測を表で得る")
    parser.add_argument("--api_key", type=str, help="ALPHAGENOME_API_KEY環境変数でもok")
    parser.add_argument("--vcf_file", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/scores.txt"))
    return parser.parse_args()

args = parse_args()
api_key = get_api_key(args.api_key)

# ----------------------------------------------------------------------
# @title Install AlphaGenome

# @markdown Run this cell to install AlphaGenome.
#from IPython.display import clear_output
#! pip install alphagenome
#clear_output()

# ----------------------------------------------------------------------
# @title Setup and imports.

from io import StringIO
from alphagenome import colab_utils
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
#from google.colab import data_table, files
import pandas as pd
from tqdm import tqdm

#data_table.enable_dataframe_formatter()

# Load the model.
#dna_model = dna_client.create(colab_utils.get_api_key())
dna_model = dna_client.create( api_key )

# ----------------------------------------------------------------------
# @title Score batch of variants.

# Load VCF file containing variants.
#vcf_file = 'placeholder'  # @param
#
# We provide an example list of variants to illustrate.
#vcf_file = """variant_id\tCHROM\tPOS\tREF\tALT
#chr3_58394738_A_T_b38\tchr3\t58394738\tA\tT
#chr8_28520_G_C_b38\tchr8\t28520\tG\tC
#chr16_636337_G_A_b38\tchr16\t636337\tG\tA
#chr16_1135446_G_T_b38\tchr16\t1135446\tG\tT
#"""
#
#vcf = pd.read_csv(StringIO(vcf_file), sep='\t')
vcf = read_vcf_as_df_4_alphagenome( args.vcf_file )

required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
for column in required_columns:
  if column not in vcf.columns:
    raise ValueError(f'VCF file is missing required column: {column}.')

organism = 'human'  # @param ["human", "mouse"] {type:"string"}

# @markdown Specify length of sequence around variants to predict:
sequence_length = '1MB'  # @param ["2KB", "16KB", "100KB", "500KB", "1MB"] { type:"string" }
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
    f'SEQUENCE_LENGTH_{sequence_length}'
]

# @markdown Specify which scorers to use to score your variants:
score_rna_seq = True  # @param { type: "boolean"}
score_cage = True  # @param { type: "boolean" }
score_procap = True  # @param { type: "boolean" }
score_atac = True  # @param { type: "boolean" }
score_dnase = True  # @param { type: "boolean" }
score_chip_histone = True  # @param { type: "boolean" }
score_chip_tf = True  # @param { type: "boolean" }
score_polyadenylation = True  # @param { type: "boolean" }
score_splice_sites = True  # @param { type: "boolean" }
score_splice_site_usage = True  # @param { type: "boolean" }
score_splice_junctions = True  # @param { type: "boolean" }

# @markdown Other settings:
#download_predictions = False  # @param { type: "boolean" }

# Parse organism specification.
organism_map = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]

# Parse scorer specification.
scorer_selections = {
    'rna_seq': score_rna_seq,
    'cage': score_cage,
    'procap': score_procap,
    'atac': score_atac,
    'dnase': score_dnase,
    'chip_histone': score_chip_histone,
    'chip_tf': score_chip_tf,
    'polyadenylation': score_polyadenylation,
    'splice_sites': score_splice_sites,
    'splice_site_usage': score_splice_site_usage,
    'splice_junctions': score_splice_junctions,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]

# Remove any scorers or output types that are not supported for the chosen organism.
unsupported_scorers = [
    scorer
    for scorer in selected_scorers
    if (
        organism.value
        not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    )
    | (
        (scorer.requested_output == dna_client.OutputType.PROCAP)
        & (organism == dna_client.Organism.MUS_MUSCULUS)
    )
]
if len(unsupported_scorers) > 0:
  print(
      f'Excluding {unsupported_scorers} scorers as they are not supported for'
      f' {organism}.'
  )
  for unsupported_scorer in unsupported_scorers:
    selected_scorers.remove(unsupported_scorer)


# Score variants in the VCF file.
results = []

for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf)):
  variant = genome.Variant(
      chromosome=str(vcf_row.CHROM),
      position=int(vcf_row.POS),
      reference_bases=vcf_row.REF,
      alternate_bases=vcf_row.ALT,
      name=vcf_row.variant_id,
  )
  interval = variant.reference_interval.resize(sequence_length)

  variant_scores = dna_model.score_variant(
      interval=interval,
      variant=variant,
      variant_scorers=selected_scorers,
      organism=organism,
  )
  results.append(variant_scores)

df_scores = variant_scorers.tidy_scores(results)

#if download_predictions:
#  df_scores.to_csv('variant_scores.csv', index=False)
#  files.download('variant_scores.csv')

df_scores.to_csv(args.out, sep='\t')

## ----------------------------------------------------------------------
## Examine just the effects of the variants on T-cells.
#columns = [c for c in df_scores.columns if c != 'ontology_curie']
#df_scores[(df_scores['ontology_curie'] == 'CL:0000084')][columns]
