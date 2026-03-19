#!/usr/bin/env python3
#
# ATAC: cells
# CAGE: assays x strands x cells
# CHIP_HISTONE: mods x cells
# CHIP_TF: tfs x cells
# DNASE: cells
# PROCAP: strands x cells
# RNA_SEQ: genes x strands x assays x cells (x gtex)
# SPLICE_JUNCTIONS: genes x track_name (x gtex)
# SPLICE_SITES: genes (x track_name x strands)
# SPLICE_SITE_USAGE: genes x cells x assay (x strands x gtex)
#
#
# cells: ontology_curie
# strands: track_strand
# assay: Assay title
# mods: histone_mark
# tfs: transcription_factor
# genes: gene_id or gene_name
# gtex: gtex_tissue
#
# 1. ファイルパスが書かれたリストを読み込む
# 2. 各ファイルを読み込んでリストに格納
# 3. 縦に結合
import pandas as pd
import gzip
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--in", required=True, help="適当な数のAGアノテーション結果表のリスト、全てを調べると大きすぎるので")
args = parser.parse_args()


# 1
with open(args.in, 'r') as f:
    filepaths = [line.strip() for line in f if line.strip()]


# 2
dfs = []
for path in filepaths:
    print(f"Reading: {path}")
    temp_df = pd.read_csv(path, sep='\t', low_memory=False)
    dfs.append(temp_df)


# 3
df = pd.concat(dfs, axis=0, ignore_index=True, sort=False)


# 4
len( combined_df.query('output_type == "ATAC"')[['variant_id', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "ATAC"')[['variant_id', 'ontology_curie']] )

len( combined_df.query('output_type == "CAGE"')[['variant_id', 'Assay title', 'track_strand', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "CAGE"')[['variant_id', 'Assay title', 'track_strand', 'ontology_curie']] )

len( combined_df.query('output_type == "CHIP_HISTONE"')[['variant_id', 'histone_mark', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "CHIP_HISTONE"')[['variant_id', 'histone_mark', 'ontology_curie']] )

len( combined_df.query('output_type == "CHIP_TF"')[['variant_id', 'transcription_factor', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "CHIP_TF"')[['variant_id', 'transcription_factor', 'ontology_curie']] )

len( combined_df.query('output_type == "DNASE"')[['variant_id', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "DNASE"')[['variant_id', 'ontology_curie']] )

len( combined_df.query('output_type == "PROCAP"')[['variant_id', 'track_strand', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "PROCAP"')[['variant_id', 'track_strand', 'ontology_curie']] )

len( combined_df.query('output_type == "RNA_SEQ"')[['variant_id', 'gene_id', 'gene_name', 'Assay title', 'track_strand', 'ontology_curie']].drop_duplicates() )
len( combined_df.query('output_type == "RNA_SEQ"')[['variant_id', 'gene_id', 'gene_name', 'Assay title', 'track_strand', 'ontology_curie']] )

len( combined_df.query('output_type == "RNA_SEQ"')[['variant_id', 'gene_id', 'gene_name', 'Assay title', 'track_strand', 'ontology_curie', 'gtex_tissue']].drop_duplicates() )
len( combined_df.query('output_type == "RNA_SEQ"')[['variant_id', 'gene_id', 'gene_name', 'Assay title', 'track_strand', 'ontology_curie', 'gtex_tissue']] )

len( combined_df.query('output_type == "SPLICE_JUNCTIONS"')[['variant_id', 'gene_id', 'gene_name', 'track_name']].drop_duplicates() )
len( combined_df.query('output_type == "SPLICE_JUNCTIONS"')[['variant_id', 'gene_id', 'gene_name', 'track_name']] )

len( combined_df.query('output_type == "SPLICE_JUNCTIONS"')[['variant_id', 'gene_id', 'gene_name', 'track_name', 'gtex_tissue']].drop_duplicates() )
len( combined_df.query('output_type == "SPLICE_JUNCTIONS"')[['variant_id', 'gene_id', 'gene_name', 'track_name', 'gtex_tissue']] )

len( combined_df.query('output_type == "SPLICE_SITES"')[['variant_id', 'gene_id', 'gene_name']].drop_duplicates() )
len( combined_df.query('output_type == "SPLICE_SITES"')[['variant_id', 'gene_id', 'gene_name']] )

len( combined_df.query('output_type == "SPLICE_SITES"')[['variant_id', 'gene_id', 'gene_name', 'track_name', 'track_strand']].drop_duplicates() )
len( combined_df.query('output_type == "SPLICE_SITES"')[['variant_id', 'gene_id', 'gene_name', 'track_name', 'track_strand']] )

len( combined_df.query('output_type == "SPLICE_SITE_USAGE"')[['variant_id', 'gene_id', 'gene_name', 'ontology_curie', 'Assay title']].drop_duplicates() )
len( combined_df.query('output_type == "SPLICE_SITE_USAGE"')[['variant_id', 'gene_id', 'gene_name', 'ontology_curie', 'Assay title']] )

len( combined_df.query('output_type == "SPLICE_SITE_USAGE"')[['variant_id', 'gene_id', 'gene_name', 'ontology_curie', 'Assay title', 'gtex_tissue', 'track_strand']].drop_duplicates() )
len( combined_df.query('output_type == "SPLICE_SITE_USAGE"')[['variant_id', 'gene_id', 'gene_name', 'ontology_curie', 'Assay title', 'gtex_tissue', 'track_strand']] )
