#!/usr/bin/env python3
# 参考：https://www.alphagenomedocs.com/colabs/visualization_modality_tour.html

# To do：
# - 引き算した絵を出せるようにする（https://www.alphagenomedocs.com/colabs/example_analysis_workflow.html）
# - 未対応のモダリティの図示
# -- ChIP-TFの図示
# -- contact mapの図示
# - プロットするデータをフィルターする引数を増やす、今はontology_curieのみ

# --- 描画パラメーター定義 ---
TSS_ALPHA = 0.5
TSS_COLOR = 'blue'
TSS_MARGIN = 1000

BED_ALPHA = 0.2
BED_LABEL_ANGLE = 90

REF_ALT_COLORS = {'REF': 'black', 'ALT': 'red'}
VARIANT_ALPHA = 0.8
DESPINE = True
DESPINE_KEEP_BOTTOM = False
TITLE = ''

TX_HEIGHT = 0.6 # ?
FIG_WIDTH = 7
FIG_HEIGHT_SCALE = 3.0
HSPACE = 0.3

import sys
import argparse
from pathlib import Path
import logging

# ロギング設定
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# 自作モジュールを追加
github_path = Path.home() / "github"
sys.path.append(str(github_path))
from utils.others import get_api_key

# alphagenome 関連
from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome, track_data, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# --- 引数パーサー設定 ---
parser = argparse.ArgumentParser(description="一つの変異の効果をプロットする")

parser.add_argument("--api_key", type=str, help="ALPHAGENOME_API_KEY 環境変数でも可")
parser.add_argument("--variant", type=str, required=True, help="効果を予測したい変異、例: chr22:36201698:A>C")
parser.add_argument("--curie_list", type=Path, required=True, help="効果を予測するオミクスプロファイルのフィルター条件（細胞種）")
parser.add_argument("--out", type=Path, default=Path("results/01/out.png"))

# Transcript filter
parser.add_argument("--tx_longest", action='store_true', help="Transcriptモデルのうち、最も長いものに絞る")
parser.add_argument("--tx_protein_coding", action='store_true', help="Transcriptモデルのうち、protein-codingのものに絞る")
parser.add_argument("--tx_high_support", action='store_true', help="Transcriptモデルのうち、GENCODE基準でhigh supportなものに絞る")

# Annotations
parser.add_argument("--anno_variant", action='store_true', help="Variantの位置をハイライトする")
parser.add_argument("--anno_bed", type=Path, help="指定したBED領域をハイライトする")
parser.add_argument("--anno_tss", action='store_true', help="TSSの位置をハイライトする")

# プロット位置グループ
center_group = parser.add_mutually_exclusive_group(required=True)
center_group.add_argument("--gene_centered_on", type=str, help="指定した遺伝子を中心にプロットする")
center_group.add_argument("--centered_on_variant", action='store_true', help="バリアントを中心にプロットする")
center_group.add_argument("--plt_interval", type=str, help="指定した領域をプロットする、例: chr1:100-1000")

# 幅の指定グループ
width_group = parser.add_mutually_exclusive_group()
width_group.add_argument("--plt_bp_width", type=int, help="プロット幅を指定する（中心はそのまま）")
width_group.add_argument("--plt_bp_margin", type=int, help="プロット幅を広げる（中心はそのまま）")

args = parser.parse_args()


# --- モデルとデータの読み込み ---
api_key = get_api_key(args.api_key)
dna_model = dna_client.create(api_key)
output_metadata = dna_model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

# GTF読み込み
gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
gtf = pd.read_feather(gtf_url)

# 変異情報の取得
variant = genome.Variant.from_str(args.variant)

# CURIEリスト読み込み
with args.curie_list.open() as f:
    ontology_terms = [line.strip() for line in f if line.strip()]


# --- プロットするインターバルを準備する（色んな指定の仕方ある）---
if args.gene_centered_on:
  plt_interval = gene_annotation.get_gene_interval( gtf, gene_symbol=args.gene_centered_on )
elif args.centered_on_variant:
  plt_interval = variant.reference_interval
elif args.plt_interval:
  chr, sta_end = args.plt_interval.split(":")
  sta, end = map(int, sta_end.split("-"))
  plt_interval = genome.Interval( chr, sta, end )
else:
  print( "エラー：プロットする領域を指定する" )

if args.plt_bp_width:
  plt_interval.resize_inplace( args.plt_bp_width )
elif args.plt_bp_margin:
  plt_interval.resize_inplace( plt_interval.width + args.plt_bp_margin )


# --- プロットするtxを選ぶ ---
gtf_tx_filtered = gtf.copy()
if args.tx_protein_coding:
    gtf_tx_filtered = gene_annotation.filter_protein_coding(gtf_tx_filtered)
if args.tx_high_support:
    gtf_tx_filtered = gene_annotation.filter_transcript_support_level(gtf_tx_filtered, ['1'])
if args.tx_longest:
    gtf_tx_filtered = gene_annotation.filter_to_longest_transcript(gtf_tx_filtered)

tx_extractor = transcript.TranscriptExtractor(gtf_tx_filtered)
txs_4_plot = tx_extractor.extract(plt_interval)


# --- アノテーションのリストを作る ---
annos = []

if args.anno_variant:
    annos.append(plot_components.VariantAnnotation( [variant], alpha = VARIANT_ALPHA ))

if args.anno_bed:
    import pybedtools
    bed = pybedtools.BedTool(str(args.anno_bed))
    bed_intervals = [genome.Interval(l.chrom, l.start, l.end, l.strand) for l in bed]
    bed_labels = [l.name for l in bed]
    annos.append(plot_components.IntervalAnnotation(
        bed_intervals, alpha=BED_ALPHA, labels=bed_labels, label_angle=BED_LABEL_ANGLE))

if args.anno_tss:
    gtf_tss = gene_annotation.extract_tss(gtf_tx_filtered)
    tss_intervals = [
        genome.Interval(
            chromosome=row.Chromosome,
            start=row.Start - TSS_MARGIN,
            end=row.End + TSS_MARGIN,
            name=row.gene_name,
        ) for _, row in gtf_tss.iterrows()
    ]
    annos.append(plot_components.IntervalAnnotation(
        tss_intervals, alpha=TSS_ALPHA, colors=TSS_COLOR))


# --- バリアントの影響を予測する ---
output = dna_model.predict_variant(
    interval = plt_interval.resize( dna_client.SEQUENCE_LENGTH_1MB ),
    variant = variant, 
    requested_outputs = {
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.CAGE,
        dna_client.OutputType.DNASE,
        dna_client.OutputType.ATAC,
        dna_client.OutputType.PROCAP,
        dna_client.OutputType.SPLICE_SITES,
        dna_client.OutputType.SPLICE_SITE_USAGE,
        dna_client.OutputType.SPLICE_JUNCTIONS,
        dna_client.OutputType.CHIP_HISTONE,
       # dna_client.OutputType.CHIP_TF,
       # dna_client.OutputType.CONTACT_MAPS,
    },
    ontology_terms = ontology_terms,
)


# --- データ（実験）を絞る---
# output.reference.atac.metadata を色んな条件でフィルターしてindexを得て、
# それを使って.select_tracks_by_indexすればよいが、
# name列がモダリティのよってちがうしAssay title列は空白で使うのが難しいし、
# ref/alt x 全モダリティにしないといけなくて面倒
# 各データに一意な名前が付けられるとか改善を期待して今は何もしない


# --- プロットする ---
plot = plot_components.plot(
    [
        plot_components.TranscriptAnnotation( txs_4_plot, fig_height = TX_HEIGHT ),
	## junction tracks.
        #plot_components.Sashimi(
        #    output.reference.splice_junctions,
        #    ylabel_template='Ref {biosample_name} ({strand})\n{name}',
        #),
        #plot_components.Sashimi(
        #    output.alternate.splice_junctions,
        #    ylabel_template='Alt {biosample_name} ({strand})\n{name}',
        #),
	# RNA-seq tracks.
        plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.rna_seq, 'ALT': output.alternate.rna_seq },
            colors=REF_ALT_COLORS,
            ylabel_template='RNA-seq:{biosample_name} ({strand})\n{name}',
        ),
	## splice site usage tracks.
        #plot_components.OverlaidTracks(
        #    tdata={ 'REF': output.reference.splice_site_usage, 'ALT': output.alternate.splice_site_usage },
        #    colors=REF_ALT_COLORS,
        #    ylabel_template='SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}',
        #),
        ## CAGE track.
        #plot_components.OverlaidTracks(
        #    tdata={ 'REF': output.reference.cage, 'ALT': output.alternate.cage },
        #    colors=REF_ALT_COLORS,
        #    ylabel_template='CAGE: {biosample_name} ({strand})\n{name}',
        #),
        ## DNase
	#plot_components.OverlaidTracks(
        #    tdata={ 'REF': output.reference.dnase, 'ALT': output.alternate.dnase },
        #    colors=REF_ALT_COLORS,
        #    ylabel_template='DNASE: {biosample_name} ({strand})\n{name}',
        #),
        ## ATAC
	#plot_components.OverlaidTracks(
        #    tdata={ 'REF': output.reference.atac, 'ALT': output.alternate.atac },
        #    colors=REF_ALT_COLORS,
        #    ylabel_template='ATAC: {biosample_name} ({strand})\n{name}',
        #),
        ## splice sites
        #plot_components.OverlaidTracks(
        #    tdata={ 'REF': output.reference.splice_sites, 'ALT': output.alternate.splice_sites },
        #    colors=REF_ALT_COLORS,
        #    ylabel_template='SPLICE SITES: {name} ({strand})',
        #),
        ## ChIP-histone
        #plot_components.Tracks(
        #    tdata=reordered_chip_histone,
        #    ylabel_template=(
        #        'CHIP HISTONE: {biosample_name} ({strand})\n{histone_mark}'
        #    ),
        #    filled=True,
        #    track_colors=track_colors,
        #),
        ## contact maps
	#plot_components.ContactMaps(
        #    tdata=output.contact_maps,
        #    ylabel_template='{biosample_name}\n{name}',
        #    cmap='autumn_r',
        #    vmax=1.0,
        #),
    ],
    annotations = annos,
    interval = plt_interval,
    despine = DESPINE,
    despine_keep_bottom = DESPINE_KEEP_BOTTOM,
    title = TITLE,
    fig_width = FIG_WIDTH,
    fig_height_scale = FIG_HEIGHT_SCALE,
    hspace = HSPACE,
)

plot.figure.savefig(args.out, dpi=300, bbox_inches='tight')
