#!/usr/bin/env python3
# 参考：https://www.alphagenomedocs.com/colabs/visualization_modality_tour.html

# To do：引き算した絵を出せるようにする（https://www.alphagenomedocs.com/colabs/example_analysis_workflow.html）

# --- 描画パラメーター定義 ---
TSS_ALPHA = 0.5
TSS_COLOR = 'blue'
TSS_MARGIN = 1000
BED_ALPHA = 0.2
BED_LABEL_ANGLE = 90
TX_HEIGHT = 0.1
REF_ALT_COLORS = {'REF': 'dimgrey', 'ALT': 'red'}
VARIANT_ALPHA = 0.8
DESPINE = True
DESPINE_KEEP_BOTTOM = False
TITLE = ''
FIG_WIDTH = 35


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
parser.add_argument("--variant", type=str, required=True, help="例: chr22:36201698:A>C")
parser.add_argument("--curie_list", type=Path, required=True)
parser.add_argument("--out", type=Path, default=Path("results/plot.png"))

# Transcript filter
parser.add_argument("--tx_longest", action='store_true')
parser.add_argument("--tx_protein_coding", action='store_true')
parser.add_argument("--tx_high_support", action='store_true')

# Annotations
parser.add_argument("--anno_variant", action='store_true')
parser.add_argument("--anno_bed", type=Path)
parser.add_argument("--anno_tss", action='store_true')

# プロット位置グループ
center_group = parser.add_mutually_exclusive_group(required=True)
center_group.add_argument("--gene_centered_on", type=str)
center_group.add_argument("--centered_on_variant", action='store_true')
center_group.add_argument("--plt_interval", type=str, help="例: chr1:100-1000")

# 幅の指定グループ
width_group = parser.add_mutually_exclusive_group()
width_group.add_argument("--plt_bp_width", type=int)
width_group.add_argument("--plt_bp_margin", type=int)

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


# --- プロットする ---
plot = plot_components.plot(
    [
        plot_components.TranscriptAnnotation( txs_4_plot, fig_height = TX_HEIGHT ),
	# junction tracks.
        plot_components.Sashimi(
            output.reference.splice_junctions,
            ylabel_template='Ref {biosample_name} ({strand})\n{name}',
        ),
        plot_components.Sashimi(
            output.alternate.splice_junctions,
            ylabel_template='Alt {biosample_name} ({strand})\n{name}',
        ),
	# RNA-seq tracks.
        plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.rna_seq, 'ALT': output.alternate.rna_seq },
            colors=REF_ALT_COLORS,
            ylabel_template='RNA-seq:{biosample_name} ({strand})\n{name}',
        ),
	# splice site usage tracks.
        plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.splice_site_usage, 'ALT': output.alternate.splice_site_usage },
            colors=REF_ALT_COLORS,
            ylabel_template='SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}',
        ),
        # CAGE track.
        plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.cage, 'ALT': output.alternate.cage },
            colors=REF_ALT_COLORS,
            ylabel_template='CAGE: {biosample_name} ({strand})\n{name}',
        ),
        # DNase
	plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.dnase, 'ALT': output.alternate.dnase },
            colors=REF_ALT_COLORS,
            ylabel_template='DNASE: {biosample_name} ({strand})\n{name}',
        ),
        # ATAC
	plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.atac, 'ALT': output.alternate.atac },
            colors=REF_ALT_COLORS,
            ylabel_template='ATAC: {biosample_name} ({strand})\n{name}',
        ),
        # splice sites
        plot_components.OverlaidTracks(
            tdata={ 'REF': output.reference.splice_sites, 'ALT': output.alternate.splice_sites },
            colors=REF_ALT_COLORS,
            ylabel_template='SPLICE SITES: {name} ({strand})',
        ),
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
)

plot.figure.savefig(args.out, dpi=300, bbox_inches='tight')



	## RNA-seq tracks.
        #plot_components.Tracks(
        #    tdata=output.rna_seq,
        #    ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
        #    shared_y_scale = True,
        #),


#output = dna_model.predict_interval(
#    interval=interval,
#    requested_outputs={
#        dna_client.OutputType.RNA_SEQ,
#        dna_client.OutputType.CAGE,
#    },
#    ontology_terms=ontology_terms,
#)


## -------------------------------------------------------
## Splicing
## List of IDs corresponding to various intestinal tissues.
#ontology_terms = [
#    'UBERON:0001157',
#    'UBERON:0001159',
#]
#
## Make predictions for splicing outputs and RNA_SEQ.
#output = dna_model.predict_interval(
#    interval=interval,
#    requested_outputs={
#        dna_client.OutputType.RNA_SEQ,
#        dna_client.OutputType.SPLICE_SITES,
#        dna_client.OutputType.SPLICE_SITE_USAGE,
#        dna_client.OutputType.SPLICE_JUNCTIONS,
#    },
#    ontology_terms=ontology_terms,
#)
#
## Build plot.
## Since APOL4 is on the negative DNA strand, we use `filter_negative_strand` to
## consider only negative stranded splice predictions.
#plot = plot_components.plot(
#    [
#        plot_components.TranscriptAnnotation(longest_transcripts),
#        plot_components.Tracks(
#            tdata=output.splice_sites.filter_to_negative_strand(),
#            ylabel_template='SPLICE SITES: {name} ({strand})',
#        ),
#    ],
#    interval=apol4_interval,
#    title='Predicted splicing effects for colon tissue',
#)
#
## Make predictions for REF and ALT alleles.
#output = dna_model.predict_variant(
#    interval=interval,
#    variant=variant,
#    requested_outputs={
#        dna_client.OutputType.RNA_SEQ,
#        dna_client.OutputType.SPLICE_SITES,
#        dna_client.OutputType.SPLICE_SITE_USAGE,
#        dna_client.OutputType.SPLICE_JUNCTIONS,
#    },
#    ontology_terms=ontology_terms,
#)
#
## Get all transcripts, not just the longest one per gene.
#transcripts = transcript_extractor.extract(interval)
#
#ref_output = output.reference
#alt_output = output.alternate
#
## Build plot.
#plot = plot_components.plot(
#    [
#        plot_components.TranscriptAnnotation(transcripts),
#        plot_components.Sashimi(
#            ref_output.splice_junctions
#            .filter_to_strand('-')
#            .filter_by_tissue('Colon_Transverse'),
#            ylabel_template='Reference {biosample_name} ({strand})\n{name}',
#        ),
#        plot_components.Sashimi(
#            alt_output.splice_junctions
#            .filter_to_strand('-')
#            .filter_by_tissue('Colon_Transverse'),
#            ylabel_template='Alternate {biosample_name} ({strand})\n{name}',
#        ),
#        plot_components.OverlaidTracks(
#            tdata={
#                'REF': ref_output.rna_seq.filter_to_nonpositive_strand(),
#                'ALT': alt_output.rna_seq.filter_to_nonpositive_strand(),
#            },
#            colors=ref_alt_colors,
#            ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
#        ),
#        plot_components.OverlaidTracks(
#            tdata={
#                'REF': ref_output.splice_sites.filter_to_nonpositive_strand(),
#                'ALT': alt_output.splice_sites.filter_to_nonpositive_strand(),
#            },
#            colors=ref_alt_colors,
#            ylabel_template='SPLICE SITES: {name} ({strand})',
#        ),
#        plot_components.OverlaidTracks(
#            tdata={
#                'REF': (
#                    ref_output.splice_site_usage.filter_to_nonpositive_strand()
#                ),
#                'ALT': (
#                    alt_output.splice_site_usage.filter_to_nonpositive_strand()
#                ),
#            },
#            colors=ref_alt_colors,
#            ylabel_template=(
#                'SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}'
#            ),
#        ),
#    ],
#    interval=apol4_interval,
#    annotations=[plot_components.VariantAnnotation([variant])],
#    title='Predicted REF vs. ALT effects of variant in colon tissue',
#)
#
#
#
## -------------------------------------------------------
## ChIP-Histone
#
#output_metadata.chip_histone[
#    output_metadata.chip_histone['biosample_name'].str.contains('colon')
#]
#
## List of IDs corresponding to various colon tissues in `CHIP_HISTONE` output.
#ontology_terms = [
#    'UBERON:0000317',
#    'UBERON:0001155',
#    'UBERON:0001157',
#    'UBERON:0001159',
#]
#
## Make predictions.
#output = dna_model.predict_interval(
#    interval=interval,
#    requested_outputs={dna_client.OutputType.CHIP_HISTONE},
#    ontology_terms=ontology_terms,
#)
#
#gtf_tss = gene_annotation.extract_tss(gtf_longest_transcript)
#
#tss_as_intervals = [
#    genome.Interval(
#        chromosome=row.Chromosome,
#        start=row.Start,
#        end=row.End + 1000,  # Add extra 1Kb so the TSSs are visible.
#        name=row.gene_name,
#    )
#    for _, row in gtf_tss.iterrows()
#]
#
#reordered_chip_histone = output.chip_histone.select_tracks_by_index(
#    output.chip_histone.metadata.sort_values('histone_mark').index
#)
#
#histone_to_color = {
#    'H3K27AC': '#e41a1c',
#    'H3K36ME3': '#ff7f00',
#    'H3K4ME1': '#377eb8',
#    'H3K4ME3': '#984ea3',
#    'H3K9AC': '#4daf4a',
#    'H3K27ME3': '#ffc0cb',
#}
#
#track_colors = (
#    reordered_chip_histone.metadata['histone_mark']
#    .map(lambda x: histone_to_color.get(x.upper(), '#000000'))
#    .values
#)
#
## Build plot.
#plot = plot_components.plot(
#    [
#        plot_components.TranscriptAnnotation(longest_transcripts),
#        plot_components.Tracks(
#            tdata=reordered_chip_histone,
#            ylabel_template=(
#                'CHIP HISTONE: {biosample_name} ({strand})\n{histone_mark}'
#            ),
#            filled=True,
#            track_colors=track_colors,
#        ),
#    ],
#    interval=interval,
#    annotations=[
#        plot_components.IntervalAnnotation(
#            tss_as_intervals, alpha=0.5, colors='blue'
#        )
#    ],
#    despine_keep_bottom=True,
#    title='Predicted histone modification markers in colon tissue',
#)
#
#
## -------------------------------------------------------
## ChIP-TF
#ontology_terms = [
#    'UBERON:0001159',  # Sigmoid colon.
#    'UBERON:0001157',  # Transverse colon.
#    'EFO:0002067',  # K562.
#    'EFO:0001187',  # HepG2.
#]
#
#output = dna_model.predict_interval(
#    interval=interval,
#    requested_outputs={dna_client.OutputType.CHIP_TF},
#    ontology_terms=ontology_terms,
#)
#
#ontology_terms = [
#    'EFO:0002067',  # K562.
#    'EFO:0001187',  # HepG2.
#]
#
#output_chip_tf = output.chip_tf.filter_tracks(
#    (output.chip_tf.metadata['ontology_curie'].isin(ontology_terms)).values
#)
#len(output_chip_tf.metadata)
#
#max_predictions = output_chip_tf.metadata[
#    ['ontology_curie', 'biosample_name', 'transcription_factor']
#].copy()
#
#max_predictions.loc[:, 'max_prediction'] = output_chip_tf.values.max(axis=0)
#max_predictions.sort_values('max_prediction', ascending=False).reset_index(
#    drop=True
#)
#
#print(f'Number of tracks before filtering: {len(output_chip_tf.metadata)}')
#
#output_filtered = output_chip_tf.filter_tracks(
#    output_chip_tf.values.max(axis=0) > 8000
#)
#print(f'Number of tracks after filtering: {len(output_filtered.metadata)}')
#
## Build plot.
#plot_components.plot(
#    components=[
#        plot_components.TranscriptAnnotation(longest_transcripts),
#        plot_components.Tracks(
#            tdata=output_filtered,
#            ylabel_template=(
#                'CHIP TF: {biosample_name} ({strand})\n{transcription_factor}'
#            ),
#            filled=True,
#        ),
#    ],
#    interval=interval,
#    title='Predicted TF-binding in K562 and HepG2 cell-lines.',
#    despine_keep_bottom=True,
#    annotations=[
#        plot_components.IntervalAnnotation(
#            tss_as_intervals, alpha=0.3, colors='blue'
#        )
#    ],
#)

