#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(Gviz)
  library(GenomicFeatures)
  library(AnnotationDbi)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="TSV input file"),
  make_option(c("-g", "--gtf"), type="character", help="GTF annotation file"),
  make_option(c("-o", "--output"), type="character", help="Output image file (PDF)"),
  make_option(c("-s", "--score_type"), type="character", default="raw_score",
              help="Which score to plot: raw_score or quantile_score [default: %default]"),
  make_option(c("-t", "--txdb_db"), type="character", default=NULL,
              help="Path to saved TxDb SQLite DB file (optional). If not provided, create TxDb from GTF."),
  make_option(c("--chr"), type="character", default=NULL, help="Chromosome to plot"),
  make_option(c("--start"), type="integer", default=NULL, help="Start position"),
  make_option(c("--end"), type="integer", default=NULL, help="End position")
)

opt <- parse_args(OptionParser(option_list=option_list))

df <- read_tsv(opt$input)

df <- df %>%
  mutate(
    chrom = str_extract(variant_id, "^chr[0-9XYM]+"),
    pos = as.numeric(str_extract(variant_id, "(?<=:)[0-9]+"))
  )

# パターン１
chrom <- unique(df$chrom)
start <- min(df$pos) - 1000
end <- max(df$pos) + 1000

# パターン２
if (!is.null(opt$chr) && !is.null(opt$start) && !is.null(opt$end)) {
  chrom <- opt$chr
  start <- opt$start
  end <- opt$end
} else {
  chrom <- unique(df$chrom)
  start <- min(df$pos) - 1000
  end <- max(df$pos) + 1000
}


if (!is.null(opt$txdb_db) && file.exists(opt$txdb_db)) {
  message("Loading TxDb object from ", opt$txdb_db)
  txdb <- loadDb(opt$txdb_db)
} else {
  message("Creating TxDb object from GTF...")
  txdb <- txdbmaker::makeTxDbFromGFF(opt$gtf, format="gtf")
  if (!is.null(opt$txdb_db)) {
    message("Saving TxDb object to ", opt$txdb_db)
    saveDb(txdb, opt$txdb_db)
  }
}

gene_track <- GeneRegionTrack(txdb, chromosome=chrom, start=start, end=end,
                              name="Gene Model", transcriptAnnotation="symbol")

# 軸トラック
axis_track <- GenomeAxisTrack()

# 各output_typeごとのtrackを作成
score_column <- opt$score_type
tracks <- list(axis_track, gene_track)

# 例：トラックごとに交互の背景色を指定
bg_colors <- rep(c("#f0f0f0", "#ffffff"), length.out = length(tracks))

for (ot in unique(df$output_type)) {
  sub <- df %>% filter(output_type == ot)
  track <- DataTrack(
    start = sub$pos,
    end = sub$pos,
    chromosome = chrom,
    genome = "hg38",  # 必要に応じて変更
    name = ot,
    type = c("p", "g"),
    data = sub[[score_column]],
    col = rgb(0, 0, 0.55, alpha = 0.1),
    #col = "darkblue",
    cex = 0.8,
    baseline = 0,
    lty.grid = 2,             # 破線のグリッド
    col.grid = "gray80",      # グリッドの色
    col.axis = "black",       # 軸の色
    col.border.title = "gray" # タイトル枠
  )
  tracks <- append(tracks, list(track))
}

# プロット
#pdf(opt$output, width=12, height=5 + length(tracks)*0.3)
png(opt$output, width=1200, height=500 + length(tracks)*100, res=150)
plotTracks(
  tracks, 
  from=start, 
  to=end, 
  chromosome=chrom,
  background.panel = "#FFFEDB", # bg_colors,  # パネルごとの背景
  background.title = "lightgray",  # タイトルの背景
  col.axis = "black",
  col.grid = "gray80",
  lty.grid = 2,
  lwd.grid = 0.5,
  showGrid = TRUE,
  frame = TRUE,
  col.frame = "gray80"
)
dev.off()
