#!/usr/bin/env Rscript

# To do:
# - 元になった実験データをダウンロードして、並べて可視化する
# - （可能なら）ドットの色を、値の正負で分ける

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(Gviz)
  library(GenomicFeatures)
  library(AnnotationDbi)
})

option_list <- list(
  make_option(c("--input"), type="character"),
  make_option(c("--gtf"), type="character", help="--gtfか--txdb_dbを指定する"),
  make_option(c("--output"), type="character"),
  make_option(c("--score_type"), type="character", default="raw_score", help="raw_score or quantile_score"),
  make_option(c("--txdb_db"), type="character", default=NULL, help="--gtfか--txdb_dbを指定する"),
  make_option(c("--chr"), type="character", default=NULL, help="プロット領域"),
  make_option(c("--start"), type="integer", default=NULL, help="プロット領域"),
  make_option(c("--end"), type="integer", default=NULL, help="プロット領域")
)

opt <- parse_args(OptionParser(option_list=option_list))

# AGの予測結果を読み込む
df <- read_tsv(opt$input)

df <- df %>%
  mutate(
    chrom = str_extract(variant_id, "^chr[0-9XYM]+"),
    pos = as.numeric(str_extract(variant_id, "(?<=:)[0-9]+"))
  )

# quantile_scoreを対数にしてわかりやすくする
# （quantile_scoreが1や-1になることはないと思うが念のため保険をかけつつ）
score_column <- opt$score_type
eps <- 1e-10

df <- df %>%
  mutate(
    quantile_score = ifelse(
      quantile_score >= 0,
      -log10(pmax(1 - quantile_score, eps)),
      log10(pmax(1 + quantile_score, eps))
    )
  )

# プロット領域を決める、引数で指定するか、scoresファイルから自動で抜き出すか
if (!is.null(opt$chr) && !is.null(opt$start) && !is.null(opt$end)) {
  chrom <- opt$chr
  start <- opt$start
  end <- opt$end
} else {
  chrom <- unique(df$chrom)
  start <- min(df$pos) - 1000
  end <- max(df$pos) + 1000
}

# GTFを引数で指定してTxDbにするのだが、TxDbを引数で指定すれば省略できる
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
tracks <- list(axis_track, gene_track)

# トラックごとに背景色を指定
bg_colors <- rep(c("#f0f0f0", "#ffffff"), length.out = length(tracks))

# 各モダリティのTrackを作る
for (ot in unique(df$output_type)) {
  sub <- df %>% filter(output_type == ot)
  
  # Y軸の範囲を対称に設定（最小でも-0.1 ~ 0.1）
  max_abs <- max(abs(sub[[score_column]]), na.rm = TRUE)
  max_abs <- max(max_abs, 0.1)  
  y_range <- c(-1.05*max_abs, 1.05*max_abs)

  track <- DataTrack(
    start = sub$pos,
    end = sub$pos,
    chromosome = chrom,
    genome = "hg38",  
    name = ot,
    type = c("p", "g"),
    data = sub[[score_column]],
    #col = rgb(0, 0, 0.55, alpha = 0.3),
    col = "darkblue",
    cex = 0.8,
    baseline = 0,
    lty.grid = 2,              # 破線のグリッド
    col.grid = "red",       # グリッドの色
    col.axis = "black",        # 軸の色
    col.border.title = "gray80", # タイトル枠
    ylim = y_range             # Y軸を0を中心に対称に
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
  #background.panel = "#FFFEDB",  # パネルごとの背景
  #background.panel =  bg_colors,  # パネルごとの背景
  background.title = "lightgray",  # タイトルの背景
  col.axis = "black",
  col.grid = "gray80",
  col.border.title = "gray80",
  lty.grid = 2,
  lwd.grid = 0.5,
  showGrid = TRUE,
  frame = TRUE,
  col.frame = "black"
)
dev.off()
