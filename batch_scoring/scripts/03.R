#!/usr/bin/env Rscript

# To do:
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
  make_option(c("--don_acc_max"), action="store_true", default=FALSE, help="SPLICE_SITESのdoc & accスコアの最大値を示す"),
  make_option(c("--chr"), type="character", default=NULL, help="プロット領域"),
  make_option(c("--start"), type="integer", default=NULL, help="プロット領域"),
  make_option(c("--end"), type="integer", default=NULL, help="プロット領域"),
  make_option(c("--width"), type="integer", default=3000),
  make_option(c("--height"), type="integer", default=2000)
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

# SPLICE_SITESをaccとdonで最大値をとる
if( opt$don_acc_max ){
  df <- df %>%
    # SPLICE_SITES とそれ以外に分ける
    group_by(output_type) %>%
    group_split() %>%
    lapply(function(sub_df) {
      if (unique(sub_df$output_type) == "SPLICE_SITES") {
        # グループ化キーを動的に決定
        exclude_cols <- c("Unnamed: 0", "track_name", "raw_score", "quantile_score")
        group_cols <- setdiff(colnames(sub_df), exclude_cols)

        sub_df %>%
          #group_by(variant_id, gene_id, gene_name, gene_strand, track_strand) %>%
          group_by(across(all_of(group_cols))) %>%
          summarise(
            n_rows = n(),
            quantile_score = max(quantile_score),
            raw_score = max(raw_score),
            .groups = "drop"
          ) %>%
          # 行数が2か確認（acc & don）
          { if (any(.$n_rows != 2)) stop("Some groups do not have 2 rows!") else . } %>%
          dplyr::select(-n_rows)
      } else {
        sub_df
      }
    }) %>%
    bind_rows()
}


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
                              name="Gene Model") #, transcriptAnnotation="symbol")

# 軸トラック
axis_track <- GenomeAxisTrack()

# 各output_typeごとのtrackを作成
tracks <- list(axis_track, gene_track)

# トラックごとに背景色を指定
bg_colors <- rep(c("#f0f0f0", "#ffffff"), length.out = length(tracks))

group_cols <- c(
  "gene_name",
  "gene_strand",
  "output_type",
  "track_name",
  "track_strand",
  "Assay title",
  "ontology_curie"
)

# 各モダリティのTrackを作る
for (sub in df %>% group_by(across(all_of(group_cols))) %>% group_split()) {
  # セミコロンで全情報を連結
  meta <- sub %>% distinct(across(all_of(group_cols)))

  track_label <- meta %>%
    unlist() %>%
    paste(collapse = "; ")

  # 値の最大絶対値を取得
  # output_typeに応じてY軸範囲を変更
  max_abs <- max(abs(sub[[score_column]]), na.rm = TRUE)
  if (meta$output_type[[1]] == "SPLICE_SITES") {
    max_abs <- max(max_abs, 0.2)
    y_range <- c(0, 1.05 * max_abs)
  } else {
    max_abs <- max(max_abs, 0.1)
    y_range <- c(-1.05 * max_abs, 1.05 * max_abs)
  }

  track <- DataTrack(
    start = sub$pos,
    end = sub$pos,
    chromosome = chrom,
    genome = "hg38",
    name = meta$output_type[[1]], #track_label,
    type = c("p", "g"),
    data = sub[[score_column]],
    #col = rgb(0, 0, 0.55, alpha = 0.3),
    col = "blue", # "darkblue",
    col.baseline = "black",
    cex = 0.8,
    baseline = 0,
    #lty.grid = 2,              # 破線のグリッド
    #col.grid = "red",       # グリッドの色
    #col.axis = "black",        # 軸の色
    #col.border.title = "gray80", # タイトル枠
    ylim = y_range             # Y軸を0を中心に対称に
  )
  tracks <- append(tracks, list(track))
}


# プロット
#pdf(opt$output, width=12, height=5 + length(tracks)*0.3)
png(opt$output, width=opt$width, height=opt$height + length(tracks)*100, res=800)
plotTracks(
  tracks,
  from=start,
  to=end,
  chromosome=chrom,
  #background.panel = "#FFFEDB",  # パネルごとの背景
  #background.panel =  bg_colors,  # パネルごとの背景
  #background.title = "lightgray",  # タイトルの背景
  col.axis = "black",
  col.grid = "gray80",
  col.border.title = "gray80",
  lty.grid = 2,
  lwd.grid = 0.5,
  showGrid = TRUE,
  frame = TRUE,
  col.frame = "gray"
)
dev.off()
