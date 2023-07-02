#' Integrate Differential Expression Results and Peak Annotation Results of ChIP-seq and ATAC-seq.
#'
#' @param de.res Data frame contains all genes of differential expression analysis.
#' @param chip.peak.res Dataframe contains peak annotation results of ChIP-seq data.
#' @param atac.peak.res Dataframe contains peak annotation results of ATAC-seq data.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes. Default: 1.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres}).
#' @param merge.key The columns used for merging, chosen from geneId, ENSEMBL, SYMBOL. Default: geneId.
#'
#' @return Dataframe contains integrated results. The results will contain eleven categories: UP, DOWN, ChIP, ATAC, UPbChIP, DOWNbChIP
#' UPbATAC, DOWNbATAC, ChIPbATAC, UPbPeak, DOWNbPeak.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select distinct mutate case_when mutate_at vars
#' @importFrom tibble rownames_to_column
#' @importFrom purrr set_names
#' @importFrom tidyr drop_na
#'
#' @examples
#' library(DEbPeak)
#' library(DESeq2)
#' # ChIP-Seq data
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.profile <- PeakProfile(peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
#' peak.anno <- AnnoPeak(
#'   peak.df = peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # RNA-Seq data
#' count.file <- system.file("extdata", "debchip_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "debchip_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' # create DESeqDataSet object
#' dds <- DESeq2::DESeqDataSetFromMatrix(
#'   countData = count.matrix, colData = meta.info,
#'   design = ~condition
#' )
#' # set control level
#' dds$condition <- relevel(dds$condition, ref = "NF")
#' # conduct differential expressed genes analysis
#' dds <- DESeq(dds)
#' # extract results
#' dds.results <- results(dds, contrast = c("condition", "RX", "NF"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' # Integrated with RNA-Seq
#' debchip.res <- DEbPeak(
#'   de.res = dds.results.ordered, peak.res = peak.anno[["df"]],
#'   peak.anno.key = "Promoter", merge.key = "SYMBOL"
#' )
#' # ATAC-seq data
#' atac.peak.file <- system.file("extdata", "debatac_peaks.bed", package = "DEbPeak")
#' atac.peak.df <- GetConsensusPeak(peak.file = atac.peak.file)
#' atac.peak.profile <- PeakProfile(atac.peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
#' atac.peak.anno <- AnnoPeak(
#'   peak.df = atac.peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' debpeak.res <- DEbCA(
#'   de.res = dds.results.ordered, chip.peak.res = peak.anno[["df"]],
#'   atac.peak.res = atac.peak.anno[["df"]],
#'   peak.anno.key = "Promoter", merge.key = "SYMBOL"
#' )
DEbCA <- function(de.res, chip.peak.res, atac.peak.res, peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                  signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, label.key = NULL, merge.key = c("geneId", "ENSEMBL", "SYMBOL")) {
  # check parameters
  peak.anno.key <- match.arg(arg = peak.anno.key)
  merge.key <- match.arg(arg = merge.key)

  # parepare peak results
  peak.alt.columns <- c("ENSEMBL", "SYMBOL", "GENENAME")
  ## process chip results
  chip.peak.alt.valid <- intersect(colnames(chip.peak.res), peak.alt.columns)
  chip.peak.df <- chip.peak.res %>%
    dplyr::filter(anno == peak.anno.key) %>%
    dplyr::select(c("geneId", "annotation", "anno", "seqnames", "start", "end", chip.peak.alt.valid)) %>%
    dplyr::distinct(geneId, anno, .keep_all = TRUE)
  ## process atac results
  atac.peak.alt.valid <- intersect(colnames(atac.peak.res), peak.alt.columns)
  atac.peak.df <- atac.peak.res %>%
    dplyr::filter(anno == peak.anno.key) %>%
    dplyr::mutate(Peak = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>%
    dplyr::select(c("geneId", "Peak", "annotation", "anno", atac.peak.alt.valid)) %>%
    dplyr::distinct(geneId, anno, .keep_all = TRUE)

  # Integrate RNA-seq and ChIP-seq results
  debchip.res <- DEbPeak(
    de.res = de.res, peak.res = chip.peak.df, peak.anno.key = peak.anno.key, signif = signif, signif.threshold = signif.threshold,
    l2fc.threshold = l2fc.threshold, label.key = label.key, merge.key = merge.key
  )
  # Integrate above results with ATAC-seq results
  # change debchip results name
  debchip.res[debchip.res$Type == "UPbPeak", "Type"] <- "UPbChIP"
  debchip.res[debchip.res$Type == "DOWNbPeak", "Type"] <- "DOWNbChIP"
  debchip.res[debchip.res$Type == "Peak", "Type"] <- "ChIP"
  colnames(debchip.res) <- gsub(pattern = "Type", replacement = "Type1", x = colnames(debchip.res))
  de.peak <- merge(debchip.res, atac.peak.df, by = merge.key, all = TRUE, suffixes = c("_ChIP", "_ATAC"))
  # eleven categories: UP, DOWN, ChIP, ATAC, UPbChIP, DOWNbChIP, UPbATAC, DOWNbATAC, ChIPbATAC, UPbPeak, DOWNbPeak
  de.peak <- de.peak %>%
    dplyr::mutate(Type = dplyr::case_when(
      is.na(annotation_ATAC) & Type1 == "UP" ~ "UP",
      is.na(annotation_ATAC) & Type1 == "DOWN" ~ "DOWN",
      is.na(annotation_ATAC) & Type1 == "ChIP" ~ "ChIP",
      is.na(annotation_ATAC) & Type1 == "UPbChIP" ~ "UPbChIP",
      is.na(annotation_ATAC) & Type1 == "DOWNbChIP" ~ "DOWNbChIP",
      !is.na(annotation_ATAC) & is.na(Type1) ~ "ATAC",
      !is.na(annotation_ATAC) & Type1 == "UP" ~ "UPbATAC",
      !is.na(annotation_ATAC) & Type1 == "DOWN" ~ "DOWNbATAC",
      !is.na(annotation_ATAC) & Type1 == "ChIP" ~ "ChIPbATAC",
      !is.na(annotation_ATAC) & Type1 == "UPbChIP" ~ "UPbPeak",
      !is.na(annotation_ATAC) & Type1 == "DOWNbChIP" ~ "DOWNbPeak"
    ))
  # process gene name
  de.peak <- de.peak %>% dplyr::mutate(geneId = ifelse(is.na(geneId_ChIP), geneId_ATAC, geneId_ChIP))
  return(de.peak)
}
