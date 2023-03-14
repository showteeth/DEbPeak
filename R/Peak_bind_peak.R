#' Integrate Two Peak Annotation/Differential Analysis Results.
#'
#' @param peak1.res Peak1 dataframe contains all peak annotation (\code{peak.mode} is consenus) or
#' differential analysis results of peak-related data (\code{peak.mode} is diff).
#' @param peak2.res Peak2 dataframe contains all peak annotation (\code{peak.mode} is consenus) or
#' differential analysis results of peak-related data (\code{peak.mode} is diff).
#' @param peak.mode The source of peak results, choose from consenus (peak annotation) and diff (differential analysis).
#' Default: consenus.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All".
#' Used when \code{peak.mode} is consenus. Default: "Promoter".
#' @param peak1.signif Used when \code{peak.mode} is diff. Significance criterion for peak-associated results for Peak1. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param peak1.signif.threshold Used when \code{peak.mode} is diff. Significance threshold for peak-associated results to get differentially accessible/binding peaks for Peak1. Default: 0.05.
#' @param peak1.l2fc.threshold Used when \code{peak.mode} is diff. Log2 fold change threshold for peak-associated results to get differentially accessible/binding peaks for Peak1. Default: 1.
#' @param peak2.signif Used when \code{peak.mode} is diff. Significance criterion for peak-associated results for Peak2. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param peak2.signif.threshold Used when \code{peak.mode} is diff. Significance threshold for peak-associated results to get differentially accessible/binding peaks for Peak2. Default: 0.05.
#' @param peak2.l2fc.threshold Used when \code{peak.mode} is diff. Log2 fold change threshold for peak-associated results to get differentially accessible/binding peaks for Peak2. Default: 1.
#'
#' @return Dataframe contains integration results. When \code{peak.mode} is 'diff', the 'Type' column contains "Down_Up", "Up_Up", "Down_Down", "Up_Down", "Peak1_Up", "Peak1_Down", "Peak2_Up", "Peak2_Down".
#' When \code{peak.mode} is 'consenus', the 'Type' column contains "Common", "Peak1", "Peak2".
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at vars filter
#' @importFrom purrr set_names
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' library(DEbPeak)
#' # ChIP-seq data
#' chip.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' chip.df <- GetConsensusPeak(peak.file = chip.file)
#' chip.anno <- AnnoPeak(
#'   peak.df = chip.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # ATAC-seq data
#' atac.file <- system.file("extdata", "debatac_peaks.bed", package = "DEbPeak")
#' atac.df <- GetConsensusPeak(peak.file = atac.file)
#' atac.anno <- AnnoPeak(
#'   peak.df = atac.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # integrate
#' chip.atac <- PeakbPeak(peak1.res = chip.anno$df, peak2.res = atac.anno$df, peak.mode = "consenus", peak.anno.key = "Promoter")
PeakbPeak <- function(peak1.res, peak2.res, peak.mode = c("consensus", "diff"),
                      peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                      peak1.signif = "padj", peak1.signif.threshold = 0.05, peak1.l2fc.threshold = 1,
                      peak2.signif = "padj", peak2.signif.threshold = 0.05, peak2.l2fc.threshold = 1) {
  # check parameters
  peak.mode <- match.arg(arg = peak.mode)
  peak.anno.key <- match.arg(arg = peak.anno.key)

  if (peak.mode == "consensus") {
    # prepare peak1
    peak1.res.used <- peak1.res %>%
      dplyr::filter(anno == peak.anno.key) %>%
      dplyr::mutate(Region = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>%
      dplyr::mutate(AnnoType = dplyr::case_when(
        anno == "Promoter" ~ "P",
        anno == "5' UTR" ~ "5U",
        anno == "3' UTR" ~ "3U",
        anno == "Exon" ~ "E",
        anno == "Intron" ~ "I",
        anno == "Downstream" ~ "D",
        anno == "Distal Intergenic" ~ "DI"
      )) %>%
      dplyr::mutate(Feature = paste(Region, SYMBOL, AnnoType, sep = "|")) %>%
      dplyr::select(c("Feature", "SYMBOL")) %>%
      purrr::set_names(c("Feature", "Gene"))
    colnames(peak1.res.used) <- gsub(pattern = "^", replacement = "P1_", x = colnames(peak1.res.used))
    # prepare peak2
    peak2.res.used <- peak2.res %>%
      dplyr::filter(anno == peak.anno.key) %>%
      dplyr::mutate(Region = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>%
      dplyr::mutate(AnnoType = dplyr::case_when(
        anno == "Promoter" ~ "P",
        anno == "5' UTR" ~ "5U",
        anno == "3' UTR" ~ "3U",
        anno == "Exon" ~ "E",
        anno == "Intron" ~ "I",
        anno == "Downstream" ~ "D",
        anno == "Distal Intergenic" ~ "DI"
      )) %>%
      dplyr::mutate(Feature = paste(Region, SYMBOL, AnnoType, sep = "|")) %>%
      dplyr::select(c("Feature", "SYMBOL")) %>%
      purrr::set_names(c("Feature", "Gene"))
    colnames(peak2.res.used) <- gsub(pattern = "^", replacement = "P2_", x = colnames(peak2.res.used))
    # merge results
    peak.peak <- merge(peak1.res.used, peak2.res.used, by.x = "P1_Gene", by.y = "P2_Gene", all = TRUE)
    # anno the results
    peak.peak <- peak.peak %>% dplyr::mutate(Type = dplyr::case_when(
      !is.na(P1_Feature) & !is.na(P2_Feature) ~ "Common",
      !is.na(P1_Feature) & is.na(P2_Feature) ~ "Peak1",
      is.na(P1_Feature) & !is.na(P2_Feature) ~ "Peak2"
    ))
    peak.peak$Type <- factor(peak.peak$Type, levels = c("Common", "Peak1", "Peak2"))
  } else if (peak.mode == "diff") {
    # prepare peak1 dataframe
    peak1.df <- PrepareDEPlot(
      deres = peak1.res, signif = peak1.signif, signif.threshold = peak1.signif.threshold,
      l2fc.threshold = peak1.l2fc.threshold, label.key = NULL
    )
    colnames(peak1.df) <- gsub(pattern = "^Gene", replacement = "Feature", x = colnames(peak1.df))
    # prepare peak2 dataframe
    peak2.df <- PrepareDEPlot(
      deres = peak2.res, signif = peak2.signif, signif.threshold = peak2.signif.threshold,
      l2fc.threshold = peak2.l2fc.threshold, label.key = NULL
    )
    colnames(peak2.df) <- gsub(pattern = "^Gene", replacement = "Feature", x = colnames(peak2.df))
    # filter with peak.anno.key
    if (peak.anno.key == "All") {
      # peak 1
      peak1.df <- peak1.df
      # peak 2
      peak2.df <- peak2.df
    } else {
      # extract peak annotated region
      peak1.df$PeakRegion <- gsub(pattern = ".*\\|.*\\|(.*)", replacement = "\\1", x = peak1.df$Feature)
      peak2.df$PeakRegion <- gsub(pattern = ".*\\|.*\\|(.*)", replacement = "\\1", x = peak2.df$Feature)
      anno.key.named <- c("P", "5U", "3U", "E", "I", "D", "DI")
      names(anno.key.named) <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
      peak1.df <- peak1.df[peak1.df$PeakRegion == anno.key.named[peak.anno.key], ]
      peak1.df$PeakRegion <- NULL
      peak2.df <- peak2.df[peak2.df$PeakRegion == anno.key.named[peak.anno.key], ]
      peak2.df$PeakRegion <- NULL
    }
    # extract de results
    peak1.deg.df <- peak1.df %>% dplyr::filter(regulation != "Not_regulated")
    # add gene info
    peak1.deg.df$Gene <- gsub(pattern = ".*\\|(.*)\\|.*", replacement = "\\1", x = peak1.deg.df$Feature)
    colnames(peak1.deg.df) <- gsub(pattern = "^", replacement = "P1_", x = colnames(peak1.deg.df))
    peak2.deg.df <- peak2.df %>% dplyr::filter(regulation != "Not_regulated")
    # add gene info
    peak2.deg.df$Gene <- gsub(pattern = ".*\\|(.*)\\|.*", replacement = "\\1", x = peak2.deg.df$Feature)
    colnames(peak2.deg.df) <- gsub(pattern = "^", replacement = "P2_", x = colnames(peak2.deg.df))
    # merge results
    peak.peak <- merge(peak1.deg.df, peak2.deg.df, by.x = "P1_Gene", by.y = "P2_Gene", all = TRUE)
    # anno the results
    peak.peak <- peak.peak %>% dplyr::mutate(Type = dplyr::case_when(
      P1_regulation == "Up_regulated" & P2_regulation == "Up_regulated" ~ "Up_Up",
      P1_regulation == "Up_regulated" & P2_regulation == "Down_regulated" ~ "Up_Down",
      P1_regulation == "Up_regulated" & is.na(P2_regulation) ~ "Peak1_Up",
      P1_regulation == "Down_regulated" & P2_regulation == "Up_regulated" ~ "Down_Up",
      P1_regulation == "Down_regulated" & P2_regulation == "Down_regulated" ~ "Down_Down",
      P1_regulation == "Down_regulated" & is.na(P2_regulation) ~ "Peak1_Down",
      is.na(P1_regulation) & P2_regulation == "Up_regulated" ~ "Peak2_Up",
      is.na(P1_regulation) & P2_regulation == "Down_regulated" ~ "Peak2_Down"
    ))
    peak.peak$Type <- factor(peak.peak$Type, levels = c(
      "Down_Up", "Up_Up", "Down_Down", "Up_Down",
      "Peak1_Up", "Peak1_Down", "Peak2_Up", "Peak2_Down"
    ))
  }
  return(peak.peak)
}


#' GO Enrichment on Two Peak Annotation/Differential Analysis Integration Results.
#'
#' @param peak.peak Dataframe contains integrated results of two peak annotation/differential analysis.
#' @param peak.fe.key The key type of integrated results ("Type" column of \code{peak.peak}) to perform functional enrichment.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param go.type GO enrichment type, chosen from ALL, BP, MF, CC. Default: ALL.
#' @param enrich.pvalue Cutoff value of pvalue. Default: 0.05.
#' @param enrich.qvalue Cutoff value of qvalue. Default: 0.05.
#' @param species Species used, chosen from "Human","Mouse","Rat","Fly","Arabidopsis","Yeast","Zebrafish","Worm","Bovine","Pig","Chicken","Rhesus",
#' "Canine","Xenopus","Anopheles","Chimp","E coli strain Sakai","Myxococcus xanthus DK 1622". Default: "Human".
#' @param padj.method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.
#' @param show.term Number of enrichment term to show. Default: 15.
#' @param str.width Length of enrichment term in plot. Default: 30.
#' @param plot.resolution Resolution of plot. Default: 300.
#' @param plot.width The width of plot. Default: 7.
#' @param plot.height The height of plot. Default: 9.
#' @param save Logical value, whether to save all results. Default: TRUE.
#'
#' @return If \code{save} is TRUE, return NULL (all results are in \code{out.folder}), else retutn result dataframe.
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr filter group_by top_n mutate arrange desc select pull
#' @importFrom magrittr %>%
#' @import clusterProfiler
#' @importFrom enrichplot dotplot
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom BiocManager install
#' @export
#'
#' @examples
#' library(DEbPeak)
#' # ChIP-seq data
#' chip.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' chip.df <- GetConsensusPeak(peak.file = chip.file)
#' chip.anno <- AnnoPeak(
#'   peak.df = chip.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # ATAC-seq data
#' atac.file <- system.file("extdata", "debatac_peaks.bed", package = "DEbPeak")
#' atac.df <- GetConsensusPeak(peak.file = atac.file)
#' atac.anno <- AnnoPeak(
#'   peak.df = atac.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # integrate
#' chip.atac <- PeakbPeak(peak1.res = chip.anno$df, peak2.res = atac.anno$df, peak.mode = "consenus", peak.anno.key = "Promoter")
#' # functional enrichment
#' chip.atac.fe <- PeakbPeakFE(
#'   peak.peak = chip.atac, peak.fe.key = "Common", gene.type = "SYMBOL",
#'   go.type = "BP", species = "Mouse", save = FALSE
#' )
PeakbPeakFE <- function(peak.peak, peak.fe.key, out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                        go.type = c("ALL", "BP", "MF", "CC"), enrich.pvalue = 0.05, enrich.qvalue = 0.05, species = c(
                          "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                          "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                        ),
                        padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                        show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  gene.type <- match.arg(arg = gene.type)
  go.type <- match.arg(arg = go.type)
  species <- match.arg(arg = species)
  padj.method <- match.arg(arg = padj.method)

  # get genes
  if (!peak.fe.key %in% as.character(unique(peak.peak$Type))) {
    stop(paste0(
      "Please provide valid functional enrichment key, choose from: ",
      paste(as.character(unique(peak.peak$Type)), collapse = ", ")
    ))
  }
  peak.genes <- peak.peak %>%
    dplyr::filter(Type == peak.fe.key) %>%
    dplyr::pull(P1_Gene)

  # prepare org db
  spe.anno <- GetSpeciesAnno(species)
  org.db <- spe.anno[["OrgDb"]]
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org_db : ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }

  # GO enrichment
  if (save) {
    SingleFE(
      genes = peak.genes, out.folder = out.folder, regulation = peak.fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(NULL)
  } else {
    debde.go.results <- SingleFE(
      genes = peak.genes, out.folder = out.folder, regulation = peak.fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(debde.go.results)
  }
}
