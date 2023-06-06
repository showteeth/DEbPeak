#' Integrate Differential Expression Results with Peak Annotation/Differential Expression Results.
#'
#' @param de.res Dataframe contains all genes of differential expression analysis.
#' @param peak.mode The source of peak results, choose from consensus (peak annotation) and diff (differential analysis).
#' Default: consensus.
#' @param peak.res Dataframe contains all peak annotation (\code{peak.mode} is consensus) or
#' differential analysis results of peak-related data (\code{peak.mode} is diff).
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All".
#' Used when \code{peak.mode} is consensus. Default: "Promoter".
#' @param signif Significance criterion for RNA-seq results. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold for RNA-seq to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold for RNA-seq to get differentially expressed genes. Default: 1.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{de.res}).
#' @param peak.signif Significance criterion for peak-associated results. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param peak.signif.threshold Significance threshold for peak-associated results to get differentially accessible/binding peaks. Default: 0.05.
#' @param peak.l2fc.threshold Log2 fold change threshold for peak-associated results to get differentially accessible/binding peaks. Default: 1.
#' @param merge.key The columns used for merging, chosen from geneId (ENTREZID), ENSEMBL, SYMBOL. Default: geneId.
#' @param org.db Organism database. Used when \code{peak.mode} is diff and \code{merge.key} is not "SYMBOL".
#' For peak-associated differential expression results, only support merging on "SYMBOL". Default: org.Mm.eg.db.
#'
#' @return Dataframe contains integrated results. The results will contain five categories: UPbPeak, DOWNbPeak, UP, DOWN, Peak when \code{peak.mode} is consensus
#' and Down_Up, Up_Up, Down_Down, Up_Down when \code{peak.mode} is diff. RNA in front.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select distinct mutate case_when mutate_at vars arrange desc
#' @importFrom purrr set_names
#' @importFrom tidyr drop_na
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom clusterProfiler bitr
#' @importFrom rlang .data
#'
#' @examples
#' library(DEbPeak)
#' library(DESeq2)
#' library(openxlsx)
#' ### consensus mode
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
#' ### diff mode
#' # ATAC-Seq data
#' peak.matrix.file <- system.file("extdata", "RA_ATAC_count.txt", package = "DEbPeak")
#' peak.meta.file <- system.file("extdata", "RA_ATAC_meta.txt", package = "DEbPeak")
#' peak.matrix <- read.table(file = peak.matrix.file, header = TRUE, sep = "\t")
#' peak.meta <- read.table(file = peak.meta.file, header = TRUE)
#' dds.peak <- DESeq2::DESeqDataSetFromMatrix(
#'   countData = peak.matrix, colData = peak.meta,
#'   design = ~condition
#' )
#' # set control level
#' dds.peak$condition <- relevel(dds.peak$condition, ref = "WT")
#' # conduct differential expressed genes analysis
#' dds.peak <- DESeq(dds.peak)
#' # extract results
#' dds.peak.results <- results(dds.peak, contrast = c("condition", "KO", "WT"))
#' dds.peak.results.ordered <- dds.peak.results[order(dds.peak.results$log2FoldChange, decreasing = TRUE), ]
#' # RNA-seq results
#' # rna.diff.file <- system.file("extdata", "RA_RNA_diff.xlsx", package = "DEbPeak")
#' rna.diff.file <- system.file("extdata", "RA_RNA_diff.txt", package = "DEbPeak")
#' # rna.diff <- openxlsx::read.xlsx(xlsxFile = rna.diff.file, rowNames = T, check.names = FALSE)
#' rna.diff <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' # integrate differential expression results of RNA-seq and ATAC-seq
#' debatac.res <- DEbPeak(
#'   de.res = rna.diff, peak.res = dds.peak.results.ordered, peak.mode = "diff", peak.anno.key = "All", l2fc.threshold = 0,
#'   peak.l2fc.threshold = 0, org.db = "org.Mm.eg.db", merge.key = "SYMBOL"
#' )
DEbPeak <- function(de.res, peak.res, peak.mode = c("consensus", "diff"),
                    peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                    signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, label.key = NULL,
                    peak.signif = "padj", peak.signif.threshold = 0.05, peak.l2fc.threshold = 1,
                    merge.key = c("geneId", "ENSEMBL", "SYMBOL"), org.db = "org.Mm.eg.db") {

  # check parameters
  peak.mode <- match.arg(arg = peak.mode)
  peak.anno.key <- match.arg(arg = peak.anno.key)
  merge.key <- match.arg(arg = merge.key)

  if (peak.mode == "consensus") {
    # process RNA de results
    de.df <- PrepareDEPlot(
      deres = de.res, signif = signif, signif.threshold = signif.threshold,
      l2fc.threshold = l2fc.threshold, label.key = label.key
    )
    # remove gene version information
    de.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = de.df$Gene)
    # get deg
    deg.df <- de.df %>% dplyr::filter(regulation != "Not_regulated")
    # parepare peak results
    peak.alt.columns <- c("ENSEMBL", "SYMBOL", "GENENAME")
    peak.alt.valid <- intersect(colnames(peak.res), peak.alt.columns)
    peak.df <- peak.res %>%
      dplyr::filter(anno == peak.anno.key) %>%
      dplyr::select(c("geneId", "annotation", "anno", peak.alt.valid)) %>%
      dplyr::distinct(geneId, anno, .keep_all = TRUE)

    # get all genes used
    if (!merge.key %in% colnames(peak.res)) {
      stop("The merge.key you provided is not valid!")
    }
    all.gene.used <- union(peak.df[, merge.key], deg.df$Gene)
    de.df.used <- de.df[de.df$Gene %in% all.gene.used, ]

    # merge DE and Peak results
    de.peak <- merge(peak.df, de.df.used, by.x = merge.key, by.y = "Gene", all = TRUE)
    # five categories: UPbPeak, DOWNbPeak, UP, DOWN, Peak
    de.peak <- de.peak %>%
      dplyr::mutate(Type = dplyr::case_when(
        is.na(annotation) & regulation == "Up_regulated" ~ "UP",
        is.na(annotation) & regulation == "Down_regulated" ~ "DOWN",
        !is.na(annotation) & (regulation == "Not_regulated" | is.na(regulation)) ~ "Peak",
        !is.na(annotation) & regulation == "Up_regulated" ~ "UPbPeak",
        !is.na(annotation) & regulation == "Down_regulated" ~ "DOWNbPeak"
      ))
  } else if (peak.mode == "diff") {
    if (length(intersect(colnames(de.res), c("ENSEMBL", "ENTREZID", "SYMBOL"))) == 0) {
      warning("In diff mode, the rownames of de table should be SYMBOL or SYMBOL in columns!
              You can use IDConversion to perfrom gene ID conversion!")
      # process RNA de results
      de.df <- PrepareDEPlot(
        deres = de.res, signif = signif, signif.threshold = signif.threshold,
        l2fc.threshold = l2fc.threshold, label.key = NULL
      )
    } else if ("SYMBOL" %in% colnames(de.res)) {
      # process RNA de results
      de.df <- PrepareDEPlot(
        deres = de.res, signif = signif, signif.threshold = signif.threshold,
        l2fc.threshold = l2fc.threshold, label.key = "SYMBOL"
      )
      # arrange by SYMBOL
      de.df <- de.df %>%
        dplyr::arrange(desc(abs(log2FoldChange))) %>%
        dplyr::filter(!is.na(SYMBOL)) %>%
        dplyr::distinct(SYMBOL, .keep_all = TRUE)
      colnames(de.df) <- gsub(pattern = "^Gene", replacement = "GeneID", x = colnames(de.df))
      colnames(de.df) <- gsub(pattern = "^SYMBOL", replacement = "Gene", x = colnames(de.df))
    } else {
      # process RNA de results
      de.df <- PrepareDEPlot(
        deres = de.res, signif = signif, signif.threshold = signif.threshold,
        l2fc.threshold = l2fc.threshold, label.key = NULL
      )
    }
    # remove gene version information
    de.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = de.df$Gene)
    # ID conversion for peak data
    peak.res <- IDConversionPeak(deres = peak.res, org.db = org.db, sort.key = "log2FoldChange")
    # prepare diff peak results
    peak.de.df <- PrepareDEPlot(
      deres = peak.res, signif = peak.signif, signif.threshold = peak.signif.threshold,
      l2fc.threshold = peak.l2fc.threshold, label.key = "SYMBOL"
    )
    # filter with peak.anno.key
    if (peak.anno.key == "All") {
      peak.de.df <- peak.de.df
    } else {
      peak.de.df$PeakRegion <- gsub(pattern = ".*\\|.*\\|(.*)", replacement = "\\1", x = peak.de.df$Gene)
      anno.key.named <- c("P", "5U", "3U", "E", "I", "D", "DI")
      names(anno.key.named) <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
      peak.de.df <- peak.de.df[peak.de.df$PeakRegion == anno.key.named[peak.anno.key], ]
      peak.de.df$PeakRegion <- NULL
    }
    # extract de results
    peak.deg.df <- peak.de.df %>% dplyr::filter(regulation != "Not_regulated")
    colnames(peak.deg.df) <- gsub(pattern = "^", replacement = "Peak_", x = colnames(peak.deg.df))
    # make sure Gene column are SYMBOL
    if (merge.key != "SYMBOL") {
      if (merge.key == "geneId") {
        merge.key <- "ENTREZID"
      }
      # gene ID conversion
      de.df <- IDConversion_internal(de.df = de.df, gene.type = merge.key, org.db = org.db, sort.key = "log2FoldChange")
      colnames(de.df) <- gsub(pattern = "^Gene$", replacement = merge.key, x = colnames(de.df))
      colnames(de.df) <- gsub(pattern = "^SYMBOL$", replacement = "Gene", x = colnames(de.df))
      colnames(de.df) <- gsub(pattern = "^ENTREZID$", replacement = "geneId", x = colnames(de.df))
    } else {
      # gene ID conversion
      de.df <- IDConversion_internal(de.df = de.df, gene.type = "SYMBOL", org.db = org.db, sort.key = "log2FoldChange")
      colnames(de.df) <- gsub(pattern = "^ENTREZID$", replacement = "geneId", x = colnames(de.df))
    }
    # get deg
    deg.df <- de.df %>% dplyr::filter(regulation != "Not_regulated")
    colnames(deg.df) <- gsub(pattern = "^", replacement = "RNA_", x = colnames(deg.df))

    # merge results
    de.peak <- merge(peak.deg.df, deg.df, by.x = "Peak_SYMBOL", by.y = "RNA_Gene", all = TRUE)
    # anno the results
    de.peak <- de.peak %>% dplyr::mutate(Type = dplyr::case_when(
      RNA_regulation == "Up_regulated" & Peak_regulation == "Up_regulated" ~ "Up_Up",
      RNA_regulation == "Up_regulated" & Peak_regulation == "Down_regulated" ~ "Up_Down",
      RNA_regulation == "Up_regulated" & is.na(Peak_regulation) ~ "RNAUp",
      RNA_regulation == "Down_regulated" & Peak_regulation == "Up_regulated" ~ "Down_Up",
      RNA_regulation == "Down_regulated" & Peak_regulation == "Down_regulated" ~ "Down_Down",
      RNA_regulation == "Down_regulated" & is.na(Peak_regulation) ~ "RNADown",
      is.na(RNA_regulation) & Peak_regulation == "Up_regulated" ~ "PeakUp",
      is.na(RNA_regulation) & Peak_regulation == "Down_regulated" ~ "PeakDown"
    ))
    de.peak$Type <- factor(de.peak$Type, levels = c(
      "Down_Up", "Up_Up", "Down_Down", "Up_Down",
      "RNAUp", "RNADown", "PeakUp", "PeakDown"
    ))
    colnames(de.peak) <- gsub(pattern = "^RNA_geneId$", replacement = "geneId", x = colnames(de.peak))
  }
  return(de.peak)
}

# prepare venn plot dataframe
PrepareVenn <- function(key, named.vec) {
  if (key %in% names(named.vec)) {
    if (named.vec[key] >= 1) {
      return(paste0(key, 1:named.vec[key]))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' Create Integrated Summary Plot.
#'
#' @param de.peak Dataframe contains integrated results.
#' @param peak.type The source of peaks, chosen from ATAC, ChIP and Peak (ChIP and ATAC). Default: ChIP.
#' @param peak.mode The source of peak results, choose from consensus (peak annotation) and diff (differential expression analysis).
#' Default: consensus.
#' @param gene.col Column of \code{inte.res} contains genes. Same as \code{merge.key} in \code{\link{DEbPeak}}.
#' @param ... Parameters for \code{\link{ggvenn}}.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble deframe
#' @import ggvenn
#' @export
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
#' # DE and ChIP venn plot
#' debchip.plot <- PlotDEbPeak(
#'   de.peak = debchip.res, peak.type = "ChIP", gene.col = "SYMBOL",
#'   show_percentage = FALSE
#' )
PlotDEbPeak <- function(de.peak, peak.type = c("ChIP", "ATAC", "Peak"), peak.mode = c("consensus", "diff"),
                        gene.col = c("geneId", "ENSEMBL", "SYMBOL"), ...) {
  # check parameters
  peak.type <- match.arg(arg = peak.type)

  # get summary results
  # distinct Type and gene, show gene number instead of peak number
  if (peak.type %in% c("ChIP", "ATAC")) {
    if (peak.mode == "consensus") {
      de.peak <- de.peak %>%
        dplyr::distinct(.data[[gene.col]], .data[["Type"]], .keep_all = TRUE)
    } else if (peak.mode == "diff") {
      de.peak <- de.peak %>%
        dplyr::distinct(Peak_SYMBOL, Type, .keep_all = TRUE)
    }
  } else if (peak.type == "Peak") {
    de.peak <- de.peak %>%
      dplyr::distinct(.data[[gene.col]], .data[["Type"]], .keep_all = TRUE)
  }
  type.summary <- table(de.peak$Type) %>%
    as.data.frame() %>%
    tibble::deframe()
  if (peak.type == "Peak") {
    # create number vector
    ChIP.vec <- c(
      PrepareVenn("ChIP", type.summary), PrepareVenn("ChIPbATAC", type.summary),
      PrepareVenn("UPbChIP", type.summary), PrepareVenn("DOWNbChIP", type.summary),
      PrepareVenn("UPbPeak", type.summary), PrepareVenn("DOWNbPeak", type.summary)
    )
    ATAC.vec <- c(
      PrepareVenn("ATAC", type.summary), PrepareVenn("ChIPbATAC", type.summary),
      PrepareVenn("UPbATAC", type.summary), PrepareVenn("DOWNbATAC", type.summary),
      PrepareVenn("UPbPeak", type.summary), PrepareVenn("DOWNbPeak", type.summary)
    )
    UP.vec <- c(
      PrepareVenn("UP", type.summary), PrepareVenn("UPbChIP", type.summary),
      PrepareVenn("UPbATAC", type.summary), PrepareVenn("UPbPeak", type.summary)
    )
    DOWN.vec <- c(
      PrepareVenn("DOWN", type.summary), PrepareVenn("DOWNbChIP", type.summary),
      PrepareVenn("DOWNbATAC", type.summary), PrepareVenn("DOWNbPeak", type.summary)
    )
    # create plot
    plot.list <- list()
    plot.list[["UP"]] <- UP.vec
    plot.list[["DOWN"]] <- DOWN.vec
    plot.list[["ChIP"]] <- ChIP.vec
    plot.list[["ATAC"]] <- ATAC.vec
    plot <- ggvenn::ggvenn(plot.list, ...)
  } else if (peak.type %in% c("ChIP", "ATAC")) {
    if (peak.mode == "consensus") {
      # create number vector
      Peak.vec <- c(
        PrepareVenn("Peak", type.summary), PrepareVenn("UPbPeak", type.summary),
        PrepareVenn("DOWNbPeak", type.summary)
      )
      UP.vec <- c(PrepareVenn("UP", type.summary), PrepareVenn("UPbPeak", type.summary))
      DOWN.vec <- c(PrepareVenn("DOWN", type.summary), PrepareVenn("DOWNbPeak", type.summary))
      # create plot list
      plot.list <- list()
      plot.list[[peak.type]] <- Peak.vec
      plot.list[["UP"]] <- UP.vec
      plot.list[["DOWN"]] <- DOWN.vec
    } else if (peak.mode == "diff") {
      # create number vector
      # peak up related
      Peak.up.vec <- c(
        PrepareVenn("Down_Up", type.summary), PrepareVenn("Up_Up", type.summary),
        PrepareVenn("PeakUp", type.summary)
      )
      # peak down related
      Peak.down.vec <- c(
        PrepareVenn("Down_Down", type.summary), PrepareVenn("Up_Down", type.summary),
        PrepareVenn("PeakDown", type.summary)
      )
      # RNA up related
      RNA.up.vec <- c(
        PrepareVenn("Up_Up", type.summary), PrepareVenn("Up_Down", type.summary),
        PrepareVenn("RNAUp", type.summary)
      )
      # RNA down related
      RNA.down.vec <- c(
        PrepareVenn("Down_Up", type.summary), PrepareVenn("Down_Down", type.summary),
        PrepareVenn("RNADown", type.summary)
      )
      # create plot list
      plot.list <- list()
      plot.list[[paste(peak.type, "UP", sep = " ")]] <- Peak.up.vec
      plot.list[[paste(peak.type, "DOWN", sep = " ")]] <- Peak.down.vec
      plot.list[["RNA UP"]] <- RNA.up.vec
      plot.list[["RNA DOWN"]] <- RNA.down.vec
    }
    plot <- ggvenn::ggvenn(plot.list, ...)
  }
  return(plot)
}

#' Create Venn Diagram for Two Differential Analysis Integration Results.
#'
#' @param inte.res Integration results, can be output of \code{DEbPeak}, \code{PeakbPeak}, \code{DEbDE}.
#' @param inte.type The integration type, choose from "DEbDE", "PeakbPeak", "DEbPeak". Default: "DEbPeak".
#' @param peak.type Used when \code{inte.type} is "DEbPeak". The source of peaks, chosen from ATAC, ChIP and Peak (ChIP and ATAC). Default: ChIP.
#' @param peak.mode Used when \code{inte.type} is "DEbPeak" or "PeakbPeak". The source of peak results, choose from consensus (peak annotation) and diff (differential expression analysis).
#' Default: consensus.
#' @param gene.col Used when \code{inte.type} is "DEbPeak" and \code{peak.type} is Peak or \code{peak.mode} is "consensus" in \code{peak.type} ATAC or ChIP.
#' Column of \code{inte.res} contains genes. Same as \code{merge.key} in \code{\link{DEbPeak}}.
#' @param ... Parameters for \code{\link{ggvenn}}.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble deframe
#' @importFrom dplyr distinct
#' @importFrom rlang .data
#' @import ggvenn
#' @export
#'
#' @examples
#' library(DEbPeak)
#' #### RNA-seq and RNA-seq
#' rna.diff.file <- system.file("extdata", "RA_RNA_diff.txt", package = "DEbPeak")
#' de1.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' de2.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' # use same file as example
#' de.de <- DEbDE(de1.res = de1.res, de2.res = de2.res, de1.l2fc.threshold = 0.5, de2.l2fc.threshold = 1)
#' de.de.venn <- InteVenn(inte.res = de.de, inte.type = "DEbDE", show_percentage = FALSE)
#' #### peak-related and peak-related
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
#' chip.atac <- PeakbPeak(peak1.res = chip.anno$df, peak2.res = atac.anno$df, peak.mode = "consensus", peak.anno.key = "Promoter")
#' # functional enrichment
#' chip.atac.venn <- InteVenn(inte.res = chip.atac, inte.type = "PeakbPeak", peak.mode = "consensus", show_percentage = FALSE)
#' #### RNA-seq and peak-related
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
#' debatac.res.venn <- InteVenn(
#'   inte.res = debatac.res, inte.type = "DEbPeak", peak.mode = "diff",
#'   peak.type = "ATAC", show_percentage = FALSE
#' )
InteVenn <- function(inte.res, inte.type = c("DEbPeak", "PeakbPeak", "DEbDE"), peak.type = c("ChIP", "ATAC", "Peak"),
                     peak.mode = c("consensus", "diff"), gene.col = c("geneId", "ENSEMBL", "SYMBOL"), ...) {
  # check parameters
  inte.type <- match.arg(arg = inte.type)
  peak.type <- match.arg(arg = peak.type)
  gene.col <- match.arg(arg = gene.col)
  peak.mode <- match.arg(arg = peak.mode)

  # get summary results
  if (inte.type == "DEbPeak") {
    # distinct Type and gene, show gene number instead of peak number
    if (peak.type %in% c("ChIP", "ATAC")) {
      if (peak.mode == "consensus") {
        inte.res <- inte.res %>%
          dplyr::distinct(.data[[gene.col]], .data[["Type"]], .keep_all = TRUE)
      } else if (peak.mode == "diff") {
        inte.res <- inte.res %>%
          dplyr::distinct(Peak_SYMBOL, Type, .keep_all = TRUE)
      }
    } else if (peak.type == "Peak") {
      inte.res <- inte.res %>%
        dplyr::distinct(.data[[gene.col]], .data[["Type"]], .keep_all = TRUE)
    }
  } else if (inte.type == "PeakbPeak") {
    inte.res <- inte.res %>%
      dplyr::distinct(P1_Gene, Type, .keep_all = TRUE)
  } else if (inte.type == "DEbDE") {
    inte.res <- inte.res %>%
      dplyr::distinct(DE1_Gene, Type, .keep_all = TRUE)
  }
  type.summary <- table(inte.res$Type) %>%
    as.data.frame() %>%
    tibble::deframe()

  # create venn plot
  if (inte.type == "DEbPeak") {
    if (peak.type == "Peak") {
      # create number vector
      ChIP.vec <- c(
        PrepareVenn("ChIP", type.summary), PrepareVenn("ChIPbATAC", type.summary),
        PrepareVenn("UPbChIP", type.summary), PrepareVenn("DOWNbChIP", type.summary),
        PrepareVenn("UPbPeak", type.summary), PrepareVenn("DOWNbPeak", type.summary)
      )
      ATAC.vec <- c(
        PrepareVenn("ATAC", type.summary), PrepareVenn("ChIPbATAC", type.summary),
        PrepareVenn("UPbATAC", type.summary), PrepareVenn("DOWNbATAC", type.summary),
        PrepareVenn("UPbPeak", type.summary), PrepareVenn("DOWNbPeak", type.summary)
      )
      UP.vec <- c(
        PrepareVenn("UP", type.summary), PrepareVenn("UPbChIP", type.summary),
        PrepareVenn("UPbATAC", type.summary), PrepareVenn("UPbPeak", type.summary)
      )
      DOWN.vec <- c(
        PrepareVenn("DOWN", type.summary), PrepareVenn("DOWNbChIP", type.summary),
        PrepareVenn("DOWNbATAC", type.summary), PrepareVenn("DOWNbPeak", type.summary)
      )
      # create plot
      plot.list <- list()
      plot.list[["UP"]] <- UP.vec
      plot.list[["DOWN"]] <- DOWN.vec
      plot.list[["ChIP"]] <- ChIP.vec
      plot.list[["ATAC"]] <- ATAC.vec
      plot <- ggvenn::ggvenn(plot.list, ...)
    } else if (peak.type %in% c("ChIP", "ATAC")) {
      if (peak.mode == "consensus") {
        # create number vector
        Peak.vec <- c(
          PrepareVenn("Peak", type.summary), PrepareVenn("UPbPeak", type.summary),
          PrepareVenn("DOWNbPeak", type.summary)
        )
        UP.vec <- c(PrepareVenn("UP", type.summary), PrepareVenn("UPbPeak", type.summary))
        DOWN.vec <- c(PrepareVenn("DOWN", type.summary), PrepareVenn("DOWNbPeak", type.summary))
        # create plot list
        plot.list <- list()
        plot.list[[peak.type]] <- Peak.vec
        plot.list[["UP"]] <- UP.vec
        plot.list[["DOWN"]] <- DOWN.vec
      } else if (peak.mode == "diff") {
        # create number vector
        # peak up related
        Peak.up.vec <- c(
          PrepareVenn("Down_Up", type.summary), PrepareVenn("Up_Up", type.summary),
          PrepareVenn("PeakUp", type.summary)
        )
        # peak down related
        Peak.down.vec <- c(
          PrepareVenn("Down_Down", type.summary), PrepareVenn("Up_Down", type.summary),
          PrepareVenn("PeakDown", type.summary)
        )
        # RNA up related
        RNA.up.vec <- c(
          PrepareVenn("Up_Up", type.summary), PrepareVenn("Up_Down", type.summary),
          PrepareVenn("RNAUp", type.summary)
        )
        # RNA down related
        RNA.down.vec <- c(
          PrepareVenn("Down_Up", type.summary), PrepareVenn("Down_Down", type.summary),
          PrepareVenn("RNADown", type.summary)
        )
        # create plot list
        plot.list <- list()
        plot.list[[paste(peak.type, "UP", sep = " ")]] <- Peak.up.vec
        plot.list[[paste(peak.type, "DOWN", sep = " ")]] <- Peak.down.vec
        plot.list[["RNA UP"]] <- RNA.up.vec
        plot.list[["RNA DOWN"]] <- RNA.down.vec
      }
      plot <- ggvenn::ggvenn(plot.list, ...)
    }
  } else if (inte.type == "PeakbPeak") {
    if (peak.mode == "consensus") {
      peak1.vec <- c(
        PrepareVenn("Peak1", type.summary), PrepareVenn("Common", type.summary)
      )
      peak2.vec <- c(
        PrepareVenn("Peak2", type.summary), PrepareVenn("Common", type.summary)
      )
      # create plot list
      plot.list <- list()
      plot.list[["Peak1"]] <- peak1.vec
      plot.list[["Peak2"]] <- peak2.vec
    } else if (peak.mode == "diff") {
      # create number vector
      # Peak1 up related
      Peak1.up.vec <- c(
        PrepareVenn("Up_Up", type.summary), PrepareVenn("Up_Down", type.summary),
        PrepareVenn("Peak1_Up", type.summary)
      )
      # Peak1 down related
      Peak1.down.vec <- c(
        PrepareVenn("Down_Up", type.summary), PrepareVenn("Down_Down", type.summary),
        PrepareVenn("Peak1_Down", type.summary)
      )
      # Peak2 up related
      Peak2.up.vec <- c(
        PrepareVenn("Down_Up", type.summary), PrepareVenn("Up_Up", type.summary),
        PrepareVenn("Peak2_Up", type.summary)
      )
      # Peak2 down related
      Peak2.down.vec <- c(
        PrepareVenn("Down_Down", type.summary), PrepareVenn("Up_Down", type.summary),
        PrepareVenn("Peak2_Down", type.summary)
      )
      # create plot list
      plot.list <- list()
      plot.list[["Peak1 UP"]] <- Peak1.up.vec
      plot.list[["Peak1 DOWN"]] <- Peak1.down.vec
      plot.list[["Peak2 UP"]] <- Peak2.up.vec
      plot.list[["Peak2 DOWN"]] <- Peak2.down.vec
    }
    plot <- ggvenn::ggvenn(plot.list, ...)
  } else if (inte.type == "DEbDE") {
    # create number vector
    # DE1 up related
    DE1.up.vec <- c(
      PrepareVenn("Up_Up", type.summary), PrepareVenn("Up_Down", type.summary),
      PrepareVenn("DE1_Up", type.summary)
    )
    # DE1 down related
    DE1.down.vec <- c(
      PrepareVenn("Down_Up", type.summary), PrepareVenn("Down_Down", type.summary),
      PrepareVenn("DE1_Down", type.summary)
    )
    # DE2 up related
    DE2.up.vec <- c(
      PrepareVenn("Down_Up", type.summary), PrepareVenn("Up_Up", type.summary),
      PrepareVenn("DE2_Up", type.summary)
    )
    # DE2 down related
    DE2.down.vec <- c(
      PrepareVenn("Down_Down", type.summary), PrepareVenn("Up_Down", type.summary),
      PrepareVenn("DE2_Down", type.summary)
    )
    # create plot list
    plot.list <- list()
    plot.list[["DE1 UP"]] <- DE1.up.vec
    plot.list[["DE1 DOWN"]] <- DE1.down.vec
    plot.list[["DE2 UP"]] <- DE2.up.vec
    plot.list[["DE2 DOWN"]] <- DE2.down.vec
    plot <- ggvenn::ggvenn(plot.list, ...)
  }
  return(plot)
}

#' Create Quadrant Diagram for Two Differential Analysis Integration Results.
#'
#' @param inte.res Integration results, can be output of \code{DEbPeak}, \code{PeakbPeak}, \code{DEbDE}.
#' @param inte.type The integration type, choose from "DEbDE", "PeakbPeak", "DEbPeak". Default: "DEbPeak".
#' @param point.alpha Opacity of a geom. Default: 0.6.
#' @param point.size.vec Point size for regular and(or) labeled points (specified by \code{label.df} or \code{label.num}).
#' Default: 2 for regular and 4 for labeled points.
#' @param f1.l2fc.threshold Log2 fold change threshold for file1 to get differential analysis results. Default: 0.
#' @param f2.l2fc.threshold Log2 fold change threshold for file2 to get differential analysis results. Default: 0.
#' @param linetype Threshold linetype. Default: 2.
#' @param point.color.vec Point color for Down_Up, Up_Up, Down_Down, Up_Down.
#' Default: grey for Down_Up, red for Up_Up, blue for Down_Down and grey for Up_Down.
#' @param legend.pos Legend position. Default: top.
#' @param show.corr Logical value, whether to add pearson correlation coefficient and its significance. Default: TRUE.
#' @param label.num Gene number to label, according to RNA_log2FoldChange and Peak_log2FoldChange. When \code{label.df} is set NULL,
#' use this to determine genes to label. Default: NULL.
#' @param label.df Label data frame, at least contains Gene column. Default: NULL(use \code{label.num}).
#' When provided, the second column should not be in \code{de.peak}.
#' @param label.color Color for labels. Default: NULL (black).
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange desc
#' @importFrom utils head
#' @importFrom ggpubr stat_cor
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' @import ggplot2
#' @export
#'
#' @examples
#' # RNA-seq and RNA-seq
#' # de.de.label.df = data.frame(Gene = c("ENSDARG00000007396", "ENSDARG00000010729", "ENSDARG00000002194", "ENSDARG00000002587"))
#' # de.de.quad = InteDiffQuad(inte.res = de.de, inte.type = "DEbDE", f1.l2fc.threshold = 0.6, f2.l2fc.threshold = 0.6,
#' #                           show.corr = FALSE, label.df = de.de.label.df)
#' # Peak and Peak
#' # peak.peak.label.df = data.frame(Gene = c("Aak1", "Adam19", "Oaf", "Nsg2"))
#' # peak.peak.quad = InteDiffQuad(inte.res = atac.atac, inte.type = "PeakbPeak", f1.l2fc.threshold = 0,
#' #                               f2.l2fc.threshold = 0.5, show.corr = TRUE, label.df = peak.peak.label.df)
#' # RNA-seq and Peak
#' # de.peak.label.df = data.frame(Gene = c("Ccl3", "Ccl5", "Cd28", "Cx3cr1", "Prdm1", "Tcf7", "Slamf6", "Id3", "Cxcr5"))
#' # de.peak.quad = InteDiffQuad(inte.res = debatac.res, inte.type = "DEbPeak", f1.l2fc.threshold = 0,
#' #                             f2.l2fc.threshold = 0, show.corr = TRUE, label.df = de.peak.label.df)
InteDiffQuad <- function(inte.res, inte.type = c("DEbPeak", "PeakbPeak", "DEbDE"), point.alpha = 0.6, point.size.vec = c(2, 4),
                         f1.l2fc.threshold = 0, f2.l2fc.threshold = 0, linetype = 2,
                         point.color.vec = c("grey", "red", "blue", "grey"), legend.pos = "top", show.corr = TRUE,
                         label.num = NULL, label.df = NULL, label.color = NULL) {
  # check parameters
  inte.type <- match.arg(arg = inte.type)

  # prepare integrate dataframe: filter out single regulation
  if (inte.type == "DEbDE") {
    inte.res <- inte.res %>%
      dplyr::filter(!is.na(DE1_regulation)) %>%
      dplyr::filter(!is.na(DE2_regulation))
    ax.x <- "DE1_log2FoldChange"
    ax.y <- "DE2_log2FoldChange"
    merge.key <- "DE1_Gene"
  } else if (inte.type == "PeakbPeak") {
    inte.res <- inte.res %>%
      dplyr::filter(!is.na(P1_regulation)) %>%
      dplyr::filter(!is.na(P2_regulation))
    ax.x <- "P1_log2FoldChange"
    ax.y <- "P2_log2FoldChange"
    merge.key <- "P1_Gene"
  } else if (inte.type == "DEbPeak") {
    inte.res <- inte.res %>%
      dplyr::filter(!is.na(Peak_regulation)) %>%
      dplyr::filter(!is.na(RNA_regulation))
    ax.x <- "RNA_log2FoldChange"
    ax.y <- "Peak_log2FoldChange"
    merge.key <- "Peak_SYMBOL"
  }
  inte.res$Type <- droplevels(inte.res$Type)

  # prepare point size
  if (length(point.size.vec) == 1) {
    point.size.vec <- rep(point.size.vec, 2)
  }
  # prepare basic plot
  p <- ggplot() +
    geom_point(
      data = inte.res,
      aes_string(x = ax.x, y = ax.y, color = "Type"),
      size = point.size.vec[1], alpha = point.alpha
    ) +
    scale_color_manual(values = point.color.vec) +
    geom_vline(xintercept = c(-f1.l2fc.threshold, f1.l2fc.threshold), lty = linetype) +
    geom_hline(yintercept = c(-f2.l2fc.threshold, f2.l2fc.threshold), lty = linetype) +
    theme_classic(base_size = 14) +
    theme(legend.position = legend.pos) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1), add = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1), add = c(0, 0)))
  # add correlation
  if (show.corr) {
    p <- p + stat_cor(
      data = inte.res,
      aes_string(x = ax.x, y = ax.y), method = "pearson"
    )
  }
  # add label
  if (is.null(label.df)) {
    if (!is.null(label.num)) {
      label.data.up <- inte.res %>%
        dplyr::arrange(desc(.data[[ax.x]]), desc(.data[[ax.y]])) %>%
        head(label.num)
      label.data.down <- inte.res %>%
        dplyr::arrange(.data[[ax.x]], .data[[ax.y]]) %>%
        head(label.num)
      label.data <- as.data.frame(rbind(label.data.up, label.data.down))
    } else {
      return(p)
    }
  } else if (nrow(label.df) >= 1) {
    label.data <- merge(inte.res, label.df, by.x = merge.key, by.y = "Gene") %>% as.data.frame()
  }
  if (is.null(label.color)) {
    p <- p +
      geom_point(
        data = label.data,
        aes_string(x = ax.x, y = ax.y),
        size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes_string(x = ax.x, y = ax.y, label = merge.key)
      )
  } else {
    p <- p +
      geom_point(
        data = label.data,
        aes_string(x = ax.x, y = ax.y),
        color = label.color, size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes_string(x = ax.x, y = ax.y, label = merge.key)
      )
  }
  return(p)
}

#' GO Enrichment on Integrated Results.
#'
#' @param inte.res Integration results, can be output of \code{DEbPeak}, \code{PeakbPeak}, \code{DEbDE}.
#' @param fe.key The key type of integrated results ("Type" column of \code{inte.res}) to perform functional enrichment.
#' @param inte.type The integration type, choose from "DEbDE", "PeakbPeak", "DEbPeak". Default: "DEbPeak".
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param gene.type Gene name type (if \code{inte.res} is from \code{DEbPeak}, this should be ENTREZID; if \code{inte.res} is from \code{PeakbPeak},
#' this should be Gene name type of P1_Gene; if \code{inte.res} is from \code{DEbDE}, this should be Gene name type of DE1_Gene).
#' Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
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
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select case_when mutate_at vars filter
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @importFrom tidyr separate
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' library(DEbPeak)
#' #### RNA-seq and RNA-seq
#' rna.diff.file <- system.file("extdata", "RA_RNA_diff.txt", package = "DEbPeak")
#' de1.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' de2.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' # use same file as example
#' de.de <- DEbDE(de1.res = de1.res, de2.res = de2.res, de1.l2fc.threshold = 0.5, de2.l2fc.threshold = 1)
#' de.de.fe <- InteFE(
#'   inte.res = de.de, fe.key = "Down_Down", inte.type = "DEbDE", gene.type = "SYMBOL", go.type = "BP",
#'   species = "Mouse", save = F
#' )
#' #### peak-related and peak-related
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
#' chip.atac <- PeakbPeak(peak1.res = chip.anno$df, peak2.res = atac.anno$df, peak.mode = "consensus", peak.anno.key = "Promoter")
#' # functional enrichment
#' chip.atac.fe <- InteFE(
#'   inte.res = chip.atac, fe.key = "Common", inte.type = "PeakbPeak", gene.type = "SYMBOL",
#'   go.type = "BP", species = "Mouse", save = FALSE
#' )
#' #### RNA-seq and peak-related
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
#' # functional enrichment on UPbPeak genes
#' upbpeak.fe.results <- InteFE(
#'   inte.res = debchip.res, fe.key = "UPbPeak", inte.type = "DEbPeak", gene.type = "ENTREZID",
#'   species = "Mouse", save = FALSE
#' )
InteFE <- function(inte.res, fe.key, inte.type = c("DEbPeak", "PeakbPeak", "DEbDE"), out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                   go.type = c("ALL", "BP", "MF", "CC"), enrich.pvalue = 0.05, enrich.qvalue = 0.05, species = c(
                     "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                     "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                   ),
                   padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                   show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  inte.type <- match.arg(arg = inte.type)
  gene.type <- match.arg(arg = gene.type)
  go.type <- match.arg(arg = go.type)
  species <- match.arg(arg = species)
  padj.method <- match.arg(arg = padj.method)

  # valid parameter value
  if (!fe.key %in% as.character(unique(inte.res$Type))) {
    stop(paste0(
      "Please provide valid functional enrichment key, choose from: ",
      paste(as.character(unique(inte.res$Type)), collapse = ", ")
    ))
  }

  # prepare genes
  if (inte.type == "DEbDE") {
    inte.genes <- inte.res %>%
      dplyr::filter(Type == fe.key) %>%
      dplyr::pull(DE1_Gene)
  } else if (inte.type == "PeakbPeak") {
    inte.genes <- inte.res %>%
      dplyr::filter(Type == fe.key) %>%
      dplyr::pull(P1_Gene)
  } else if (inte.type == "DEbPeak") {
    inte.genes <- inte.res %>%
      dplyr::filter(Type == fe.key) %>%
      dplyr::pull(geneId)
    # check gene type
    if (gene.type != "ENTREZID") {
      warning("To perform GO enrichment, the gene type should be ENTREZID!")
    }
    gene.type <- "ENTREZID"
  }

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
      genes = inte.genes, out.folder = out.folder, regulation = fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(NULL)
  } else {
    inte.go.results <- SingleFE(
      genes = inte.genes, out.folder = out.folder, regulation = fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(inte.go.results)
  }
}


#' GO Enrichment on Integrated Results.
#'
#' @param de.peak Dataframe contains integrated results.
#' @param peak.fe.key The key type of integrated results ("Type" column of \code{de.peak}) to perform functional enrichment.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
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
#' # functional enrichment on UPbPeak genes
#' upbpeak.fe.results <- DEbPeakFE(
#'   de.peak = debchip.res, peak.fe.key = "UPbPeak",
#'   species = "Mouse", save = FALSE
#' )
#' # functional enrichment on DOWNbPeak genes
#' downbpeak.fe.results <- DEbPeakFE(
#'   de.peak = debchip.res, peak.fe.key = "DOWNbPeak",
#'   species = "Mouse", save = FALSE
#' )
DEbPeakFE <- function(de.peak, peak.fe.key, out.folder = NULL,
                      go.type = c("ALL", "BP", "MF", "CC"), enrich.pvalue = 0.05, enrich.qvalue = 0.05, species = c(
                        "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                        "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                      ),
                      padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                      show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  go.type <- match.arg(arg = go.type)
  species <- match.arg(arg = species)
  padj.method <- match.arg(arg = padj.method)

  # get genes
  if (!peak.fe.key %in% as.character(unique(de.peak$Type))) {
    stop(paste0(
      "Please provide valid functional enrichment key, choose from: ",
      paste(as.character(unique(de.peak$Type)), collapse = ", ")
    ))
  }
  inte.genes <- de.peak %>%
    dplyr::filter(Type == peak.fe.key) %>%
    dplyr::pull(geneId)
  # UPbPeak.genes <- de.peak %>%
  #   dplyr::filter(Type == "UPbPeak") %>%
  #   dplyr::pull(geneId)
  # DOWNbPeak.genes <- de.peak %>%
  #   dplyr::filter(Type == "DOWNbPeak") %>%
  #   dplyr::pull(geneId)

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

  # regulation string
  # up.reg.str <- paste0("UPb", peak.type)
  # down.reg.str <- paste0("DOWNb", peak.type)

  if (save) {
    SingleFE(
      genes = inte.genes, out.folder = out.folder, regulation = peak.fe.key, gene.type = "ENTREZID", enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # # for positive
    # SingleFE(
    #   genes = UPbPeak.genes, out.folder = out.folder, regulation = up.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
    #   enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
    #   show.term = show.term, str.width = str.width, save = save
    # )
    # # for negative
    # SingleFE(
    #   genes = DOWNbPeak.genes, out.folder = out.folder, regulation = down.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
    #   enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
    #   show.term = show.term, str.width = str.width, save = save
    # )
    return(NULL)
  } else {
    # DEbPeak.go.results <- list()
    # # for positive
    # UPbPeak.go <- SingleFE(
    #   genes = UPbPeak.genes, out.folder = out.folder, regulation = up.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
    #   enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
    #   show.term = show.term, str.width = str.width, save = save
    # )
    # # for negative
    # DOWNbPeak.go <- SingleFE(
    #   genes = DOWNbPeak.genes, out.folder = out.folder, regulation = down.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
    #   enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
    #   show.term = show.term, str.width = str.width, save = save
    # )
    # DEbPeak.go.results[[up.reg.str]] <- UPbPeak.go
    # DEbPeak.go.results[[down.reg.str]] <- DOWNbPeak.go
    DEbPeak.go.results <- SingleFE(
      genes = inte.genes, out.folder = out.folder, regulation = peak.fe.key, gene.type = "ENTREZID", enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(DEbPeak.go.results)
  }
}
