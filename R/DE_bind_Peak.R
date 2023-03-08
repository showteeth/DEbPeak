#' Integrate Differential Expression Results and Peak Annotation/Differential Expression Results.
#'
#' @param de.res Data frame contains all genes of differential expression analysis.
#' @param peak.mode The source of peak results, choose from consenus (peak annotation) and diff (differential expression analysis).
#' Default: consenus.
#' @param peak.res Dataframe contains all peak annotation (\code{peak.mode} is consenus) or
#' differential analysis results of peak-related data (\code{peak.mode} is diff).
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All".
#' Used when \code{peak.mode} is consenus. Default: "Promoter".
#' @param signif Significance criterion for RNA-seq results. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold for RNA-seq to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold for RNA-seq to get differentially expressed genes. Default: 1.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres}).
#' @param peak.signif Significance criterion for peak-associated results. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param peak.signif.threshold Significance threshold for peak-associated results to get differentially accessible/binding peaks. Default: 0.05.
#' @param peak.l2fc.threshold Log2 fold change threshold for peak-associated results to get differentially accessible/binding peaks. Default: 1.
#' @param merge.key The columns used for merging, chosen from geneId (ENTREZID), ENSEMBL, SYMBOL. Default: geneId.
#' @param org.db Organism database. Used when \code{peak.mode} is diff and \code{merge.key} is not "SYMBOL".
#' For peak-associated differential expression results, only support merging on "SYMBOL". Default: org.Mm.eg.db.
#'
#' @return Dataframe contains integrated results. The results will contain five categories: UPbPeak, DOWNbPeak, UP, DOWN, Peak when \code{peak.mode} is consenus
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
#' rna.diff.file <- system.file("extdata", "RA_RNA_diff.txt", package = "DEbPeak")
#' rna.diff <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' # integrate differential expression results of RNA-seq and ATAC-seq
#' debatac.res <- DEbPeak(
#'   de.res = rna.diff, peak.res = dds.peak.results.ordered, peak.mode = "diff", peak.anno.key = "All", l2fc.threshold = 0,
#'   peak.l2fc.threshold = 0, org.db = "org.Mm.eg.db", merge.key = "SYMBOL"
#' )
DEbPeak <- function(de.res, peak.res, peak.mode = c("consenus", "diff"),
                    peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                    signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, label.key = NULL,
                    peak.signif = "padj", peak.signif.threshold = 0.05, peak.l2fc.threshold = 1,
                    merge.key = c("geneId", "ENSEMBL", "SYMBOL"), org.db = "org.Mm.eg.db") {

  # check parameters
  peak.mode <- match.arg(arg = peak.mode)
  peak.anno.key <- match.arg(arg = peak.anno.key)
  merge.key <- match.arg(arg = merge.key)

  # process RNA de results
  de.df <- PrepareDEPlot(
    deres = de.res, signif = signif, signif.threshold = signif.threshold,
    l2fc.threshold = l2fc.threshold, label.key = label.key
  )
  # remove gene version information
  de.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = de.df$Gene)

  if (peak.mode == "consenus") {
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
    return(paste0(key, 1:named.vec[key]))
  } else {
    return(NULL)
  }
}

#' Create Integrated Summary Plot.
#'
#' @param de.peak Dataframe contains integrated results.
#' @param peak.type The source of peaks, chosen from ATAC, ChIP and Peak (ChIP and ATAC). Default: ChIP.
#' @param peak.mode The source of peak results, choose from consenus (peak annotation) and diff (differential expression analysis).
#' Default: consenus.
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
#' debchip.plot <- PlotDEbPeak(de.peak = debchip.res, peak.type = "ChIP", show_percentage = FALSE)
PlotDEbPeak <- function(de.peak, peak.type = c("ChIP", "ATAC", "Peak"), peak.mode = c("consenus", "diff"), ...) {
  # check parameters
  peak.type <- match.arg(arg = peak.type)

  # get summary results
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
    if (peak.mode == "consenus") {
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

#' Create Quadrant Diagram for Differential Expression Analysis of RNA-seq and Peak-related Data.
#'
#' @param de.peak Dataframe contains integrated results of differential analysis of RNA-seq and Peak-related data.
#' @param peak.type Peak data type, choose from ChIP, ATAC. Default: ChIP.
#' @param point.alpha Opacity of a geom. Default: 0.6.
#' @param point.size.vec Point size for regular(DE or non-DE) and(or) labeled points.
#' Default: 2 for regular and 4 for labeled points.
#' @param rna.l2fc.threshold Log2 fold change threshold for RNA-seq to get differentially expressed results. Default: 0.
#' @param peak.l2fc.threshold Log2 fold change threshold for peak-related data to get differentially accessible/binding peaks. Default: 0.
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
#' @import ggplot2
#' @export
#'
#' @examples
#' # # prepare integrated differential expression results of RNA-seq and ATAC-seq
#' # # label dataframe
#' # label.df = data.frame(Gene = c("Ccl3", "Ccl5", "Cd28", "Cx3cr1", "Prdm1", "Tcf7", "Slamf6", "Id3", "Cxcr5"))
#' # # create plot
#' # quad.plot = InteDiffQuad(de.peak = debatac.res, peak.type = "ATAC", label.df = label.df, label.color = "green")
InteDiffQuad <- function(de.peak, peak.type = c("ChIP", "ATAC"), point.alpha = 0.6, point.size.vec = c(2, 4), rna.l2fc.threshold = 0,
                         peak.l2fc.threshold = 0, linetype = 2, point.color.vec = c("grey", "red", "blue", "grey"),
                         legend.pos = "top", show.corr = TRUE, label.num = NULL, label.df = NULL, label.color = NULL) {
  # check parameters
  peak.type <- match.arg(arg = peak.type)

  # prepare integrate dataframe: filter out single regulation
  de.peak <- de.peak %>%
    dplyr::filter(!is.na(Peak_regulation)) %>%
    dplyr::filter(!is.na(RNA_regulation))
  # prepare point size
  if (length(point.size.vec) == 1) {
    point.size.vec <- rep(point.size.vec, 2)
  }
  # prepare basic plot
  p <- ggplot() +
    geom_point(
      data = de.peak,
      aes(x = RNA_log2FoldChange, y = Peak_log2FoldChange, color = Type),
      size = point.size.vec[1], alpha = point.alpha
    ) +
    scale_color_manual(values = point.color.vec) +
    geom_vline(xintercept = c(-rna.l2fc.threshold, rna.l2fc.threshold), lty = linetype) +
    geom_hline(yintercept = c(-peak.l2fc.threshold, peak.l2fc.threshold), lty = linetype) +
    theme_classic(base_size = 14) +
    theme(legend.position = legend.pos) +
    labs(x = "RNA log2Foldchange", y = paste(peak.type, "log2Foldchange", sep = " "), color = "") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1), add = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1), add = c(0, 0)))
  # add correlation
  if (show.corr) {
    p <- p + stat_cor(
      data = de.peak,
      aes(x = RNA_log2FoldChange, y = Peak_log2FoldChange), method = "pearson"
    )
  }
  # add label
  if (is.null(label.df)) {
    if (!is.null(label.num)) {
      label.data.up <- de.peak %>%
        dplyr::arrange(desc(RNA_log2FoldChange), desc(Peak_log2FoldChange)) %>%
        head(label.num)
      label.data.down <- de.peak %>%
        dplyr::arrange(RNA_log2FoldChange, Peak_log2FoldChange) %>%
        head(label.num)
      label.data <- as.data.frame(rbind(label.data.up, label.data.down))
    } else {
      return(p)
    }
  } else if (nrow(label.df) >= 1) {
    label.data <- merge(de.peak, label.df, by.x = "Peak_SYMBOL", by.y = "Gene") %>% as.data.frame()
  }
  if (is.null(label.color)) {
    p <- p +
      geom_point(
        data = label.data,
        aes(x = RNA_log2FoldChange, y = Peak_log2FoldChange),
        size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(x = RNA_log2FoldChange, y = Peak_log2FoldChange, label = Peak_SYMBOL)
      )
  } else {
    p <- p +
      geom_point(
        data = label.data,
        aes(x = RNA_log2FoldChange, y = Peak_log2FoldChange),
        color = label.color, size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(x = RNA_log2FoldChange, y = Peak_log2FoldChange, label = Peak_SYMBOL)
      )
  }
  return(p)
}

#' GO Enrichment on Integrated Results.
#'
#' @param de.peak Dataframe contains integrated results.
#' @param peak.fe.key The key type of integrated results ("Type" column of \code{de.peak}) to perform functional enrichment.
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
#'   de.peak = debchip.res, peak.fe.key = "UPbPeak", gene.type = "ENTREZID",
#'   species = "Mouse", save = FALSE
#' )
#' # functional enrichment on DOWNbPeak genes
#' downbpeak.fe.results <- DEbPeakFE(
#'   de.peak = debchip.res, peak.fe.key = "DOWNbPeak", gene.type = "ENTREZID",
#'   species = "Mouse", save = FALSE
#' )
DEbPeakFE <- function(de.peak, peak.fe.key, out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
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
      genes = inte.genes, out.folder = out.folder, regulation = peak.fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
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
      genes = inte.genes, out.folder = out.folder, regulation = peak.fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(DEbPeak.go.results)
  }
}
