#' Extract Differential Analysis Results.
#'
#' @param deres Data frame contains all genes/peaks.
#' @param data.type Input data type, choose from RNA, ChIP, ATAC. Default: RNA.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#'
#' @return A data frame.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate case_when mutate_at vars filter
#' @importFrom tidyr drop_na separate
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' dds.degs <- ExtractDA(
#'   deres = dds.results.ordered, signif = "padj", signif.threshold = 0.05,
#'   l2fc.threshold = 1
#' )
ExtractDA <- function(deres, data.type = c("RNA", "ChIP", "ATAC"),
                      peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                      signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1) {
  # check parameters
  data.type <- match.arg(arg = data.type)

  # make sure the input is dataframe
  deres <- as.data.frame(deres)
  # prepare for different data type
  if (data.type != "RNA") {
    # seperate the data
    deres <- deres %>%
      tibble::rownames_to_column(var = "Gene") %>%
      dplyr::filter(!is.na(Gene)) %>%
      dplyr::mutate(name = Gene) %>%
      tidyr::separate(col = Gene, into = c("Peak", "PeakGene", "Annotation"), sep = "\\|") %>%
      tidyr::separate(col = Peak, into = c("Chr", "region"), sep = ":") %>%
      tidyr::separate(col = region, into = c("Start", "End"), sep = "-")
    if (peak.anno.key == "All") {
      deres <- deres
    } else {
      # change annotation name
      anno.key.named <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
      names(anno.key.named) <- c("P", "5U", "3U", "E", "I", "D", "DI")
      deres$Annotation <- anno.key.named[deres$Annotation]
      # filter annotation
      deres <- deres[deres$Annotation == peak.anno.key, ]
    }
    # change row names
    rownames(deres) <- deres$name
    deres$name <- NULL
  }
  # extract the results
  if ("FDR" %in% colnames(deres)) {
    message("Differential expression analysis with edgeR!")
    if (!signif %in% colnames(deres)) {
      stop(paste0(signif, " is not valid!"))
    }
    de.df <- deres %>% tibble::rownames_to_column(var = "Gene")
    # deal with pvalue or padj is zero, min(pvalue)*0.1 or min(padj)*0.1, https://github.com/kevinblighe/EnhancedVolcano/issues/43
    signif.min <- min(de.df[, signif][de.df[, signif] > 0], na.rm = TRUE)
    de.df <- de.df %>% dplyr::mutate(!!signif := ifelse(.data[[signif]] == 0,
      signif.min * 0.1, .data[[signif]]
    ))
    # classfy results
    de.df <- de.df %>%
      dplyr::mutate(regulation = dplyr::case_when(
        logFC > l2fc.threshold & .data[[signif]] < signif.threshold ~ "Up_regulated",
        logFC < -l2fc.threshold & .data[[signif]] < signif.threshold ~ "Down_regulated",
        abs(logFC) <= l2fc.threshold | .data[[signif]] >= signif.threshold ~ "Not_regulated"
      )) %>%
      dplyr::mutate_at(vars(regulation), as.factor) %>%
      tidyr::drop_na(regulation)
    de.df[[signif]] <- -log10(de.df[[signif]])
  } else if ("padj" %in% colnames(deres)) {
    message("Differential expression analysis with DESeq2!")
    if (!signif %in% colnames(deres)) {
      stop(paste0(signif, " is not valid!"))
    }
    de.df <- deres %>% tibble::rownames_to_column(var = "Gene")
    # deal with pvalue or padj is zero, min(pvalue)*0.1 or min(padj)*0.1, https://github.com/kevinblighe/EnhancedVolcano/issues/43
    signif.min <- min(de.df[, signif][de.df[, signif] > 0], na.rm = TRUE)
    de.df <- de.df %>% dplyr::mutate(!!signif := ifelse(.data[[signif]] == 0,
      signif.min * 0.1, .data[[signif]]
    ))
    # classfy results
    de.df <- de.df %>%
      dplyr::mutate(regulation = dplyr::case_when(
        log2FoldChange > l2fc.threshold & .data[[signif]] < signif.threshold ~ "Up_regulated",
        log2FoldChange < -l2fc.threshold & .data[[signif]] < signif.threshold ~ "Down_regulated",
        abs(log2FoldChange) <= l2fc.threshold | .data[[signif]] >= signif.threshold ~ "Not_regulated"
      )) %>%
      dplyr::mutate_at(vars(regulation), as.factor) %>%
      tidyr::drop_na(regulation)
    de.df[[signif]] <- -log10(de.df[[signif]])
  }
  # modify colnames
  if (data.type != "RNA") {
    colnames(de.df) <- gsub(pattern = "^Gene$", replacement = "Peak", x = colnames(de.df))
  }
  return(de.df)
}

# classify up, down and not differentially expressed
PrepareDEPlot <- function(deres, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, label.key = NULL) {
  deres <- as.data.frame(deres)
  if ("FDR" %in% colnames(deres)) {
    message("Differential expression analysis with edgeR!")
    if (!signif %in% colnames(deres)) {
      stop(paste0(signif, " is not valid!"))
    }
    # modify edgeR results
    if (is.null(label.key)) {
      de.df <- deres %>%
        tibble::rownames_to_column(var = "Gene") %>%
        dplyr::select(c("Gene", "logFC", "logCPM", `signif`)) %>%
        purrr::set_names(c("Gene", "log2FoldChange", "abundance", "signif"))
    } else {
      de.df <- deres %>%
        tibble::rownames_to_column(var = "Gene") %>%
        dplyr::select(c("Gene", "logFC", "logCPM", `signif`, `label.key`)) %>%
        purrr::set_names(c("Gene", "log2FoldChange", "abundance", "signif", label.key))
    }
  } else if ("padj" %in% colnames(deres)) {
    message("Differential expression analysis with DESeq2!")
    if (!signif %in% colnames(deres)) {
      stop(paste0(signif, " is not valid!"))
    }
    if (is.null(label.key)) {
      de.df <- deres %>%
        tibble::rownames_to_column(var = "Gene") %>%
        dplyr::select(c("Gene", "log2FoldChange", "baseMean", `signif`)) %>%
        purrr::set_names(c("Gene", "log2FoldChange", "abundance", "signif"))
    } else {
      de.df <- deres %>%
        tibble::rownames_to_column(var = "Gene") %>%
        dplyr::select(c("Gene", "log2FoldChange", "baseMean", `signif`, `label.key`)) %>%
        purrr::set_names(c("Gene", "log2FoldChange", "abundance", "signif", label.key))
    }
  }
  # deal with pvalue or padj is zero, min(pvalue)*0.1 or min(padj)*0.1, https://github.com/kevinblighe/EnhancedVolcano/issues/43
  signif.min <- min(de.df$signif[de.df$signif > 0])
  de.df <- de.df %>% dplyr::mutate(signif = ifelse(signif == 0, signif.min * 0.1, signif))
  # filter na
  de.df <- de.df %>%
    tidyr::drop_na(signif) %>%
    tidyr::drop_na(log2FoldChange)
  # anno the results
  de.df <- de.df %>%
    dplyr::mutate(regulation = dplyr::case_when(
      log2FoldChange > l2fc.threshold & signif < signif.threshold ~ "Up_regulated",
      log2FoldChange < -l2fc.threshold & signif < signif.threshold ~ "Down_regulated",
      abs(log2FoldChange) <= l2fc.threshold | signif >= signif.threshold ~ "Not_regulated"
    )) %>%
    dplyr::mutate_at(vars(regulation), as.factor) %>%
    tidyr::drop_na(regulation)
  de.df$signif <- -log10(de.df$signif)

  return(de.df)
}

#' VolcanoPlot for Differential Analysis Results.
#'
#' @param deres Data frame contains all genes/peaks.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param point.alpha Opacity of a geom. Default: 0.6.
#' @param point.size.vec Point size for regular (all points) and(or) labeled points.
#' Default: 2 for regular and 4 for labeled points.
#' @param linetype Threshold linetype. Default: 2.
#' @param point.color.vec Point color for Down, Not, Up regulated results.
#' Default: red for Up, grey for Not and blue for Down.
#' @param legend.pos Legend position. Default: top.
#' @param label.num Gene/Peak number to label, choose according to log2FoldChange.
#' When \code{label.df} is set NULL, use this to determine genes to label. Default: NULL.
#' @param label.df Label data frame, at least contains Gene column. Default: NULL(use \code{label.num}).
#' When provided, the second column should not be in \code{deres}.
#' @param label.key Which column to use as label.
#' Default: NULL (use rownames of \code{deres} or Gene column of \code{gene.df}).
#' @param label.color Color vector for labels. Default: NULL.
#' @param tick.trans Scale transformations for y axis. Default: sqrt.
#' @param ticklabel.break Y axis tick label break. Default: NULL.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at filter top_n vars
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' VolcanoPlot(dds.results.ordered, signif = "pvalue", l2fc.threshold = 0.3, label.num = 2, point.alpha = 0.8, label.color = c("purple", "green"), tick.trans = NULL)
VolcanoPlot <- function(deres, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1,
                        point.alpha = 0.6, point.size.vec = c(2, 4), linetype = 2, point.color.vec = c("blue", "grey", "red"), legend.pos = "top",
                        label.num = NULL, label.df = NULL, label.key = NULL, label.color = NULL,
                        tick.trans = NULL, ticklabel.break = NULL) {
  # preapare DE dataframe
  de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = label.key)
  up.de.df <- de.df %>% dplyr::filter(regulation == "Up_regulated")
  down.de.df <- de.df %>% dplyr::filter(regulation == "Down_regulated")

  if (length(point.size.vec) == 1) {
    point.size.vec <- rep(point.size.vec, 2)
  }
  p <- ggplot() +
    geom_point(
      data = de.df,
      aes(log2FoldChange, signif, color = regulation),
      size = point.size.vec[1], alpha = point.alpha
    ) +
    scale_color_manual(values = point.color.vec) +
    geom_vline(xintercept = c(-l2fc.threshold, l2fc.threshold), lty = linetype) +
    geom_hline(yintercept = -log10(signif.threshold), lty = linetype) +
    theme_classic(base_size = 14) +
    theme(legend.position = legend.pos) +
    labs(x = "log2FoldChange", y = paste0("-log10(", signif, ")"), color = "")
  if (!is.null(tick.trans)) {
    if (!is.null(ticklabel.break)) {
      p <- p +
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.1), add = c(0, 0)),
          trans = tick.trans, breaks = ticklabel.break
        )
    } else {
      p <- p +
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.1), add = c(0, 0)),
          trans = tick.trans
        )
    }
  } else {
    p <- p +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = c(0, 0)))
  }
  if (is.null(label.df)) {
    if (!is.null(label.num)) {
      label.data.up <- up.de.df %>% dplyr::top_n(n = label.num, wt = log2FoldChange)
      label.data.down <- down.de.df %>% dplyr::top_n(n = label.num, wt = -log2FoldChange)
      if (is.null(label.key)) {
        label.key <- "Label"
        label.data.up[, label.key] <- label.data.up$Gene
        label.data.down[, label.key] <- label.data.down$Gene
      }
      label.data <- as.data.frame(rbind(label.data.up, label.data.down))
    } else {
      return(p)
    }
  } else if (nrow(label.df) >= 1) {
    if (is.null(label.key)) {
      label.key <- "Label"
      label.df[, label.key] <- label.df$Gene
    }
    label.data.up <- merge(up.de.df, label.df, by = "Gene") %>% as.data.frame()
    label.data.down <- merge(down.de.df, label.df, by = "Gene") %>% as.data.frame()
    label.data <- as.data.frame(rbind(label.data.up, label.data.down))
  }
  if (is.null(label.color)) {
    p <- p +
      geom_point(
        data = label.data,
        aes(log2FoldChange, signif),
        size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(log2FoldChange, signif),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 1) {
    p <- p +
      geom_point(
        data = label.data,
        aes(log2FoldChange, signif),
        color = label.color, size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(log2FoldChange, signif),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 2) {
    p <- p +
      geom_point(
        data = label.data.up,
        aes(log2FoldChange, signif),
        color = label.color[1], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.up,
        aes(log2FoldChange, signif),
        label = label.data.up[, label.key]
      ) +
      geom_point(
        data = label.data.down,
        aes(log2FoldChange, signif),
        color = label.color[2], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.down,
        aes(log2FoldChange, signif),
        label = label.data.down[, label.key]
      )
  }
  return(p)
}

#' ScatterPlot for Differential Analysis Results.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param deres Data frame contains all genes/peaks.
#' @param group.key Sample group information. When set NULL, select first column of metadata. Default: NULL.
#' @param ref.group Reference group name. When set NULL, select first element of groups. Default: NULL.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Default: 10.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param point.alpha Opacity of a geom. Default: 0.6.
#' @param point.size.vec Point size for regular (all points) and(or) labeled points. Default: 2 for regular and 4 for labeled points.
#' @param linetype Diagonal linetype. Default: 2.
#' @param point.color.vec Point color for Down, Not, Up regulated genes peaks. Default: red for Up, grey for Not and blue for Down.
#' @param legend.pos Legend position. Default: top.
#' @param label.num Gene/Peak number to label, choose according to log2FoldChange. When \code{label.df} is set NULL, use this to determine genes to label. Default: NULL.
#' @param label.df Label data frame, at least contains Gene column. Default: NULL(use \code{label.num}). When provided,
#' the second column should not be in \code{deres}.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres} or Gene column of \code{gene.df}).
#' @param label.color Color vector for labels. Default: NULL.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at pull filter top_n vars
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom purrr set_names
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' ScatterPlot(deobj = dds, deres = dds.results.ordered, group.key = "condition", ref.group = "WT", signif = "pvalue", l2fc.threshold = 0.3, label.num = 2, point.alpha = 0.8, label.color = c("purple", "green"))
ScatterPlot <- function(deobj, deres, group.key = NULL, ref.group = NULL, base = 10,
                        signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1,
                        point.alpha = 0.6, point.size.vec = c(2, 4), linetype = 2, point.color.vec = c("blue", "grey", "red"), legend.pos = "top",
                        label.num = NULL, label.df = NULL, label.key = NULL, label.color = NULL) {
  # identify analysis method
  if (class(deobj) == "DGEList") {
    counts.matrix <- edgeR::cpm(deobj, normalized.lib.sizes = TRUE)
    metadata <- as.data.frame(deobj$samples) %>% tibble::rownames_to_column(var = "Sample")
    group.key <- "group"
  } else if (class(deobj) == "DESeqDataSet") {
    if (!"sizeFactor" %in% colnames(SummarizedExperiment::colData(deobj))) {
      deobj <- DESeq2::estimateSizeFactors(deobj)
    }
    counts.matrix <- DESeq2::counts(deobj, normalized = TRUE)
    metadata <- as.data.frame(SummarizedExperiment::colData(deobj)) %>% tibble::rownames_to_column(var = "Sample")
  } else {
    stop("Input object is either DESeq2 or edgeR results!")
  }
  # get group key
  if (is.null(group.key)) {
    group.key <- colnames(metadata)[2]
  } else if (!group.key %in% colnames(metadata)) {
    stop(paste0("group.key you provided is not in ", colnames(metadata)))
  }
  # get ref condition
  groups <- unique(as.character(metadata[, group.key]))
  if (is.null(ref.group)) {
    ref.group <- groups[1]
  } else if (!ref.group %in% groups) {
    stop(paste0("ref.group you provided is not in ", groups))
  }
  # get mean matrix
  treat.group <- setdiff(groups, ref.group)
  ref.samples <- metadata[metadata[, group.key] == ref.group, ] %>% dplyr::pull(Sample)
  treat.samples <- metadata[metadata[, group.key] == treat.group, ] %>% dplyr::pull(Sample)
  ref.matrix <- counts.matrix[, ref.samples]
  treat.matrix <- counts.matrix[, treat.samples]
  ref.mean.matrix <- apply(ref.matrix, 1, mean) %>%
    as.data.frame() %>%
    purrr::set_names(c("Group1"))
  treat.mean.matrix <- apply(treat.matrix, 1, mean) %>%
    as.data.frame() %>%
    purrr::set_names(c("Group2"))
  mean.matrix <- merge(ref.mean.matrix, treat.mean.matrix, by = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "Row.names")
  # get log transformed matrix
  log.mean.matrix <- round(log(mean.matrix + 1, base = base), digits = 4)
  # preapare DE dataframe
  de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = label.key)
  # merge info
  de.df <- merge(de.df, log.mean.matrix, by.x = "Gene", by.y = 0)

  # create basic plot
  if (length(point.size.vec) == 1) {
    point.size.vec <- rep(point.size.vec, 2)
  }
  p <- ggplot() +
    geom_point(
      data = de.df,
      aes(Group1, Group2, color = regulation),
      size = point.size.vec[1], alpha = point.alpha
    ) +
    scale_color_manual(values = point.color.vec) +
    theme_classic(base_size = 14) +
    theme(panel.background = element_rect(color = "black")) +
    annotation_logticks(base = base) +
    labs(
      x = paste0("Log", base, "(Mean of Normalized Counts (", ref.group, "))"),
      y = paste0("Log", base, "(Mean of Normalized Counts (", treat.group, "))"), color = ""
    ) +
    theme(legend.position = legend.pos) +
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = linetype, size = 0.5)
  # get up and down regulated data
  up.de.df <- de.df %>% dplyr::filter(regulation == "Up_regulated")
  down.de.df <- de.df %>% dplyr::filter(regulation == "Down_regulated")
  if (is.null(label.df)) {
    if (!is.null(label.num)) {
      label.data.up <- up.de.df %>% dplyr::top_n(n = label.num, wt = log2FoldChange)
      label.data.down <- down.de.df %>% dplyr::top_n(n = label.num, wt = -log2FoldChange)
      if (is.null(label.key)) {
        label.key <- "Label"
        label.data.up[, label.key] <- label.data.up$Gene
        label.data.down[, label.key] <- label.data.down$Gene
      }
      label.data <- as.data.frame(rbind(label.data.up, label.data.down))
    } else {
      return(p)
    }
  } else if (nrow(label.df) >= 1) {
    if (is.null(label.key)) {
      label.key <- "Label"
      label.df[, label.key] <- label.df$Gene
    }
    label.data.up <- merge(up.de.df, label.df, by = "Gene") %>% as.data.frame()
    label.data.down <- merge(down.de.df, label.df, by = "Gene") %>% as.data.frame()
    label.data <- as.data.frame(rbind(label.data.up, label.data.down))
  }
  if (is.null(label.color)) {
    p <- p +
      geom_point(
        data = label.data,
        aes(Group1, Group2),
        size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(Group1, Group2),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 1) {
    p <- p +
      geom_point(
        data = label.data,
        aes(Group1, Group2),
        color = label.color, size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(Group1, Group2),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 2) {
    p <- p +
      geom_point(
        data = label.data.up,
        aes(Group1, Group2),
        color = label.color[1], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.up,
        aes(Group1, Group2),
        label = label.data.up[, label.key]
      ) +
      geom_point(
        data = label.data.down,
        aes(Group1, Group2),
        color = label.color[2], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.down,
        aes(Group1, Group2),
        label = label.data.down[, label.key]
      )
  }
  return(p)
}

#' MA-plot for Differential Analysis Results.
#'
#' @param deres Data frame contains all genes/peaks.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param point.alpha Opacity of a geom. Default: 0.6.
#' @param point.size.vec Point size for regular (all points) and(or) labeled points. Default: 2 for regular and 4 for labeled points.
#' @param linetype Threshold linetype. Default: 2.
#' @param point.color.vec Point color for Down, Not, Up regulated genes or peaks. Default: red for Up, grey for Not and blue for Down.
#' @param legend.pos Legend position. Default: top.
#' @param label.num Gene/Peak number to label, choose according to log2FoldChange. When \code{label.df} is set NULL, use this to determine genes to label. Default: NULL.
#' @param label.df Label data frame, at least contains Gene column. Default: NULL(use \code{label.num}). When provided,
#' the second column should not be in \code{deres}.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres} or Gene column of \code{gene.df}).
#' @param label.color Color vector for labels. Default: NULL.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at filter top_n vars
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' MAPlot(dds.results.ordered, signif = "pvalue", l2fc.threshold = 0.3, label.num = 2, point.alpha = 0.8, label.color = c("purple", "green"))
MAPlot <- function(deres, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1,
                   point.alpha = 0.6, point.size.vec = c(2, 4), linetype = 2, point.color.vec = c("blue", "grey", "red"), legend.pos = "top",
                   label.num = NULL, label.df = NULL, label.key = NULL, label.color = NULL) {
  # preapare DE dataframe
  de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = label.key)
  # determine analysis method
  if ("FDR" %in% colnames(deres)) {
    x.label <- "Average log CPM"
  } else if ("padj" %in% colnames(deres)) {
    x.label <- "Log2(mean of normalized counts)"
    de.df$abundance <- round(log2(de.df$abundance + 1), digits = 4)
  }
  up.de.df <- de.df %>% dplyr::filter(regulation == "Up_regulated")
  down.de.df <- de.df %>% dplyr::filter(regulation == "Down_regulated")
  # create plot
  if (length(point.size.vec) == 1) {
    point.size.vec <- rep(point.size.vec, 2)
  }
  p <- ggplot() +
    geom_point(
      data = de.df,
      aes(abundance, log2FoldChange, color = regulation),
      size = point.size.vec[1], alpha = point.alpha
    ) +
    scale_color_manual(values = point.color.vec) +
    geom_hline(yintercept = c(-l2fc.threshold, l2fc.threshold), lty = linetype) +
    theme_classic(base_size = 14) +
    theme(legend.position = legend.pos) +
    labs(x = x.label, color = "") +
    theme(panel.background = element_rect(color = "black"))
  if (is.null(label.df)) {
    if (!is.null(label.num)) {
      label.data.up <- up.de.df %>% dplyr::top_n(n = label.num, wt = log2FoldChange)
      label.data.down <- down.de.df %>% dplyr::top_n(n = label.num, wt = -log2FoldChange)
      if (is.null(label.key)) {
        label.key <- "Label"
        label.data.up[, label.key] <- label.data.up$Gene
        label.data.down[, label.key] <- label.data.down$Gene
      }
      label.data <- as.data.frame(rbind(label.data.up, label.data.down))
    } else {
      return(p)
    }
  } else if (nrow(label.df) >= 1) {
    if (is.null(label.key)) {
      label.key <- "Label"
      label.df[, label.key] <- label.df$Gene
    }
    label.data.up <- merge(up.de.df, label.df, by = "Gene") %>% as.data.frame()
    label.data.down <- merge(down.de.df, label.df, by = "Gene") %>% as.data.frame()
    label.data <- as.data.frame(rbind(label.data.up, label.data.down))
  }
  if (is.null(label.color)) {
    p <- p +
      geom_point(
        data = label.data,
        aes(abundance, log2FoldChange),
        size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(abundance, log2FoldChange),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 1) {
    p <- p +
      geom_point(
        data = label.data,
        aes(abundance, log2FoldChange),
        color = label.color, size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(abundance, log2FoldChange),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 2) {
    p <- p +
      geom_point(
        data = label.data.up,
        aes(abundance, log2FoldChange),
        color = label.color[1], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.up,
        aes(abundance, log2FoldChange),
        label = label.data.up[, label.key]
      ) +
      geom_point(
        data = label.data.down,
        aes(abundance, log2FoldChange),
        color = label.color[2], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.down,
        aes(abundance, log2FoldChange),
        label = label.data.down[, label.key]
      )
  }
  return(p)
}

#' Rank plot for Differential Analysis Results.
#'
#' @param deres Data frame contains all genes/peaks.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param point.alpha Opacity of a geom. Default: 0.6.
#' @param point.size.vec Point size for regular (all points) and(or) labeled points. Default: 2 for regular and 4 for labeled points.
#' @param linetype Threshold linetype. Default: 2.
#' @param point.color.vec Point color for Down, Up regulated genes or peaks. Default: red for Up and blue for Down.
#' @param legend.pos Legend position. Default: top.
#' @param label.num Gene/Peak number to label, choose according to log2FoldChange. When \code{label.df} is set NULL, use this to determine genes to label. Default: NULL.
#' @param label.df Label data frame, at least contains Gene column. Default: NULL(use \code{label.num}). When provided,
#' the second column should not be in \code{deres}.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres} or Gene column of \code{gene.df}).
#' @param label.color Color vector for labels. Default: NULL.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at filter arrange desc top_n vars
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' RankPlot(dds.results.ordered, signif = "pvalue", l2fc.threshold = 0.3, label.num = 2, point.alpha = 0.8, point.color.vec = c("purple", "green"))
RankPlot <- function(deres, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1,
                     point.alpha = 0.6, point.size.vec = c(2, 4), linetype = 2, point.color.vec = c("blue", "red"), legend.pos = "top",
                     label.num = NULL, label.df = NULL, label.key = NULL, label.color = NULL) {
  # preapare DE dataframe
  de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = label.key)
  up.de.df <- de.df %>% dplyr::filter(regulation == "Up_regulated")
  down.de.df <- de.df %>% dplyr::filter(regulation == "Down_regulated")
  deg.df <- as.data.frame(rbind(up.de.df, down.de.df)) %>% dplyr::arrange(dplyr::desc(log2FoldChange))
  deg.df$Order <- 1:nrow(deg.df)
  up.deg.df <- deg.df %>% dplyr::filter(regulation == "Up_regulated")
  down.deg.df <- deg.df %>% dplyr::filter(regulation == "Down_regulated")

  # create plot
  if (length(point.size.vec) == 1) {
    point.size.vec <- rep(point.size.vec, 2)
  }
  p <- ggplot() +
    geom_point(
      data = deg.df,
      aes(Order, log2FoldChange, color = regulation),
      size = point.size.vec[1], alpha = point.alpha
    ) +
    scale_color_manual(values = point.color.vec) +
    theme_classic(base_size = 14) +
    theme(legend.position = legend.pos) +
    labs(x = "Gene", color = "") +
    theme(panel.background = element_rect(color = "black"))
  if (is.null(label.df)) {
    if (!is.null(label.num)) {
      label.data.up <- up.deg.df %>% dplyr::top_n(n = label.num, wt = log2FoldChange)
      label.data.down <- down.deg.df %>% dplyr::top_n(n = label.num, wt = -log2FoldChange)
      if (is.null(label.key)) {
        label.key <- "Label"
        label.data.up[, label.key] <- label.data.up$Gene
        label.data.down[, label.key] <- label.data.down$Gene
      }
      label.data <- as.data.frame(rbind(label.data.up, label.data.down))
    } else {
      return(p)
    }
  } else if (nrow(label.df) >= 1) {
    if (is.null(label.key)) {
      label.key <- "Label"
      label.df[, label.key] <- label.df$Gene
    }
    label.data.up <- merge(up.deg.df, label.df, by = "Gene") %>% as.data.frame()
    label.data.down <- merge(down.deg.df, label.df, by = "Gene") %>% as.data.frame()
    label.data <- as.data.frame(rbind(label.data.up, label.data.down))
  }
  if (is.null(label.color)) {
    p <- p +
      geom_point(
        data = label.data,
        aes(Order, log2FoldChange),
        size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(Order, log2FoldChange),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 1) {
    p <- p +
      geom_point(
        data = label.data,
        aes(Order, log2FoldChange),
        color = label.color, size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data,
        aes(Order, log2FoldChange),
        label = label.data[, label.key]
      )
  } else if (length(label.color) == 2) {
    p <- p +
      geom_point(
        data = label.data.up,
        aes(Order, log2FoldChange),
        color = label.color[1], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.up,
        aes(Order, log2FoldChange),
        label = label.data.up[, label.key]
      ) +
      geom_point(
        data = label.data.down,
        aes(Order, log2FoldChange),
        color = label.color[2], size = point.size.vec[2]
      ) +
      geom_text_repel(
        data = label.data.down,
        aes(Order, log2FoldChange),
        label = label.data.down[, label.key]
      )
  }
  return(p)
}

#' Gene Expresion or Peak Accessibility/Binding Plot.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param deres Data frame contains all genes/peaks.
#' @param group.key Sample group information. When set NULL, select first column of metadata. Default: NULL.
#' @param ref.group Reference group name. When set NULL, select first element of groups. Default: NULL.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Default: 10.
#' @param fill.color Color for box,
#' @param fill.alpha Opacity of a geom. Default: 0.6.
#' @param gene.num Gene/Peak number to plot, choose according to log2FoldChange. When \code{gene.df} is set NULL, use this to determine genes/peak to plot. Default: NULL.
#' @param gene.df Gene data frame, at least contains Gene column. Default: NULL. When set NULL, use \code{gene.num}.
#' @param label.key Column name in \code{gene.df} or \code{deres} to use as gene plot title. Default: NULL. When set NULL, use Gene column.
#' @param plot.col Column number of final plot. Default: 2.
#' @param scales Scales same as \code{\link{facet_wrap}}.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param legend.pos Legend position. Default: top.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select mutate case_when mutate_at filter top_n vars
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @import ggplot2
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' GenePlot(deobj = dds, deres = dds.results.ordered, group.key = "condition", ref.group = "WT", fill.color = c("red", "blue"), fill.alpha = 0.8, gene.num = 2, signif = "pvalue", l2fc.threshold = 0.3)
GenePlot <- function(deobj, deres, group.key = NULL, ref.group = NULL, base = 10, fill.color = c("blue", "red"), fill.alpha = 0.6,
                     gene.df = NULL, label.key = NULL, gene.num = 2, plot.col = 2, scales = "free",
                     signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, legend.pos = "top") {
  # identify analysis method
  if (class(deobj) == "DGEList") {
    counts.matrix <- edgeR::cpm(deobj, normalized.lib.sizes = TRUE)
    metadata <- as.data.frame(deobj$samples) %>% tibble::rownames_to_column(var = "Sample")
    group.key <- "group"
  } else if (class(deobj) == "DESeqDataSet") {
    if (!"sizeFactor" %in% colnames(SummarizedExperiment::colData(deobj))) {
      deobj <- DESeq2::estimateSizeFactors(deobj)
    }
    counts.matrix <- DESeq2::counts(deobj, normalized = TRUE)
    metadata <- as.data.frame(SummarizedExperiment::colData(deobj)) %>% tibble::rownames_to_column(var = "Sample")
  } else {
    stop("Input object is either DESeq2 or edgeR results!")
  }
  # get group key
  if (is.null(group.key)) {
    group.key <- colnames(metadata)[2]
  } else if (!group.key %in% colnames(metadata)) {
    stop(paste0("group.key you provided is not in ", colnames(metadata)))
  }
  # get ref and treat condition
  groups <- unique(as.character(metadata[, group.key]))
  if (is.null(ref.group)) {
    ref.group <- groups[1]
  } else if (!ref.group %in% groups) {
    stop(paste0("ref.group you provided is not in ", groups))
  }
  treat.group <- setdiff(groups, ref.group)
  # get log transformed matrix
  log.counts.matrix <- round(log(counts.matrix + 1, base = base), digits = 4)
  if (is.null(gene.df)) {
    if (is.null(gene.num)) {
      stop("Please specify gene or gene.num!")
    } else {
      # preapare DE dataframe
      de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = label.key)
      up.plot.df <- de.df %>%
        dplyr::filter(regulation == "Up_regulated") %>%
        dplyr::top_n(n = gene.num, wt = log2FoldChange)
      down.plot.df <- de.df %>%
        dplyr::filter(regulation == "Down_regulated") %>%
        dplyr::top_n(n = gene.num, wt = -log2FoldChange)
      deg.plot.df <- as.data.frame(rbind(up.plot.df, down.plot.df))
      if (is.null(label.key)) {
        label.key <- "Label"
        gene.df <- deg.plot.df["Gene"]
        gene.df[, label.key] <- gene.df$Gene
      } else {
        gene.df <- deg.plot.df[c("Gene", label.key)]
      }
    }
  } else if (nrow(gene.df) >= 1) {
    if (is.null(label.key)) {
      label.key <- "Label"
      gene.df[, label.key] <- gene.df$Gene
    }
  }
  # get plot matrix
  gene.matrix <- merge(log.counts.matrix, gene.df, by.y = "Gene", by.x = 0) %>%
    tibble::column_to_rownames(var = "Row.names") %>%
    reshape2::melt(id.vars = label.key, variable.name = "Sample", value.name = "Expression")
  # merge with metadata
  use.meta <- metadata[c("Sample", group.key)]
  colnames(use.meta) <- c("Sample", "Condition")
  plot.df <- merge(gene.matrix, use.meta, by = "Sample") %>%
    as.data.frame() %>%
    purrr::set_names(c("Sample", "Label", "Expression", "Condition"))
  plot.df$Label <- factor(plot.df$Label, levels = as.character(gene.df[, label.key]))
  plot.df$Condition <- factor(plot.df$Condition, levels = c(ref.group, treat.group))
  # create basic plot
  p <- ggplot(plot.df, aes(x = Label, y = Expression, fill = Condition)) +
    geom_boxplot(alpha = fill.alpha) +
    theme_classic(base_size = 14) +
    theme(
      panel.background = element_rect(color = "black"), legend.position = legend.pos,
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = fill.color) +
    labs(
      x = "Gene",
      y = paste0("Log", base, "(Normalized Counts)"), fill = ""
    ) +
    facet_wrap(~Label, scales = scales, ncol = plot.col)
  return(p)
}


#' Heatmap for Differential Analysis Results.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param deres Data frame contains all genes/peaks.
#' @param group.key Sample group information. When set NULL, select first column of metadata. Default: NULL.
#' @param ref.group Reference group name. When set NULL, select first element of groups. Default: NULL.
#' @param group.color Color for different sample group. Default: blue for \code{ref.group}, and red for the other group.
#' @param gene.df Gene data frame, at least contains Gene column. Default: NULL. When provided,
#' the second column should not be in \code{deres}.
#' @param label.key Which column to use as label. Default: NULL. When set NULL, use rownames of \code{deres} or Gene column of \code{gene.df}.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param exp.range Z-score range to plot. Default: c(-2,2).
#' @param exp.color Color map used for heatmap. Default: c("green","black","red").
#' @param heatmap.height The height of whole heatmap. Default: 20cm.
#' @param heatmap.width The width of whole heatmap. Default: 20cm.
#' @param col.gap Gap between column slices. Default: 2mm.
#' @param legend.height The height of legend. Default: 5cm.
#' @param link.height The height of the segments. Default: 4mm.
#'
#' @return A \code{\link{Heatmap-class}}.
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select mutate case_when mutate_at pull filter arrange desc vars
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' DEHeatmap(deobj = dds, deres = dds.results.ordered, group.key = "condition", ref.group = "WT", signif = "pvalue", l2fc.threshold = 0.3)
DEHeatmap <- function(deobj, deres, group.key = NULL, ref.group = NULL, group.color = c("blue", "red"), gene.df = NULL, label.key = NULL,
                      signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1,
                      exp.range = c(-2, 2), exp.color = c("green", "black", "red"), heatmap.height = 20, heatmap.width = 20,
                      col.gap = 2, legend.height = 5, link.height = 4) {
  # identify analysis method
  if (class(deobj) == "DGEList") {
    counts.matrix <- edgeR::cpm(deobj, normalized.lib.sizes = TRUE)
    metadata <- as.data.frame(deobj$samples) %>% tibble::rownames_to_column(var = "Sample")
    group.key <- "group"
  } else if (class(deobj) == "DESeqDataSet") {
    if (!"sizeFactor" %in% colnames(SummarizedExperiment::colData(deobj))) {
      deobj <- DESeq2::estimateSizeFactors(deobj)
    }
    counts.matrix <- DESeq2::counts(deobj, normalized = TRUE)
    metadata <- as.data.frame(SummarizedExperiment::colData(deobj)) %>% tibble::rownames_to_column(var = "Sample")
  } else {
    stop("Input object is either DESeq2 or edgeR results!")
  }
  # get group key
  if (is.null(group.key)) {
    group.key <- colnames(metadata)[2]
  } else if (!group.key %in% colnames(metadata)) {
    stop(paste0("group.key you provided is not in ", colnames(metadata)))
  }
  # get ref condition
  groups <- unique(as.character(metadata[, group.key]))
  if (is.null(ref.group)) {
    ref.group <- groups[1]
  } else if (!ref.group %in% groups) {
    stop(paste0("ref.group you provided is not in ", groups))
  }
  treat.group <- setdiff(groups, ref.group)
  ref.samples <- metadata[metadata[, group.key] == ref.group, ] %>% dplyr::pull(Sample)
  treat.samples <- metadata[metadata[, group.key] == treat.group, ] %>% dplyr::pull(Sample)

  # preapare DE dataframe
  de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = label.key)
  deg.df <- de.df %>% dplyr::filter(regulation != "Not_regulated")
  # get plot genes
  if (is.null(gene.df)) {
    selected.df <- deg.df
    if (is.null(label.key)) {
      label.key <- "Label"
      selected.df[, label.key] <- selected.df$Gene
      selected.df <- selected.df[c("Gene", "log2FoldChange", label.key)]
      mark.genes <- NULL
    }
  } else {
    selected.df <- merge(deg.df, gene.df, all = T) %>% as.data.frame()
    label.selected.df <- merge(deg.df, gene.df, all.y = T) %>% as.data.frame()
    if (is.null(label.key)) {
      label.key <- "Label"
      selected.df <- selected.df[c("Gene", "log2FoldChange")]
      selected.df[, label.key] <- selected.df$Gene
      mark.genes <- gene.df$Gene
    } else {
      selected.df <- selected.df[c("Gene", "log2FoldChange", label.key)]
      selected.df[, label.key] <- ifelse(is.na(selected.df[, label.key]), selected.df[, "Gene"], selected.df[, label.key])
      mark.genes <- label.selected.df[, label.key]
    }
  }
  # get normalized counts and order
  selected.matrix <- merge(selected.df, counts.matrix, by.x = "Gene", by.y = 0) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
    dplyr::select(-c("Gene", "log2FoldChange")) %>%
    tibble::column_to_rownames(var = label.key)
  # log transform and score
  log.selected.matrix <- t(round(log2(selected.matrix + 1), digits = 4))
  scale.selected.matrix <- scale(log.selected.matrix)
  # adjust sample order
  scale.selected.matrix <- scale.selected.matrix[c(ref.samples, treat.samples), ]
  # adjust value range
  if (is.null(exp.range)) {
    exp.range <- c(-2, 2)
  } else if (length(exp.range) == 1) {
    exp.range <- c(-exp.range, exp.range)
  }
  scale.selected.matrix[scale.selected.matrix < exp.range[1]] <- exp.range[1]
  scale.selected.matrix[scale.selected.matrix > exp.range[2]] <- exp.range[2]
  plot.matrix <- t(scale.selected.matrix)
  # prepare annotation
  split <- data.frame(Condition = rep(c(ref.group, treat.group), times = c(length(ref.samples), length(treat.samples))))
  split$Condition <- factor(split$Condition, levels = c(ref.group, treat.group))
  if (is.null(group.color)) {
    group.color <- c("blue", "red")
  } else if (length(group.color) < 2) {
    stop("Length of group.color is less than group number.")
  }
  names(group.color) <- c(ref.group, treat.group)
  col.anno <- ComplexHeatmap::HeatmapAnnotation(df = split, show_legend = F, col = list(Condition = group.color))
  # prepare color
  if (is.null(exp.color)) {
    exp.color <- c("green", "black", "red")
  } else if (length(exp.color) < 3) {
    stop("Length of exp.color should be 3.")
  }
  exp.color.map <- circlize::colorRamp2(
    breaks = c(min(plot.matrix), 0, max(plot.matrix)),
    colors = exp.color
  )
  if (is.null(mark.genes)) {
    ht <- ComplexHeatmap::Heatmap(plot.matrix,
      show_row_names = FALSE,
      show_column_dend = FALSE, show_row_dend = FALSE,
      name = "Z-score", col = exp.color.map,
      row_order = rownames(plot.matrix), column_order = colnames(plot.matrix), row_names_side = "right",
      heatmap_width = unit(heatmap.width, "cm"), heatmap_height = unit(heatmap.height, "cm"),
      column_split = split, column_title = NULL, column_gap = unit(col.gap, "mm"),
      top_annotation = col.anno,
      heatmap_legend_param = list(
        color_bar = "continuous", legend_height = unit(legend.height, "cm"),
        title_position = "leftcenter-rot"
      )
    )
  } else {
    # mark gene position and annotation
    mark.pos <- match(mark.genes, rownames(plot.matrix))
    row.anno <- ComplexHeatmap::rowAnnotation(link = anno_mark(
      at = mark.pos,
      labels = mark.genes,
      link_height = unit(link.height, "mm")
    ))
    # create heatmap
    ht <- ComplexHeatmap::Heatmap(plot.matrix,
      show_row_names = FALSE,
      show_column_dend = FALSE, show_row_dend = FALSE,
      name = "Z-score", col = exp.color.map,
      row_order = rownames(plot.matrix), column_order = colnames(plot.matrix), row_names_side = "right",
      heatmap_width = unit(heatmap.width, "cm"), heatmap_height = unit(heatmap.height, "cm"),
      column_split = split, column_title = NULL, column_gap = unit(col.gap, "mm"),
      right_annotation = row.anno, top_annotation = col.anno,
      heatmap_legend_param = list(
        color_bar = "continuous", legend_height = unit(legend.height, "cm"),
        title_position = "leftcenter-rot"
      )
    )
  }
  return(ht)
}

#' Stat Genomic Regions of Differential Peaks with Pie Plot.
#'
#' @param deres Data frame contains all peaks.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially accessible/binding peaks. Default: 1.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate case_when mutate_at vars filter
#' @importFrom tidyr drop_na separate
#' @importFrom rlang :=
#' @importFrom ggpie ggpie
#' @export
#'
DiffPeakPie <- function(deres, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1) {
  # extract differential analysis results
  de.res <- ExtractDA(
    deres = deres, data.type = "ATAC", peak.anno.key = "All",
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold
  )
  # change annotation information
  anno.key.named <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
  names(anno.key.named) <- c("P", "5U", "3U", "E", "I", "D", "DI")
  de.res$Annotation <- anno.key.named[de.res$Annotation]
  # filter not regulated
  de.degs <- de.res[de.res$regulation != "Not_regulated", ]
  # create plot
  diff.peak.pie <-
    ggpie::ggpie(
      data = de.degs, group_key = "Annotation", count_type = "full",
      label_info = "ratio", label_type = "horizon", label_split = NULL,
      label_size = 4, label_pos = "in", label_threshold = 10
    )
  return(diff.peak.pie)
}
