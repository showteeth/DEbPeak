# modify from countsbio.plot in NOISeq
CPMPlot <- function(dat, samples = NULL, toplot = "global",
                    cat.colors = NULL, legend.col = NULL, legend.inset = -0.2, ...) {

  ## Preparing data
  if (is.null(samples)) {
    if (NCOL(dat$result) == 1) {
      samples <- 1
    } else {
      samples <- 1:NCOL(dat$result)
    }
  }
  if (is.numeric(toplot)) toplot <- names(dat$summary)[toplot]
  if (is.numeric(samples) && !is.null(colnames(dat$result))) samples <- colnames(dat$result)[samples]
  datos <- dat$summary[[toplot]]
  mytotal <- as.numeric(datos[, "total"])
  datos <- as.matrix(datos[, samples])
  rownames(datos) <- as.character(dat$summary[[toplot]][, 1])

  # get color
  if (is.null(cat.colors)) {
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    miscolores <- getPalette(nrow(datos))
  }

  # creat plot
  if ((exists("ylab") && !is.character(ylab)) || !exists("ylab")) ylab <- ""
  graphics::barplot(as.numeric(datos[1, ]),
    col = miscolores[1], las = 2, main = "", ylab = "", density = 70,
    ylim = c(0, 100), cex.axis = 0.8, names.arg = "", ...
  )
  for (i in 2:(length(mytotal) - 2)) {
    graphics::barplot(as.numeric(datos[i, ]),
      col = miscolores[i], las = 2, main = "", ylab = "", add = TRUE,
      density = 70, ylim = c(0, 100), cex.axis = 0.8, names.arg = "", ...
    )
  }
  bp <- graphics::barplot(as.numeric(datos[(length(mytotal) - 1), ]),
    col = miscolores[(length(mytotal) - 1)], las = 2,
    add = TRUE, names.arg = colnames(datos), cex.axis = 0.8,
    density = 70, ylim = c(0, 100), cex.names = 0.8, ...
  )
  for (j in 1:(length(mytotal) - 1)) graphics::abline(h = mytotal[j], col = miscolores[j], lwd = 2)

  # add text
  if (length(samples) <= 10) {
    graphics::mtext(side = 3, text = datos["depth", ], adj = 0.5, at = bp, cex = 0.8)
  } else {
    graphics::mtext(side = 3, text = datos["depth", ], at = bp, cex = 0.7, las = 2)
  }
  # add legend
  if (is.null(legend.col)) {
    legend.col <- length(mytotal) - 1
  }
  graphics::legend("bottom", rownames(datos)[-length(mytotal)],
    fill = miscolores,
    density = 70, bty = "n", ncol = legend.col, inset = legend.inset, xpd = TRUE
  )
}

#' Count QC plot.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param group.key Sample group information. When set NULL, select first column of metadata. Default: NULL.
#' @param type QC plot type, chosen from saturation and cpm. Default: saturation.
#' @param min.count A feature is considered to be detected if the corresponding number of read counts is > \code{min.count}. By default, \code{min.count} = 0. . This parameter is used by type "saturation".
#' @param ndepth Number of different sequencing depths to be simulated and plotted apart from the real depth. Default: 10. This parameter is only used by type "saturation".
#' @param cat.colors Color used for different CPM Threshold groups or samples. Default: NULL (auto-selection).
#' @param ... Parameters for \code{\link{dat}}.
#'
#' @return Count QC plot
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr arrange select
#' @importFrom scales comma
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom NOISeq readData dat
#' @importFrom DESeq2 counts
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
#' CountQC(deobj = dds, group.key = "condition", type = "cpm")
#' CountQC(deobj = dds, group.key = "condition", type = "saturation")
CountQC <- function(deobj, group.key = NULL, type = c("saturation", "cpm"), min.count = 0, ndepth = 10,
                    cat.colors = NULL, ...) {
  # check parameters, default：saturation
  type <- match.arg(arg = type)

  # identify analysis method
  if (class(deobj) == "DGEList") {
    message("Differential expression analysis with edgeR!")
    counts.matrix <- deobj$counts
    metadata <- as.data.frame(deobj$samples) %>% tibble::rownames_to_column(var = "Sample")
    group.key <- "group"
  } else if (class(deobj) == "DESeqDataSet") {
    message("Differential expression analysis with DESeq2!")
    counts.matrix <- DESeq2::counts(deobj, normalized = FALSE)
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
  group.meta <- metadata[group.key]
  # create ExpressionSet
  data <- NOISeq::readData(data = counts.matrix, factors = group.meta)
  if (type == "saturation") {
    # get saturation info
    saturation.info <- NOISeq::dat(data, k = min.count, ndepth = ndepth, type = "saturation", ...)
    # prepare the data
    sample.saturation.df <- saturation.info@dat$saturation$global %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Rank") %>%
      reshape2::melt(id = "Rank") %>%
      purrr::set_names(c("Rank", "Sample", "Genes"))
    sample.depth.df <- saturation.info@dat$depth %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Rank") %>%
      reshape2::melt(id = "Rank") %>%
      purrr::set_names(c("Rank", "Sample", "Depth"))
    sample.depth.df$Depth <- sample.depth.df$Depth / 1e6
    # merge data
    saturation.plot.data <- merge(sample.saturation.df, sample.depth.df, by = c("Rank", "Sample")) %>%
      as.data.frame()
    saturation.plot.data$Rank <- as.numeric(saturation.plot.data$Rank)
    saturation.plot.data <- saturation.plot.data %>%
      dplyr::arrange(Rank, Sample)
    # real data
    sample.real.df <- saturation.info@dat$real$global %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Sample") %>%
      purrr::set_names(c("Sample", "RealDepth", "RealGenes"))
    sample.real.df$RealDepth <- sample.real.df$RealDepth
    sample.real.df$MinGenes <- min(saturation.plot.data$Genes)
    # prepare color
    sample.num <- length(unique(saturation.plot.data$Sample))
    if (is.null(cat.colors)) {
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      saturation.colors <- getPalette(sample.num)
    } else if (length(cat.colors) < sample.num) {
      warning("The number of colors you provided with cat.colors is smaller than sample number, use auto-selection!")
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      saturation.colors <- getPalette(sample.num)
    } else {
      saturation.colors <- cat.colors
    }
    # create plot
    saturation.plot <- ggplot() +
      geom_line(data = saturation.plot.data, aes_string(x = "Depth", y = "Genes", color = "Sample")) +
      geom_point(data = saturation.plot.data, aes_string(x = "Depth", y = "Genes", color = "Sample"), pch = 21) +
      theme_classic() +
      theme(
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        legend.title = element_blank()
      ) +
      geom_segment(data = sample.real.df, aes_string(x = "RealDepth", y = "MinGenes", xend = "RealDepth", yend = "RealGenes", color = "Sample"), show.legend = F) +
      scale_y_continuous(labels = scales::comma) +
      scale_color_manual(values = saturation.colors) +
      labs(x = "Sequencing depth (million reads)", y = "Number of detected features")
    return(saturation.plot)
  } else if (type == "cpm") {
    countsbio.data <- NOISeq::dat(data, factor = NULL, type = "countsbio", ...)
    # prepare plot.data
    plot.data.total <- countsbio.data@dat$summary$global
    plot.data.bar <- plot.data.total[plot.data.total$global != "depth", colnames(plot.data.total) != "total"] %>%
      reshape2::melt(id = "global")
    colnames(plot.data.bar) <- c("CPM", "Sample", "Percentage")
    plot.data.line <- plot.data.total[plot.data.total$global != "depth", c("global", "total")]
    colnames(plot.data.line) <- c("CPM", "Total")
    plot.data.depth <- plot.data.total[plot.data.total$global == "depth", colnames(plot.data.total) != "total"] %>%
      dplyr::select(-"global") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Sample")
    colnames(plot.data.depth) <- c("Sample", "Depth")
    # # prepare color
    if (is.null(cat.colors)) {
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      cpm.colors <- getPalette(5)
      names(cpm.colors) <- rev(unique(plot.data.bar$CPM))
    } else if (length(cat.colors) < 5) {
      warning("The number of colors you provided with cat.colors is smaller than 5, use auto-selection!")
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      cpm.colors <- getPalette(5)
      names(cpm.colors) <- rev(unique(plot.data.bar$CPM))
    } else {
      cpm.colors <- cat.colors
    }
    # order
    plot.data.bar$CPM <- factor(plot.data.bar$CPM, levels = unique(plot.data.bar$CPM))
    cpm.plot <- ggplot() +
      geom_bar(data = plot.data.bar, aes_string(x = "Sample", y = "Percentage", fill = "CPM"), stat = "identity", position = "identity") +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      geom_hline(data = plot.data.line, aes_string(yintercept = "Total", color = "CPM"), show.legend = F) +
      theme_classic() +
      theme(
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        legend.title = element_blank(),
        axis.title.x = element_blank()
      ) +
      scale_fill_manual(values = cpm.colors) +
      scale_color_manual(values = cpm.colors) +
      geom_text(data = plot.data.depth, aes_string(x = "Sample", y = "98", label = "Depth"), size = 3, vjust = "inward", hjust = "inward") +
      labs(y = "Percentage (%)")
    return(cpm.plot)
  }
}

# CountQC <- function(deobj, group.key = NULL, type = c("saturation", "cpm"), min.count = 0, ndepth = 10,
#                     cat.colors = NULL, legend.col = NULL, legend.inset = -0.2, ...) {
#   # check parameters, default：saturation
#   type <- match.arg(arg = type)
#
#   # identify analysis method
#   if (class(deobj) == "DGEList") {
#     message("Differential expression analysis with edgeR!")
#     counts.matrix <- deobj$counts
#     metadata <- as.data.frame(deobj$samples) %>% tibble::rownames_to_column(var = "Sample")
#     group.key <- "group"
#   } else if (class(deobj) == "DESeqDataSet") {
#     message("Differential expression analysis with DESeq2!")
#     counts.matrix <- DESeq2::counts(deobj, normalized = FALSE)
#     metadata <- as.data.frame(SummarizedExperiment::colData(deobj)) %>% tibble::rownames_to_column(var = "Sample")
#   } else {
#     stop("Input object is either DESeq2 or edgeR results!")
#   }
#   # get group key
#   if (is.null(group.key)) {
#     group.key <- colnames(metadata)[2]
#   } else if (!group.key %in% colnames(metadata)) {
#     stop(paste0("group.key you provided is not in ", colnames(metadata)))
#   }
#   group.meta <- metadata[group.key]
#   # create ExpressionSet
#   data <- NOISeq::readData(data = counts.matrix, factors = group.meta)
#   if (type == "saturation") {
#     # get saturation info
#     saturation.info <- NOISeq::dat(data, k = min.count, ndepth = ndepth, type = "saturation", ...)
#     NOISeq::explo.plot(saturation.info,
#       toplot = 1, samples = 1:ncol(data),
#       yleftlim = NULL, yrightlim = NULL, col.main = "white"
#     )
#     graphics::title("Sequencing Saturation")
#   } else if (type == "cpm") {
#     countsbio.data <- NOISeq::dat(data, factor = NULL, type = "countsbio", ...)
#     CPMPlot(countsbio.data@dat,
#       samples = NULL, toplot = 1,
#       cat.colors = cat.colors, legend.col = legend.col, legend.inset = legend.inset
#     )
#     graphics::title("CPM Threshold", ylab = "Percentage (%)")
#   }
# }
