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
#' @param cat.colors Color used for different CPM Threshold group. Default: NULL.
#' @param legend.col Column number of legend. Default: NULL(length of CPM Threshold group).
#' @param legend.inset Inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword. Default: -0.2.
#' @param ... Parameters for \code{\link{dat}}.
#'
#' @return Count QC plot
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics barplot abline mtext legend title
#' @importFrom RColorBrewer brewer.pal
#' @importFrom NOISeq readData dat explo.plot
#' @importFrom DESeq2 counts
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbChIP)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbChIP")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbChIP")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' CountQC(deobj = dds, group.key = "condition", type = "cpm")
#' CountQC(deobj = dds, group.key = "condition", type = "saturation")
CountQC <- function(deobj, group.key = NULL, type = c("saturation", "cpm"), min.count = 0, ndepth = 10,
                    cat.colors = NULL, legend.col = NULL, legend.inset = -0.2, ...) {
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
    NOISeq::explo.plot(saturation.info,
      toplot = 1, samples = 1:ncol(data),
      yleftlim = NULL, yrightlim = NULL, col.main = "white"
    )
    graphics::title("Sequencing Saturation")
  } else if (type == "cpm") {
    countsbio.data <- NOISeq::dat(data, factor = NULL, type = "countsbio", ...)
    CPMPlot(countsbio.data@dat,
      samples = NULL, toplot = 1,
      cat.colors = cat.colors, legend.col = legend.col, legend.inset = legend.inset
    )
    graphics::title("CPM Threshold", ylab = "Percentage (%)")
  }
}
