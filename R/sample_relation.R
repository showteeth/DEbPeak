#' Sample to sample relationship.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param gene Genes used to create sample relationship. Default: NULL.
#' @param remove.sample Sample(s) to remove. Default: NULL.
#' @param transform.method Data transformation methods, chosen from rlog, vst and ntd. Default: rlog.
#' @param anno.key Group to annotate samples. When set NULL, select first column of metadata. Default: NULL.
#' @param ... Parameters for \code{\link{pheatmap}}.
#'
#' @return Plot of sample distances and sample pearson correlation.
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats dist cor
#' @importFrom ComplexHeatmap pheatmap draw
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplotify as.ggplot
#' @importFrom grid grid.grabExpr
#' @importFrom cowplot plot_grid
#' @importFrom DESeq2 rlog vst normTransform
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
#' SampleRelation(deobj = dds, transform.method = "rlog", anno.key = "condition")
SampleRelation <- function(deobj, gene = NULL, remove.sample = NULL, transform.method = c("rlog", "vst", "ntd"),
                           anno.key = NULL, ...) {
  # check parameter
  transform.method <- match.arg(arg = transform.method)

  # identify analysis method
  if (class(deobj) == "DGEList") {
    message("Differential expression analysis with edgeR!")
    # convert to DESeqDataSet
    deobj <- DEFormats::as.DESeqDataSet(deobj)
  } else if (class(deobj) == "DESeqDataSet") {
    message("Differential expression analysis with DESeq2!")
  } else {
    stop("Input object is either DESeq2 or edgeR results!")
  }
  # remove sample
  if (!is.null(remove.sample)) {
    valid.sample <- setdiff(colnames(deobj), remove.sample)
    deobj <- deobj[, valid.sample]
  }
  # get metadata
  metadata <- as.data.frame(SummarizedExperiment::colData(deobj))
  # get transformed data
  transform.data <- DataTrans(deobj, transform.method)
  trans.data <- SummarizedExperiment::assay(transform.data)
  # select genes data
  if (!is.null(gene)) {
    trans.data <- trans.data[gene, ]
  }
  # sample dist and sample correlation
  sample.dist <- stats::dist(t(trans.data))
  sample.cor <- stats::cor(trans.data, method = "pearson")
  # create annotation dataframe
  if (is.null(anno.key)) {
    anno.key <- colnames(metadata)[1]
  } else if (!anno.key %in% colnames(metadata)) {
    stop("anno.key you provided is not in ", colnames(metadata))
  }
  anno.col <- metadata[anno.key]
  anno.color.list <- list()
  anno.color.vec <- c("red", "blue")
  names(anno.color.vec) <- as.character(unique(metadata[, anno.key]))
  anno.color.list[[anno.key]] <- anno.color.vec
  # sample dist heatmap
  plot.list <- list()
  sample.dist.plot <- ComplexHeatmap::pheatmap(as.matrix(sample.dist),
    color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))(100),
    clustering_distance_rows = sample.dist, clustering_distance_cols = sample.dist, name = " ",
    annotation_col = anno.col, main = "Heatmap of the sample distances", annotation_colors = anno.color.list
  )
  plot.list[["dist"]] <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(sample.dist.plot)))
  # sample pearson correlation heatmap
  sample.pcc.plot <- ComplexHeatmap::pheatmap(as.matrix(sample.cor),
    annotation_col = anno.col, name = " ",
    main = "Heatmap of the sample pearson correlation", annotation_colors = anno.color.list
  )
  plot.list[["corr"]] <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(sample.pcc.plot)))
  # arrange plot
  sample.relation.plot <- cowplot::plot_grid(plotlist = plot.list, ncol = 2)
  return(sample.relation.plot)
}
