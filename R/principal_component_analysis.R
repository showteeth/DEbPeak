# data transformation
DataTrans <- function(obj, transform.method) {
  if (transform.method == "rlog") {
    trans.data <- DESeq2::rlog(obj)
  } else if (transform.method == "vst") {
    trans.data <- DESeq2::vst(obj)
  } else if (transform.method == "ntd") {
    trans.data <- DESeq2::normTransform(obj)
  }
  return(trans.data)
}

PreparePCA <- function(deobj, var.genes = NULL, remove.sample = NULL, transform.method = c("rlog", "vst", "ntd")) {
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
  transform.data <- DataTrans(deobj, transform.method)
  # prepare metadata
  metadata <- SummarizedExperiment::colData(deobj)
  # prepare data for PCA
  trans.data <- SummarizedExperiment::assay(transform.data)
  td.rv <- matrixStats::rowVars(trans.data)
  if (is.null(var.genes)) {
    message("Use all genes for PCA!")
    select <- order(td.rv, decreasing = TRUE)
  } else {
    message("Use Top: ", var.genes, " genes for PCA!")
    var.genes <- as.integer(var.genes)
    select <- order(td.rv, decreasing = TRUE)[seq_len(min(var.genes, length(td.rv)))]
  }
  # get final transformed data
  trans.data.used <- t(trans.data[select, ])
  trans_list <- list(td = trans.data.used, metadata = metadata)
  return(trans_list)
}

#' Calculate PCA.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param var.genes Select genes with larger variance for PCA analysis. Default value is NULL, using all genes.
#' @param remove.sample Sample(s) to remove. Default: NULL.
#' @param transform.method Data transformation methods, chosen from rlog, vst and ntd. Default: rlog.
#'
#' @return List suitable for PCAtools.
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom SummarizedExperiment colData assay
#' @importFrom matrixStats rowVars
#' @importFrom DESeq2 rlog vst normTransform
#' @importFrom stats prcomp
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
#' pca_res <- PCA(deobj = dds, transform.method = "rlog")
PCA <- function(deobj, var.genes = NULL, remove.sample = NULL, transform.method = c("rlog", "vst", "ntd")) {
  # check parameters, default：rlog
  transform.method <- match.arg(arg = transform.method)
  # calculate PCA
  prepare.list <- PreparePCA(deobj = deobj, var.genes = var.genes, remove.sample = remove.sample, transform.method = transform.method)
  trans.data.used <- prepare.list[["td"]]
  metadata <- prepare.list[["metadata"]]
  pca <- stats::prcomp(trans.data.used)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
  # preapre results suitable for PCAtools
  # modified from https://github.com/kevinblighe/PCAtools/blob/master/R/pca.R
  pcaobj <- list(
    rotated = data.frame(pca$x),
    loadings = data.frame(pca$rotation),
    variance = percentVar,
    sdev = pca$sdev,
    metadata = metadata,
    xvars = colnames(trans.data.used),
    yvars = rownames(trans.data.used),
    components = colnames(pca$x)
  )
  rownames(pcaobj$rotated) <- pcaobj$yvars
  rownames(pcaobj$loadings) <- pcaobj$xvars
  names(pcaobj$variance) <- pcaobj$components
  return(pcaobj)
}

#' Detect outlier with robust PCA
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param var.genes Select genes with larger variance for PCA analysis. Default: all genes.
#' @param remove.sample Sample(s) to remove. Default: NULL.
#' @param transform.method Data transformation methods, chosen from rlog, vst and ntd. Default: rlog.
#' @param rpca.method Robust PCA method, chosen from \code{\link{PcaGrid}}, \code{\link{PcaHubert}}. Default: PcaGrid.
#' @param k Number of principal components to compute, for \code{\link{PcaGrid}}, \code{\link{PcaHubert}}. Default: 2.
#' @param ... Parameter for \code{\link{PcaGrid}}, \code{\link{PcaHubert}}.
#'
#' @return A ggplot2 object
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom SummarizedExperiment colData assay
#' @importFrom matrixStats rowVars
#' @importFrom DESeq2 rlog vst normTransform
#' @importFrom rrcov PcaGrid PcaHubert
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
#' outlier.res <- OutlierDetection(deobj = dds, var.genes = NULL, transform.method = "rlog")
#' outlier.res$outlier
#' outlier.res$plot
OutlierDetection <- function(deobj, var.genes = NULL, remove.sample = NULL, transform.method = c("rlog", "vst", "ntd"),
                             rpca.method = c("PcaGrid", "PcaHubert"), k = 2, ...) {
  # check parameters, default：rlog
  transform.method <- match.arg(arg = transform.method)
  # check rPCA method used
  rpca.method <- match.arg(arg = rpca.method)
  # prepare PCA
  prepare.list <- PreparePCA(deobj = deobj, var.genes = var.genes, remove.sample = remove.sample, transform.method = transform.method)
  trans.data.used <- prepare.list[["td"]]
  # run PCA
  if (rpca.method == "PcaGrid") {
    pca <- PcaGrid(trans.data.used, k = k, ...)
  } else if (rpca.method == "PcaHubert") {
    pca <- PcaHubert(trans.data.used, k = k, ...)
  }
  # get outiler
  outliers.flag <- pca@flag
  outliers <- names(outliers.flag[which(outliers.flag == FALSE)])
  plot.df <- data.frame(score_distance = pca@sd, orthogonal_distance = pca@od)
  plot.df$SampleName <- rownames(plot.df)
  if (length(outliers) > 0) {
    outliers.df <- data.frame(SampleName = outliers, Type = "outlier")
    message("Detecting ", length(outliers), " outlier(s): ", paste(outliers, collapse = ","))
    plot.df <- as.data.frame(merge(plot.df, outliers.df, by = "SampleName", all.x = TRUE))
  } else {
    plot.df$Type <- NA
  }
  plot.df$Type <- as.character(plot.df$Type)
  plot.df$Type <- ifelse(is.na(plot.df$Type), "normal", plot.df$Type)
  rpca.plot <- ggplot(data = plot.df, aes(x = score_distance, y = orthogonal_distance, color = Type)) +
    geom_point() +
    geom_hline(yintercept = pca@cutoff.od, color = "red") +
    geom_vline(xintercept = pca@cutoff.sd, color = "red") +
    theme_bw() +
    geom_text_repel(aes(label = SampleName), show.legend = FALSE) +
    labs(x = "Score distance", y = "Orthogonal distance", title = paste0("Outlier detection with ", rpca.method)) +
    theme(plot.title = element_text(hjust = 0.5))
  rpca.res <- list(outlier = outliers, plot = rpca.plot)
  return(rpca.res)
}


#' PCA related functions used in quality control.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param var.genes Select genes with larger variance for PCA analysis. Default: all genes.
#' @param remove.sample Sample(s) to remove. Default: NULL.
#' @param transform.method Data transformation methods, chosen from rlog, vst and ntd. Default: rlog.
#' @param batch Batch column to conduct batch correction. Default value is NULL, do not conduct batch correction.
#' @param raw.deobj Object created by DESeq2 or edgeR before filtering with counts. Default: NULL.
#' @param min.count A feature is considered to be detected if the corresponding number of read counts is > \code{min.count}. By default, \code{min.count} = 10.
#' @param colby Group information to color samples. Default: NULL.
#' @param legend.pos Position of legend ('top', 'bottom', 'left', 'right', 'none'). Default: bottom.
#' @param outlier.detection Logical value. If TRUE, conduct outlier detection with robust PCA. Default: TRUE.
#' @param rpca.method Robust PCA method, chosen from \code{\link{PcaGrid}}, \code{\link{PcaHubert}}. Default: PcaGrid.
#' @param k Number of principal components to compute, for \code{\link{PcaGrid}}, \code{\link{PcaHubert}}. Default: 2.
#' @param ... Parameter for \code{\link{PcaGrid}}, \code{\link{PcaHubert}}.
#'
#' @return List contains all step plots, final PCA results, final deobj.
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom SummarizedExperiment colData assay
#' @importFrom matrixStats rowVars
#' @importFrom DESeq2 rlog vst normTransform counts DESeqDataSetFromMatrix
#' @importFrom stats prcomp
#' @importFrom PCAtools biplot
#' @importFrom sva ComBat_seq
#' @importFrom rrcov PcaGrid PcaHubert
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
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
#' pcs_res <- QCPCA(
#'   deobj = dds, var.genes = NULL, remove.sample = NULL, transform.method = "rlog",
#'   outlier.detection = TRUE, rpca.method = "PcaGrid"
#' )
QCPCA <- function(deobj, var.genes = NULL, remove.sample = NULL, transform.method = c("rlog", "vst", "ntd"),
                  batch = NULL, raw.deobj = NULL, min.count = 10, colby = NULL, legend.pos = c("bottom", "top", "left", "right", "none"),
                  outlier.detection = TRUE, rpca.method = c("PcaGrid", "PcaHubert"), k = 2, ...) {
  # check parameters
  transform.method <- match.arg(arg = transform.method)
  legend.pos <- match.arg(arg = legend.pos)
  rpca.method <- match.arg(arg = rpca.method)

  message("Conduct PCA analysis!")
  plot.list <- list()
  raw.pca <- PCA(deobj = deobj, var.genes = var.genes, remove.sample = remove.sample, transform.method = transform.method)
  metadata <- raw.pca$metadata
  raw.pca.plot <- PCAtools::biplot(raw.pca, colby = colby, legendPosition = legend.pos)
  raw.pca.plot <- raw.pca.plot + labs(title = "PCA before processing") +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)))
  plot.list[["raw"]] <- raw.pca.plot
  pca.res <- raw.pca
  if (!is.null(batch)) {
    if (is.null(raw.deobj)) {
      stop("Batch correction needs dds before filtering!")
    }
    if (batch %in% colnames(raw.pca$metadata)) {
      message("Conduct batch correction!")
      # identify analysis method
      if (class(raw.deobj) == "DGEList") {
        # convert to DESeqDataSet
        raw.deobj <- DEFormats::as.DESeqDataSet(raw.deobj)
      } else if (class(raw.deobj) != "DESeqDataSet") {
        stop("Input object is either DESeq2 or edgeR results!")
      }
      # remove sample
      if (!is.null(remove.sample)) {
        valid.sample <- setdiff(colnames(raw.deobj), remove.sample)
        raw.deobj <- raw.deobj[, valid.sample]
      }
      # get raw counts for batch correction
      raw.counts <- DESeq2::counts(raw.deobj, normalized = FALSE)
      design <- raw.deobj@design
      adjusted.count <- sva::ComBat_seq(as.matrix(raw.counts), batch = as.character(metadata[, batch]), group = NULL)
      # create dds with corrected counts
      deobj <- DESeq2::DESeqDataSetFromMatrix(
        countData = adjusted.count,
        colData = metadata,
        design = design
      )
      # filter dds
      keep.genes <- rowSums(DESeq2::counts(deobj, normalized = FALSE)) >= min.count
      deobj <- deobj[keep.genes, ]
      # after batch correction
      batch.pca <- PCA(deobj = deobj, var.genes = var.genes, remove.sample = remove.sample, transform.method = transform.method)
      batch.pca.plot <- PCAtools::biplot(batch.pca, colby = colby, legendPosition = legend.pos)
      batch.pca.plot <- batch.pca.plot + labs(title = "PCA after batch correction") +
        theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)))
      plot.list[["batch"]] <- batch.pca.plot
      pca.res <- batch.pca
    } else {
      stop("batch you provided is not in ", colnames(raw.pca$metadata))
    }
  }
  if (outlier.detection) {
    message("Conduct outlier detection!")
    outlier.res <- OutlierDetection(
      deobj = deobj, var.genes = var.genes, remove.sample = remove.sample,
      rpca.method = rpca.method, transform.method = transform.method, k = k, ...
    )
    samples.to.remove <- outlier.res$outlier
    outlier.plot <- outlier.res$plot
    plot.list[["rPCA"]] <- outlier.plot
    # after outliter removal
    outlier.pca <- PCA(deobj = deobj, var.genes = var.genes, remove.sample = samples.to.remove, transform.method = transform.method)
    outlier.pca.plot <- PCAtools::biplot(outlier.pca, colby = colby, legendPosition = legend.pos)
    outlier.pca.plot <- outlier.pca.plot + labs(title = "PCA after outliter removal") +
      theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 0, 0)))
    plot.list[["outlier"]] <- outlier.pca.plot
    pca.res <- outlier.pca
    # get dds
    valid.sample <- setdiff(colnames(deobj), samples.to.remove)
    deobj <- deobj[, valid.sample]
  }
  # final plot
  final.pca.plot <- cowplot::plot_grid(plotlist = plot.list)
  final.res <- list(plot = final.pca.plot, deobj = deobj, pca = pca.res)
  return(final.res)
}

#' Generated PCA baisc plots, including screen plot, biplot and pairs plot.
#'
#' @param pca PCA results of \code{\link{PCA}}.
#' @param x The principal component to display on the x axis. Default: PC1.
#' @param y The principal component to display on the y axis. Default: PC2.
#' @param explain.threshold The threshold of explained variance. Default: 90.
#' @param loading.num Select loading gene number based on absolute ordered variable loading for each PC in the biplot. Default: 5.
#' @param loading.label.size Size of loading label. Default: 3.
#' @param loading.label.color Color of loading label. Default: red.
#' @param colby Group information to color samples. Default: NULL.
#' @param pair.pc The principal components to be included in the plot. Default: 5.
#' @param pair.label.size Size of p rincipal component label. Default: 14.
#' @param legend.pos Position of legend ('top', 'bottom', 'left', 'right', 'none'). Default: right.
#'
#' @return A ggplot2 object.
#' @importFrom PCAtools screeplot biplot pairsplot
#' @importFrom cowplot plot_grid
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
#' pca_res <- PCA(deobj = dds, transform.method = "rlog")
#' PCABasic(pca_res, colby = "condition")
PCABasic <- function(pca, x = "PC1", y = "PC2", explain.threshold = 90, loading.num = 5, loading.label.size = 3, loading.label.color = "red",
                     colby = NULL, pair.pc = 5, pair.label.size = 14, legend.pos = c("right", "bottom", "top", "left", "none")) {
  # check parameters
  legend.pos <- match.arg(arg = legend.pos)

  plot.list <- list()
  # screen plot
  screen.plot <- PCAtools::screeplot(pca, hline = explain.threshold)
  plot.list[["screen"]] <- screen.plot
  # biplot
  bi.plot <- PCAtools::biplot(pca,
    x = x, y = y, colby = colby, legendPosition = legend.pos, showLoadings = TRUE, ntopLoadings = loading.num,
    sizeLoadingsNames = loading.label.size, colLoadingsNames = loading.label.color
  )
  plot.list[["biplot"]] <- bi.plot
  # pair plot
  pair_components <- paste0("PC", 1:pair.pc)
  pair.plot <- PCAtools::pairsplot(pca, colby = colby, trianglelabSize = pair.label.size, components = pair_components)
  plot.list[["pairs"]] <- pair.plot
  return(plot.list)
}

#' Export selected PC genes.
#'
#' @param pca PCA results of \code{\link{PCA}}.
#' @param data.type Input data type, choose from RNA, ChIP, ATAC. Default: RNA.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param pc Specify PC to export genes. Default: 1:5.
#' @param gene.num Gene number to export for every PC. Default: 10.
#'
#' @return Dataframe contains PC, Gene and Type (Positive or Negative).
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr filter group_by top_n mutate arrange desc select
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
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
#' pca_res <- PCA(deobj = dds, transform.method = "rlog")
#' ExportPCGenes(pca = pca_res, pc = 1:2, gene.num = 200)
ExportPCGenes <- function(pca, data.type = c("RNA", "ChIP", "ATAC"),
                          peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                          pc = 1:5, gene.num = 10) {

  # check parameters
  data.type <- match.arg(arg = data.type)
  peak.anno.key <- match.arg(arg = peak.anno.key)

  # prepare pca for different data type
  anno.key.named <- c("P", "5U", "3U", "E", "I", "D", "DI")
  names(anno.key.named) <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
  if (data.type != "RNA") {
    pca.loading <- pca$loadings %>%
      tibble::rownames_to_column(var = "feature") %>%
      tidyr::separate(feature, sep = "\\|", into = c("peak", "Gene", "Type"))
    if (peak.anno.key == "All") {
      pca.loading <- pca.loading
    } else {
      pca.loading <- pca.loading[pca.loading$Type == anno.key.named[peak.anno.key], ]
    }
    pca.loading <- pca.loading %>%
      dplyr::mutate(name = paste(peak, Gene, Type, sep = "|")) %>%
      dplyr::select(-c(peak, Gene, Type))
    rownames(pca.loading) <- pca.loading$name
    pca.loading$name <- NULL
  } else {
    pca.loading <- as.data.frame(pca$loadings)
  }

  # # get pca loading info
  # pca.loading <- as.data.frame(pca$loadings)
  # get used PC
  pc.use <- paste0("PC", pc)
  pca.loading.used <- pca.loading[pc.use] %>%
    tibble::rownames_to_column(var = "Gene") %>%
    reshape2::melt(id = "Gene") %>%
    purrr::set_names(c("Gene", "PC", "Loadding"))
  pc_positive_top <- pca.loading.used %>%
    dplyr::filter(Loadding > 0) %>%
    dplyr::group_by(PC) %>%
    dplyr::top_n(gene.num) %>%
    dplyr::mutate(Type = "Positive") %>%
    dplyr::arrange(PC, dplyr::desc(Loadding))
  pc_negative_top <- pca.loading.used %>%
    dplyr::filter(Loadding < 0) %>%
    dplyr::group_by(PC) %>%
    dplyr::top_n(-gene.num) %>%
    dplyr::mutate(Type = "Negative") %>%
    dplyr::arrange(PC, Loadding)
  # combine pos and neg values
  pc_top_gene <- rbind(pc_positive_top, pc_negative_top) %>%
    as.data.frame() %>%
    dplyr::arrange(PC) %>%
    dplyr::select(c(PC, Gene, Loadding, Type))
  return(pc_top_gene)
}

# bar plot for loading information
LoadingBar <- function(pca, data.type = c("RNA", "ChIP", "ATAC"),
                       peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                       pc = 1:5, gene.num = 10, ncol = 2) {
  loading.info <- ExportPCGenes(
    pca = pca, data.type = data.type, peak.anno.key = peak.anno.key,
    pc = pc, gene.num = gene.num
  )
  # get used PC
  pc.use <- paste0("PC", pc)
  plot.list <- list()
  for (pcu in pc.use) {
    pcu.loadding <- loading.info[loading.info$PC == pcu, ] %>% dplyr::arrange(dplyr::desc(Loadding))
    pcu.loadding$Gene <- factor(pcu.loadding$Gene, levels = unique(pcu.loadding$Gene))
    plot.list[[pcu]] <- ggplot(pcu.loadding, aes(x = Gene, y = Loadding, fill = Type)) +
      geom_bar(stat = "identity", show.legend = FALSE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(title = paste0("Lodding info of ", pcu)) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  loading.bar.plot <- cowplot::plot_grid(plotlist = plot.list, ncol = ncol)
  return(loading.bar.plot)
}

# heatmap for loadding
LoadingHeat <- function(deobj, pca, data.type = c("RNA", "ChIP", "ATAC"),
                        peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                        pc = 1:5, gene.num = 10, ncol = 2) {
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
  # size factor or not
  metadata <- SummarizedExperiment::colData(deobj)
  if (!"sizeFactor" %in% colnames(metadata)) {
    deobj <- DESeq2::estimateSizeFactors(deobj)
  }
  # get normalized data
  normalized.counts <- DESeq2::counts(deobj, normalized = TRUE)
  norm.matrix <- t(as.matrix(log2(normalized.counts + 1)))

  loading.info <- ExportPCGenes(
    pca = pca, data.type = data.type, peak.anno.key = peak.anno.key,
    pc = pc, gene.num = gene.num
  )
  # get used PC
  pc.use <- paste0("PC", pc)
  plot.list <- list()
  for (pcu in pc.use) {
    pcu.loadding <- loading.info[loading.info$PC == pcu, ] %>% dplyr::arrange(desc(Loadding))
    pcu.norm.matrix <- norm.matrix[, pcu.loadding$Gene]
    pcu.norm.matrix.scale <- scale(pcu.norm.matrix)
    pcu.norm.matrix.scale[pcu.norm.matrix.scale > 2] <- 2
    pcu.norm.matrix.scale[pcu.norm.matrix.scale <- 2] <- -2
    pcu.heat <- ComplexHeatmap::pheatmap(pcu.norm.matrix.scale,
      cluster_cols = FALSE, legend = FALSE, name = " ",
      main = paste0("Loading heatmap of ", pcu)
    )
    plot.list[[pcu]] <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(pcu.heat)))
  }
  loading.heat.plot <- cowplot::plot_grid(plotlist = plot.list, ncol = ncol)
  return(loading.heat.plot)
}

#' PCA loading plot.
#'
#' @param pca PCA results of \code{\link{PCA}}.
#' @param deobj Object created by DESeq2 or edgeR.
#' @param type loading plot type, chosen from bar, heat. Default: bar.
#' @param data.type Input data type, choose from RNA, ChIP, ATAC. Default: RNA.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param pc Specify PC to export genes. Default: 1:5.
#' @param gene.num Gene number to export for every PC. Default: 10.
#' @param ncol Column of final plots. Default: 2.
#'
#' @return Loading plot.
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr filter group_by top_n mutate arrange desc select
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom cowplot plot_grid
#' @importFrom SummarizedExperiment colData assay
#' @importFrom DESeq2 counts estimateSizeFactors
#' @importFrom ComplexHeatmap pheatmap draw
#' @importFrom ggplotify as.ggplot
#' @importFrom grid grid.grabExpr
#' @importFrom cowplot plot_grid
#' @importFrom tidyr separate
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
#' pca_res <- PCA(deobj = dds, transform.method = "rlog")
#' LoadingPlot(pca = pca_res, type = "bar")
#' LoadingPlot(pca = pca_res, deobj = dds, type = "heat")
LoadingPlot <- function(pca, deobj = NULL, type = c("bar", "heat"),
                        data.type = c("RNA", "ChIP", "ATAC"),
                        peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                        pc = 1:5, gene.num = 10, ncol = 2) {
  # check plot type
  type <- match.arg(arg = type)
  data.type <- match.arg(arg = data.type)
  peak.anno.key <- match.arg(arg = peak.anno.key)

  if (type == "heat") {
    if (is.null(deobj)) {
      stop("Please provide deobj when create loading heatmap!")
    } else {
      loading.plot <- LoadingHeat(
        deobj = deobj, pca = pca, data.type = data.type, peak.anno.key = peak.anno.key,
        pc = pc, gene.num = gene.num, ncol = ncol
      )
    }
  } else if (type == "bar") {
    loading.plot <- LoadingBar(
      pca = pca, data.type = data.type, peak.anno.key = peak.anno.key,
      pc = pc, gene.num = gene.num, ncol = ncol
    )
  }
  return(loading.plot)
}

#' GO enrichment on PC loading genes.
#'
#' @param pca PCA results of \code{\link{PCA}}.
#' @param pc Selected PC. Default: 1.
#' @param gene.num Gene number to export for every PC. Default: 200.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param data.type Input data type, choose from RNA, ChIP, ATAC. Default: RNA.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param go.type GO enrichment type, chosen from ALL, BP, MF, CC. Default: ALL.
#' @param enrich.pvalue Cutoff value of pvalue. Default: 0.05.
#' @param enrich.qvalue Cutoff value of qvalue. Default: 0.05.
#' @param org.db Organism database. Default: org.Mm.eg.db.
#' @param padj.method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.
#' @param show.term Number of enrichment term to show. Default: 15.
#' @param str.width Length of enrichment term in plot. Default: 30.
#' @param plot.resolution Resolution of plot. Default: 300.
#' @param plot.width The width of plot. Default: 7.
#' @param plot.height The height of plot. Default: 9.
#' @param save Logical value, whether to save all results. Default: TRUE.
#'
#' @return If \code{save} is TRUE, return NULL (all results are in \code{out.folder}), else retutn list contains all results.
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr filter group_by top_n mutate arrange desc select
#' @importFrom magrittr %>%
#' @import clusterProfiler
#' @importFrom enrichplot dotplot
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom BiocManager install
#' @importFrom tidyr separate
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
#' pca_res <- PCA(deobj = dds, transform.method = "rlog")
#' LoadingGO(pca_res, gene.type = "ENSEMBL", go.type = "BP", padj.method = "BH", save = TRUE)
LoadingGO <- function(pca, pc = 1, gene.num = 200, out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), data.type = c("RNA", "ChIP", "ATAC"),
                      peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                      go.type = c("ALL", "BP", "MF", "CC"), enrich.pvalue = 0.05, enrich.qvalue = 0.05, org.db = "org.Mm.eg.db",
                      padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                      show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  gene.type <- match.arg(arg = gene.type)
  data.type <- match.arg(arg = data.type)
  peak.anno.key <- match.arg(arg = peak.anno.key)
  go.type <- match.arg(arg = go.type)
  padj.method <- match.arg(arg = padj.method)

  # get positive and negative genes
  pc.loading.genes <- ExportPCGenes(
    pca = pca, data.type = data.type, peak.anno.key = peak.anno.key,
    pc = pc, gene.num = gene.num
  )
  # prepare genes for ChIP-seq and ATAC-seq
  if (data.type != "RNA") {
    pc.loading.genes$Gene <- gsub(pattern = ".*\\|(.*)\\|.*", replacement = "\\1", x = pc.loading.genes$Gene)
  }
  pc.loading.genes.positive <- pc.loading.genes[pc.loading.genes$Type == "Positive", "Gene"]
  pc.loading.genes.negative <- pc.loading.genes[pc.loading.genes$Type == "Negative", "Gene"]

  # prepare org db
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org_db : ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (save) {
    # for positive
    SingleFE(
      genes = pc.loading.genes.positive, out.folder = out.folder, regulation = "Positive", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # for negative
    SingleFE(
      genes = pc.loading.genes.negative, out.folder = out.folder, regulation = "Negative", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(NULL)
  } else {
    loading.go.results <- list()
    # for positive
    loading.positive.go <- SingleFE(
      genes = pc.loading.genes.positive, out.folder = out.folder, regulation = "Positive", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # for negative
    loading.negative.go <- SingleFE(
      genes = pc.loading.genes.negative, out.folder = out.folder, regulation = "Negative", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    loading.go.results[["Positive"]] <- loading.positive.go
    loading.go.results[["Negative"]] <- loading.negative.go
    return(loading.go.results)
  }
}

#' 3D PCA plot.
#'
#' @param pca PCA results of \code{\link{PCA}}.
#' @param x The principal component to display on the x axis. Default: PC1.
#' @param y The principal component to display on the y axis. Default: PC2.
#' @param z The principal component to display on the z axis. Default: PC3.
#' @param color.key Group to color the points. When set NULL, select first column of metadata. Default: NULL.
#' @param color Colors used. When set NULL, select color from Set1 of \code{RColorBrewer} package. Default: NULL.
#' @param cex Point size. Default: 1.5.
#' @param col Point edge color. Default: black.
#' @param ticktype Axis ticktype, chosen from: detailed and simple. Default: detailed.
#' @param bty The type of the box, the default draws only the back panels,
#' chosen from: f, b, b2, g, bl, bl2, u, n. Default: g.
#' @param box Should the bounding box for the surface be displayed. Default: TRUE.
#' @param theta Angles defining the viewing direction. \code{theta} gives the azimuthal direction. Default: 140.
#' @param phi Angles defining the viewing direction. \code{phi} gives the colatitude direction. Default: 20.
#' @param d a value which can be used to vary the strength of the perspective transformation.
#' Values of d greater than 1 will lessen the perspective effect and
#' values less and 1 will exaggerate it. Default: 3.
#' @param colkey \code{\link{colkey}}. Default: FALSE.
#' @param label.size Label font size. Default: 0.5.
#' @param legend.pos Legend position. Default: bottom.
#' @param legend.horiz Logical; if TRUE, set the legend horizontally rather than
#' vertically (specifying horiz overrides the ncol specification). Default: TRUE.
#' @param legend.inset Inset distance(s) from the margins as a fraction of the plot region when
#' legend is placed by keyword. Default: -0.1.
#' @param legend.bg The background color for the legend box.
#' Note that this is only used if bty != "n". Default: white.
#' @param legend.bty The type of box to be drawn around the legend,
#' chosen from 0 and n. Default: n.
#' @param legend.ncol The number of columns in which to set the legend items.
#' Default: NULL (set the legend horizontally).
#' @param main Title of the plot. Default: NULL.
#' @param ... Parameters for \code{\link{scatter3D}}.
#'
#' @return 3D PCA plot.
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plot3D scatter3D text3D
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
#' pca_res <- PCA(deobj = dds, transform.method = "rlog")
#' PCA3D(pca = pca_res, color.key = "condition", main = "3D PCA")
PCA3D <- function(pca, x = "PC1", y = "PC2", z = "PC3", color.key = NULL,
                  color = NULL, cex = 1.5, col = "black", ticktype = "detailed", bty = "g",
                  box = TRUE, theta = 140, phi = 20, d = 3, colkey = FALSE, label.size = 0.5,
                  legend.pos = "bottom", legend.horiz = TRUE, legend.inset = -0.1,
                  legend.bg = "white", legend.bty = "n", legend.ncol = NULL, main = NULL, ...) {
  # get plot info and meta info
  pca.rotated <- pca$rotated
  pca.meta <- pca$metadata
  if (is.null(color.key)) {
    color.key <- colnames(pca.meta)[1]
  } else if (!color.key %in% colnames(pca.meta)) {
    stop("color.key you provided is not in ", colnames(pca.meta))
  }
  # get plot colors
  if (is.null(color)) {
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    color.len <- length(unique(as.character(pca.meta[, color.key])))
    color.used <- getPalette(color.len)
    point.colors <- color.used[as.integer(as.factor(pca.meta[, color.key]))]
  } else {
    color.len <- length(unique(as.character(pca.meta[, color.key])))
    if (color.len != length(color)) {
      stop("Color provided does not match ", color.key)
    }
    color.used <- color
    point.colors <- color.used[as.integer(as.factor(pca.meta[, color.key]))]
  }
  # create 3D plot
  plot3D::scatter3D(
    x = pca.rotated[, x], y = pca.rotated[, y], z = pca.rotated[, z],
    pch = 21, cex = cex, col = col, bg = point.colors,
    xlab = x, ylab = y, zlab = z,
    ticktype = ticktype, bty = bty, box = box,
    theta = theta, phi = phi, d = d, colkey = colkey, main = main, ...
  )
  plot3D::text3D(
    x = pca.rotated[, x], y = pca.rotated[, y], z = pca.rotated[, z], labels = rownames(pca.meta),
    add = TRUE, colkey = FALSE, cex = label.size
  )
  if (is.null(legend.ncol)) {
    graphics::legend(
      x = legend.pos, legend = levels(as.factor(pca.meta[, color.key])), pch = 21,
      pt.bg = color.used, bg = legend.bg, bty = legend.bty,
      inset = legend.inset, xpd = TRUE, horiz = legend.horiz
    )
  } else {
    graphics::legend(
      x = legend.pos, legend = levels(as.factor(pca.meta[, color.key])), pch = 21,
      pt.bg = color.used, bg = legend.bg, bty = legend.bty,
      inset = legend.inset, xpd = TRUE, ncol = legend.ncol
    )
  }
}
