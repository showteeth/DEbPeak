#' Conduct Differential Axpression Analysis with DESeq2.
#'
#' @param counts.folder Folder contains all sample's count file. Count file should be SampleName.txt.
#' @param count.matrix.file File contains count matrix, if provided, use this instead of \code{counts.folder}. Default: NULL.
#' @param meta.file File contains sample metadata.
#' @param group.key Column in \code{meta.file} that represents sample group information. Default: NULL.
#' @param count.type The source of count file, chosen from htseq-count, featurecounts. Default: htseq-count.
#' @param min.count A feature is considered to be detected if the corresponding number of read counts is > \code{min.count}. By default, \code{min.count} = 10.
#' @param ref.group Reference group name. When set NULL, select first element of groups. Default: NULL.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param data.type Input data type, choose from RNA, ChIP, ATAC. Default: RNA.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param qc.ndepth Number of different sequencing depths to be simulated and plotted apart from the real depth. Default: 10. This parameter is only used by type "saturation".
#' @param transform.method Data transformation methods, chosen from rlog, vst and ntd. Default: rlog.
#' @param var.genes Select genes with larger variance for PCA analysis. Default: all genes.
#' @param batch Batch column to conduct batch correction. Default value is NULL, do not conduct batch correction.
#' @param outlier.detection Logical value. If TRUE, conduct outlier detection with robust PCA.
#' @param rpca.method robust PCA method, chosen from \code{\link{PcaGrid}}, \code{\link{PcaHubert}}. Default: PcaGrid.
#' @param k number of principal components to compute, for \code{\link{PcaGrid}}, \code{\link{PcaHubert}}. Default: 2.
#' @param pca.x The principal component to display on the x axis. Default: PC1.
#' @param pca.y The principal component to display on the y axis. Default: PC2.
#' @param pca.z The principal component to display on the z axis. Default: PC3.
#' @param loding.pc Specify PC to create loding plot. Default: 1:5.
#' @param loading.gene.num Specify gene number of PC to create loding plot. Default: 10.
#' @param loading.ncol The columns of loading bar or heatmap. Default: 2.
#' @param enrich.loading.pc Specify PC to conduct enrichment analysis. Default: 1:5.
#' @param enrich.loading.gene Specify gene number of PC to conduct enrichment analysis. Default: 200.
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param enrich.type Enrichment type, chosen from ALL, GO, KEGG. Default: ALL.
#' @param go.type GO enrichment type, chosen from ALL, BP, MF, CC. Default: ALL.
#' @param enrich.pvalue Cutoff value of pvalue. Default: 0.05.
#' @param enrich.qvalue Cutoff value of qvalue. Default: 0.05.
#' @param org.db Organism database. Default: org.Mm.eg.db.
#' @param organism Supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'. Default: mmu.
#' @param padj.method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.
#' @param show.term Number of enrichment term to show. Default: 15.
#' @param str.width Length of enrichment term in plot. Default: 30.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes. Default: 1.
#' @param gene.map Use data frame instead of \code{org.db} to conduct gene id conversion, first column should be Gene. Default: NULL.
#' @param gtf.file Gene annotation file used to get gene length, used if \code{norm.type=="RPKM"} or \code{norm.type=="TPM"}. Default: NULL.
#' @param norm.type Normalization method, chosen from DESeq2, TMM, CPM, RPKM, TPM. Default: DESeq2.
#' @param log.counts Logical value, if TRUE, export log2(normalized.counts + 1), else export normalized.counts. Default: TRUE.
#' @param deg.label.df Label data frame, at least contains Gene column. Default: NULL. When set NULL, use \code{deg.label.num}. When provided,
#' the second column should not be SYMBOL, ENSEMBL, ENTREZID.
#' @param deg.label.key Which column to use as label. Default: NULL (use Gene column of \code{deg.label.df}).
#' @param deg.label.num Gene number to label, choose according to log2FoldChange. When \code{deg.label.df} is set NULL, use this to determine genes to label. Default: NULL.
#' @param deg.label.color Color vector for labels. Default: NULL.
#' @param fe.gene.key Column name to conduct enrichment analysis. Default: NULL.
#' @param gmt.file Gene Matrix Transposed file format.
#' @param gene.sets Gene sets information, containing two columns: gs_name, entrez_gene. Default: NULL.
#' @param minGSSize Minimal size of each geneSet for analyzing. Default: 10.
#' @param maxGSSize Maximal size of genes annotated for testing. Default: 500.
#' @param gsea.pvalue Cutoff value of pvalue. Default: 0.05.
#'
#' @return NULL
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames deframe
#' @importFrom SummarizedExperiment colData assay
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom NOISeq readData dat
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom stats dist cor prcomp as.formula
#' @importFrom ggplotify as.ggplot
#' @importFrom grid grid.grabExpr
#' @importFrom cowplot plot_grid
#' @importFrom DESeq2 rlog vst normTransform counts estimateSizeFactors
#' @importFrom matrixStats rowVars
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr filter group_by top_n mutate arrange desc select case_when mutate_at pull
#' @importFrom scales comma
#' @importFrom rrcov PcaGrid PcaHubert
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom PCAtools biplot
#' @importFrom sva ComBat_seq
#' @importFrom rrcov PcaGrid PcaHubert
#' @importFrom plot3D scatter3D text3D
#' @importFrom edgeR cpm
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @import clusterProfiler
#' @importFrom enrichplot dotplot gseaplot2
#' @importFrom stringr str_wrap
#' @importFrom BiocManager install
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' # library(DESeq2)
#' # library(DEbPeak)
#' # count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' # meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' # gmt.file <- system.file("extdata", "m5.go.bp.v2022.1.Mm.entrez.gmt", package = "DEbPeak")
#' # ConductDESeq2(count.matrix.file = count.file, meta.file = meta.file, group.key = "condition",
#' #               count.type = "htseq-count", ref.group = "WT", signif = "pvalue", l2fc.threshold = 0.3, gmt.file = gmt.file)
ConductDESeq2 <- function(counts.folder, count.matrix.file = NULL, meta.file, group.key = NULL,
                          count.type = c("htseq-count", "featurecounts"), min.count = 10, ref.group = NULL,
                          out.folder = NULL, data.type = c("RNA", "ChIP", "ATAC"),
                          peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                          qc.ndepth = 10, transform.method = c("rlog", "vst", "ntd"),
                          var.genes = NULL, batch = NULL, outlier.detection = T, rpca.method = c("PcaGrid", "PcaHubert"), k = 2,
                          pca.x = "PC1", pca.y = "PC2", pca.z = "PC3", loding.pc = 1:5, loading.gene.num = 10, loading.ncol = 2, enrich.loading.pc = 1:5, enrich.loading.gene = 200,
                          gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), enrich.type = c("ALL", "GO", "KEGG"), go.type = c("ALL", "BP", "MF", "CC"),
                          enrich.pvalue = 0.05, enrich.qvalue = 0.05, org.db = "org.Mm.eg.db", organism = "mmu", padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                          show.term = 15, str.width = 30, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, gene.map = NULL, gtf.file = NULL, norm.type = c("DESeq2", "TMM", "CPM", "RPKM", "TPM"), log.counts = TRUE,
                          deg.label.df = NULL, deg.label.key = NULL, deg.label.num = 2, deg.label.color = NULL,
                          fe.gene.key = NULL, gmt.file, gene.sets = NULL, minGSSize = 10, maxGSSize = 500, gsea.pvalue = 0.05) {
  # check parameters
  count.type <- match.arg(arg = count.type)
  data.type <- match.arg(arg = data.type)
  peak.anno.key <- match.arg(arg = peak.anno.key)
  transform.method <- match.arg(arg = transform.method)
  rpca.method <- match.arg(arg = rpca.method)
  gene.type <- match.arg(arg = gene.type)
  enrich.type <- match.arg(arg = enrich.type)
  go.type <- match.arg(arg = go.type)
  padj.method <- match.arg(arg = padj.method)
  norm.type <- match.arg(arg = norm.type)

  # create DESeqDataSet section
  meta.info <- read.table(file = meta.file, header = TRUE)
  if (is.null(group.key)) {
    group.key <- colnames(meta.info)[1]
  } else if (!group.key %in% colnames(meta.info)) {
    stop(paste0("group.key you provided is not in ", colnames(meta.info)))
  }
  if (!is.null(count.matrix.file)) {
    message("Create DESeqDataSet from count matrix!")
    count.matrix <- read.table(file = count.matrix.file, header = TRUE, sep = "\t")
    # create raw dds
    raw.dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count.matrix,
      colData = meta.info,
      design = stats::as.formula(paste0("~", group.key))
    )
  } else {
    count.files <- list.files(path = counts.folder, full.names = TRUE)
    sample.names <- gsub(pattern = "\\.txt", replacement = "", basename(count.files))
    if (count.type == "htseq-count") {
      message("Create DESeqDataSet from count obtained from htseq-count!")
      # create matrix
      count.list <- lapply(as.character(count.files), function(fn) read.table(fn, fill = TRUE, stringsAsFactors = FALSE))
      if (!all(sapply(count.list, function(a) all(a$V1 == count.list[[1]]$V1)))) {
        stop("Gene IDs (first column) differ between files.")
      }
      count.list.pure <- sapply(count.list, function(a) a[, ncol(a)])
      colnames(count.list.pure) <- sample.names
      rownames(count.list.pure) <- count.list[[1]]$V1
      # remove special terms
      special.names <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
      special.rows <- (substr(rownames(count.list.pure), 1, 1) == "_") | rownames(count.list.pure) %in% special.names
      count.list.pure <- count.list.pure[!special.rows, , drop = FALSE]
    } else if (count.type == "featurecounts") {
      message("Create DESeqDataSet from count obtained from featureCounts!")
      # create matrix
      count.list <- lapply(as.character(count.files), function(fn) {
        read.table(fn, fill = TRUE, comment.char = "#", header = TRUE, stringsAsFactors = FALSE)
      })
      if (!all(sapply(count.list, function(a) all(a[, 1] == count.list[[1]][, 1])))) {
        stop("Gene IDs (first column) differ between files.")
      }
      count.list.pure <- sapply(count.list, function(a) a[, ncol(a)])
      colnames(count.list.pure) <- sample.names
      rownames(count.list.pure) <- count.list[[1]][, 1]
    }
    raw.dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count.list.pure,
      colData = meta.info,
      design = stats::as.formula(paste0("~", group.key))
    )
  }
  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # conduct Count quality Control
  # on raw dds
  message("Conduct Quality Control on Counts!")
  grDevices::pdf(file.path(out.folder, "CountQC_CPM.pdf"), width = 10, height = 6)
  CountQC(deobj = raw.dds, group.key = group.key, type = "cpm")
  grDevices::dev.off()
  grDevices::pdf(file.path(out.folder, "CountQC_saturation.pdf"), width = 10, height = 6)
  CountQC(deobj = raw.dds, group.key = group.key, type = "saturation", ndepth = qc.ndepth, min.count = min.count)
  grDevices::dev.off()

  # filter with counts and save results to dds
  message("Raw gene number: ", nrow(raw.dds))
  keep.genes <- rowSums(DESeq2::counts(raw.dds, normalized = FALSE)) >= min.count
  dds <- raw.dds[keep.genes, ]
  message("Gene number after filtered with counts bigger than ", min.count, " : ", nrow(dds))

  # get sample relation
  message("Conduct Quality Control on Samples!")
  grDevices::pdf(file.path(out.folder, "SampleQC_dist_pcc.pdf"), width = 20, height = 10)
  SampleRelation(deobj = dds, transform.method = transform.method, anno.key = group.key)
  grDevices::dev.off()

  # conduct PCA analysis
  message("Conduct PCA analysis!")
  pca.results <- QCPCA(
    deobj = dds, var.genes = var.genes, colby = group.key, transform.method = transform.method, batch = batch, raw.deobj = raw.dds, min.count = min.count,
    outlier.detection = outlier.detection, rpca.method = rpca.method, k = k
  )
  dds <- pca.results$deobj
  pca.results.folder <- file.path(out.folder, "PCA")
  dir.create(pca.results.folder, showWarnings = FALSE, recursive = TRUE)
  setwd(pca.results.folder)

  grDevices::pdf("PCA_overview.pdf", width = 20, height = 15)
  pca.results$plot
  grDevices::dev.off()
  pca.info <- pca.results$pca
  grDevices::pdf("PCA_screen_plot.pdf", width = 10, height = 6)
  PCAtools::screeplot(pca.info, axisLabSize = 18, titleLabSize = 22)
  grDevices::dev.off()
  grDevices::pdf("PCA_biplot.pdf", width = 8, height = 8)
  PCAtools::biplot(pca.info, x = pca.x, y = pca.y, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
  grDevices::dev.off()
  grDevices::pdf("PCA_pairs_plot.pdf", width = 8, height = 8)
  PCAtools::pairsplot(pca.info, colby = group.key)
  grDevices::dev.off()
  grDevices::pdf("PCA_loading_bar.pdf", width = 10, height = 6)
  LoadingPlot(pca = pca.info, type = "bar", pc = loding.pc, gene.num = loading.gene.num)
  grDevices::dev.off()
  grDevices::pdf("PCA_loading_heat.pdf", width = 10, height = 6)
  LoadingPlot(pca = pca.info, deobj = dds, type = "heat", pc = loding.pc, gene.num = loading.gene.num)
  grDevices::dev.off()
  grDevices::pdf("PCA_3DPCA.pdf", width = 8, height = 8)
  PCA3D(pca = pca.info, x = pca.x, y = pca.y, z = pca.z, color.key = group.key, main = "3D PCA")
  grDevices::dev.off()

  # prepare org db
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org_db : ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # PC loading gene enrichment analysis
  message("Conduct Functional Enrichment on PC Loading genes!")
  pc.loading.gene.df <- ExportPCGenes(pca.info,
    data.type = data.type, peak.anno.key = peak.anno.key,
    pc = enrich.loading.pc, gene.num = enrich.loading.gene
  )
  # prepare genes for ChIP-seq and ATAC-seq
  if (data.type != "RNA") {
    pc.loading.gene.df$Gene <- gsub(pattern = ".*\\|(.*)\\|.*", replacement = "\\1", x = pc.loading.gene.df$Gene)
  }
  for (pc in paste0("PC", enrich.loading.pc)) {
    pc.gene.df <- pc.loading.gene.df[pc.loading.gene.df["PC"] == pc, ]
    # on PC positive genes
    pc.positive.genes <- pc.gene.df[pc.gene.df["Type"] == "Positive", "Gene"]
    pc.positive.folder <- file.path(pca.results.folder, pc, "Positive")
    dir.create(pc.positive.folder, showWarnings = FALSE, recursive = TRUE)
    SingleFE(
      genes = pc.positive.genes, out.folder = pc.positive.folder, regulation = paste0(pc, "_Positive"),
      gene.type = gene.type, enrich.type = enrich.type, go.type = go.type, enrich.pvalue = enrich.pvalue,
      enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method,
      show.term = show.term, str.width = str.width
    )
    # on PC negative genes
    pc.negative.genes <- pc.gene.df[pc.gene.df["Type"] == "Negative", "Gene"]
    pc.negative.folder <- file.path(pca.results.folder, pc, "Negative")
    dir.create(pc.negative.folder, showWarnings = FALSE, recursive = TRUE)
    SingleFE(
      genes = pc.negative.genes, out.folder = pc.negative.folder, regulation = paste0(pc, "_Negative"),
      gene.type = gene.type, enrich.type = enrich.type, go.type = go.type, enrich.pvalue = enrich.pvalue,
      enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method,
      show.term = show.term, str.width = str.width
    )
  }
  setwd(out.folder)

  message("Conduct Differential Expression Analysis!")
  # get ref condition
  groups <- unique(as.character(meta.info[, group.key]))
  if (is.null(ref.group)) {
    ref.group <- groups[1]
  } else if (!ref.group %in% groups) {
    stop(paste0("ref.group you provided is not in ", groups))
  }
  treat.group <- setdiff(groups, ref.group)
  # conduct differential expression analysis
  dds <- DESeq2::DESeq(dds)
  # extract results
  dds.results <- DESeq2::results(dds, contrast = c(group.key, treat.group, ref.group))
  dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]

  message("Conduct Gene ID Conversion!")
  if (data.type != "RNA") {
    dds.results.ordered <- IDConversionPeak(deres = dds.results.ordered, org.db = org.db, gene.map = gene.map)
    dds.results.selected <- as.data.frame(dds.results.ordered) %>% tibble::rownames_to_column(var = "Peak")
  } else {
    dds.results.ordered <- IDConversion(deres = dds.results.ordered, gene.type = gene.type, org.db = org.db, gene.map = gene.map)
    dds.results.selected <- as.data.frame(dds.results.ordered) %>% tibble::rownames_to_column(var = "Gene")
  }
  # remove na value in selected columns
  dds.results.selected <- dds.results.selected[!is.na(dds.results.selected[signif]), ]
  dds.results.selected <- dds.results.selected[!is.na(dds.results.selected["log2FoldChange"]), ]
  dds.results.sig <- dds.results.selected[dds.results.selected[signif] <= signif.threshold & abs(dds.results.selected["log2FoldChange"]) >= l2fc.threshold, ]
  # create output folder
  deg.results.folder <- file.path(out.folder, "DEG")
  dir.create(deg.results.folder, showWarnings = FALSE, recursive = TRUE)
  setwd(deg.results.folder)
  # save degs
  out.deg.str <- paste0(signif, signif.threshold, "_FC", l2fc.threshold, ".csv")
  utils::write.csv(as.data.frame(dds.results.sig),
    row.names = FALSE,
    file = paste("Condition", treat.group, ref.group, out.deg.str, sep = "_")
  )
  utils::write.csv(as.data.frame(dds.results.ordered),
    file = paste("Condition", treat.group, ref.group, "all.csv", sep = "_")
  )
  # save normalized counts
  normalized.counts <- NormalizedCount(deobj = dds, gtf.file = gtf.file, norm.type = norm.type)
  if (log.counts) {
    normalized.counts <- as.matrix(log2(normalized.counts + 1))
  }
  utils::write.csv(as.data.frame(normalized.counts),
    file = "normalized_counts.csv"
  )

  message("Visualize Differential Expression Analysis!")
  grDevices::pdf("DEG_VolcanoPlot.pdf", width = 8, height = 10)
  VolcanoPlot(dds.results.ordered,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold,
    label.num = deg.label.num, label.df = deg.label.df, label.key = deg.label.key, label.color = deg.label.color
  )
  grDevices::dev.off()
  grDevices::pdf("DEG_ScatterPlot.pdf", width = 8, height = 9)
  ScatterPlot(
    deobj = dds, deres = dds.results.ordered, group.key = group.key, ref.group = ref.group,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold,
    label.num = deg.label.num, label.df = deg.label.df, label.key = deg.label.key, label.color = deg.label.color
  )
  grDevices::dev.off()
  grDevices::pdf("DEG_MAPlot.pdf", width = 10, height = 8)
  MAPlot(dds.results.ordered,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold,
    label.num = deg.label.num, label.df = deg.label.df, label.key = deg.label.key, label.color = deg.label.color
  )
  grDevices::dev.off()
  grDevices::pdf("DEG_RankPlot.pdf", width = 8, height = 9)
  RankPlot(dds.results.ordered,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold,
    label.num = deg.label.num, label.df = deg.label.df, label.key = deg.label.key, label.color = deg.label.color
  )
  grDevices::dev.off()
  grDevices::pdf("DEG_GenePlot.pdf", width = 8, height = 9)
  GenePlot(
    deobj = dds, deres = dds.results.ordered, group.key = group.key, ref.group = ref.group,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold,
    gene.num = deg.label.num, gene.df = deg.label.df, label.key = deg.label.key
  )
  grDevices::dev.off()
  grDevices::pdf("DEG_Heatmap.pdf", width = 10, height = 8)
  DEHeatmap(
    deobj = dds, deres = dds.results.ordered, group.key = group.key, ref.group = ref.group,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, gene.df = deg.label.df, label.key = deg.label.key
  )
  grDevices::dev.off()
  setwd(out.folder)

  message("Conduct Functional Enrichment Analysis!")
  # create output folder
  fe.results.folder <- file.path(out.folder, "FE")
  dir.create(fe.results.folder, showWarnings = FALSE, recursive = TRUE)
  setwd(fe.results.folder)
  ConductFE(
    deres = dds.results.ordered, out.folder = fe.results.folder, data.type = data.type, peak.anno.key = peak.anno.key,
    signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold,
    gene.key = fe.gene.key, gene.type = gene.type, enrich.type = enrich.type, go.type = go.type, enrich.pvalue = enrich.pvalue,
    enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method, show.term = show.term, str.width = str.width
  )
  setwd(out.folder)
  if (data.type != "RNA") {
    message("Gene Set Enrichment Analysis is not available for ChIP-seq or ATAC-seq data!")
  } else {
    message("Conduct Gene Set Enrichment Analysis!")
    # create output folder
    gsea.results.folder <- file.path(out.folder, "GSEA")
    dir.create(gsea.results.folder, showWarnings = FALSE, recursive = TRUE)
    setwd(gsea.results.folder)
    ConductGSEA(
      deres = dds.results.ordered, gmt.file = gmt.file, gene.sets = gene.sets, out.folder = gsea.results.folder, gene.key = fe.gene.key, gene.type = gene.type,
      org.db = org.db, minGSSize = minGSSize, maxGSSize = maxGSSize, pvalue = gsea.pvalue, padj.method = padj.method
    )
  }
  setwd(out.folder)
  message("All Analysis Done!")
}
