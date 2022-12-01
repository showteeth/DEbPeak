#' Integrate Differential Expression Results and Peak Annotation Results.
#'
#' @param de.res Data frame contains all genes of differential expression analysis.
#' @param peak.res Dataframe contains all peak annotation results.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes. Default: 1.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres}).
#' @param merge.key The columns used for merging, chosen from geneId, ENSEMBL, SYMBOL. Default: geneId.
#'
#' @return Dataframe contains integrated results. The results will contain five categories: UPbPeak, DOWNbPeak, UP, DOWN, Peak.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select distinct mutate case_when mutate_at vars
#' @importFrom tibble rownames_to_column
#' @importFrom purrr set_names
#' @importFrom tidyr drop_na
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
DEbPeak <- function(de.res, peak.res, peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                    signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, label.key = NULL, merge.key = c("geneId", "ENSEMBL", "SYMBOL")) {

  # check parameters
  peak.anno.key <- match.arg(arg = peak.anno.key)
  merge.key <- match.arg(arg = merge.key)

  # parepare peak results
  peak.alt.columns <- c("ENSEMBL", "SYMBOL", "GENENAME")
  peak.alt.valid <- intersect(colnames(peak.res), peak.alt.columns)
  peak.df <- peak.res %>%
    dplyr::filter(anno == peak.anno.key) %>%
    dplyr::select(c("geneId", "annotation", "anno", peak.alt.valid)) %>%
    dplyr::distinct(geneId, anno, .keep_all = TRUE)

  # prepare DE results
  de.df <- PrepareDEPlot(
    deres = de.res, signif = signif, signif.threshold = signif.threshold,
    l2fc.threshold = l2fc.threshold, label.key = label.key
  )
  # remove gene version information
  de.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = de.df$Gene)
  deg.df <- de.df %>% dplyr::filter(regulation != "Not_regulated")
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
PlotDEbPeak <- function(de.peak, peak.type = c("ChIP", "ATAC", "Peak"), ...) {
  # check parameters
  peak.type <- match.arg(arg = peak.type)

  # get summary results
  type.summary <- table(de.peak$Type) %>%
    as.data.frame() %>%
    tibble::deframe()
  if (peak.type == "Peak") {
    #
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
    # create number vector
    Peak.vec <- c(
      PrepareVenn("Peak", type.summary), PrepareVenn("UPbPeak", type.summary),
      PrepareVenn("DOWNbPeak", type.summary)
    )
    UP.vec <- c(PrepareVenn("UP", type.summary), PrepareVenn("UPbPeak", type.summary))
    DOWN.vec <- c(PrepareVenn("DOWN", type.summary), PrepareVenn("DOWNbPeak", type.summary))
    # create plot
    plot.list <- list()
    plot.list[[peak.type]] <- Peak.vec
    plot.list[["UP"]] <- UP.vec
    plot.list[["DOWN"]] <- DOWN.vec
    plot <- ggvenn::ggvenn(plot.list, ...)
  }
  return(plot)
}

#' GO Enrichment on Integrated Results.
#'
#' @param de.peak Dataframe contains integrated results.
#' @param peak.type The source of peaks, chosen from ATAC, ChIP and Peak (ChIP and ATAC). Default: ChIP.
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
#' @return If \code{save} is TRUE, return NULL (all results are in \code{out.folder}), else retutn list contains all results.
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
#' # functional enrichment on genes
#' fe.results <- DEbPeakFE(
#'   de.peak = debchip.res, peak.type = "ChIP", gene.type = "ENTREZID",
#'   species = "Mouse", save = FALSE
#' )
DEbPeakFE <- function(de.peak, peak.type = c("ChIP", "ATAC", "Peak"), out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                      go.type = c("ALL", "BP", "MF", "CC"), enrich.pvalue = 0.05, enrich.qvalue = 0.05, species = c(
                        "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                        "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                      ),
                      padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                      show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  peak.type <- match.arg(arg = peak.type)
  gene.type <- match.arg(arg = gene.type)
  go.type <- match.arg(arg = go.type)
  species <- match.arg(arg = species)
  padj.method <- match.arg(arg = padj.method)

  # get UPbPeak and DOWNbPeak genes
  UPbPeak.genes <- de.peak %>%
    dplyr::filter(Type == "UPbPeak") %>%
    dplyr::pull(geneId)
  DOWNbPeak.genes <- de.peak %>%
    dplyr::filter(Type == "DOWNbPeak") %>%
    dplyr::pull(geneId)

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
  up.reg.str <- paste0("UPb", peak.type)
  down.reg.str <- paste0("DOWNb", peak.type)

  if (save) {
    # for positive
    SingleFE(
      genes = UPbPeak.genes, out.folder = out.folder, regulation = up.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # for negative
    SingleFE(
      genes = DOWNbPeak.genes, out.folder = out.folder, regulation = down.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(NULL)
  } else {
    DEbPeak.go.results <- list()
    # for positive
    UPbPeak.go <- SingleFE(
      genes = UPbPeak.genes, out.folder = out.folder, regulation = up.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # for negative
    DOWNbPeak.go <- SingleFE(
      genes = DOWNbPeak.genes, out.folder = out.folder, regulation = down.reg.str, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    DEbPeak.go.results[[up.reg.str]] <- UPbPeak.go
    DEbPeak.go.results[[down.reg.str]] <- DOWNbPeak.go
    return(DEbPeak.go.results)
  }
}
