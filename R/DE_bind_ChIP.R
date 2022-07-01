#' Integrate Differential Expression Results and Peak Annotation Results.
#'
#' @param de.res Data frame contains all genes of differential expression analysis.
#' @param chip.res Dataframe contains all peak annotation results.
#' @param chip.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threashold Significance threashold to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threashold Log2 fold change threashold to get differentially expressed genes. Default: 1.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{deres}).
#' @param merge.key The columns used for merging, chosen from geneId, ENSEMBL, SYMBOL. Default: geneId.
#'
#' @return Dataframe contains integrated results. The results will contain five categories: UPbChIP, DOWNbChIP, UP, DOWN, ChIP.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select distinct mutate case_when mutate_at vars
#' @importFrom tibble rownames_to_column
#' @importFrom purrr set_names
#' @importFrom tidyr drop_na
#' @examples
#' library(DEbChIP)
#' library(DESeq2)
#' # ChIP-Seq data
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbChIP")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.profile <- PeakProfile(peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
#' peak.anno <- AnnoPeak(
#'   peak.df = peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # RNA-Seq data
#' count.file <- system.file("extdata", "debchip_count.txt", package = "DEbChIP")
#' meta.file <- system.file("extdata", "debchip_meta.txt", package = "DEbChIP")
#' count.matrix <- read.table(file = count.file, header = T, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = T)
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
#' debchip.res <- DEbChIP(
#'   de.res = dds.results.ordered, chip.res = peak.anno.df,
#'   chip.anno.key = "Promoter", merge.key = "SYMBOL"
#' )
DEbChIP <- function(de.res, chip.res, chip.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                    signif = "padj", signif.threashold = 0.05, l2fc.threashold = 1, label.key = NULL, merge.key = c("geneId", "ENSEMBL", "SYMBOL")) {

  # check parameters
  chip.anno.key <- match.arg(arg = chip.anno.key)
  merge.key <- match.arg(arg = merge.key)

  # parepare ChIP results
  chip.alt.columns <- c("ENSEMBL", "SYMBOL", "GENENAME")
  chip.alt.valid <- intersect(colnames(chip.res), chip.alt.columns)
  chip.df <- chip.res %>%
    dplyr::filter(anno == chip.anno.key) %>%
    dplyr::select(c("geneId", "annotation", "anno", chip.alt.valid)) %>%
    dplyr::distinct(geneId, anno, .keep_all = TRUE)

  # prepare DE results
  de.df <- PrepareDEPlot(
    deres = de.res, signif = signif, signif.threashold = signif.threashold,
    l2fc.threashold = l2fc.threashold, label.key = label.key
  )
  # remove gene version information
  de.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = de.df$Gene)
  deg.df <- de.df %>% dplyr::filter(regulation != "Not_regulated")
  # get all genes used
  if (!merge.key %in% colnames(chip.res)) {
    stop("The merge.key you provided is not valid!")
  }
  all.gene.used <- union(chip.df[, merge.key], deg.df$Gene)
  de.df.used <- de.df[de.df$Gene %in% all.gene.used, ]

  # merge DE and ChIP results
  de.chip <- merge(chip.df, de.df.used, by.x = merge.key, by.y = "Gene", all = TRUE)
  # five categories: UPbChIP, DOWNbChIP, UP, DOWN, ChIP
  de.chip <- de.chip %>%
    dplyr::mutate(Type = dplyr::case_when(
      is.na(annotation) & regulation == "Up_regulated" ~ "UP",
      is.na(annotation) & regulation == "Down_regulated" ~ "DOWN",
      !is.na(annotation) & (regulation == "Not_regulated" | is.na(regulation)) ~ "ChIP",
      !is.na(annotation) & regulation == "Up_regulated" ~ "UPbChIP",
      !is.na(annotation) & regulation == "Down_regulated" ~ "DOWNbChIP"
    ))
  return(de.chip)
}

#' Create Integrated Summary Plot.
#'
#' @param de.chip Dataframe contains integrated results.
#' @param ... Parameters for \code{\link{ggvenn}}.
#'
#' @return A ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom tibble deframe
#' @import ggvenn
#' @export
#'
#' @examples
#' library(DEbChIP)
#' library(DESeq2)
#' # ChIP-Seq data
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbChIP")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.profile <- PeakProfile(peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
#' peak.anno <- AnnoPeak(
#'   peak.df = peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # RNA-Seq data
#' count.file <- system.file("extdata", "debchip_count.txt", package = "DEbChIP")
#' meta.file <- system.file("extdata", "debchip_meta.txt", package = "DEbChIP")
#' count.matrix <- read.table(file = count.file, header = T, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = T)
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
#' debchip.res <- DEbChIP(
#'   de.res = dds.results.ordered, chip.res = peak.anno.df,
#'   chip.anno.key = "Promoter", merge.key = "SYMBOL"
#' )
#' # DE and ChIP venn plot
#' debchip.plot <- PlotDEbChIP(debchip.res, show_percentage = FALSE)
PlotDEbChIP <- function(de.chip, ...) {
  # get summary results
  type.summary <- table(de.chip$Type) %>%
    as.data.frame() %>%
    tibble::deframe()
  # create number vector
  ChIP.vec <- c(
    paste0("ChIP", 1:type.summary["ChIP"]), paste0("UPbChIP", 1:type.summary["UPbChIP"]),
    paste0("DOWNbChIP", 1:type.summary["DOWNbChIP"])
  )
  UP.vec <- c(paste0("UP", 1:type.summary["UP"]), paste0("UPbChIP", 1:type.summary["UPbChIP"]))
  DOWN.vec <- c(paste0("DOWN", 1:type.summary["DOWN"]), paste0("DOWNbChIP", 1:type.summary["DOWNbChIP"]))
  # create plot
  plot <- ggvenn::ggvenn(list(ChIP = ChIP.vec, UP = UP.vec, DOWN = DOWN.vec), ...)
  return(plot)
}
# PlotDEbChIP(de.chip,show_percentage=FALSE)


#' GO Enrichment on Integrated Results.
#'
#' @param de.chip Dataframe contains integrated results.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
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
#' library(DEbChIP)
#' library(DESeq2)
#' # ChIP-Seq data
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbChIP")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.profile <- PeakProfile(peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
#' peak.anno <- AnnoPeak(
#'   peak.df = peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' # RNA-Seq data
#' count.file <- system.file("extdata", "debchip_count.txt", package = "DEbChIP")
#' meta.file <- system.file("extdata", "debchip_meta.txt", package = "DEbChIP")
#' count.matrix <- read.table(file = count.file, header = T, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = T)
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
#' debchip.res <- DEbChIP(
#'   de.res = dds.results.ordered, chip.res = peak.anno.df,
#'   chip.anno.key = "Promoter", merge.key = "SYMBOL"
#' )
#' # functional enrichment on genes
#' fe.results <- DEbChIPFE(de.chip = debchip.res, gene.type = "ENTREZID", species = "Mouse", save = F)
DEbChIPFE <- function(de.chip, out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), go.type = c("ALL", "BP", "MF", "CC"),
                      enrich.pvalue = 0.05, enrich.qvalue = 0.05, species = c(
                        "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                        "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                      ),
                      padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                      show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = T) {
  # check parameter
  gene.type <- match.arg(arg = gene.type)
  go.type <- match.arg(arg = go.type)
  species <- match.arg(arg = species)
  padj.method <- match.arg(arg = padj.method)

  # get UPbChIP and DOWNbChIP genes
  UPbChIP.genes <- de.chip %>%
    dplyr::filter(Type == "UPbChIP") %>%
    dplyr::pull(geneId)
  DOWNbChIP.genes <- de.chip %>%
    dplyr::filter(Type == "DOWNbChIP") %>%
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
  if (save) {
    # for positive
    SingleFE(
      genes = UPbChIP.genes, out.folder = out.folder, regulation = "UPbChIP", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # for negative
    SingleFE(
      genes = DOWNbChIP.genes, out.folder = out.folder, regulation = "DOWNbChIP", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(NULL)
  } else {
    DEbChIP.go.results <- list()
    # for positive
    UPbChIP.go <- SingleFE(
      genes = UPbChIP.genes, out.folder = out.folder, regulation = "UPbChIP", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    # for negative
    DOWNbChIP.go <- SingleFE(
      genes = DOWNbChIP.genes, out.folder = out.folder, regulation = "DOWNbChIP", gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    DEbChIP.go.results[["UPbChIP"]] <- UPbChIP.go
    DEbChIP.go.results[["DOWNbChIP"]] <- DOWNbChIP.go
    return(DEbChIP.go.results)
  }
}
