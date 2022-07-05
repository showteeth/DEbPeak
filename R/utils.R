#' Gene ID Conversion.
#'
#' @param deres Data frame contains all genes.
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param org.db Organism database. Default: org.Mm.eg.db.
#' @param gene.map Use data frame instead of \code{org.db} to conduct gene id conversion, first column should be Gene. Default: NULL.
#' @param sort.key Column to sort the results. Default: log2FoldChange.
#'
#' @return Data frame contains gene id conversion results.
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr distinct mutate select arrange desc
#' @importFrom clusterProfiler bitr
#' @export
#'
#' @examples
#' library(DEbChIP)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbChIP")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' IDConversion(count.matrix, gene.type = "ENSEMBL", sort.key = NULL)
IDConversion <- function(deres, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), org.db = "org.Mm.eg.db", gene.map = NULL, sort.key = "log2FoldChange") {
  # check parameters
  gene.type <- match.arg(arg = gene.type)

  # library
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  de.df <- as.data.frame(deres) %>% tibble::rownames_to_column(var = "Gene")

  if (!is.null(gene.map)) {
    final.df <- merge(de.df, gene.map, by = "Gene", all.x = TRUE) %>%
      as.data.frame() %>%
      dplyr::distinct(Gene, .keep_all = TRUE) %>%
      tibble::column_to_rownames(var = "Gene")
  } else {
    de.df <- de.df %>% dplyr::mutate(Gene2Convert = gsub(pattern = "\\.[0-9]*$", replacement = "", x = as.character(Gene)))
    convert.types <- setdiff(c("ENSEMBL", "ENTREZID", "SYMBOL"), gene.type)
    convert.df <- clusterProfiler::bitr(de.df[, "Gene2Convert"],
      fromType = gene.type,
      toType = convert.types,
      OrgDb = get(org.db), drop = FALSE
    )
    final.df <- merge(de.df, convert.df, by.x = "Gene2Convert", by.y = gene.type, all.x = TRUE) %>%
      as.data.frame() %>%
      dplyr::distinct(Gene2Convert, .keep_all = TRUE) %>%
      dplyr::select(-Gene2Convert) %>%
      tibble::column_to_rownames(var = "Gene")
  }
  if (!is.null(sort.key)) {
    final.df <- final.df %>% dplyr::arrange(dplyr::desc(sort.key))
  }
  return(final.df)
}

# get gene length from gene annotation file
GetLength <- function(gtf.file) {
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.file, format = "gtf")
  exons.list <- GenomicFeatures::exonsBy(txdb, by = "gene")
  reduce.exons.list <- GenomicRanges::reduce(exons.list)
  gene.length <- sum(GenomicRanges::width(reduce.exons.list))
  gene.length.df <- as.data.frame(gene.length) %>%
    tibble::rownames_to_column(var = "Gene") %>%
    purrr::set_names(c("Gene", "Length"))
  # remove gene version information
  gene.length.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = as.character(gene.length.df$Gene))
  return(gene.length.df)
}

# convert RPKM to TPM
RPKM2TPM <- function(RPKM) {
  exp(log(RPKM) - log(sum(RPKM)) + log(1e6))
}

#' Perform Counts Normalization.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param gtf.file Gene annotation file used to get gene length, used if \code{norm.type=="RPKM"} or \code{norm.type=="TPM"}. Default: NULL.
#' @param norm.type Normalization method, chosen from DESeq2, TMM, CPM, RPKM, TPM. Default: DESeq2.
#'
#' @return Matrix contains normalized counts.
#' @importFrom DEFormats as.DESeqDataSet as.DGEList
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom edgeR calcNormFactors cpm rpkm
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr mutate select
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicRanges reduce width
#' @importFrom purrr set_names
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
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' NormalizedCount(dds, norm.type = "CPM")
NormalizedCount <- function(deobj, gtf.file = NULL, norm.type = c("DESeq2", "TMM", "CPM", "RPKM", "TPM")) {
  # check parameters
  norm.type <- match.arg(arg = norm.type)

  # identify analysis method
  if (class(deobj) == "DGEList") {
    if (norm.type == "DESeq2") {
      # convert to DESeqDataSet
      deobj <- DEFormats::as.DESeqDataSet(deobj)
      if (!"sizeFactor" %in% colnames(SummarizedExperiment::colData(deobj))) {
        deobj <- DESeq2::estimateSizeFactors(deobj)
      }
      normalized.matrix <- DESeq2::counts(deobj, normalized = TRUE)
      return(normalized.matrix)
    }
  } else if (class(deobj) == "DESeqDataSet") {
    if (norm.type == "DESeq2") {
      if (!"sizeFactor" %in% colnames(SummarizedExperiment::colData(deobj))) {
        deobj <- DESeq2::estimateSizeFactors(deobj)
      }
      normalized.matrix <- DESeq2::counts(deobj, normalized = TRUE)
      return(normalized.matrix)
    } else {
      # convert to DGEList
      deobj <- DEFormats::as.DGEList(deobj)
    }
  } else {
    stop("Input object is either DESeq2 or edgeR results!")
  }
  # use edgeR to calculate TMM, CPM and RPKM
  if (norm.type == "TMM") {
    if (all(deobj$samples$norm.factors == 1)) {
      deobj <- edgeR::calcNormFactors(deobj)
    }
    normalized.matrix <- edgeR::cpm(deobj, normalized.lib.sizes = TRUE)
  } else if (norm.type == "CPM") {
    normalized.matrix <- edgeR::cpm(deobj$counts)
  } else if (norm.type %in% c("RPKM", "TPM")) {
    if (is.null(gtf.file)) {
      stop("Calculate RPKM/TPM needs gtf file!")
    }
    gene.length <- GetLength(gtf.file)
    raw.counts <- deobj$counts %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Gene") %>%
      dplyr::mutate(Gene = gsub(pattern = "\\.[0-9]*$", replacement = "", x = as.character(Gene)))
    raw.counts.length <- merge(raw.counts, gene.length, by = "Gene")
    raw.counts.length.vec <- as.numeric(raw.counts.length$Length)
    raw.counts.new <- raw.counts.length %>%
      dplyr::select(-Length) %>%
      tibble::column_to_rownames(var = "Gene") %>%
      as.matrix()
    normalized.matrix <- edgeR::rpkm(y = raw.counts.new, gene.length = raw.counts.length.vec)
    if (norm.type == "TPM") {
      normalized.matrix <- RPKM2TPM(normalized.matrix)
    }
  }
  return(normalized.matrix)
}
