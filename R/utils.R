# internal ID conversion function
IDConversion_internal <- function(de.df, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                                  org.db = "org.Mm.eg.db", gene.map = NULL, sort.key = "log2FoldChange") {
  # library
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  if (!is.null(gene.map)) {
    # remove duplicated
    gene.map <- gene.map %>% dplyr::distinct(.data[["Gene"]], .keep_all = TRUE)
    final.df <- merge(de.df, gene.map, by = "Gene", all.x = TRUE) %>%
      as.data.frame()
  } else {
    de.df <- de.df %>% dplyr::mutate(Gene2Convert = gsub(pattern = "\\.[0-9]*$", replacement = "", x = as.character(Gene)))
    convert.types <- setdiff(c("ENSEMBL", "ENTREZID", "SYMBOL"), gene.type)
    convert.df <- clusterProfiler::bitr(de.df[, "Gene2Convert"],
      fromType = gene.type,
      toType = convert.types,
      OrgDb = get(org.db), drop = FALSE
    )
    # remove duplicated
    convert.df <- convert.df %>% dplyr::distinct(.data[[gene.type]], .keep_all = TRUE)
    final.df <- merge(de.df, convert.df, by.x = "Gene2Convert", by.y = gene.type, all.x = TRUE) %>%
      as.data.frame() %>%
      dplyr::select(-Gene2Convert)
  }
  if (!is.null(sort.key)) {
    final.df <- final.df %>% dplyr::arrange(dplyr::desc(sort.key))
  }
  return(final.df)
}

#' Gene ID Conversion.
#'
#' @param deres Data frame contains all genes.
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param org.db Organism database. Default: org.Mm.eg.db.
#' @param gene.map Use data frame instead of \code{org.db} to conduct gene id conversion, first column should be Gene. Default: NULL.
#' @param sort.key Column to sort the results. Default: log2FoldChange.
#'
#' @return Data frame contains gene id conversion results.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr distinct mutate select arrange desc
#' @importFrom clusterProfiler bitr
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' IDConversion(count.matrix, gene.type = "ENSEMBL", sort.key = NULL)
IDConversion <- function(deres, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                         org.db = "org.Mm.eg.db", gene.map = NULL, sort.key = "log2FoldChange") {
  # check parameters
  gene.type <- match.arg(arg = gene.type)
  # create dea dataframe
  de.df <- as.data.frame(deres) %>% tibble::rownames_to_column(var = "Gene")
  # ID conversion
  final.df <- IDConversion_internal(
    de.df = de.df, gene.type = gene.type, org.db = org.db,
    gene.map = gene.map, sort.key = sort.key
  )
  final.df <- final.df %>% tibble::column_to_rownames(var = "Gene")
  return(final.df)
}


IDConversionPeak <- function(deres, org.db = "org.Mm.eg.db", gene.map = NULL, sort.key = "log2FoldChange") {
  # prepare dataframe
  de.df <- as.data.frame(deres) %>%
    tibble::rownames_to_column(var = "FullGene") %>%
    dplyr::mutate(Gene = gsub(pattern = ".*\\|(.*)\\|.*", replacement = "\\1", x = FullGene))
  # ID conversion
  final.df <- IDConversion_internal(
    de.df = de.df, gene.type = "SYMBOL", org.db = org.db,
    gene.map = gene.map, sort.key = sort.key
  )
  final.df <- final.df %>% tibble::column_to_rownames(var = "FullGene")
  # change Gene column name to SYMBOL
  colnames(final.df) <- gsub(pattern = "^Gene$", replacement = "SYMBOL", colnames(final.df))
  return(final.df)
}


#' Get Gene Length from GTF.
#'
#' @param gtf.file GTF file.
#' @param save Logical value, whether to save the dataframe.Default: FALSE.
#' @param out.file File to save gene length, used when \code{save} is TRUE. Default: NULL.
#'
#' @return A dataframe or NULL.
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicRanges reduce width
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom purrr set_names
#' @importFrom utils write.table
#' @export
#'
#' @examples
#' # GetGeneLength(gtf.file, save = TRUE)
GetGeneLength <- function(gtf.file, save = FALSE, out.file = NULL) {
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.file, format = "gtf")
  exons.list <- GenomicFeatures::exonsBy(txdb, by = "gene")
  reduce.exons.list <- GenomicRanges::reduce(exons.list)
  gene.length <- sum(GenomicRanges::width(reduce.exons.list))
  gene.length.df <- as.data.frame(gene.length) %>%
    tibble::rownames_to_column(var = "Gene") %>%
    purrr::set_names(c("Gene", "Length"))
  # remove gene version information
  gene.length.df$Gene <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = as.character(gene.length.df$Gene))
  if (save) {
    if (is.null(out.file)) {
      utils::write.table(x = gene.length.df, file = "DEbPeak_geneLength.txt", sep = "\t", quote = FALSE, row.names = F)
    } else {
      utils::write.table(x = gene.length.df, file = out.file, sep = "\t", quote = FALSE, row.names = F)
    }
  } else {
    return(gene.length.df)
  }
}

# convert RPKM to TPM
RPKM2TPM <- function(RPKM) {
  exp(log(RPKM) - log(sum(RPKM)) + log(1e6))
}

#' Perform Counts Normalization.
#'
#' @param deobj Object created by DESeq2 or edgeR.
#' @param gtf.file Gene annotation file used to get gene length, used if \code{norm.type=="RPKM"} or \code{norm.type=="TPM"}. Default: NULL.
#' @param gene.length.file Gene length file, tab-separated, should contains "Gene" and "Length" columns.
#' Used if \code{norm.type=="RPKM"} or \code{norm.type=="TPM"}. Default: NULL.
#' @param norm.type Normalization method, chosen from DESeq2, TMM, CPM, RPKM, TPM. Default: DESeq2.
#'
#' @return Matrix contains normalized counts.
#' @importFrom magrittr %>%
#' @importFrom DEFormats as.DESeqDataSet as.DGEList
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom edgeR calcNormFactors cpm rpkm
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr mutate select
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicRanges reduce width
#' @importFrom purrr set_names
#' @importFrom utils read.table
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
#' NormalizedCount(dds, norm.type = "CPM")
NormalizedCount <- function(deobj, gtf.file = NULL, gene.length.file = NULL,
                            norm.type = c("DESeq2", "TMM", "CPM", "RPKM", "TPM")) {
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
  # use edgeR to calculate TMM, CPM, TPM and RPKM
  if (norm.type == "TMM") {
    if (all(deobj$samples$norm.factors == 1)) {
      deobj <- edgeR::calcNormFactors(deobj)
    }
    normalized.matrix <- edgeR::cpm(deobj, normalized.lib.sizes = TRUE)
  } else if (norm.type == "CPM") {
    normalized.matrix <- edgeR::cpm(deobj$counts)
  } else if (norm.type %in% c("RPKM", "TPM")) {
    if (!is.null(gtf.file)) {
      gene.length <- GetGeneLength(gtf.file)
    } else if (!is.null(gene.length.file)) {
      gene.length <- utils::read.table(file = gene.length.file, header = TRUE, sep = "\t")
      if (!all(c("Gene", "Length") %in% colnames(gene.length))) {
        stop("gene length file should contains 'Gene' and 'Length' columns!")
      }
    } else {
      stop("Calculate RPKM/TPM needs gtf file or gene length file!")
    }
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
