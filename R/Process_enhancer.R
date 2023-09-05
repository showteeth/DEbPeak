#' Get Genes Near Differential Peaks.
#'
#' @param de.res Dataframe contains all differential analysis results of peak-related data.
#' @param signif Significance criterion for RNA-seq results. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold for RNA-seq to get differentially expressed genes. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold for RNA-seq to get differentially expressed genes. Default: 1.
#' @param label.key Which column to use as label. Default: NULL (use rownames of \code{de.res}).
#' @param species Species used, chosen from "Human","Mouse","Rat","Fly","Arabidopsis",
#' "Yeast","Zebrafish","Worm","Bovine","Pig","Chicken","Rhesus",
#' "Canine","Xenopus","Anopheles","Chimp","E coli strain Sakai","Myxococcus xanthus DK 1622". Default: "Human".
#' @param seq.style The style of sequence, chosen from UCSC, NCBI, Ensembl, None. This should be compatible with the genome and gtf file you used
#' to generate count matrix and peak files. Default: "UCSC".
#' @param gtf.file GTF file used to create TxDb object. Useful when specie you used is not available in \code{species}. Default: NULL.
#' @param dis.threshold Distance threshold. Default: 250000.
#' @param n.cores The number of cores to be used for this job. Default:1.
#'
#' @return Dataframe contains differential peaks and their near genes.
#' @export
#' @importFrom GenomicFeatures makeTxDbFromGFF genes
#' @importFrom magrittr %>%
#' @importFrom BiocManager install
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select mutate case_when mutate_at vars filter
#' @importFrom purrr set_names
#' @importFrom tidyr drop_na separate
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle keepSeqlevels seqnames
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges distance
#' @importFrom BiocParallel register MulticoreParam bplapply
#'
#' @examples
#' # # multi-core (faster)
#' # enchancer.genes = ProcessEnhancer(de.res = all.enhancer.res, signif = "pvalue", signif.threshold = 0.05,
#' #                                   l2fc.threshold = 0, label.key = NULL, species = "Mouse", n.cores = 5)
#' # # single-core (slower)
#' # enchancer.genes2 = ProcessEnhancer(de.res = all.enhancer.res, signif = "pvalue", signif.threshold = 0.05,
#' #                                    l2fc.threshold = 0, label.key = NULL, species = "Mouse", n.cores = 1)
ProcessEnhancer <- function(de.res, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1, label.key = NULL,
                            species = c(
                              "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                              "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                            ), seq.style = "Ensembl", gtf.file = NULL,
                            dis.threshold = 250000, n.cores = 1) {
  # check parameters
  species <- match.arg(arg = species)

  # get orgdb and txdb
  spe.anno <- GetSpeciesAnno(species)
  # prepare txdb
  if (!is.null(gtf.file)) {
    message("Create txdb from gtf file!")
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.file)
    # annotation
    txdb.obj <- txdb
  } else {
    txdb <- spe.anno[["txdb"]]
    # library txdb
    if (!require(txdb, quietly = TRUE, character.only = TRUE)) {
      message("Install txdb: ", txdb)
      BiocManager::install(txdb)
    }
    suppressWarnings(suppressMessages(library(txdb, character.only = TRUE)))
    # annotation
    txdb.obj <- get(txdb)
  }

  # get diff enhancer
  de.df <- PrepareDEPlot(
    deres = de.res, signif = signif, signif.threshold = signif.threshold,
    l2fc.threshold = l2fc.threshold, label.key = label.key
  )
  # get diff enhancer
  diff.res <- de.df %>%
    dplyr::filter(regulation != "Not_regulated") %>%
    dplyr::mutate(Range = gsub(pattern = "([^\\|]*)\\|.*", replacement = "\\1", x = Gene))
  # extract range
  diff.res.range <- diff.res["Range"] %>%
    dplyr::mutate(SRange = Range) %>%
    tidyr::separate(col = SRange, into = c("chr", "region"), sep = ":") %>%
    tidyr::separate(col = region, into = c("start", "end"), sep = "-")

  # change seq level
  if (seq.style != "None") {
    GenomeInfoDb::seqlevelsStyle(txdb.obj) <- seq.style
    diff.res.seqs <- unique(as.character(diff.res.range$chr))
    # filter seqs
    txdb.seqs <- GenomeInfoDb::seqlevels(txdb.obj)
    common.seqs <- intersect(diff.res.seqs, txdb.seqs)
    diff.res.range <- diff.res.range[diff.res.range$chr %in% common.seqs, ]
    # to avoid Unable to find an inherited method for function 'NSBS' for signature '"SortedByQueryHits"'
    txdb.obj <- GenomeInfoDb::keepSeqlevels(x = txdb.obj, value = common.seqs)
  }

  # get gene region
  gene.gr <- GenomicFeatures::genes(txdb.obj)
  gene.df <- as.data.frame(gene.gr)

  # get near genes
  if (is.null(n.cores) || n.cores == 1) {
    diff.res.range.list <- lapply(1:nrow(diff.res.range), function(x) {
      Region2Gene(x.df = diff.res.range[x, ], gene.gr = gene.gr, dis.threshold = dis.threshold)
    })
  } else {
    # register
    BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
    diff.res.range.sl <- split(diff.res.range, diff.res.range$Range)
    diff.res.range.list <- BiocParallel::bplapply(diff.res.range.sl, BPPARAM = BiocParallel::MulticoreParam(), FUN = function(x) {
      Region2Gene(x.df = x, gene.gr = gene.gr, dis.threshold = dis.threshold)
    })
  }
  diff.res.range.df <- do.call(rbind, diff.res.range.list)
  diff.res.range.df <- diff.res.range.df %>%
    dplyr::mutate(Range = paste(chr, paste(start, end, sep = "-"), sep = ":")) %>%
    dplyr::select(-c(chr, start, end))

  # get final results
  diff.res.gene <- merge(diff.res, diff.res.range.df, by = "Range", all.x = TRUE) %>%
    as.data.frame() %>%
    dplyr::select(-Range) %>%
    tibble::column_to_rownames(var = "Gene")

  return(diff.res.gene)
}

# get near gene of peak
Region2Gene <- function(x.df, gene.gr, dis.threshold) {
  # convert df to granges
  x.gr <- GenomicRanges::makeGRangesFromDataFrame(x.df)
  # distance
  dis <- IRanges::distance(x.gr, gene.gr)
  valid.dis <- which(dis < dis.threshold)
  # convert gene granges to df
  gene.df <- as.data.frame(gene.gr)
  if (length(valid.dis) > 0) {
    x.dis <- gene.df[valid.dis, ]
    x.dis$dist <- dis[valid.dis]
    x.dis$gene_id <- paste(x.dis$gene_id, x.dis$dist, sep = "|")
    x.df$Near_Gene <- paste0(x.dis$gene_id, collapse = ", ")
  } else {
    x.df$Near_Gene <- ""
  }
  return(x.df)
}

# get gene name type
CheckGeneName <- function(genes, orgdb) {
  gene.types <- c("ENTREZID", "SYMBOL", "ENSEMBL")
  # library
  suppressWarnings(suppressMessages(library(orgdb, character.only = TRUE)))
  orgdb.obj <- get(orgdb)
  # get gene types
  all.keys <- AnnotationDbi::keytypes(orgdb.obj)
  # check gene types
  valid.gene.types <- intersect(all.keys, gene.types)
  if (length(valid.gene.types) >= 1) {
    # get matched gene name number across gene name type
    valid.gene.len <- sapply(valid.gene.types, function(x) {
      length(intersect(genes, AnnotationDbi::keys(orgdb.obj, keytype = x)))
    })
    # get max matched gene name number
    valid.gene.len.max <- valid.gene.len[which.max(valid.gene.len)]
    # get final gene name type
    if (valid.gene.len.max == 0) {
      stop("No valid gene name type detected in: ", paste0(valid.gene.types, collapse = ", "))
    } else {
      gene.type <- names(valid.gene.len.max)
      return(gene.type)
    }
  } else {
    stop("ENTREZID, SYMBOL, and ENSEMBL are not in provided annotation: ", orgdb)
  }
}
