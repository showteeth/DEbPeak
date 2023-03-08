#' Prepare Count Matrix and Sample Metadata for Peak-related data.
#'
#' @param meta.file Sample metadata contains peak related information (eg: sample, peakPath, bamPath, condition).
#' @param min.overlap Only include peaks in at least this many peaksets in the main binding matrix. Default: 2.
#' Parameter of \code{\link{dba}}.
#' @param submits If the value is greater than zero, all consensus peaks will be re-centered around a consensus summit,
#' with the value of summits indicating how many base pairs to include upstream and downstream of the
#' summit (so all consensus peaks will be of the same width, namely \code{2 * summits + 1}).
#' Default: 200. Parameter of \code{\link{dba.count}}.
#' @param use.summarizeOverlaps Logical value, indicating that \code{\link{summarizeOverlaps}} should be used for counting
#' instead of the built-in counting code. Default: TRUE. Parameter of \code{\link{dba.count}}.
#' @param blacklist Species-specific abnormal regions to be removed. Choose from "DBA_BLACKLIST_HG19", "DBA_BLACKLIST_HG38", "DBA_BLACKLIST_GRCH37",
#' "DBA_BLACKLIST_GRCH38", "DBA_BLACKLIST_MM9", "DBA_BLACKLIST_MM10", "DBA_BLACKLIST_CE10", "DBA_BLACKLIST_CE11", "DBA_BLACKLIST_DM3", "DBA_BLACKLIST_DM6",
#' TRUE (auto-detection genome), a GRanges object containing the blacklisted regions. Default:TRUE. Parameter of \code{\link{dba.blacklist}}.
#' @param sub.control Logical value, whether Control read counts are subtracted for each site in each sample.
#' Default: TRUE. Parameter of \code{\link{dba.count}}.
#' @param used.cols Used columns used to create sample metadata. If specified, sampleID should be placed first. Default: c("SampleID", "Condition").
#' @param out.folder Output folder to save created count matrix and sample metadata. Default: NULL (current working directory).
#' @param species Species used, chosen from "Human","Mouse","Rat","Fly","Arabidopsis","Yeast","Zebrafish","Worm","Bovine","Pig","Chicken","Rhesus",
#' "Canine","Xenopus","Anopheles","Chimp","E coli strain Sakai","Myxococcus xanthus DK 1622". Default: "Human".
#' @param seq.style The style of sequence, chosen from UCSC, NCBI, Ensembl, None. This should be compatible with the genome and gtf file you used
#' to generate count matrix and peak files. Default: "UCSC".
#' @param gtf.file GTF file used to create TxDb object. Useful when specie you used is not available in \code{species}. Default: NULL.
#' @param up.dist The upstream distance from the TSS. Default: 3000bp.
#' @param down.dist The downstream distance from the TSS. Default: 3000bp.
#' @param ... Parameters for \code{\link{annotatePeak}}.
#'
#' @return A dataframe contains count matrix, peak annotation and sample metadata (if provided \code{used.cols}).
#' @importFrom utils read.table write.table
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select case_when
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom BiocManager install
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle keepSeqlevels seqnames
#' @importFrom ggplotify as.ggplot
#' @importFrom cowplot plot_grid
#' @import ChIPseeker
#' @import ggupset
#' @import DiffBind
#' @import parallel
#'
#' @export
#'
#' @examples
#' # library(DEbPeak)
#' # library(DESeq2)
#' # metadata file contains peak and bam information
#' # meta.file = 'path/to/metadata'
#' # PeakMatrix(meta.file = meta.file, species = "Human",  seq.style = "UCSC",
#' #            up.dist = 20000, down.dist = 20000)
PeakMatrix <- function(meta.file, min.overlap = 2, submits = 200, use.summarizeOverlaps = TRUE,
                       blacklist = TRUE, sub.control = TRUE, used.cols = c("SampleID", "Condition"),
                       out.folder = NULL, species = c(
                         "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                         "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                       ),
                       seq.style = c("UCSC", "NCBI", "Ensembl", "None"), gtf.file = NULL, up.dist = 3000, down.dist = 3000, ...) {
  # load sample meta info
  sample.info <- read.table(meta.file, sep = "\t", header = TRUE, check.names = FALSE)
  # check columns
  required.cols <- c("SampleID", "bamReads", "bamControl", "Peaks", "PeakCaller")
  if (!all(intersect(colnames(sample.info), required.cols) == required.cols)) {
    stop("The metadata should contain SampleID, bamReads, bamControl, Peaks, PeakCaller columns!")
  }
  # create dba object
  dbaObj <- DiffBind::dba(sampleSheet = sample.info, minOverlap = min.overlap)
  # count
  suppressWarnings(suppressMessages(library(parallel)))
  dbaObj <- DiffBind::dba.count(dbaObj, bUseSummarizeOverlaps = use.summarizeOverlaps)
  # remove blacklist region
  if (!is.null(blacklist)) {
    if (isTRUE(blacklist)) {
      dbaObj <- DiffBind::dba.blacklist(dbaObj, blacklist = blacklist, greylist = FALSE)
    } else if (blacklist %in% c(
      "DBA_BLACKLIST_HG19", "DBA_BLACKLIST_HG38", "DBA_BLACKLIST_GRCH37",
      "DBA_BLACKLIST_GRCH38", "DBA_BLACKLIST_MM9", "DBA_BLACKLIST_MM10",
      "DBA_BLACKLIST_CE10", "DBA_BLACKLIST_CE11", "DBA_BLACKLIST_DM3", "DBA_BLACKLIST_DM6"
    )) {
      dbaObj <- DiffBind::dba.blacklist(dbaObj, blacklist = get(blacklist), greylist = FALSE)
    } else {
      stop("Please provide valid blacklist!")
    }
  }
  # extract count
  if (sub.control) {
    dbaObj <- DiffBind::dba.count(dbaObj, peaks = NULL, score = DBA_SCORE_READS_MINUS)
  } else {
    dbaObj <- DiffBind::dba.count(dbaObj, peaks = NULL, score = DBA_SCORE_READS)
  }
  count.df <- DiffBind::dba.peakset(dbaObj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
  # prepare peak dataframe for peak annotation
  peak.df <- count.df[c("CHR", "START", "END")]
  colnames(peak.df) <- c("chr", "start", "stop")
  peak.df$name <- paste(peak.df$chr, paste(peak.df$start, peak.df$stop, sep = "-"), sep = ":")
  # peak annotation
  peak.anno <- AnnoPeak(
    peak.df = peak.df, species = species, seq.style = seq.style,
    gtf.file = gtf.file, up.dist = up.dist, down.dist = down.dist, ...
  )
  peak.anno.df <- peak.anno$df
  # change anno type
  peak.anno.df <- peak.anno.df %>%
    dplyr::mutate(AnnoType = dplyr::case_when(
      anno == "Promoter" ~ "P",
      anno == "5' UTR" ~ "5U",
      anno == "3' UTR" ~ "3U",
      anno == "Exon" ~ "E",
      anno == "Intron" ~ "I",
      anno == "Downstream" ~ "D",
      anno == "Distal Intergenic" ~ "DI"
    )) %>%
    dplyr::mutate(FullName = paste(name, SYMBOL, AnnoType, sep = "|")) %>%
    dplyr::select(c("name", "FullName"))
  # merge count columns
  count.df <- count.df %>%
    dplyr::mutate(name = paste(CHR, paste(START, END, sep = "-"), sep = ":")) %>%
    dplyr::select(-c(CHR, START, END))
  # merge paek annotation and count
  count.df <- merge(count.df, peak.anno.df, by = "name", all.x = TRUE)
  # filter no peak annotation info, usually abnormal chromosome, eg: chr1_GL456211_random
  count.df <- count.df[!is.na(count.df$FullName), ]
  rownames(count.df) <- count.df$FullName
  count.df$name <- NULL
  count.df$FullName <- NULL
  # save data
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # write count matrix
  write.table(
    x = count.df, file = file.path(out.folder, "consensus_peak_matrix.txt"),
    sep = "\t", quote = FALSE, col.names = TRUE
  )
  # write peak annotation results
  write.table(
    x = peak.anno$df, file = file.path(out.folder, "consensus_peak_anno.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  # create metadata
  if (!is.null(used.cols)) {
    valid.used.cols <- intersect(colnames(sample.info), used.cols)
    if (!all(valid.used.cols == used.cols)) {
      stop("The metadata does not contain all columns specified by used.cols!")
    }
    sample.meta <- sample.info[used.cols]
    rownames(sample.meta) <- sample.meta[, 1]
    sample.meta[, 1] <- NULL
    # write sample meta
    write.table(
      x = sample.meta, file = file.path(out.folder, "peak_metadata.txt"),
      sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
    )
    return(list(count = count.df, meta = sample.meta, peak_anno = peak.anno$df))
  } else {
    return(list(count = count.df, peak_anno = peak.anno$df))
  }
}
