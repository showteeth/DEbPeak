#' Prepare Count Matrix and Sample Metadata for Peak data (ChIP-seq/ATAC-seq).
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
#'
#' @return NULL
#' @importFrom utils read.table write.table
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @import DiffBind
#'
#' @export
#'
#' @examples
#' # library(DEbPeak)
#' # library(DESeq2)
#' # metadata file contains peak and bam information
#' # meta.file = 'path/to/metadata'
#' # PeakMatrix(meta.file = meta.file, blacklist = )
PeakMatrix <- function(meta.file, min.overlap = 2, submits = 200, use.summarizeOverlaps = TRUE,
                       blacklist = TRUE, sub.control = TRUE, used.cols = c("SampleID", "Condition"), out.folder = NULL) {
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
  # merge columns
  count.df <- count.df %>%
    dplyr::mutate(feature = paste(CHR, paste(START, END, sep = "-"), sep = ":")) %>%
    dplyr::select(-c(CHR, START, END))
  rownames(count.df) <- count.df$feature
  count.df$feature <- NULL
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
  }
  # save data
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # write count matrix
  write.table(
    x = count.df, file = file.path(out.folder, "consensus_peak_matrix.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
}
