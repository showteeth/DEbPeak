#' Find Motif of Peaks from RNA-seq and Peak-related Data Integrated Results.
#'
#' @param inte.res DEGs and peak annotation integration results.
#' @param peak.motif.key The key type of integrated results ("Type" column of \code{de.peak}) to find motif.
#' @param homer.motif.path The path to Homer's \code{findMotifsGenome.pl}. Default: NULL.
#' @param genome Parameter for \code{findMotifsGenome.pl}, can be genome FASTA files or pre-build genome info. Default: mm10.
#' @param out.folder The output folder. Default: NULL (current directory).
#' @param other.paras Parameter for \code{findMotifsGenome.pl}, can be '-len 8,10,12 -size -100,50 -S 25'. Default: NULL.
#'
#' @return List contains known motif results of UPbATAC and DOWNbATAC.
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr separate
#' @export
#'
FindMotif <- function(inte.res, peak.motif.key, peak.mode = c("consensus", "diff"),
                      homer.motif.path = NULL, genome = "mm10", out.folder = NULL, other.paras = NULL) {
  # check parameters
  peak.mode <- match.arg(arg = peak.mode)

  # get homer motif find path
  if (is.null(homer.motif.path)) {
    # detect homer motif find path
    homer.motif.path <- Sys.which("findMotifsGenome.pl")
    if (homer.motif.path == "") {
      stop("Can not find findMotifsGenome.pl automatically, please specify the path!")
    }
  } else {
    homer.motif.path <- homer.motif.path
  }
  export.path <- paste0('export PATH="$PATH:', dirname(homer.motif.path), '";')

  # set output folder
  if (is.null(out.folder)) {
    out.folder <- file.path(getwd(), "Motif")
  } else {
    out.folder <- file.path(out.folder, "Motif")
  }

  # prepare peak dataframe
  # get genes
  # inte.genes <- inte.res[inte.res$Type == peak.motif.key, gene.key]
  used.inte <- inte.res[inte.res$Type == peak.motif.key, ]
  # upbchip.genes <- inte.res[inte.res$Type == "UPbPeak", gene.key]
  # downbchip.genes <- inte.res[inte.res$Type == "DOWNbPeak", gene.key]

  # find motif
  if (nrow(used.inte) >= 1) {
    # get peak
    if (peak.mode == "consensus") {
      inte.peak <- inte.res %>%
        dplyr::filter(Type == peak.motif.key) %>%
        dplyr::filter(!is.na(Peak)) %>%
        dplyr::select(Peak) %>%
        dplyr::mutate(name = Peak) %>%
        tidyr::separate(col = Peak, into = c("seqnames", "region"), sep = ":") %>%
        tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>%
        dplyr::select(c("seqnames", "start", "end", "name"))
      inte.peak$score <- 1
    } else if (peak.mode == "diff") {
      inte.peak <- inte.res %>%
        dplyr::filter(Type == peak.motif.key) %>%
        dplyr::filter(!is.na(Peak_Gene)) %>%
        dplyr::select(Peak_Gene) %>%
        dplyr::mutate(name = Peak_Gene) %>%
        tidyr::separate(col = Peak_Gene, into = c("Peak", "Peak_Gene", "Peak_Type"), sep = "\\|") %>%
        tidyr::separate(col = Peak, into = c("seqnames", "region"), sep = ":") %>%
        tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>%
        dplyr::select(c("seqnames", "start", "end", "name"))
      inte.peak$score <- 1
    }
    inte.peak$strand <- "+"
    # write bed
    peak.bed.file <- tempfile(pattern = peak.motif.key, fileext = ".bed")
    write.table(
      x = inte.peak, file = peak.bed.file, quote = FALSE, sep = "\t",
      col.names = FALSE, row.names = FALSE
    )
    # command used
    peak.out.folder <- file.path(out.folder, peak.motif.key)
    peak.homer.motif.cmd <- paste(homer.motif.path, peak.bed.file, genome, peak.out.folder, other.paras)
    # add path
    peak.homer.motif.cmd <- paste0(export.path, peak.homer.motif.cmd)

    # run command
    message(paste("Calling findMotifsGenome.pl on", peak.motif.key, ":", peak.homer.motif.cmd))
    peak.homer.motif.status <- system(peak.homer.motif.cmd, intern = TRUE)
    peak.homer.motif.status.code <- attr(peak.homer.motif.status, "status")
    if (!is.null(peak.homer.motif.status.code)) {
      stop("Run findMotifsGenome.pl error!")
    }
    # read peak known motif
    peak.known.motif <- data.table::fread(file = file.path(peak.out.folder, "knownResults.txt"), sep = "\t") %>%
      as.data.frame()
  } else {
    peak.known.motif <- data.frame()
  }

  return(peak.known.motif)
  # # extract peaks
  # ## up bind peak
  # if (length(upbchip.genes) >= 1) {
  #   # get peak
  #   upbchip.peak <- peak.anno.res[
  #     peak.anno.res[[gene.key]] %in% upbchip.genes,
  #     c("seqnames", "start", "end", "name", "score")
  #   ]
  #   upbchip.peak$strand <- "+"
  #   # write bed
  #   up.bed.file <- tempfile(pattern = "UP", fileext = ".bed")
  #   write.table(
  #     x = upbchip.peak, file = up.bed.file, quote = FALSE, sep = "\t",
  #     col.names = FALSE, row.names = FALSE
  #   )
  #   # command used
  #   up.out.folder <- file.path(out.folder, "UPbATAC")
  #   up.homer.motif.cmd <- paste(homer.motif.path, up.bed.file, genome, up.out.folder, other.paras)
  #   # add path
  #   up.homer.motif.cmd <- paste0(export.path, up.homer.motif.cmd)
  #
  #   # run command
  #   message(paste("Calling findMotifsGenome.pl on UPbATAC: ", up.homer.motif.cmd))
  #   up.homer.motif.status <- system(up.homer.motif.cmd, intern = TRUE)
  #   up.homer.motif.status.code <- attr(up.homer.motif.status, "status")
  #   if (!is.null(up.homer.motif.status.code)) {
  #     stop("Run findMotifsGenome.pl error!")
  #   }
  #   # read up known motif
  #   up.known.motif <- data.table::fread(file = file.path(up.out.folder, "knownResults.txt"), sep = "\t")
  # } else {
  #   up.known.motif <- data.frame()
  # }
  # ## down bind peak
  # if (length(downbchip.genes) >= 1) {
  #   downbchip.peak <- peak.anno.res[
  #     peak.anno.res[[gene.key]] %in% downbchip.genes,
  #     c("seqnames", "start", "end", "name", "score")
  #   ]
  #   downbchip.peak$strand <- "+"
  #   # write bed
  #   down.bed.file <- tempfile(pattern = "DOWN", fileext = ".bed")
  #   write.table(
  #     x = downbchip.peak, file = down.bed.file, quote = FALSE, sep = "\t",
  #     col.names = FALSE, row.names = FALSE
  #   )
  #   # command used
  #   down.out.folder <- file.path(out.folder, "DOWNbATAC")
  #   down.homer.motif.cmd <- paste(homer.motif.path, down.bed.file, genome, down.out.folder, other.paras)
  #   # add path
  #   down.homer.motif.cmd <- paste0(export.path, down.homer.motif.cmd)
  #
  #   # run command
  #   message(paste("Calling findMotifsGenome.pl on DOWNbATAC: ", down.homer.motif.cmd))
  #   down.homer.motif.status <- system(down.homer.motif.cmd, intern = TRUE)
  #   down.homer.motif.status.code <- attr(down.homer.motif.status, "status")
  #   if (!is.null(down.homer.motif.status.code)) {
  #     stop("Run findMotifsGenome.pl error!")
  #   }
  #   # read down known motif
  #   down.known.motif <- data.table::fread(file = file.path(down.out.folder, "knownResults.txt"), sep = "\t")
  # } else {
  #   down.known.motif <- data.frame()
  # }
  # # final results
  # final.res <- list(UPbATAC = up.known.motif, DOWNbATAC = down.known.motif)
  # return(final.res)
}
