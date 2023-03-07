#' Motif Enrichment for Differentially Accessible/Binding Peaks.
#'
#' @param de.res Differential analysis results of peak-related data.
#' @param reg.type The regulation type, choose from Up_regulated, Not_regulated, Down_regulated.
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed results. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed results. Default: 0.
#' @param homer.motif.path The path to Homer's \code{findMotifsGenome.pl}. Default: NULL.
#' @param genome Parameter for \code{findMotifsGenome.pl}, can be genome FASTA files or pre-build genome info. Default: mm10.
#' @param out.folder The output folder. Default: NULL (current directory).
#' @param other.paras Parameter for \code{findMotifsGenome.pl}, can be '-len 8,10,12 -size -100,50 -S 25'. Default: NULL.
#'
#' @return
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select case_when mutate_at vars filter
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @importFrom tidyr separate
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' # dds.peak.results.ordered is differential analysis of peak-related data.
#' # peak.motif = MotifEnrich(de.res = dds.peak.results.ordered, reg.type = "Up_regulated",
#' #                          homer.motif.path = '~/anaconda3/bin/findMotifsGenome.pl',
#' #                          out.folder = "/path/to/out/folder", other.paras = "-len 8,10,12 -size -100,50 -S 25")
MotifEnrich <- function(de.res, reg.type, signif = "padj", signif.threshold = 0.05, l2fc.threshold = 0,
                        homer.motif.path = NULL, genome = "mm10", out.folder = NULL, other.paras = NULL) {

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
  # prepare used dataframe
  de.df <- PrepareDEPlot(
    deres = de.res, signif = signif, signif.threshold = signif.threshold,
    l2fc.threshold = l2fc.threshold, label.key = NULL
  )
  reg.de.df <- de.df[de.df$regulation == reg.type, ]
  # run
  if (nrow(reg.de.df) >= 1) {
    # prepare bed file
    reg.de.peak <- reg.de.df %>%
      dplyr::filter(!is.na(Gene)) %>%
      dplyr::select(Gene) %>%
      dplyr::mutate(name = Gene) %>%
      tidyr::separate(col = Gene, into = c("Peak", "Peak_Gene", "Peak_Type"), sep = "\\|") %>%
      tidyr::separate(col = Peak, into = c("seqnames", "region"), sep = ":") %>%
      tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>%
      dplyr::select(c("seqnames", "start", "end", "name"))
    reg.de.peak$score <- 1
    # write bed
    reg.bed.file <- tempfile(pattern = reg.type, fileext = ".bed")
    write.table(
      x = reg.de.peak, file = reg.bed.file, quote = FALSE, sep = "\t",
      col.names = FALSE, row.names = FALSE
    )
    # command used
    reg.out.folder <- file.path(out.folder, reg.type)
    reg.homer.motif.cmd <- paste(homer.motif.path, reg.bed.file, genome, reg.out.folder, other.paras)
    # add path
    reg.homer.motif.cmd <- paste0(export.path, reg.homer.motif.cmd)
    # run command
    message(paste("Calling findMotifsGenome.pl on", reg.type, ":", reg.homer.motif.cmd))
    reg.homer.motif.status <- system(reg.homer.motif.cmd, intern = TRUE)
    reg.homer.motif.status.code <- attr(reg.homer.motif.status, "status")
    if (!is.null(reg.homer.motif.status.code)) {
      stop("Run findMotifsGenome.pl error!")
    }
    # read reg known motif
    reg.known.motif <- data.table::fread(file = file.path(reg.out.folder, "knownResults.txt"), sep = "\t") %>%
      as.data.frame()
  } else {
    reg.known.motif <- data.frame()
  }
  return(reg.known.motif)
}
