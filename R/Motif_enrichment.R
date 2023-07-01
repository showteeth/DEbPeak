#' Motif Enrichment for Differentially Accessible/Binding Peaks.
#'
#' @param de.res Differential analysis results of peak-related data.
#' @param reg.type The regulation type, choose from Up_regulated, Not_regulated, Down_regulated.
#' @param signif Significance criterion to get differential peaks. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differential peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differential peaks. Default: 0.
#' @param homer.motif.path The path to Homer's \code{findMotifsGenome.pl}. Default: NULL.
#' @param genome Parameter for \code{findMotifsGenome.pl}, can be genome FASTA files or pre-build genome info. Default: mm10.
#' @param out.folder The output folder. Default: NULL (current directory).
#' @param other.paras Parameter for \code{findMotifsGenome.pl}, can be '-len 8,10,12 -size -100,50 -S 25'. Default: NULL.
#'
#' @return
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select case_when mutate_at vars filter
#' @importFrom tidyr drop_na separate
#' @importFrom purrr set_names
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

#' de novo Motif Discovery with STREME.
#'
#' @param peak.df Dataframe contains peaks, should contains chr, start, stop columns. Default: NULL.
#' @param diff.peak Dataframe contains differential peaks, the rownames should be peak annotations results. Default: NULL.
#' @param inte.res The Dataframe contains integrated results. \code{peak.df}, \code{diff.peak} and \code{inte.res}
#' cannot be empty at the same time. Dafault: NULL.
#' @param peak.motif.key The key type of integrated results ("Type" column of \code{inte.res}) or
#' regulation type of differential peaks (Up_regulated, Not_regulated, Down_regulated) to perfrom de novo motif discovery.
#' Used when \code{inte.res} is not NULL. Dafault: NULL.
#' @param genome The genome fasta path.
#' @param signif Significance criterion to get differential peaks. Used when \code{diff.peak} is not NULL.
#' For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differential peaks. Used when \code{diff.peak} is not NULL. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differential peaks. Used when \code{diff.peak} is not NULL. Default: 0.
#' @param peak.mode The source of \code{inte.res} results, choose from consensus (consensus mode) and diff (differential analysis).
#' Default: consensus.
#' @param samtools.path The path to samtools (the version should >=1.9). Default: NULL (conduct automatic detection).
#' @param streme.path The path to MEME's streme. Default: NULL (conduct automatic detection).
#' @param streme.paras Parameter for MEME's streme.
#' Default: "--nmotifs 50 --minw 8 --maxw 15 --thresh 0.05 --align center".
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param show.html Logical value, whether to show the results with pop-up window. Default: TRUE.
#'
#' @return NULL.
#' @importFrom magrittr %>%
#' @importFrom utils write.table browseURL
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select case_when mutate_at vars filter
#' @importFrom tidyr drop_na separate
#' @importFrom purrr set_names
#' @export
#'
#' @examples
#' # library(DEbPeak)
#' # # de novo motif discovery with peak sets
#' # peak.file = system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' # peak.df = GetConsensusPeak(peak.file = peak.file)
#' # MotifDiscovery(peak.df = peak.df, genome = '/path/to/genome.fa',
#' #                streme.path = "/path/to/streme", samtools.path = "/path/to/samtools", out.folder = "/path/to/output")
#' # # de novo motif discovery with differential peaks
#' # MotifDiscovery(diff.peak = dds.peak.results.ordered, peak.motif.key = "Up_regulated",
#' #                genome = '/path/to/genome.fa', streme.path = "/path/to/streme",
#' #                samtools.path = "/path/to/samtools", out.folder = "/path/to/output")
#' # # de novo motif discovery with integrated results
#' # MotifDiscovery(inte.res = debatac.res, peak.motif.key = "Up_Up",
#' #                genome = '/path/to/genome.fa', peak.mode = "diff", streme.path = "/path/to/streme",
#' #                samtools.path = "/path/to/samtools", out.folder = "/path/to/output")
MotifDiscovery <- function(peak.df = NULL, diff.peak = NULL, inte.res = NULL, peak.motif.key = NULL, genome,
                           signif = "padj", signif.threshold = 0.05, l2fc.threshold = 0,
                           peak.mode = c("consensus", "diff"), samtools.path = NULL, streme.path = NULL,
                           streme.paras = "--nmotifs 50 --minw 8 --maxw 15 --thresh 0.05 --align center",
                           out.folder = NULL, show.html = TRUE) {
  # check prameters
  peak.mode <- match.arg(arg = peak.mode)
  # get samtools path
  if (is.null(samtools.path)) {
    # detect samtools path
    samtools.path <- Sys.which("samtools")
    if (samtools.path == "") {
      stop("Can not find samtools automatically, please specify the path!")
    }
  } else {
    samtools.path <- samtools.path
  }
  # check samtools version
  samtools.faidx.run <- system(paste0(samtools.path, " faidx -h"), intern = TRUE)
  samtools.faidx.run.code <- attr(samtools.faidx.run, "status")
  if (!is.null(samtools.faidx.run.code)) {
    stop("Run samtools version error!")
  }
  if (!grepl(pattern = "region-file", x = paste0(samtools.faidx.run, collapse = " "))) {
    stop("Please update you samtools to support --region-file parameter (faidx)!")
  }
  # get meme path
  if (is.null(streme.path)) {
    # detect streme path
    streme.path <- Sys.which("streme")
    if (streme.path == "") {
      stop("Can not find streme automatically, please specify the path!")
    }
  } else {
    streme.path <- streme.path
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # prepare peak
  if (is.null(peak.df) && is.null(inte.res) && is.null(diff.peak)) {
    stop("Please provide one of peak.df, diff.peak and inte.res!")
  }
  if (!is.null(peak.df)) {
    message("Run motif discovery on consensus peaks!")
    # check peak.df
    used.cols <- c("chr", "start", "stop")
    if (!all(used.cols %in% colnames(peak.df))) {
      stop("The peak.df you provided must contains chr, start and stop columns, please check!")
    }
    used.peak <- peak.df %>%
      dplyr::mutate(Peak = paste(chr, paste(start, stop, sep = "-"), sep = ":")) %>%
      dplyr::select(Peak)
  } else if (!is.null(diff.peak)) {
    message("Run motif discovery on differential peaks!")
    # prepare used dataframe
    diff.peak.df <- PrepareDEPlot(
      deres = diff.peak, signif = signif, signif.threshold = signif.threshold,
      l2fc.threshold = l2fc.threshold, label.key = NULL
    )
    used.peak <- diff.peak.df %>%
      dplyr::filter(regulation == peak.motif.key) %>%
      dplyr::select(Peak_Gene = Gene) %>%
      dplyr::filter(!is.na(Peak_Gene)) %>%
      tidyr::separate(col = Peak_Gene, into = c("Peak", "Peak_Gene", "Peak_Type"), sep = "\\|") %>%
      dplyr::select(Peak)
  } else if (!is.null(inte.res)) {
    if (peak.mode == "consensus") {
      message("Run motif discovery on DEbPeak (consensus mode)!")
      used.peak <- inte.res %>%
        dplyr::filter(Type == peak.motif.key) %>%
        dplyr::filter(!is.na(Peak)) %>%
        dplyr::select(Peak)
    } else if (peak.mode == "diff") {
      message("Run motif discovery on DEbPeak (diff mode)!")
      used.peak <- inte.res %>%
        dplyr::filter(Type == peak.motif.key) %>%
        dplyr::filter(!is.na(Peak_Gene)) %>%
        dplyr::select(Peak_Gene) %>%
        tidyr::separate(col = Peak_Gene, into = c("Peak", "Peak_Gene", "Peak_Type"), sep = "\\|") %>%
        dplyr::select(Peak)
    }
  }
  # write region file
  peak.region.file <- file.path(out.folder, paste0("DEbPeak_", peak.motif.key, ".txt"))
  write.table(x = used.peak, file = peak.region.file, quote = FALSE, col.names = FALSE, row.names = FALSE)
  # extract region fasta
  peak.fasta.file <- file.path(out.folder, paste0("DEbPeak_", peak.motif.key, ".fa"))
  # run samtools command
  samtools.faidx.cmd <- paste(samtools.path, "faidx", genome, "--region-file", peak.region.file, ">", peak.fasta.file)
  message(paste("Calling samtools on", peak.motif.key, ":", samtools.faidx.cmd))
  samtools.faidx.cmd.status <- system(samtools.faidx.cmd, intern = TRUE)
  samtools.faidx.cmd.status.code <- attr(samtools.faidx.cmd.status, "status")
  if (!is.null(samtools.faidx.cmd.status.code)) {
    stop("Run samtools faidx error!")
  }
  # run streme
  streme.cmd <- paste(streme.path, "--p", peak.fasta.file, "--dna --oc", file.path(out.folder, "streme_denovo"), streme.paras)
  message(paste("Calling streme on", peak.motif.key, ":", streme.cmd))
  streme.cmd.status <- system(streme.cmd, intern = TRUE)
  streme.cmd.status.code <- attr(streme.cmd.status, "status")
  if (!is.null(streme.cmd.status.code)) {
    stop("Run streme error!")
  }
  if (show.html) {
    # open result html
    browseURL(file.path(out.folder, "streme_denovo", "streme.html"))
  }
}


#' Map Motifs against a Motif Database.
#'
#' @param motif.folder Folder contains \code{\link{MotifDiscovery}} results.
#' @param known.motif Motif database, should be in meme format.
#' @param tomtom.path The path to MEME's tomtom. Default: NULL (conduct automatic detection).
#' @param tomtom.paras Parameter for MEME's tomtom. Default: "-no-ssc -min-overlap 5 -dist pearson -evalue -thresh 10.0".
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param show.html Logical value, whether to show the results with pop-up window. Default: TRUE.
#'
#' @return Dataframe contains motif mapping results.
#' @export
#'
#' @examples
#' # library(DEbPeak)
#' # motif.cmp = MotifCompare(motif.folder = "/path/to/MotifDiscovery/out", known.motif = "/path/to/database.meme",
#' #                          tomtom.path = '/path/to/tomtom', out.folder = "/path/to/output")
MotifCompare <- function(motif.folder, known.motif, tomtom.path = NULL,
                         tomtom.paras = "-no-ssc -min-overlap 5 -dist pearson -evalue -thresh 10.0",
                         out.folder = NULL, show.html = TRUE) {
  # get tomtom path
  if (is.null(tomtom.path)) {
    # detect tomtom path
    tomtom.path <- Sys.which("tomtom")
    if (tomtom.path == "") {
      stop("Can not find tomtom automatically, please specify the path!")
    }
  } else {
    tomtom.path <- tomtom.path
  }
  # detect motif files
  if (file.exists(file.path(motif.folder, "streme_denovo", "streme.txt"))) {
    motif.file <- file.path(motif.folder, "streme_denovo", "streme.txt")
  } else if (file.exists(file.path(motif.folder, "streme.txt"))) {
    motif.file <- file.path(motif.folder, "streme.txt")
  } else {
    stop("There is no streme.txt detected!")
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # run tomtom
  tomtom.cmd <- paste(tomtom.path, "-oc", file.path(out.folder, "tomtom"), tomtom.paras, motif.file, known.motif)
  message(paste("Calling tomtom:", tomtom.cmd))
  tomtom.cmd.status <- system(tomtom.cmd, intern = TRUE)
  tomtom.cmd.status.code <- attr(tomtom.cmd.status, "status")
  if (!is.null(tomtom.cmd.status.code)) {
    stop("Run tomtom error!")
  }
  tomtom.res <- utils::read.table(file = file.path(out.folder, "tomtom", "tomtom.tsv"), header = TRUE)
  if (show.html) {
    # open result html
    browseURL(file.path(out.folder, "tomtom", "tomtom.html"))
  }
  return(tomtom.res)
}
