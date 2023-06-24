# get txdb and orgdb for species
GetSpeciesAnno <- function(species) {
  spe_names <- c(
    "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
    "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
  )
  # Genome wide annotation
  spe_db <- c(
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.At.tair.db", "org.Sc.sgd.db", "org.Dr.eg.db", "org.Ce.eg.db", "org.Bt.eg.db", "org.Ss.eg.db",
    "org.Gg.eg.db", "org.Mmu.eg.db", "org.Cf.eg.db", "org.Xl.eg.db", "org.Ag.eg.db", "org.Pt.eg.db", "org.EcSakai.eg.db", "org.Mxanthus.db"
  )
  # KEGG organism
  spe_organ <- c("hsa", "mmu", "rno", "dme", "ath", "sce", "dre", "cel", "bta", "ssc", "gga", "mcc", "cfa", "xla", "aga", "ptr", "ecs", "mxa")
  # txdb
  spe_txdb <- c(
    "TxDb.Hsapiens.UCSC.hg38.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Rnorvegicus.UCSC.rn6.refGene",
    "TxDb.Dmelanogaster.UCSC.dm6.ensGene", "TxDb.Athaliana.BioMart.plantsmart28", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
    "TxDb.Drerio.UCSC.danRer11.refGene", "TxDb.Celegans.UCSC.ce11.refGene", "TxDb.Btaurus.UCSC.bosTau9.refGene",
    "TxDb.Sscrofa.UCSC.susScr11.refGene", "TxDb.Ggallus.UCSC.galGal6.refGene", "TxDb.Mmulatta.UCSC.rheMac10.refGene",
    "TxDb.Cfamiliaris.UCSC.canFam3.refGene", "", "", "TxDb.Ptroglodytes.UCSC.panTro6.refGene", "", ""
  )
  names(spe_db) <- spe_names
  names(spe_organ) <- spe_names
  names(spe_txdb) <- spe_names
  species_info <- list(
    OrgDb = as.character(spe_db[species]),
    Organism = as.character(spe_organ[species]),
    txdb = as.character(spe_txdb[species])
  )
  return(species_info)
}

#' Create ChIP Peak Binding Profile.
#'
#' @param peak.df Dataframe contains all consensus peaks.
#' @param species Species used, chosen from "Human","Mouse","Rat","Fly","Arabidopsis",
#' "Yeast","Zebrafish","Worm","Bovine","Pig","Chicken","Rhesus",
#' "Canine","Xenopus","Anopheles","Chimp","E coli strain Sakai","Myxococcus xanthus DK 1622". Default: "Human".
#' @param gtf.file GTF file used to create TxDb object. Useful when specie you used is not available in \code{species}. Default: NULL.
#' @param weight.col The column name of weight. Parameter of \code{\link{peakHeatmap}},
#' \code{\link{plotAvgProf2}}, \code{\link{plotPeakProf2}}. Default: NULL.
#' @param up.dist The upstream distance from the TSS. Default: 3000bp.
#' @param down.dist The downstream distance from the TSS. Default: 3000bp.
#' @param color The color used for heatmap. Parameter of \code{\link{peakHeatmap}}. Default: red.
#' @param conf The confidence interval. Parameter of \code{\link{plotAvgProf2}}, \code{\link{plotPeakProf2}}. Default: 0.95.
#' @param up.rel The percentage of total length of body regions used as upstream distance from the TSS.
#' Parameter of \code{\link{plotPeakProf2}}. When set NULL, use \code{up.dist}. Default: 0.2.
#' @param down.rel The percentage of total length of body regions used as downstream distance from the TSS.
#' Parameter of \code{\link{plotPeakProf2}}. When set NULL, use \code{down.dist}. Default: 0.2.
#' @param by One of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR'. Parameter of \code{\link{plotPeakProf2}}. Default: 'gene'.
#' @param region.type One of "start_site", "end_site", "body". Parameter of \code{\link{plotPeakProf2}}. Default: 'start_site'.
#' @param nbin The amount of nbines. Parameter of \code{\link{plotPeakProf2}}. Default: NULL.
#'
#' @return List contains all profile plots.
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom BiocManager install
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom ggplotify as.ggplot
#' @importFrom cowplot plot_grid
#' @import ChIPseeker
#' @export
#'
#' @examples
#' library(DEbPeak)
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.profile <- PeakProfile(peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
PeakProfile <- function(peak.df, species = c(
                          "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                          "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                        ),
                        gtf.file = NULL, weight.col = NULL, up.dist = 3000, down.dist = 3000,
                        color = "red", conf = 0.95, up.rel = 0.2, down.rel = 0.2, by = c("gene", "transcript", "exon", "intron", "3UTR", "5UTR"),
                        region.type = c("start_site", "end_site", "body"), nbin = NULL) {
  # check parameters
  species <- match.arg(arg = species)
  by <- match.arg(arg = by)
  region.type <- match.arg(arg = region.type)

  # prepare txdb
  if (!is.null(gtf.file)) {
    message("Create txdb from gtf file!")
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.file)
    # annotation
    txdb.obj <- txdb
  } else {
    spe.anno <- GetSpeciesAnno(species)
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

  # preare peak information used
  # notice that the output file of MSPC is in bed format
  peak.gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = peak.df, keep.extra.columns = TRUE, ignore.strand = TRUE,
    seqnames.field = "chr", start.field = "start", end.field = "stop",
    starts.in.df.are.0based = TRUE
  )

  # create heatmap plot
  peak.heatmap <- ggplotify::as.ggplot(function() {
    peakHeatmap(
      peak = peak.gr, weightCol = weight.col,
      TxDb = txdb.obj, upstream = up.dist, downstream = down.dist,
      title = "Heatmap of Peak binding to TSS regions", color = color
    )
  })

  # create Average Profile
  avg.profile.plot <- plotAvgProf2(peak.gr,
    TxDb = txdb.obj, upstream = up.dist, downstream = down.dist,
    conf = conf, weightCol = weight.col,
    xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency"
  )
  # Profile of ChIP peaks binding to body regions
  if (is.null(up.rel) | is.null(down.rel)) {
    peak.type.profile <- plotPeakProf2(
      peak = peak.gr, upstream = up.dist, downstream = down.dist, conf = conf,
      weightCol = weight.col, by = by, type = region.type, nbin = nbin,
      TxDb = txdb.obj
    )
  } else if (up.rel <= 1 & up.rel > 0 & down.rel > 0 & down.rel <= 1) {
    peak.type.profile <- plotPeakProf2(
      peak = peak.gr, upstream = rel(up.rel), downstream = rel(down.rel), conf = conf,
      weightCol = weight.col, by = by, type = region.type, nbin = nbin,
      TxDb = txdb.obj
    )
  }

  # merge plot
  profile.plot <- cowplot::plot_grid(peak.heatmap, cowplot::plot_grid(avg.profile.plot, peak.type.profile, ncol = 1))
  plot.list <- list(
    profile.plot = profile.plot, peak.heatmap = peak.heatmap,
    avg.profile.plot = avg.profile.plot, peak.type.profile = peak.type.profile
  )
  return(plot.list)
}

#' Create Peak Annotation Pie Plot with ggplot2.
#'
#' @param peak.obj Peak annotation object, returned from \code{\link{annotatePeak}}.
#' @param ... Parameters for \code{\link{ggpie}}.
#'
#' @return A ggplot2 object.
#' @importFrom ggpie ggpie
#' @export
#'
#' @examples
#' library(DEbPeak)
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.anno <- AnnoPeak(
#'   peak.df = peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
#' PeakAnnoPie(peak.anno$anno.obj)
PeakAnnoPie <- function(peak.obj, ...) {
  freq.stat <- peak.obj@annoStat
  colnames(freq.stat) <- c("group", "count")
  freq.stat$count <- round(freq.stat$count, 2)
  pie.plot <- ggpie(
    data = freq.stat, count_type = "count", label_pos = "out",
    label_type = "horizon", label_info = "ratio", ...
  ) +
    theme(legend.title = element_blank())
  return(pie.plot)
}

#' Conduct Peak Annotation.
#'
#' @param peak.df Dataframe contains all consensus peaks.
#' @param species Species used, chosen from "Human","Mouse","Rat","Fly","Arabidopsis","Yeast","Zebrafish","Worm","Bovine","Pig","Chicken","Rhesus",
#' "Canine","Xenopus","Anopheles","Chimp","E coli strain Sakai","Myxococcus xanthus DK 1622". Default: "Human".
#' @param seq.style The style of sequence, chosen from UCSC, NCBI, Ensembl, None. This should be compatible with the genome and gtf file you used
#' to generate count matrix and peak files. Default: "UCSC".
#' @param gtf.file GTF file used to create TxDb object. Useful when specie you used is not available in \code{species}. Default: NULL.
#' @param up.dist The upstream distance from the TSS. Default: 3000bp.
#' @param down.dist The downstream distance from the TSS. Default: 3000bp.
#' @param ... Parameters for \code{\link{annotatePeak}}.
#'
#' @return List contains peak annotation dataframe and plot.
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom BiocManager install
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle keepSeqlevels seqnames
#' @importFrom ggplotify as.ggplot
#' @importFrom cowplot plot_grid
#' @import ChIPseeker
#' @import ggupset
#' @export
#'
#' @examples
#' library(DEbPeak)
#' peak.file <- system.file("extdata", "debchip_peaks.bed", package = "DEbPeak")
#' peak.df <- GetConsensusPeak(peak.file = peak.file)
#' peak.profile <- PeakProfile(peak.df, species = "Mouse", by = "gene", region.type = "body", nbin = 800)
#' peak.anno <- AnnoPeak(
#'   peak.df = peak.df, species = "Mouse",
#'   seq.style = "UCSC", up.dist = 20000, down.dist = 20000
#' )
AnnoPeak <- function(peak.df, species = c(
                       "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                       "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                     ),
                     seq.style = c("UCSC", "NCBI", "Ensembl", "None"), gtf.file = NULL, up.dist = 3000, down.dist = 3000, ...) {
  # check parameters
  species <- match.arg(arg = species)
  seq.style <- match.arg(arg = seq.style)

  # prepare orgdb and txdb
  spe.anno <- GetSpeciesAnno(species)
  org.db <- spe.anno[["OrgDb"]]
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
  # library orgdb
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org.db: ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # process sequence level
  if (seq.style != "None") {
    GenomeInfoDb::seqlevelsStyle(txdb.obj) <- seq.style
    peak.seqs <- unique(as.character(peak.df$chr))
    # filter peak
    txdb.seqs <- GenomeInfoDb::seqlevels(txdb.obj)
    common.seqs <- intersect(peak.seqs, txdb.seqs)
    peak.df <- peak.df[peak.df$chr %in% common.seqs, ]
    # to avoid Unable to find an inherited method for function 'NSBS' for signature '"SortedByQueryHits"'
    txdb.obj <- keepSeqlevels(x = txdb.obj, value = common.seqs)
  }
  # preare peak information used
  # notice that the output file of MSPC is in bed format
  peak.gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = peak.df, keep.extra.columns = TRUE, ignore.strand = TRUE,
    seqnames.field = "chr", start.field = "start", end.field = "stop",
    starts.in.df.are.0based = TRUE
  )
  peak.anno <- ChIPseeker::annotatePeak(peak.gr,
    tssRegion = c(-up.dist, down.dist),
    TxDb = txdb.obj, annoDb = org.db, ...
  )

  # anno dataframe
  anno.df <- data.frame(peak.anno@anno)
  # modify annotaion
  anno.df$anno <- gsub(pattern = "[[:space:]]+\\(.*\\)", replacement = "", x = anno.df$annotation)

  # plot
  anno.pie <- ggplotify::as.ggplot(function() plotAnnoPie(peak.anno))
  anno.upset <- upsetplot(peak.anno)
  anno.distribution <- plotDistToTSS(peak.anno, title = "Distribution of binding loci relative to TSS")
  # merge plots
  anno.plots <- cowplot::plot_grid(cowplot::plot_grid(anno.pie, anno.distribution, ncol = 2),
    anno.upset,
    ncol = 1
  )

  # get return results
  anno.res <- list(
    anno.obj = peak.anno, df = anno.df, plots = anno.plots,
    pie = anno.pie, upset = anno.upset, distribution = anno.distribution
  )
  return(anno.res)
}
