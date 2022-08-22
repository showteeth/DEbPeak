#' Find Motif of Peaks from DEGs and Peak Annotation Integration Results.
#'
#' @param inte.res DEGs and peak annotation integration results.
#' @param peak.anno.res Peak annotation results.
#' @param gene.key Gene key name, chosen from geneId, ENSEMBL, SYMBOL. Default: geneId.
#' Should be consistent with \code{merge.key} of \code{DEbChIP}.
#' @param homer.motif.path The path to Homer's \code{findMotifsGenome.pl}. Default: NULL.
#' @param genome Parameter for \code{findMotifsGenome.pl}, can be genome FASTA files or pre-build genome info. Default: mm10.
#' @param out.folder The output folder. Default: NULL (current directory).
#' @param other.paras Parameter for \code{findMotifsGenome.pl}, can be '-len 8,10,12 -size -100,50 -S 25'. Default: NULL.
#'
#' @return List contains known motif results of UPbATAC and DOWNbATAC.
#' @importFrom data.table fread
#' @export
#'
#' @examples
FindMotif=function(inte.res, peak.anno.res, gene.key = "geneId", homer.motif.path = NULL, genome='mm10',
                   out.folder=NULL, other.paras=NULL){

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
  export.path = paste0('export PATH="$PATH:', dirname(homer.motif.path),'";')

  # set output folder
  if (is.null(out.folder)) {
    out.folder <- file.path(getwd(),"Motif")
  }

  # prepare peak dataframe
  # get genes
  upbchip.genes = inte.res[inte.res$Type=="UPbChIP", gene.key]
  downbchip.genes = inte.res[inte.res$Type=="DOWNbChIP", gene.key]

  # extract peaks
  ## up bind peak
  if(length(upbchip.genes)>=1){
    # get peak
    upbchip.peak = peak.anno.res[peak.anno.res[[gene.key]] %in% upbchip.genes,
                                 c("seqnames", "start", "end", "name", "score")]
    upbchip.peak$strand = '+'
    # write bed
    up.bed.file = tempfile(pattern = "UP",fileext = ".bed")
    write.table(x = upbchip.peak, file = up.bed.file, quote = FALSE, sep="\t",
                col.names = FALSE, row.names = FALSE)
    # command used
    up.out.folder = file.path(out.folder, "UPbATAC")
    up.homer.motif.cmd <- paste(homer.motif.path, up.bed.file, genome, up.out.folder, other.paras)
    # add path
    up.homer.motif.cmd <- paste0(export.path, up.homer.motif.cmd)

    # run command
    message(paste("Calling findMotifsGenome.pl on UPbATAC: ", up.homer.motif.cmd))
    up.homer.motif.status <- system(up.homer.motif.cmd, intern = TRUE)
    up.homer.motif.status.code <- attr(up.homer.motif.status, "status")
    if (!is.null(up.homer.motif.status.code)) {
      stop("Run findMotifsGenome.pl error!")
    }
    # read up known motif
    up.known.motif = data.table::fread(file = file.path(up.out.folder,"knownResults.txt"),sep = "\t")
  }else{
    up.known.motif = "empty"
  }
  ## down bind peak
  if(length(downbchip.genes)>=1){
    downbchip.peak = peak.anno.res[peak.anno.res[[gene.key]] %in% downbchip.genes,
                                   c("seqnames", "start", "end", "name", "score")]
    downbchip.peak$strand = '+'
    # write bed
    down.bed.file = tempfile(pattern = "DOWN",fileext = ".bed")
    write.table(x = downbchip.peak, file = down.bed.file, quote = FALSE, sep="\t",
                col.names = FALSE, row.names = FALSE)
    # command used
    down.out.folder = file.path(out.folder, "DOWNbATAC")
    down.homer.motif.cmd <- paste(homer.motif.path, down.bed.file, genome, down.out.folder, other.paras)
    # add path
    down.homer.motif.cmd <- paste0(export.path, down.homer.motif.cmd)

    # run command
    message(paste("Calling findMotifsGenome.pl on DOWNbATAC: ", down.homer.motif.cmd))
    down.homer.motif.status <- system(down.homer.motif.cmd, intern = TRUE)
    down.homer.motif.status.code <- attr(down.homer.motif.status, "status")
    if (!is.null(down.homer.motif.status.code)) {
      stop("Run findMotifsGenome.pl error!")
    }
    # read down known motif
    down.known.motif = data.table::fread(file = file.path(down.out.folder,"knownResults.txt"),sep = "\t")
  }else{
    down.known.motif = "empty"
  }
  # final results
  final.res = list(UPbATAC=up.known.motif, DOWNbATAC=down.known.motif)
  return(final.res)
}


