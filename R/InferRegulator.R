#' Identify TFs from Gene Expression Reults with ChEA3.
#'
#' @param genes Input gene symbol of human or mouse, can be a set of differentially expressed genes.
#' @param species The species of input gene set. Choose from "Human" and "Mouse". Default: "Human".
#' @param projet.name The name of the project used as Query Name. Default: NULL (DEbPeak).
#' @param library The output library type, choose from "meanRank", "topRank", "GTEx", "ReMap",
#' "Enrichr", "ENCODE", "ARCHS4", "Literature". Default: "meanRank".
#'
#' @return Dataframe contains identified TFs.
#' @export
#' @importFrom BiocManager install
#' @importFrom AnnotationDbi keytypes keys
#' @importFrom httr POST content
#' @importFrom jsonlite fromJSON
#'
#' @examples
#' # library(DEbPeak)
#' # RunChEA3(genes = c("SMAD9","FOXO1","MYC","STAT1",'STAT3',"SMAD3"), species = "Human")
RunChEA3 <- function(genes, species = c("Human", "Mouse"), projet.name = NULL,
                     library = c(
                       "meanRank", "topRank", "GTEx", "ReMap",
                       "Enrichr", "ENCODE", "ARCHS4", "Literature"
                     )) {
  # check parameters
  species <- match.arg(arg = species)
  library <- match.arg(arg = library)

  # cehck gene type
  # get org.db
  spe.anno <- GetSpeciesAnno(species)
  org.db <- spe.anno[["OrgDb"]]
  # library orgdb
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org.db: ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))
  gene.type <- CheckGeneName(as.character(genes), org.db)
  if (gene.type != "SYMBOL") {
    stop("ChEA3 only support gene symbol as input!")
  }
  # prepare project name
  if (is.null(projet.name)) {
    projet.name <- "DEbPeak"
  }
  payload <- list(query_name = projet.name, gene_set = genes)
  # POST to ChEA3 server
  response <- httr::POST(
    url = "https://maayanlab.cloud/chea3/api/enrich/",
    body = payload, encode = "json"
  )
  json <- httr::content(response, "text")
  # results as list of R dataframes
  results <- jsonlite::fromJSON(json)
  # named vector
  result.vec <- c(
    "Integrated--meanRank", "Integrated--topRank", "GTEx--Coexpression", "ReMap--ChIP-seq",
    "Enrichr--Queries", "ENCODE--ChIP-seq", "ARCHS4--Coexpression", "Literature--ChIP-seq"
  )
  names(result.vec) <- c(
    "meanRank", "topRank", "GTEx", "ReMap",
    "Enrichr", "ENCODE", "ARCHS4", "Literature"
  )
  # select results
  select.result <- results[[result.vec[library]]]
  return(select.result)
}

#' Identify TFs from Gene Expression Reults with BART2.
#'
#' @param genes Input gene symbol of human or mouse, can be a set of differentially expressed genes.
#' @param species The species of input gene set. Choose from "Human" and "Mouse". Default: "Human".
#' @param bart2.path BART2 path. Default: NULL (conduct automatic detection).
#' @param projet.name The name of the project used as Query Name. Default: NULL (DEbPeak).
#'
#' @return Dataframe contains identified TFs.
#' @export
#' @importFrom BiocManager install
#' @importFrom AnnotationDbi keytypes keys
#' @importFrom utils read.table
#'
#' @examples
#' # RunBART2(genes = c("SMAD9","FOXO1","MYC","STAT1",'STAT3',"SMAD3"), species = "Human",
#' #          bart2.path = "~/anaconda3/envs/bart2/bin/bart2")
RunBART2 <- function(genes, species = c("Human", "Mouse"), bart2.path = NULL, projet.name = NULL) {
  # check parameters
  species <- match.arg(arg = species)
  # cehck gene type
  # get org.db
  spe.anno <- GetSpeciesAnno(species)
  org.db <- spe.anno[["OrgDb"]]
  # library orgdb
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org.db: ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))
  gene.type <- CheckGeneName(as.character(genes), org.db)
  if (gene.type != "SYMBOL") {
    stop("BART2 only support gene symbol as input!")
  }
  # prepare tmp folder
  tmp.folder <- tempdir()
  # write genes to tmp file
  bart2.gene.file <- tempfile(pattern = "bart2_genes_", tmpdir = tmp.folder, fileext = ".txt")
  write.table(x = genes, file = bart2.gene.file, quote = FALSE, col.names = FALSE, row.names = FALSE)
  # get bart2 path
  if (is.null(bart2.path)) {
    # detect bart2 path
    bart2.path <- Sys.which("bart2")
    if (bart2.path == "") {
      stop("Can not find bart2 automatically, please specify the path!")
    }
  } else {
    bart2.path <- bart2.path
  }
  # get species
  spe.anno <- ifelse(species == "Human", "hg38", "mm10")
  # prepare project name
  if (is.null(projet.name)) {
    projet.name <- "DEbPeak"
  }
  # prepare BART2 cmd
  out.folder <- file.path(out.folder, "bart2")
  bart2.cmd <- paste(
    bart2.path, "geneset", "-i", bart2.gene.file, "-s", spe.anno,
    "--outdir", out.folder, "-o", projet.name
  )
  # run command
  message(paste("Running BART2: ", bart2.cmd))
  bart2.status <- system(bart2.cmd, intern = TRUE)
  bart2.status.code <- attr(bart2.status, "status")
  if (!is.null(bart2.status.code)) {
    stop("Run BART2 error!")
  }
  # obtain results
  bart2.results.file <- file.path(out.folder, paste0(projet.name, "_bart_results.txt"))
  if (!file.exists(bart2.results.file)) {
    out.base <- basename(out.folder)
    all.tmp.dirs <- sort(dir(path = dirname(out.folder), pattern = out.base, full.names = TRUE))
    out.folder <- all.tmp.dirs[length(all.tmp.dirs)]
  }
  bart2.results.df <- read.table(file = bart2.results.file, sep = "\t", header = TRUE)
  return(bart2.results.df)
}
