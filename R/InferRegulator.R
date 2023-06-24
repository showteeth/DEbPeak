#' Identify TFs from Gene Expression Reults.
#'
#' @param genes Input genes, can be a set of differentially expressed genes. When \code{method} is "ChEA3" or "BART2",
#' the genes should be SYMBOL of human or mouse. When \code{method} is "TFEA.ChIP", the genes can be
#' SYMBOL/ENSEMBL/ENTREZ of human.
#' @param control.genes Input genes of human, can be a set of non-differentially expressed genes.
#' Used when \code{method} is "TFEA.ChIP".
#' @param method The method used to identify TFs, choose from "ChEA3", "BART2", "TFEA.ChIP". Default: "ChEA3".
#' @param species The species of input gene set. Choose from "Human" and "Mouse" when \code{method} is "ChEA3" or "BART2",
#' "Human" when \code{method} is "TFEA.ChIP". Default: "Human".
#' @param projet.name The name of the project used as Query Name. Default: NULL (DEbPeak).
#' @param library The output library type, choose from "meanRank", "topRank", "GTEx", "ReMap",
#' "Enrichr", "ENCODE", "ARCHS4", "Literature". Used when \code{method} is "ChEA3". Default: "meanRank".
#' @param bart2.path BART2 path. Default: NULL (conduct automatic detection).
#' @param filter.threshold The threshold used to filter the results. When \code{method} is "ChEA3", filter on Score column.
#' When \code{method} is "BART2", filter on pvalue column. When \code{method} is "TFEA.ChIP", filter on pVal column.
#'
#' @return Dataframe contains identified TFs.
#' @export
#' @importFrom BiocManager install
#' @importFrom AnnotationDbi keytypes keys
#' @importFrom httr POST content
#' @importFrom jsonlite fromJSON
#' @importFrom TFEA.ChIP GeneID2entrez contingency_matrix getCMstats rankTFs
#' @importFrom dplyr arrange desc
#'
#' @examples
#' # library(DEbPeak)
#' # InferRegulator(genes = c("SMAD9","FOXO1","MYC","STAT1",'STAT3',"SMAD3"), method = "ChEA3", species = "Human")
InferRegulator <- function(genes, control.genes, method = c("ChEA3", "BART2", "TFEA.ChIP"),
                           species = c("Human", "Mouse"), projet.name = NULL,
                           library = c(
                             "meanRank", "topRank", "GTEx", "ReMap",
                             "Enrichr", "ENCODE", "ARCHS4", "Literature"
                           ), bart2.path = NULL, filter.threshold = NULL) {
  # check parameters
  method <- match.arg(arg = method)
  species <- match.arg(arg = species)
  library <- match.arg(arg = library)

  # obtain results
  if (method == "ChEA3") {
    infer.res <- RunChEA3(genes = genes, species = species, projet.name = projet.name, library = library)
    # filter the results
    if (!is.null(filter.threshold)) {
      infer.res <- infer.res[infer.res$Score < filter.threshold, ]
    }
  } else if (method == "BART2") {
    infer.res <- RunBART2(genes = genes, species = species, projet.name = projet.name, bart2.path = bart2.path)
    # filter the results
    if (!is.null(filter.threshold)) {
      infer.res <- infer.res[infer.res$pvalue < filter.threshold, ]
    }
  } else if (method == "TFEA.ChIP") {
    if (species == "Mouse") {
      stop("TFEA.ChIP only support human! You can convert mouse genes to their human orthologs to perform this analysis!")
    }
    infer.res <- RunTFEA(genes = genes, control.genes = control.genes)
    # filter the results
    if (!is.null(filter.threshold)) {
      infer.res <- infer.res[infer.res$pVal < filter.threshold, ]
    }
  }

  return(infer.res)
}


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
#' @importFrom utils read.table
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
#' # library(DEbPeak)
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

#' Identify TFs from Gene Expression Reults with TFEA.ChIP.
#'
#' @param genes Input genes of human, can be a set of differentially expressed genes.
#' @param control.genes Input genes of human, can be a set of non-differentially expressed genes.
#'
#' @return Dataframe contains identified TFs with GSEA enrichment score.
#' @importFrom BiocManager install
#' @importFrom AnnotationDbi keytypes keys
#' @importFrom TFEA.ChIP GeneID2entrez contingency_matrix getCMstats rankTFs
#' @importFrom dplyr arrange desc
#' @export
#'
#' @examples
#' # library(DEbPeak)
#' # genes = c("EGLN3","NFYA","ALS2","MYC","ARNT" )
#' # control.genes = c("MTX1", "OSTF1", "HDAC3", "MTG2", "OPA3")
#' # RunTFEA(genes = genes, control.genes = control.genes)
RunTFEA <- function(genes, control.genes) {
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
  # input gene type
  gene.type <- CheckGeneName(as.character(genes), org.db)
  # control gene type
  control.gene.type <- CheckGeneName(as.character(control.genes), org.db)
  # gene id conversion
  if (gene.type != "ENTREZID") {
    message("Convert ", gene.type, " to ENTREZID!")
    genes.id <- TFEA.ChIP::GeneID2entrez(gene.IDs = genes, mode = "h2h")
  } else {
    genes.id <- genes
  }
  if (control.gene.type != "ENTREZID") {
    message("Convert ", control.gene.type, " to ENTREZID!")
    control.genes.id <- TFEA.ChIP::GeneID2entrez(gene.IDs = control.genes, mode = "h2h")
  } else {
    control.genes.id <- control.genes
  }
  # analysis
  contingency.matrix <- TFEA.ChIP::contingency_matrix(test_list = genes.id, control_list = control.genes.id)
  contingency.matrix.stats <- TFEA.ChIP::getCMstats(contMatrix_list = contingency.matrix)
  # rank TF
  TF.ranking <- TFEA.ChIP::rankTFs(resultsTable = contingency.matrix.stats, rankMethod = "gsea")
  # sort by ES
  TF.ranking <- TF.ranking %>% dplyr::arrange(desc(ES), arg.ES)
  return(TF.ranking)
}

#' Visualize the Identified TFs.
#'
#' @param infer.res The identidied result dataframe, can be the output of \code{\link{InferRegulator}}.
#' @param method The method used to identify TFs, choose from "ChEA3", "BART2", "TFEA.ChIP". Default: "ChEA3".
#' @param label.gene The highlight genes. Default: NULL.
#' @param point.alpha The point alpha of identify TFs. Default: 0.5.
#' @param label.color The text/point color for the highlight genes. Default: "red".
#'
#' @return A ggplot2 object.
#' @export
#' @importFrom ggplot2 aes_string ggplot geom_point theme_classic theme element_rect
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' # library(DEbPeak)
#' # infer.res = InferRegulator(genes = c("SMAD9","FOXO1","MYC","STAT1",'STAT3',"SMAD3"),
#' #                            method = "ChEA3", species = "Human")
#' # infer.plot = VizRegulator(infer.res = infer.res, method = "ChEA3")
VizRegulator <- function(infer.res, method = c("ChEA3", "BART2", "TFEA.ChIP"), label.gene = NULL,
                         point.alpha = 0.5, label.color = "red") {
  # check parameters
  method <- match.arg(arg = method)
  # prepare data and mapping
  if (method == "ChEA3") {
    mapping <- aes_string(x = "Rank", y = "Score")
    label.mapping <- aes_string(x = "Rank", y = "Score", label = "TF")
  } else if (method == "BART2") {
    infer.res$Rank <- 1:nrow(infer.res)
    mapping <- aes_string(x = "Rank", y = "max_auc")
    label.mapping <- aes_string(x = "Rank", y = "max_auc", label = "TF")
  } else if (method == "TFEA.ChIP") {
    infer.res$Rank <- 1:nrow(infer.res)
    mapping <- aes_string(x = "Rank", y = "ES")
    label.mapping <- aes_string(x = "Rank", y = "ES", label = "TF")
  }
  # create basic plot
  basic.plot <- ggplot() +
    geom_point(data = infer.res, mapping = mapping, alpha = point.alpha) +
    theme_classic(base_size = 14) +
    theme(panel.background = element_rect(color = "black"))
  # label gene
  if (!is.null(label.gene)) {
    label.data <- infer.res[infer.res$TF %in% label.gene, ]
    final.plot <- basic.plot +
      geom_point(data = label.data, mapping = mapping, color = label.color) +
      geom_text_repel(
        data = label.data,
        color = label.color,
        mapping = label.mapping,
        box.padding = 1,
        nudge_y = 1, angle = 90, direction = "y"
      )
  } else {
    final.plot <- basic.plot
  }
  return(final.plot)
}
