#' Integrate Two Differential Expression Results.
#'
#' @param de1.res DE1 dataframe contains all genes of differential expression analysis.
#' @param de2.res DE2 dataframe contains all genes of differential expression analysis.
#' @param de1.signif Significance criterion for DE1. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param de1.signif.threshold Significance threshold for DE1 to get differentially expressed genes. Default: 0.05.
#' @param de1.l2fc.threshold Log2 fold change threshold for DE1 to get differentially expressed genes. Default: 1.
#' @param de2.signif Significance criterion for DE2. For DESeq2 results, can be chosen from padj, pvalue.
#' For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param de2.signif.threshold Significance threshold for DE2 to get differentially expressed genes. Default: 0.05.
#' @param de2.l2fc.threshold Log2 fold change threshold for DE2 to get differentially expressed genes. Default: 1.
#'
#' @return Dataframe contains integration results, the 'Type' column contains "Down_Up", "Up_Up", "Down_Down", "Up_Down",
#' "DE1_Up", "DE1_Down", "DE2_Up", "DE2_Down", "Not_Not".
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at vars filter
#' @importFrom purrr set_names
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' library(DEbPeak)
#' rna.diff.file <- system.file("extdata", "RA_RNA_diff.txt", package = "DEbPeak")
#' de1.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' de2.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' # use same file as example
#' de.de <- DEbDE(de1.res = de1.res, de2.res = de2.res, de1.l2fc.threshold = 0.5, de2.l2fc.threshold = 1)
DEbDE <- function(de1.res, de2.res,
                  de1.signif = "padj", de1.signif.threshold = 0.05, de1.l2fc.threshold = 1,
                  de2.signif = "padj", de2.signif.threshold = 0.05, de2.l2fc.threshold = 1) {
  # prepare de1 dataframe
  de1.df <- PrepareDEPlot(
    deres = de1.res, signif = de1.signif, signif.threshold = de1.signif.threshold,
    l2fc.threshold = de1.l2fc.threshold, label.key = NULL
  )
  # extract de results
  de1.deg.df <- de1.df %>% dplyr::filter(regulation != "Not_regulated")
  # modify rownames
  colnames(de1.deg.df) <- gsub(pattern = "^", replacement = "DE1_", x = colnames(de1.deg.df))

  # prepare de2 dataframe
  de2.df <- PrepareDEPlot(
    deres = de2.res, signif = de2.signif, signif.threshold = de2.signif.threshold,
    l2fc.threshold = de2.l2fc.threshold, label.key = NULL
  )
  # extract de results
  de2.deg.df <- de2.df %>% dplyr::filter(regulation != "Not_regulated")
  # modify rownames
  colnames(de2.deg.df) <- gsub(pattern = "^", replacement = "DE2_", x = colnames(de2.deg.df))
  # merge
  de.de <- merge(de1.deg.df, de2.deg.df, by.x = "DE1_Gene", by.y = "DE2_Gene", all = TRUE)
  # anno the results
  de.de <- de.de %>% dplyr::mutate(Type = dplyr::case_when(
    DE1_regulation == "Up_regulated" & DE2_regulation == "Up_regulated" ~ "Up_Up",
    DE1_regulation == "Up_regulated" & DE2_regulation == "Down_regulated" ~ "Up_Down",
    DE1_regulation == "Up_regulated" & is.na(DE2_regulation) ~ "DE1_Up",
    DE1_regulation == "Down_regulated" & DE2_regulation == "Up_regulated" ~ "Down_Up",
    DE1_regulation == "Down_regulated" & DE2_regulation == "Down_regulated" ~ "Down_Down",
    DE1_regulation == "Down_regulated" & is.na(DE2_regulation) ~ "DE1_Down",
    is.na(DE1_regulation) & DE2_regulation == "Up_regulated" ~ "DE2_Up",
    is.na(DE1_regulation) & DE2_regulation == "Down_regulated" ~ "DE2_Down"
  ))
  de.de$Type <- factor(de.de$Type, levels = c(
    "Down_Up", "Up_Up", "Down_Down", "Up_Down",
    "DE1_Up", "DE1_Down", "DE2_Up", "DE2_Down"
  ))
  return(de.de)
}

#' GO Enrichment on Two Differential Expression Integration Results.
#'
#' @param de.de Dataframe contains integrated results of two differential expression analysis of RNA-seq.
#' @param de.fe.key The key type of integrated results ("Type" column of \code{de.de}) to perform functional enrichment.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param gene.type Gene name type of \code{DE1_Gene} column. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param go.type GO enrichment type, chosen from ALL, BP, MF, CC. Default: ALL.
#' @param enrich.pvalue Cutoff value of pvalue. Default: 0.05.
#' @param enrich.qvalue Cutoff value of qvalue. Default: 0.05.
#' @param species Species used, chosen from "Human","Mouse","Rat","Fly","Arabidopsis","Yeast","Zebrafish","Worm","Bovine","Pig","Chicken","Rhesus",
#' "Canine","Xenopus","Anopheles","Chimp","E coli strain Sakai","Myxococcus xanthus DK 1622". Default: "Human".
#' @param padj.method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.
#' @param show.term Number of enrichment term to show. Default: 15.
#' @param str.width Length of enrichment term in plot. Default: 30.
#' @param plot.resolution Resolution of plot. Default: 300.
#' @param plot.width The width of plot. Default: 7.
#' @param plot.height The height of plot. Default: 9.
#' @param save Logical value, whether to save all results. Default: TRUE.
#'
#' @return If \code{save} is TRUE, return NULL (all results are in \code{out.folder}), else retutn result dataframe.
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom dplyr filter group_by top_n mutate arrange desc select pull
#' @importFrom magrittr %>%
#' @import clusterProfiler
#' @importFrom enrichplot dotplot
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom BiocManager install
#' @export
#'
#' @examples
#' library(DEbPeak)
#' rna.diff.file <- system.file("extdata", "RA_RNA_diff.txt", package = "DEbPeak")
#' de1.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' de2.res <- read.table(file = rna.diff.file, header = TRUE, sep = "\t")
#' # use same file as example
#' de.de <- DEbDE(de1.res = de1.res, de2.res = de2.res, de1.l2fc.threshold = 0.5, de2.l2fc.threshold = 1)
#' de.de.fe <- DEbDEFE(de.de = de.de, de.fe.key = "Down_Down", gene.type = "SYMBOL", go.type = "BP", species = "Mouse", save = F)
DEbDEFE <- function(de.de, de.fe.key, out.folder = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                    go.type = c("ALL", "BP", "MF", "CC"), enrich.pvalue = 0.05, enrich.qvalue = 0.05, species = c(
                      "Human", "Mouse", "Rat", "Fly", "Arabidopsis", "Yeast", "Zebrafish", "Worm", "Bovine", "Pig", "Chicken", "Rhesus",
                      "Canine", "Xenopus", "Anopheles", "Chimp", "E coli strain Sakai", "Myxococcus xanthus DK 1622"
                    ),
                    padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                    show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  gene.type <- match.arg(arg = gene.type)
  go.type <- match.arg(arg = go.type)
  species <- match.arg(arg = species)
  padj.method <- match.arg(arg = padj.method)

  # get genes
  if (!de.fe.key %in% as.character(unique(de.de$Type))) {
    stop(paste0(
      "Please provide valid functional enrichment key, choose from: ",
      paste(as.character(unique(de.de$Type)), collapse = ", ")
    ))
  }
  de.genes <- de.de %>%
    dplyr::filter(Type == de.fe.key) %>%
    dplyr::pull(DE1_Gene)

  # prepare org db
  spe.anno <- GetSpeciesAnno(species)
  org.db <- spe.anno[["OrgDb"]]
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org_db : ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }

  # GO enrichment
  if (save) {
    SingleFE(
      genes = de.genes, out.folder = out.folder, regulation = de.fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(NULL)
  } else {
    debde.go.results <- SingleFE(
      genes = de.genes, out.folder = out.folder, regulation = de.fe.key, gene.type = gene.type, enrich.type = "GO", go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, padj.method = padj.method,
      show.term = show.term, str.width = str.width, save = save
    )
    return(debde.go.results)
  }
}
