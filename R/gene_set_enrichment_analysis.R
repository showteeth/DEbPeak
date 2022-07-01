#' Conduct Gene Set Enrichment Analysis (GSEA).
#'
#' @param deres Data frame contains all genes.
#' @param gmt.file Gene Matrix Transposed file format.
#' @param gene.sets Gene sets information, containing two columns: gs_name, entrez_gene. Default: NULL.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param gene.key Column name in \code{deres} to conduct analysis. Default: NULL (use rownames of \code{deres}).
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param org.db Organism database. Default: org.Mm.eg.db.
#' @param minGSSize Minimal size of each geneSet for analyzing. Default: 10.
#' @param maxGSSize Maximal size of genes annotated for testing. Default: 500.
#' @param pvalue Cutoff value of pvalue. Default: 0.05.
#' @param padj.method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.
#' @param plot.resolution Resolution of plot. Default: 300.
#' @param plot.width The width of plot. Default: 7.
#' @param plot.height The height of plot. Default: 9.
#' @param save Logical value, whether to save results. Default: TRUE.
#' @param ... Parameters for \code{\link{GSEA}}.
#'
#' @return NULL or list contains all results (\code{GSEA} is FALSE)
#' @importFrom tibble rownames_to_column deframe
#' @importFrom dplyr select
#' @importFrom purrr set_names
#' @importFrom tidyr drop_na
#' @import clusterProfiler
#' @importFrom enrichplot gseaplot2
#' @export
#'
#' @examples
#' library(airway)
#' library(msigdbr)
#' library(DEbChIP)
#' dds <- DESeqDataSet(airway, design = ~ cell + dex)
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("dex", "trt", "untrt"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' h_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::select(gs_name, entrez_gene)
#' gsea.results <- ConductGSEA(deres = dds.results.ordered, gmt.file = NULL, gene.sets = h_t2g, org.db = "org.Hs.eg.db", pvalue = 0.05, save = F)
ConductGSEA <- function(deres, gmt.file, gene.sets = NULL, out.folder = NULL, gene.key = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), org.db = "org.Mm.eg.db",
                        minGSSize = 10, maxGSSize = 500, pvalue = 0.05, padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                        plot.resolution = 300, plot.width = 10, plot.height = 6, save = T, ...) {
  # check parameter
  gene.type <- match.arg(arg = gene.type)
  padj.method <- match.arg(arg = padj.method)
  # get gene and log2foldchange info
  deres <- as.data.frame(deres)
  # which column to use
  if (is.null(gene.key)) {
    gene.key <- "Gene"
  }
  all_results <- list()
  if ("FDR" %in% colnames(deres)) {
    message("Differential expression analysis with edgeR!")
    # modify edgeR results
    de.df <- deres %>%
      tibble::rownames_to_column(var = "Gene") %>%
      dplyr::select(c(`gene.key`, "logFC")) %>%
      purrr::set_names(c(gene.key, "log2FoldChange"))
  } else if ("padj" %in% colnames(deres)) {
    message("Differential expression analysis with DESeq2!")
    de.df <- deres %>%
      tibble::rownames_to_column(var = "Gene") %>%
      dplyr::select(c(`gene.key`, "log2FoldChange")) %>%
      purrr::set_names(c(gene.key, "log2FoldChange"))
  }
  # prepare org db
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org_db : ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # map gene to ENTREZID
  de.df[, gene.key] <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = de.df[, gene.key])
  if (gene.type != "ENTREZID") {
    message("Convert ", gene.type, " to ENTREZID!")
    convert.df <- clusterProfiler::bitr(de.df[, gene.key],
      fromType = gene.type,
      toType = c("ENTREZID"),
      OrgDb = get(org.db), drop = F
    )
    de.df <- merge(de.df, convert.df, by.x = gene.key, by.y = gene.type, allx.x = T) %>%
      tidyr::drop_na() %>%
      dplyr::select(c("ENTREZID", "log2FoldChange")) %>%
      dplyr::distinct(ENTREZID, .keep_all = TRUE)
  } else {
    de.df <- de.df %>% tidyr::drop_na()
  }
  # deframe and sort
  gene.list <- de.df %>% tibble::deframe()
  gene.list.sort <- sort(gene.list, decreasing = TRUE)

  if (is.null(gene.sets)) {
    # read gmt info
    gmt.info <- clusterProfiler::read.gmt(gmt.file)
  } else {
    gmt.info <- gene.sets
  }

  # conduct GSEA analysis
  message("conduct GSEA anaysis.")
  gsea <- clusterProfiler::GSEA(gene.list.sort,
    TERM2GENE = gmt.info, minGSSize = minGSSize,
    maxGSSize = maxGSSize, pvalueCutoff = pvalue,
    pAdjustMethod = padj.method, ...
  )
  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }

  gsea.info <- as.data.frame(gsea)
  if (nrow(na.omit(gsea.info)) >= 1) {
    p1 <- enrichplot::gseaplot2(gsea, geneSetID = 1, title = gsea$Description[1])
  } else {
    p1 <- ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "N/A")) +
      xlab(NULL)
  }
  if (save) {
    # save results
    gsea.out.file <- file.path(out.folder, "GSEA_enrich_result.csv")
    write.csv(gsea.info, file = gsea.out.file, row.names = F)

    ggsave(file.path(out.folder, "GSEA_enrich_result.png"),
      plot = p1, dpi = plot.resolution,
      width = plot.width, height = plot.height
    )
    pdf(file.path(out.folder, "GSEA_enrich_result.pdf"))
    print(p1)
    dev.off()
  } else {
    all_results[["gsea"]] <- gsea
    all_results[["table"]] <- gsea.info
    all_results[["plot"]] <- p1
  }
  return(all_results)
}
