# function for GO analysis
GO_func <- function(pvalue, qvalue, GO_type, entrez_id, org_db, padj.method) {
  erich.go <- clusterProfiler::enrichGO(
    gene = entrez_id,
    OrgDb = org_db,
    ont = GO_type,
    keyType = "ENTREZID",
    pvalueCutoff = pvalue,
    qvalueCutoff = qvalue,
    pAdjustMethod = padj.method,
    readable = TRUE
  )
  return(erich.go)
}
# function for KEGG analysis
KEGG_func <- function(pvalue, qvalue, entrez_id, organism, org_db, padj.method) {
  enrich.kegg <- clusterProfiler::enrichKEGG(
    gene = entrez_id,
    organism = organism,
    keyType = "kegg",
    pvalueCutoff = pvalue,
    qvalueCutoff = qvalue,
    pAdjustMethod = padj.method
  )
  if (!is.null(enrich.kegg)) {
    enrich.kegg <- setReadable(enrich.kegg, OrgDb = org_db, keyType = "ENTREZID")
  }
  return(enrich.kegg)
}

SingleFE <- function(genes, out.folder, regulation, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), enrich.type = c("ALL", "GO", "KEGG"), go.type = c("ALL", "BP", "MF", "CC"),
                     enrich.pvalue = 0.05, enrich.qvalue = 0.05, org.db = "org.Mm.eg.db", organism = "mmu",
                     padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                     show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # convert gene types
  # remove gene version information
  genes <- gsub(pattern = "\\.[0-9]*$", replacement = "", x = genes)
  # convert gene ids
  if (gene.type != "ENTREZID") {
    message("Convert ", gene.type, " to ENTREZID!")
    convert.df <- clusterProfiler::bitr(genes,
      fromType = gene.type,
      toType = c("ENTREZID"),
      OrgDb = get(org.db)
    )
    entrez.id <- convert.df[, "ENTREZID"]
  } else {
    entrez.id <- genes
  }
  all.results <- list()
  # conduct enrichment analysis
  if (enrich.type == "ALL") {
    if (go.type == "ALL") {
      message("conduct ALL GO enrichment analysis on: ", regulation)
      GO_enrich <- GO_func(
        pvalue = enrich.pvalue, qvalue = enrich.qvalue, GO_type = "ALL",
        entrez_id = entrez.id, org_db = org.db, padj.method = padj.method
      )
      BP_enrich <- clusterProfiler::filter(GO_enrich, ONTOLOGY == "BP")
      MF_enrich <- clusterProfiler::filter(GO_enrich, ONTOLOGY == "MF")
      CC_enrich <- clusterProfiler::filter(GO_enrich, ONTOLOGY == "CC")
      # write GO enrich results
      GO_enrich_df <- stats::na.omit(as.data.frame(GO_enrich))
      if (save) {
        go_out_file <- file.path(out.folder, paste0(regulation, "_ALL_GO.csv"))
        utils::write.table(GO_enrich_df, file = go_out_file, row.names = FALSE, quote = FALSE, sep = ",")
      } else {
        all.results[["GO"]][["table"]] <- GO_enrich_df
      }
      message("conduct KEGG enrichment analysis: ", regulation)
      KEGG_enrich <- KEGG_func(
        pvalue = enrich.pvalue, qvalue = enrich.qvalue, entrez_id = entrez.id,
        organism = organism, org_db = org.db, padj.method = padj.method
      )
      result_KEGG <- stats::na.omit(as.data.frame(KEGG_enrich))
      if (save) {
        kegg_out_file <- file.path(out.folder, paste0(regulation, "_KEGG.csv"))
        utils::write.table(result_KEGG, file = kegg_out_file, row.names = FALSE, quote = FALSE, sep = ",")
      } else {
        all.results[["KEGG"]][["table"]] <- result_KEGG
      }
      if (nrow(stats::na.omit(as.data.frame(BP_enrich))) >= 1) {
        p1 <- dotplot(BP_enrich, showCategory = show.term, title = paste0("Biological Process of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p1 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (nrow(stats::na.omit(as.data.frame(MF_enrich))) >= 1) {
        p2 <- dotplot(MF_enrich, showCategory = show.term, title = paste0("Molecular Function of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p2 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (nrow(stats::na.omit(as.data.frame(CC_enrich))) >= 1) {
        p3 <- dotplot(CC_enrich, showCategory = show.term, title = paste0("Cellular Component of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p3 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (nrow(stats::na.omit(as.data.frame(KEGG_enrich))) >= 1) {
        p4 <- dotplot(KEGG_enrich, showCategory = show.term, title = paste0("KEGG Enrichment of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p4 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }

      if (save) {
        # save png
        ggsave(paste(regulation, "Biological_Process.png", sep = "_"),
          plot = p1, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        ggsave(paste(regulation, "Molecular_Function.png", sep = "_"),
          plot = p2, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        ggsave(paste(regulation, "Cellular_Component.png", sep = "_"),
          plot = p3, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        ggsave(paste(regulation, "KEGG_Enrichment.png", sep = "_"),
          plot = p4, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        # save pdf file
        grDevices::pdf(paste(regulation, "Functional_Enrichment.pdf", sep = "_"))
        print(p1)
        print(p2)
        print(p3)
        print(p4)
        grDevices::dev.off()
      } else {
        all.results[["GO"]][["plot"]] <- patchwork::wrap_plots(p1, p2, p3, ncol = 1)
        all.results[["KEGG"]][["plot"]] <- p4
      }
    } else {
      message("conduct ", go.type, " GO enrichment analysis.")
      part_go <- GO_func(
        pvalue = enrich.pvalue, qvalue = enrich.qvalue, GO_type = go.type,
        entrez_id = entrez.id, org_db = org.db, padj.method = padj.method
      )
      go_df <- stats::na.omit(as.data.frame(part_go))

      message("conduct KEGG enrichment analysis.")
      KEGG_enrich <- KEGG_func(
        pvalue = enrich.pvalue, qvalue = enrich.qvalue, entrez_id = entrez.id,
        organism = organism, org_db = org.db, padj.method = padj.method
      )
      result_KEGG <- stats::na.omit(as.data.frame(KEGG_enrich))
      if (save) {
        go_out_file <- file.path(out.folder, paste0(regulation, "_", go.type, "_GO.csv"))
        utils::write.table(go_df, file = go_out_file, row.names = FALSE, quote = FALSE, sep = ",")
        kegg_out_file <- file.path(out.folder, paste0(regulation, "_KEGG.csv"))
        utils::write.table(result_KEGG, file = kegg_out_file, row.names = FALSE, quote = FALSE, sep = ",")
      } else {
        all.results[["GO"]][["table"]] <- go_df
        all.results[["KEGG"]][["table"]] <- result_KEGG
      }

      if (go.type == "BP") {
        go_png_file <- paste(regulation, "Biological_Process.png", sep = "_")
        plot.title <- paste0("Biological Process of ", regulation, " Genes")
      } else if (go.type == "MF") {
        go_png_file <- paste(regulation, "Molecular_Function.png", sep = "_")
        plot.title <- paste0("Molecular Function of ", regulation, " Genes")
      } else if (go.type == "CC") {
        go_png_file <- paste(regulation, "Cellular_Component.png", sep = "_")
        plot.title <- paste0("Cellular Component of ", regulation, " Genes")
      }

      if (nrow(stats::na.omit(as.data.frame(part_go))) >= 1) {
        p1 <- dotplot(part_go, showCategory = show.term, title = plot.title) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p1 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (nrow(stats::na.omit(as.data.frame(KEGG_enrich))) >= 1) {
        p2 <- dotplot(KEGG_enrich, showCategory = show.term, title = paste0("KEGG Enrichment of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p2 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }

      if (save) {
        # save png
        ggsave(go_png_file, plot = p1, dpi = plot.resolution, width = plot.width, height = plot.height)
        ggsave(paste(regulation, "KEGG_Enrichment.png", sep = "_"),
          plot = p2, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        # save pdf file
        grDevices::pdf(paste(regulation, "Functional_Enrichment.pdf", sep = "_"))
        print(p1)
        print(p2)
        grDevices::dev.off()
      } else {
        all.results[["GO"]][["plot"]] <- p1
        all.results[["KEGG"]][["table"]] <- p2
      }
    }
  } else if (enrich.type == "GO") {
    if (go.type == "ALL") {
      message("conduct ALL GO enrichment analysis on: ", regulation)
      GO_enrich <- GO_func(
        pvalue = enrich.pvalue, qvalue = enrich.qvalue, GO_type = "ALL",
        entrez_id = entrez.id, org_db = org.db, padj.method = padj.method
      )
      BP_enrich <- clusterProfiler::filter(GO_enrich, ONTOLOGY == "BP")
      MF_enrich <- clusterProfiler::filter(GO_enrich, ONTOLOGY == "MF")
      CC_enrich <- clusterProfiler::filter(GO_enrich, ONTOLOGY == "CC")
      # write GO enrich results
      GO_enrich_df <- stats::na.omit(as.data.frame(GO_enrich))
      if (save) {
        go_out_file <- file.path(out.folder, paste0(regulation, "_ALL_GO.csv"))
        utils::write.table(GO_enrich_df, file = go_out_file, row.names = FALSE, quote = FALSE, sep = ",")
      } else {
        all.results[["GO"]][["table"]] <- GO_enrich_df
      }

      if (nrow(stats::na.omit(as.data.frame(BP_enrich))) >= 1) {
        p1 <- dotplot(BP_enrich, showCategory = show.term, title = paste0("Biological Process of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p1 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (nrow(stats::na.omit(as.data.frame(MF_enrich))) >= 1) {
        p2 <- dotplot(MF_enrich, showCategory = show.term, title = paste0("Molecular Function of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p2 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (nrow(stats::na.omit(as.data.frame(CC_enrich))) >= 1) {
        p3 <- dotplot(CC_enrich, showCategory = show.term, title = paste0("Cellular Component of ", regulation, " Genes")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p3 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }

      if (save) {
        # save png
        ggsave(paste(regulation, "Biological_Process.png", sep = "_"),
          plot = p1, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        ggsave(paste(regulation, "Molecular_Function.png", sep = "_"),
          plot = p2, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        ggsave(paste(regulation, "Cellular_Component.png", sep = "_"),
          plot = p3, dpi = plot.resolution, width = plot.width, height = plot.height
        )
        # save pdf file
        grDevices::pdf(paste(regulation, "Functional_Enrichment.pdf", sep = "_"))
        print(p1)
        print(p2)
        print(p3)
        grDevices::dev.off()
      } else {
        all.results[["GO"]][["plot"]] <- patchwork::wrap_plots(p1, p2, p3, ncol = 1)
      }
    } else {
      message("conduct ", go.type, " GO enrichment analysis.")
      part_go <- GO_func(
        pvalue = enrich.pvalue, qvalue = enrich.qvalue, GO_type = go.type,
        entrez_id = entrez.id, org_db = org.db, padj.method = padj.method
      )
      go_df <- stats::na.omit(as.data.frame(part_go))
      if (save) {
        go_out_file <- file.path(out.folder, paste0(regulation, "_", go.type, "_GO.csv"))
        utils::write.table(go_df, file = go_out_file, row.names = FALSE, quote = FALSE, sep = ",")
      } else {
        all.results[["GO"]][["table"]] <- go_df
      }
      if (go.type == "BP") {
        go_png_file <- paste(regulation, "Biological_Process.png", sep = "_")
        plot.title <- paste0("Biological Process of ", regulation, " Genes")
      } else if (go.type == "MF") {
        go_png_file <- paste(regulation, "Molecular_Function.png", sep = "_")
        plot.title <- paste0("Molecular Function of ", regulation, " Genes")
      } else if (go.type == "CC") {
        go_png_file <- paste(regulation, "Cellular_Component.png", sep = "_")
        plot.title <- paste0("Cellular Component of ", regulation, " Genes")
      }
      if (nrow(stats::na.omit(as.data.frame(part_go))) >= 1) {
        p1 <- dotplot(part_go, showCategory = show.term, title = plot.title) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
      } else {
        p1 <- ggplot() +
          theme_void() +
          geom_text(aes(0, 0, label = "N/A")) +
          xlab(NULL)
      }
      if (save) {
        ggsave(go_png_file, plot = p1, dpi = plot.resolution, width = plot.width, height = plot.height)
        grDevices::pdf(paste(regulation, "Functional_Enrichment.pdf", sep = "_"))
        print(p1)
        grDevices::dev.off()
      } else {
        all.results[["GO"]][["plot"]] <- p1
      }
    }
  } else if (enrich.type == "KEGG") {
    message("conduct KEGG enrichment analysis.")
    KEGG_enrich <- KEGG_func(
      pvalue = enrich.pvalue, qvalue = enrich.qvalue, entrez_id = entrez.id,
      organism = organism, org_db = org.db, padj.method = padj.method
    )
    result_KEGG <- stats::na.omit(as.data.frame(KEGG_enrich))
    if (save) {
      kegg_out_file <- file.path(out.folder, paste0(regulation, "_KEGG.csv"))
      utils::write.table(result_KEGG, file = kegg_out_file, row.names = FALSE, quote = FALSE, sep = ",")
    } else {
      all.results[["KEGG"]][["table"]] <- result_KEGG
    }
    if (nrow(stats::na.omit(as.data.frame(KEGG_enrich))) >= 1) {
      p1 <- dotplot(KEGG_enrich, showCategory = show.term, title = paste0("KEGG Enrichment of ", regulation, " Genes")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
    } else {
      p1 <- ggplot() +
        theme_void() +
        geom_text(aes(0, 0, label = "N/A")) +
        xlab(NULL)
    }
    if (save) {
      ggsave(paste(regulation, "KEGG_Enrichment.png", sep = "_"),
        plot = p1, dpi = plot.resolution, width = plot.width, height = plot.height
      )
      grDevices::pdf(paste(regulation, "Functional_Enrichment.pdf", sep = "_"))
      print(p1)
      grDevices::dev.off()
    } else {
      all.results[["KEGG"]][["plot"]] <- p1
    }
  }
  return(all.results)
}

#' Conduct Functional Enrichment Analysis.
#'
#' @param deres Data frame contains all genes.
#' @param out.folder Folder to save enrichment results. Default: wording directory.
#' @param data.type Input data type, choose from RNA, ChIP, ATAC. Default: RNA.
#' @param peak.anno.key Peak location, chosen from "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic","All". Default: "Promoter".
#' @param signif Significance criterion. For DESeq2 results, can be chosen from padj, pvalue. For edgeR results, can be chosen from FDR, PValue. Default: padj.
#' @param signif.threshold Significance threshold to get differentially expressed genes or accessible/binding peaks. Default: 0.05.
#' @param l2fc.threshold Log2 fold change threshold to get differentially expressed genes or accessible/binding peaks. Default: 1.
#' @param gene.key Column name in \code{deres} to conduct analysis. Default: NULL (use rownames of \code{deres}).
#' @param gene.type Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.
#' @param org.db Organism database. Default: org.Mm.eg.db.
#' @param enrich.type Enrichment type, chosen from ALL, GO, KEGG. Default: ALL.
#' @param go.type GO enrichment type, chosen from ALL, BP, MF, CC. Default: ALL.
#' @param enrich.pvalue Cutoff value of pvalue. Default: 0.05.
#' @param enrich.qvalue Cutoff value of qvalue. Default: 0.05.
#' @param organism Supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'. Default: mmu.
#' @param padj.method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.
#' @param show.term Number of enrichment term to show. Default: 15.
#' @param str.width Length of enrichment term in plot. Default: 30.
#' @param plot.resolution Resolution of plot. Default: 300.
#' @param plot.width The width of plot. Default: 7.
#' @param plot.height The height of plot. Default: 9.
#' @param save Logical value, whether to save all results. Default: TRUE.
#'
#' @return If \code{save} is TRUE, return NULL (all results are in \code{out.folder}), else retutn list contains all results.
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate case_when mutate_at filter
#' @importFrom tidyr drop_na
#' @importFrom purrr set_names
#' @import clusterProfiler
#' @importFrom enrichplot dotplot
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom BiocManager install
#' @importFrom patchwork wrap_plots
#' @importFrom stats na.omit
#' @importFrom utils write.table
#' @importFrom grDevices pdf dev.off
#' @export
#'
#' @examples
#' library(DESeq2)
#' library(DEbPeak)
#' count.file <- system.file("extdata", "snon_count.txt", package = "DEbPeak")
#' meta.file <- system.file("extdata", "snon_meta.txt", package = "DEbPeak")
#' count.matrix <- read.table(file = count.file, header = TRUE, sep = "\t")
#' meta.info <- read.table(file = meta.file, header = TRUE)
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix, colData = meta.info, design = ~condition)
#' keep.genes <- rowSums(DESeq2::counts(dds, normalized = FALSE)) >= 10
#' dds <- dds[keep.genes, ]
#' dds$condition <- relevel(dds$condition, ref = "WT")
#' dds <- DESeq(dds)
#' dds.results <- results(dds, contrast = c("condition", "KO", "WT"))
#' dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
#' ConductFE(deres = dds.results.ordered, signif = "pvalue", l2fc.threshold = 0.3)
ConductFE <- function(deres, out.folder = NULL, data.type = c("RNA", "ChIP", "ATAC"),
                      peak.anno.key = c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "All"),
                      signif = "padj", signif.threshold = 0.05, l2fc.threshold = 1,
                      gene.key = NULL, gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"), org.db = "org.Mm.eg.db",
                      enrich.type = c("ALL", "GO", "KEGG"), go.type = c("ALL", "BP", "MF", "CC"),
                      enrich.pvalue = 0.05, enrich.qvalue = 0.05, organism = "mmu",
                      padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                      show.term = 15, str.width = 30, plot.resolution = 300, plot.width = 7, plot.height = 9, save = TRUE) {
  # check parameter
  data.type <- match.arg(arg = data.type)
  peak.anno.key <- match.arg(arg = peak.anno.key)
  gene.type <- match.arg(arg = gene.type)
  enrich.type <- match.arg(arg = enrich.type)
  go.type <- match.arg(arg = go.type)
  padj.method <- match.arg(arg = padj.method)

  # prepare org db
  if (!require(org.db, quietly = TRUE, character.only = TRUE)) {
    message("Install org_db : ", org.db)
    BiocManager::install(org.db)
  }
  suppressWarnings(suppressMessages(library(org.db, character.only = TRUE)))

  # prepare deres
  if (data.type != "RNA") {
    # filter dataframe based on region
    anno.key.named <- c("P", "5U", "3U", "E", "I", "D", "DI")
    names(anno.key.named) <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
    deres <- as.data.frame(deres)
    deres$Type <- gsub(pattern = ".*\\|.*\\|(.*)", replacement = "\\1", x = rownames(deres))
    if (peak.anno.key == "All") {
      deres <- deres
    } else {
      deres <- deres[deres$Type == anno.key.named[peak.anno.key], ]
    }
    deres <- deres %>%
      dplyr::select(-c(Type))
    deres$SYMBOL <- gsub(pattern = ".*\\|(.*)\\|.*", replacement = "\\1", x = rownames(deres))
    gene.key <- "SYMBOL"
    gene.type <- "SYMBOL"
  }

  # preapare DE dataframe
  de.df <- PrepareDEPlot(deres = deres, signif = signif, signif.threshold = signif.threshold, l2fc.threshold = l2fc.threshold, label.key = gene.key)
  up.de.df <- de.df %>% dplyr::filter(regulation == "Up_regulated")
  down.de.df <- de.df %>% dplyr::filter(regulation == "Down_regulated")
  if (is.null(gene.key)) {
    gene.key <- "Gene"
  }
  up.degs <- as.character(up.de.df[, gene.key])
  down.degs <- as.character(down.de.df[, gene.key])

  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (save) {
    # conduct GO and KEGG analysis for UP regulated genes
    SingleFE(
      genes = up.degs, out.folder = out.folder, regulation = "UP", gene.type = gene.type, enrich.type = enrich.type, go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method,
      show.term = show.term, str.width = str.width, plot.resolution = plot.resolution, plot.width = plot.width, plot.height = plot.height, save = save
    )
    # conduct GO and KEGG analysis for DOWN regulated genes
    SingleFE(
      genes = down.degs, out.folder = out.folder, regulation = "DOWN", gene.type = gene.type, enrich.type = enrich.type, go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method,
      show.term = show.term, str.width = str.width, plot.resolution = plot.resolution, plot.width = plot.width, plot.height = plot.height, save = save
    )
    return(NULL)
  } else {
    go.results <- list()
    # conduct GO and KEGG analysis for UP regulated genes
    up.go.results <- SingleFE(
      genes = up.degs, out.folder = out.folder, regulation = "UP", gene.type = gene.type, enrich.type = enrich.type, go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method,
      show.term = show.term, str.width = str.width, plot.resolution = plot.resolution, plot.width = plot.width, plot.height = plot.height, save = save
    )
    # conduct GO and KEGG analysis for DOWN regulated genes
    down.go.results <- SingleFE(
      genes = down.degs, out.folder = out.folder, regulation = "DOWN", gene.type = gene.type, enrich.type = enrich.type, go.type = go.type,
      enrich.pvalue = enrich.pvalue, enrich.qvalue = enrich.qvalue, org.db = org.db, organism = organism, padj.method = padj.method,
      show.term = show.term, str.width = str.width, plot.resolution = plot.resolution, plot.width = plot.width, plot.height = plot.height, save = save
    )
    go.results[["UP"]] <- up.go.results
    go.results[["DOWN"]] <- down.go.results
    return(go.results)
  }
}

#' Create Plot for Functional Enrichment Analysis.
#'
#' @param fe.res Functional enrichment analysis results to plot, usually selected from raw results.
#' @param signif Significance criterion of the results, choose from p.adjust, pvalue, qvalue,
#' NULL (do not show significance information). Default: p.adjust.
#' @param show.num The number of terms to show, select in descending order of the FoldEnrichment. Default: 10.
#' @param grad.col.vec Two color gradient vector (low, high) for \code{signif}.
#' Used when \code{signif} is not NULL. Default: c("grey", "red").
#' @param fill.color Color used to color the geom when \code{signif} is NULL. Default: "red".
#' @param str.width Length of enrichment term in plot. If NULL, use raw length. Default: 30.
#' @param plot.type The type of the plot, choose from bar and dot. Default: bar.
#' @param plot.resolution Resolution of plot. Default: 300.
#' @param plot.width The width of plot. Default: 7.
#' @param plot.height The height of plot. Default: 9.
#' @param save Logical value, whether to save all results. Default: FALSE.
#'
#' @return If \code{save} is TRUE, return NULL, else return a ggplot2 object.
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc select
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom grDevices pdf dev.off
#' @export
#'
#' @examples
#' # show p.adjust
#' # EnrichPlot(fe.res = fe.res, grad.col.vec=c("grey", "red"))
#' # EnrichPlot(fe.res = fe.res, grad.col.vec=c("grey", "red"), plot.type = "dot")
#' # do not show signif
#' # EnrichPlot(fe.res = fe.res, signif = NULL, fill.color = "red")
#' # EnrichPlot(fe.res = fe.res, signif = NULL, fill.color = "red", plot.type = "dot")
EnrichPlot <- function(fe.res, signif = "p.adjust", show.num = 10, grad.col.vec = c("grey", "red"),
                       fill.color = "red", str.width = 30, plot.type = c("bar", "dot"),
                       plot.resolution = 300, plot.width = 7, plot.height = 9, save = FALSE) {
  # check parameters
  plot.type <- match.arg(arg = plot.type)
  # create dataframe
  fe.res <- as.data.frame(fe.res)
  # check dataframe
  if (nrow(fe.res) < 1) {
    warning("The functional enrichment dataframe is empty!")
    enrich.plot <- ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "N/A")) +
      xlab(NULL)
  } else {
    # calculated fold enrichment
    fe.res$BgRatio <- sapply(X = fe.res$BgRatio, FUN = function(x) {
      eval(parse(text = x))
    })
    fe.res$GeneRatio <- sapply(X = fe.res$GeneRatio, FUN = function(x) {
      eval(parse(text = x))
    })
    fe.res$FoldEnrichment <- round(fe.res$GeneRatio / fe.res$BgRatio, 4)
    fe.res <- fe.res %>% dplyr::arrange(desc(FoldEnrichment))
    # prepare show number
    if (show.num < nrow(fe.res)) {
      fe.res.used <- fe.res[1:show.num, ]
    } else {
      fe.res.used <- fe.res[1:nrow(fe.res), ]
    }
    fe.res.used$Description <- factor(fe.res.used$Description, levels = rev(fe.res.used$Description))
    if (is.null(signif)) {
      # process functional enrichment results
      fe.res.used <- fe.res.used %>% dplyr::select(c("ID", "Description", "FoldEnrichment", "Count"))
      # create plot
      if (plot.type == "bar") {
        enrich.plot <- ggplot(fe.res.used) +
          geom_bar(aes_string(x = "FoldEnrichment", y = "Description"),
            stat = "identity", width = 0.6, fill = fill.color
          ) +
          labs(x = "Fold Enrichment", y = "Enriched Terms")
      } else if (plot.type == "dot") {
        enrich.plot <- ggplot(fe.res.used) +
          geom_point(aes_string(x = "FoldEnrichment", y = "Description", size = "Count"),
            colour = fill.color
          ) +
          labs(x = "Fold Enrichment", y = "Enriched Terms")
      }
    } else {
      # process functional enrichment results
      fe.res.used <- fe.res.used %>% dplyr::select(c("ID", "Description", "FoldEnrichment", signif, "Count"))
      # calculated -log10(pvalue)
      colnames(fe.res.used) <- c("ID", "Description", "FoldEnrichment", "signif", "Count")
      fe.res.used$signif <- -1 * log10(fe.res.used$signif)
      # create plot
      if (plot.type == "bar") {
        enrich.plot <- ggplot(fe.res.used) +
          geom_bar(aes_string(x = "FoldEnrichment", y = "Description", fill = "signif"),
            stat = "identity", width = 0.6
          ) +
          scale_fill_gradient(low = grad.col.vec[1], high = grad.col.vec[2]) +
          labs(fill = paste0("-log10 (", signif, ")"), x = "Fold Enrichment", y = "Enriched Terms")
      } else if (plot.type == "dot") {
        enrich.plot <- ggplot(fe.res.used) +
          geom_point(aes_string(x = "FoldEnrichment", y = "Description", size = "Count", color = "signif")) +
          labs(x = "FoldEnrichment", y = "Enriched Terms") +
          scale_color_gradient(low = grad.col.vec[1], high = grad.col.vec[2]) +
          labs(color = paste0("-log10 (", signif, ")"), x = "Fold Enrichment", y = "Enriched Terms")
      }
    }
  }
  # chenge term length
  if (!is.null(str.width)) {
    enrich.plot <- enrich.plot +
      scale_y_discrete(labels = function(enrich) stringr::str_wrap(enrich, width = str.width))
  }
  # change theme
  enrich.plot <- enrich.plot +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  # plot
  if (save) {
    # save png
    ggsave("Enrichment_plot.png", plot = enrich.plot, dpi = plot.resolution, width = plot.width, height = plot.height)
    # save pdf
    grDevices::pdf("Enrichment_plot.pdf")
    print(enrich.plot)
    grDevices::dev.off()
    return(NULL)
  } else {
    return(enrich.plot)
  }
}
