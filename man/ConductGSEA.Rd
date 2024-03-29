% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_set_enrichment_analysis.R
\name{ConductGSEA}
\alias{ConductGSEA}
\title{Conduct Gene Set Enrichment Analysis (GSEA).}
\usage{
ConductGSEA(
  deres,
  gmt.file,
  gene.sets = NULL,
  out.folder = NULL,
  gene.key = NULL,
  gene.type = c("ENSEMBL", "ENTREZID", "SYMBOL"),
  org.db = "org.Mm.eg.db",
  minGSSize = 10,
  maxGSSize = 500,
  pvalue = 0.05,
  padj.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr",
    "none"),
  save = TRUE,
  ...
)
}
\arguments{
\item{deres}{Data frame contains all genes.}

\item{gmt.file}{Gene Matrix Transposed file format.}

\item{gene.sets}{Gene sets information, containing two columns: gs_name, entrez_gene. Default: NULL.}

\item{out.folder}{Folder to save enrichment results. Default: wording directory.}

\item{gene.key}{Column name in \code{deres} to conduct analysis. Default: NULL (use rownames of \code{deres}).}

\item{gene.type}{Gene name type. Chosen from ENSEMBL, ENTREZID,SYMBOL. Default: ENSEMBL.}

\item{org.db}{Organism database. Default: org.Mm.eg.db.}

\item{minGSSize}{Minimal size of each geneSet for analyzing. Default: 10.}

\item{maxGSSize}{Maximal size of genes annotated for testing. Default: 500.}

\item{pvalue}{Cutoff value of pvalue. Default: 0.05.}

\item{padj.method}{One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: BH.}

\item{save}{Logical value, whether to save results. Default: TRUE.}

\item{...}{Parameters for \code{\link{GSEA}}.}
}
\value{
NULL or list contains all results (\code{GSEA} is FALSE)
}
\description{
Conduct Gene Set Enrichment Analysis (GSEA).
}
\examples{
library(airway)
library(msigdbr)
library(DESeq2)
library(DEbPeak)
# load the data
data(airway)
dds <- DESeq2::DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq2::DESeq(dds)
dds.results <- results(dds, contrast = c("dex", "trt", "untrt"))
dds.results.ordered <- dds.results[order(dds.results$log2FoldChange, decreasing = TRUE), ]
h_t2g <- msigdbr(species = "Homo sapiens", category = "C5") \%>\% dplyr::select(gs_name, entrez_gene)
gsea.results <- ConductGSEA(deres = dds.results.ordered, gmt.file = NULL, gene.sets = h_t2g, org.db = "org.Hs.eg.db", pvalue = 0.05, save = FALSE)
}
