# DEbPeak 1.2.0
## New features
* Added `ProcessEnhancer` to get near genes of differential peaks.
* Updated `DEbPeak` to support enhancer-related data (use `species` instead of `org.db`. add `enhancer`, `seq.style`, `gtf.file`, `dis.threshold` and `n.cores`).
* Added `NetViz` to visualize enhancer-gene network.

## Minor changes
* Fixed bug in `ConductFE` when no KEGG enrichment results.
* Fixed bug in `ConductDESeq2` when total PCs number is smaller than 5.
* Fixed bug in `AnnoPeak` when using gtf file.
* Fixed bus in `PeakMatrix`.
* Added `count.matrix` parameter in `PeakMatrix` to support custom count matrix.
* Simplified the output of enhancer when integrating with RNA-seq.
* Fixed bug when `tibble` version > 3.1.8 (do not allow multiple columns with same name).

-------------

# DEbPeak 1.1.0
## New features
* Added `ParseGEO` to parse data from GEO.

-------------

# DEbPeak 1.0.1
## New features
* Added `DiffPeakPie` to stat genomic regions of differential peaks with pie plot.
* Added `EnrichPlot` to create bar or dot plot for functional enrichment analysis, this is useful for visualizing selected terms.

## Minor changes
* Fixed bug in `FindMotif`

-------------

# DEbPeak 1.0.0
## New features
* Added `DEbDE` and `DEbDEFE` to support integrating two differential expression results.
* Added `PeakbPeak` and `PeakbPeakFE` to support integrating peak annotation/differential expression results.
* Rewrote `InteDiffQuad` to support integration results of peak-related data and peak-related data (PeakbPeak), RNA-seq and RNA-seq (DEbDE), RNA-seq and peak-related data (DEbPeak).
* Added `InteFE` to replace `DEbDEFE`, `PeakbPeakFE` and `DEbPeakFE`.
* Added `InteVenn` to support integration results of peak-related data and peak-related data (PeakbPeak), RNA-seq and RNA-seq (DEbDE), RNA-seq and peak-related data (DEbPeak).

## Minor changes
* Fixed the way `PrepareDEPlot` filtering.
* Removed parameter `gene.type` in `DEbPeakFE`.
* Added parameter `gene.col` in `PlotDEbPeak` to show genes instead of peaks.
* Fixed bug in `PrepareVenn`.
* Fixed spelling mistakes (`consenus` -> `consensus`).

-------------

# DEbPeak 0.9.0
## New features
* Added `MotifEnrich` to support motif enrichment for differentially accessible/binding peaks.
* Updated documents.

## Minor changes
* Fixed bugs in `ExtractDA`.
* Fixed bugs in `DEbPeak`.

-------------

# DEbPeak 0.8.0
## New features
* Changed `DEbPeak`, `PlotDEbPeak`, `DEbPeakFE`, `FindMotif` to integrate differential expression analysis results of RNA-seq and ChIP-seq/ATAC-seq.
* Added `InteDiffQuad` to create quadrant diagram for differential expression analysis of RNA-seq and Peak-related data.

## Minor changes
* Added `peak.anno.key` to filter differential expression analysis results of ChIP-seq/ATAC-seq in `ConductDESeq2`.

-------------

# DEbPeak 0.7.0
## New features
* Added `PeakMatrix` to prepare peak matrix for differential expression analysis.
* Changed `ExportPCGenes`, `LoadingBar`, `LoadingHeat`, `LoadingPlot`, `LoadingGO`, `ConductDESeq2`, `ConductFE` to support ChIP-seq/ATAC-seq differential expression analysis results.

## Minor changes
* Fixed bugs in `ConductDESeq2`.

-------------

# DEbPeak 0.6.0
## New features
* Changed `CountQC` realted plots to ggplot2.
* Changed `AnnoPeak` pie plots to ggplot2.
* Exported `GetGeneLength` to get gene length. 
* Added `gene.length.file` to `NormalizedCount` to support gene length file as input.

## Minor changes
* Fixed some spelling mistakes (`threashold` -> `threshold`).

-------------

# DEbPeak 0.5.0
## New features
* Added `VisGSEA` to visualize GSEA results.

-------------

# DEbPeak 0.4.0
## New features
* Added `DEbCA` to support integrating RNA-seq, ChIP-seq and ATAC-seq.
* Modified `PlotDEbPeak` to support integrating RNA-seq, ChIP-seq and ATAC-seq.
* Updated all vignettes.

## Minor changes
* Fixed bugs in plot.

-------------

# DEbPeak 0.3.0
## New features
* Changed package name to `DEbPeak` to broaden its usage for ChIP-seq, ATAC-seq.

-------------

# DEbPeak 0.2.0
## New features
* Added `FindMotif` to find motif in interested peaks (suitable for RNA-seq and ATAC-seq integration).

-------------

# DEbPeak 0.1.0

* Added a `NEWS.md` file to track changes to the package.
