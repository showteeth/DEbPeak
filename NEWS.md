# DEbPeak 0.8.0
## New features
* Changed `DEbPeak`, `PlotDEbPeak`, `DEbPeakFE`, `FindMotif` to integrate differential expression analysis results of RNA-seq and ChIP-seq/ATAC-seq.

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
