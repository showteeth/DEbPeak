#' Parse GEO Data.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field.
#' @param supp.idx The index of supplementary files to download. This should be consistent with \code{platform}. Default: 1.
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or emoty,
#' download supplementary files automatically). Default: FALSE.
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return List contains GEO object of platform, study information, raw count matrix and metadata.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO getGEOSuppFiles gunzip
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom tools file_ext
#' @importFrom utils untar
#' @importFrom data.table fread
#' @importFrom openxlsx read.xlsx
#' @export
#'
#' @examples
#' # GSE149838.list = ParseGEO(acce = "GSE149838", platform = "GPL21626")
#' # GPL28369 is consistent with first supplementary file
#' # GSE147507.list = ParseGEO(acce = "GSE147507", platform = "GPL28369", supp.idx = 1)
#' # GSE122774.list = ParseGEO(acce = "GSE122774", platform = "GPL17021")
ParseGEO <- function(acce, platform, supp.idx = 1, down.supp = FALSE, ...) {

  # get GEO object
  pf.obj <- GEOobj(acce = acce, platform = platform, ...)
  # extract general information
  pf.info <- ExtractGEOInfo(pf.obj = pf.obj, sample.wise = FALSE)
  # extract raw counts
  pf.count <- ExtractGEOExp(pf.obj = pf.obj, acce = acce, supp.idx = supp.idx, down.supp = down.supp)
  # select meta data
  pf.meta <- ExtractGEOMeta(pf.obj = pf.obj)
  # return list
  res.list <- list(
    obj = pf.obj, exp.info = pf.info,
    count = pf.count, metadata = pf.meta
  )
  return(res.list)
}

# connect to GEO, extract GEO object, extract platform object
GEOobj <- function(acce, platform, ...) {
  # obtain GEO object
  geo.obj <- GEOquery::getGEO(GEO = acce, ...)

  # extract platform
  pfs <- sapply(geo.obj, function(x) {
    Biobase::annotation(x)
  })
  if (!platform %in% pfs) {
    stop(paste("The platform you provides is not valid!", paste(pfs, collapse = ", ")))
  }
  # extract platform data
  pf.idx <- which(pfs == platform[1])
  pf.obj <- geo.obj[[pf.idx]]
  return(pf.obj)
}

# merge cols into one
SimpleCol <- function(df, col) {
  cols <- grep(pattern = col, x = colnames(df), value = T)
  df.cols <- df[cols]
  value <- paste(c(t(df.cols)), collapse = ". ")
  return(value)
}

#' Extract GEO Study Information.
#'
#' @param pf.obj GEO object of platform.
#' @param sample.wise Logical value, whether to extract sample-wise information. Default: FALSE.
#'
#' @return A dataframe.
#'
#' @examples
#' # pf.obj = GEOobj(acce = "GSE147507", platform = "GPL28369")
#' # pf.info = ExtractGEOInfo(pf.obj)
ExtractGEOInfo <- function(pf.obj, sample.wise = FALSE) {
  # platform information
  pf.info <- Biobase::experimentData(pf.obj)
  # additional information
  pf.info.add.df <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  if (sample.wise) {
    pf.info.final <- cbind(data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds), pf.info.add.df
    )) %>%
      as.data.frame()
  } else {
    used.cols <- c("organism", "molecule", "strategy", "extract_protocol", "data_processing")
    pf.info.add.used <- pf.info.add.df[, grep(pattern = paste(used.cols, collapse = "|"), colnames(pf.info.add.df), value = T)]
    pf.info.add.sim <- apply(pf.info.add.used, 2, function(x) {
      paste(unique(x), collapse = ", ")
    }) %>%
      t() %>%
      as.data.frame()
    # final information
    pf.info.final <- data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Organism = SimpleCol(df = pf.info.add.sim, col = "organism"),
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      Molecule = SimpleCol(df = pf.info.add.sim, col = "molecule"),
      ExtractProtocol = SimpleCol(df = pf.info.add.sim, col = "extract_protocol"),
      LibraryStrategy = SimpleCol(df = pf.info.add.sim, col = "strategy"),
      DataProcessing = SimpleCol(df = pf.info.add.sim, col = "data_processing"),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      Contact = paste(pf.info@name, pf.info@contact, sep = "; "),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds)
    )
  }
  return(pf.info.final)
}

#' Extract Sample Metadata.
#'
#' @param pf.obj GEO object of platform.
#'
#' @return A dataframe.
#'
#' @examples
#' # pf.obj = GEOobj(acce = "GSE147507", platform = "GPL28369")
#' # pf.info = ExtractGEOMeta(pf.obj)
ExtractGEOMeta <- function(pf.obj) {
  # extract sample detail information
  pf.info <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  # select used basic cols
  pf.info.used <- pf.info[c("title", "geo_accession", "source_name_ch1", "description")]
  # process characteristics
  pf.info.charac <- pf.info[grep(pattern = "^characteristics", x = colnames(pf.info))]
  ## modify colnames
  pf.info.charac.colnames <-
    unique(apply(pf.info.charac, 2, function(x) {
      gsub(pattern = "(.*?): (.*)", replacement = "\\1", x = x)
    }))
  colnames(pf.info.charac) <- pf.info.charac.colnames
  ## modify values
  pf.info.charac <- apply(pf.info.charac, 2, function(x) {
    gsub(pattern = "(.*?): (.*)", replacement = "\\2", x = x)
  })
  ## final meta
  pf.meta <- cbind(pf.info.used, pf.info.charac) %>% as.data.frame()

  return(pf.meta)
}


#' Extract Raw Count Matrix from Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#'
#' @return A dataframe.
#'
#' @examples
#' # count.mat = ExtractGEOExpSupp(acce = "GSE149838")
#' # count.mat = ExtractGEOExpSupp(acce = "GSE147507")
#' # count.mat = ExtractGEOExpSupp(acce = "GSE147507", supp.idx = 2)
#' # count.mat = ExtractGEOExpSupp(acce = "GSE122774")
ExtractGEOExpSupp <- function(acce, supp.idx = 1) {
  # create tmp folder
  tmp.folder <- tempdir()
  # download supplementary file
  supp.down.log <- GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
  if (supp.idx > nrow(supp.down.log)) {
    stop("Please provide valid supplementary file index.")
  }
  supp.file.path <- rownames(supp.down.log)[supp.idx]
  # remove unused supplementary file
  unused.supp <- setdiff(rownames(supp.down.log), supp.file.path)
  unused.supp.remove <- file.remove(unused.supp)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    GEOquery::gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # unzip
    unzip.log <- sapply(
      list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "gz$"),
      function(x) {
        GEOquery::gunzip(x, overwrite = TRUE)
      }
    )
    # read files
    count.list <- lapply(
      list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE),
      function(x) {
        sample.count <- data.table::fread(file = x) %>% as.data.frame()
        colnames(sample.count) <- c("GeneName", gsub(pattern = "(GSM[0-9]*).*", replacement = "\\1", x = basename(x)))
        sample.count
      }
    )
    # create count matrix
    count.mat <- Reduce(f = function(x, y) {
      merge.mat <- merge(x, y, by = "GeneName", all = T)
    }, x = count.list)
    rownames(count.mat) <- count.mat$GeneName
    count.mat$GeneName <- NULL
  } else {
    if (file.ext %in% c("xlsx", "xls")) {
      # read excel file
      count.mat <- openxlsx::read.xlsx(xlsxFile = supp.file.path, rowNames = TRUE)
    } else if (file.ext %in% c("csv", "tsv", "txt")) {
      # read text file
      count.mat <- data.table::fread(file = supp.file.path) %>% as.data.frame()
      # the first column must be gene
      rownames(count.mat) <- count.mat[, 1]
      count.mat[, 1] <- NULL
    }
  }
  return(count.mat)
}

#' Extract Raw Count Matrix.
#'
#' @param pf.obj GEO object of platform.
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or emoty,
#' download supplementary files automatically). Default: FALSE.
#'
#' @return A dataframe.
#'
#' @examples
#' # pf.obj = GEOobj(acce = "GSE147507", platform = "GPL28369")
#' # count.mat = ExtractGEOExp(pf.obj, acce = "GSE147507")
ExtractGEOExp <- function(pf.obj, acce, supp.idx = 1, down.supp = FALSE) {
  # download supplementary files
  if (down.supp) {
    count.mat <- ExtractGEOExpSupp(acce = acce, supp.idx = supp.idx)
  } else {
    expr.mat <- Biobase::exprs(pf.obj)
    if (nrow(expr.mat) == 0) {
      message("Matrix not available! Downloading supplementary files.")
      count.mat <- ExtractGEOExpSupp(acce = acce, supp.idx = supp.idx)
    } else {
      if (all(expr.mat %% 1 == 0)) {
        count.mat <- expr.mat
      } else {
        message("Matrix contains non-integer values! Downloading supplementary files.")
        count.mat <- ExtractGEOExpSupp(acce = acce, supp.idx = supp.idx)
      }
    }
  }
  return(count.mat)
}


#' Download GEO Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. If NULL, download all supplementary files. Default: 1.
#' @param out.folder The output folder used to save downloaded results.
#'
#' @return A vector contains paths.
#' @importFrom GEOquery getGEOSuppFiles
#' @export
#'
#' @examples
#' # DownloadGEOSupp(acce = "GSE149836", supp.idx = NULL, out.folder = "/path/to/suppfiles")
DownloadGEOSupp <- function(acce, supp.idx = 1, out.folder = NULL) {
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # download supplementary file
  supp.down.log <- GEOquery::getGEOSuppFiles(GEO = acce, baseDir = out.folder)
  if (is.null(supp.idx)) {
    supp.file.path <- rownames(supp.down.log)
  } else {
    if (supp.idx > nrow(supp.down.log)) {
      stop("Please provide valid supplementary file index.")
    }
    supp.file.path <- rownames(supp.down.log)[supp.idx]
    # remove unused supplementary file
    unused.supp <- setdiff(rownames(supp.down.log), supp.file.path)
    unused.supp.remove <- file.remove(unused.supp)
  }
  return(supp.file.path)
}
