# The following discrete color palettes are from:
# * scanpy (with scanpy suffix): https://github.com/scverse/scanpy/blob/master/scanpy/plotting/palettes.py
# * seurat (with seurat suffix): https://github.com/satijalab/seurat/blob/master/R/visualization.R
# * dittoseq (with dittoseq suffix): https://github.com/dtm2451/dittoSeq/blob/master/R/dittoColors.R

#' Create Discrete Color Palette.
#'
#' @param palette Palette used, chosen from "dark2_8", "accent_8", "set2_8", "set1_9", "vega_10_scanpy",
#' "set3_12", "vega_20_scanpy", "stepped_24_seurat", "alphabet_26_seurat", "alphabet2_26_seurat",
#' "zeileis_28_scanpy", "glasbey_32_seurat", "polychrome_36_seurat", "cbf_40_dittoseq", "godsnot_102_scanpy".
#' The last number indicates the length of this palette. Default: "dark2_8".
#'
#' @return A vector of colors.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics image
#' @export
#'
#' @examples
#' DiscreteColorPalette(palette = "dark2_8")
DiscreteColorPalette <- function(palette = c(
                                   "dark2_8", "accent_8", "set2_8", "set1_9", "vega_10_scanpy",
                                   "set3_12", "vega_20_scanpy", "stepped_24_seurat", "alphabet_26_seurat", "alphabet2_26_seurat",
                                   "zeileis_28_scanpy", "glasbey_32_seurat", "polychrome_36_seurat", "cbf_40_dittoseq", "godsnot_102_scanpy"
                                 )) {
  # check parameters
  palette <- match.arg(arg = palette)

  # all available palettes
  palettes.list <- list(
    dark2_8 = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
    accent_8 = RColorBrewer::brewer.pal(n = 8, name = "Accent"),
    set2_8 = RColorBrewer::brewer.pal(n = 8, name = "Set2"),
    set1_9 = RColorBrewer::brewer.pal(n = 9, name = "Set1"),
    vega_10_scanpy = c(
      "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc",
      "#8c564b", "#e377c2", "#7f7f7f", "#b5bd61", "#17becf"
    ),
    set3_12 = RColorBrewer::brewer.pal(n = 12, name = "Set3"),
    vega_20_scanpy = c(
      "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc",
      "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8",
      "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
      "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31"
    ),
    stepped_24_seurat = c(
      "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F",
      "#B3823E", "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E",
      "#A3CC7A", "#CFE6B8", "#0F8299", "#3E9FB3", "#7ABECC",
      "#B8DEE6", "#3D0F99", "#653EB3", "#967ACC", "#C7B8E6",
      "#333333", "#666666", "#999999", "#CCCCCC"
    ),
    alphabet_26_seurat = c(
      "#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919",
      "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5",
      "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405",
      "#FFA8BB", "#426600", "#FF0010", "#5EF1F2", "#00998F",
      "#E0FF66", "#740AFF", "#990000", "#FFFF80", "#FFE100",
      "#FF5005"
    ),
    alphabet2_26_seurat = c(
      "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656",
      "#1C8356", "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F",
      "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16",
      "#F8A19F", "#90AD1C", "#F6222E", "#1CFFCE", "#2ED9FF",
      "#B10DA1", "#C075A6", "#FC1CBF", "#B00068", "#FBE426",
      "#FA0087"
    ),
    zeileis_28_scanpy = c(
      "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784",
      "#8e063b", "#4a6fe3", "#8595e1", "#b5bbe3", "#e6afb9",
      "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7",
      "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6",
      "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4", "#7f7f7f",
      "#c7c7c7", "#1CE6FF", "#336600"
    ),
    glasbey_32_seurat = c(
      "#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6",
      "#005300", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE",
      "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", "#F1085C",
      "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95",
      "#02AD24", "#C8FF00", "#886C00", "#FFB79F", "#858567",
      "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
      "#004CFF", "#F2F318"
    ),
    polychrome_36_seurat = c(
      "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32",
      "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
      "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B",
      "#C4451C", "#1C8356", "#85660D", "#B10DA1", "#FBE426",
      "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6",
      "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5",
      "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF",
      "#3B00FB"
    ),
    cbf_40_dittoseq = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
      "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
      "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
      "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
      "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57",
      "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
      "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45",
      "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D"
    ),
    godsnot_102_scanpy = c(
      "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6",
      "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601",
      "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA",
      "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101",
      "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66",
      "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459",
      "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
      "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
      "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625",
      "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98",
      "#A4E804", "#324E72"
    )
  )

  # get color palette
  color.palette <- palettes.list[[palette]]
  # show selected color palette
  image(
    1:length(color.palette), 1, as.matrix(1:length(color.palette)),
    col = color.palette,
    xlab = palette, ylab = "", xaxt = "n", yaxt = "n", bty = "n"
  )
  return(color.palette)
}

# The following continuous color palettes are from:
# viridis color maps: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
# custom palette: https://github.com/satijalab/seurat/blob/master/R/visualization.R
SeuratCustomPalette <- function(low = "white",
                                high = "red",
                                mid = NULL,
                                k = 50) {
  low <- grDevices::col2rgb(col = low) / 255
  high <- grDevices::col2rgb(col = high) / 255
  if (is.null(x = mid)) {
    r <- seq(from = low[1], to = high[1], len = k)
    g <- seq(from = low[2], to = high[2], len = k)
    b <- seq(from = low[3], to = high[3], len = k)
  } else {
    k2 <- round(x = k / 2)
    mid <- grDevices::col2rgb(col = mid) / 255
    r <- c(
      seq(from = low[1], to = mid[1], len = k2),
      seq(from = mid[1], to = high[1], len = k2)
    )
    g <- c(
      seq(from = low[2], to = mid[2], len = k2),
      seq(from = mid[2], to = high[2], len = k2)
    )
    b <- c(
      seq(from = low[3], to = mid[3], len = k2),
      seq(from = mid[3], to = high[3], len = k2)
    )
  }
  return(grDevices::rgb(red = r, green = g, blue = b))
}


#' Create Continuous Color Map.
#'
#' @param n The number of colors (≥ 1) to be in the palette.
#' @param direction Sets the order of colors in the scale. If 1, the default, colors are ordered from darkest to lightest.
#' If -1, the order of colors is reversed.
#' @param maps Existing color maps, chosen from "magma","inferno","plasma","viridis","cividis","rocket","mako","turbo".
#' Default: NULL (use \code{low}, \code{high}, \code{mid} to create).
#' @param low The low color. Default: white.
#' @param high The high color. Default: red.
#' @param mid The mid color. Default: NULL.
#' @param ... Parameters for \code{\link{viridis}}.
#'
#' @return A vector of colors.
#' @importFrom viridis viridis
#' @importFrom graphics image
#' @importFrom grDevices col2rgb rgb
#' @export
#'
#' @examples
#' # use existing color map
#' ContinuousColorMap(maps = "viridis")
#' # create color map with low, high, mid.
#' ContinuousColorMap(maps = NULL, low = "magenta", high = "yellow", mid = "black")
ContinuousColorMap <- function(n = 50, direction = 1, maps = NULL,
                               low = "white", high = "red", mid = NULL, ...) {
  # all available maps
  all.maps <- c("A", "B", "C", "D", "E", "F", "G", "H")
  names(all.maps) <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")

  # use existing color map
  if (!is.null(maps)) {
    if (maps %in% names(all.maps)) {
      color.palette <- viridis::viridis(n = n, direction = direction, option = all.maps[maps], ...)
      image(
        1:n, 1, as.matrix(1:n),
        col = viridis::viridis(n = n, direction = direction, option = all.maps[maps], ...),
        xlab = paste0("map= ", maps, " n= ", n, " direction= ", direction), ylab = "", xaxt = "n", yaxt = "n", bty = "n"
      )
    } else {
      stop("The maps you provided is not valid!")
    }
  } else {
    if (!is.null(low) & !is.null(high)) {
      color.palette <- SeuratCustomPalette(low = low, high = high, mid = mid, k = n)
      image(
        1:n, 1, as.matrix(1:n),
        col = color.palette,
        xlab = paste0("Low=", low, " High=", high, " Mid=", ifelse(is.null(mid), "null", mid), " n=", n), ylab = "", xaxt = "n", yaxt = "n", bty = "n"
      )
    } else {
      stop("Neither low or high should be NULL!")
    }
  }
  return(color.palette)
}
