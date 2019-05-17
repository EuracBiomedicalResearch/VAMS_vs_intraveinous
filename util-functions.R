## General utility functions for the analysis.

#' @title Plot of PCA results
#'
#' `plot_pca` is a simple utility function to plot the results from a PCA
#' analysis.
#'
#' @param pc the result from a principal component analysis (i.e. the result
#'     returned by `prcomp`.
#'
#' @param pch the point character. See [plot()] or [par()] for more information.
#'
#' @param col the color to be used for each data point/sample.
#'
#' @param pc_x `integer(1)` defining which principal component should be drawn
#'     on the x-axis.
#'
#' @param pc_y `integer(1)` defining the principal component to be drawn on the
#'     y-axis.
#'
#' @param main `character(1)` with the optional title of the plot.
#'
#' @param labels `character` with length equal to the number of samples. If
#'     provided, these will be displayed instead of data points.
#'
#' @param ... additional arguments to be passed to the [points()] or [text()]
#'     calls (if `labels = NULL` or not).
#'
#' @author Johannes Rainer
#'
#' @md
plot_pca <- function(pc, pch = 16, col = "#000000", pc_x = 1, pc_y = 2, 
                     main = "", labels = NULL, ...) {
    pcSummary <- summary(pc)
    plot(pc$x[, pc_x], pc$x[, pc_y], pch = NA, main = main,
         xlab = paste0("PC", pc_x, ": ",
                       format(pcSummary$importance[2, pc_x] * 100, 
                              digits = 3), " % variance"),
         ylab = paste0("PC", pc_y, ": ",
                       format(pcSummary$importance[2, pc_y] * 100, 
                              digits = 3), " % variance"))
    grid()
    if (!is.null(labels)) 
        text(pc$x[, pc_x], pc$x[, pc_y], labels = labels, col = col, 
             ...)
    else points(pc$x[, pc_x], pc$x[, pc_y], pch = pch, col = col, 
                ...)
}

#' @title Calculate relative standard deviations
#'
#' `rsd` and `rowRsd` are convenience functions to calculate the relative
#' standard deviation of a numerical vector or for rows of a numerical matrix,
#' respectively.
#'
#' @param x for `rsd` a `numeric`, for `rowRsd` a numerical `matrix`.
#'
#' @param na.rm `logical(1)` whether missing values (`NA`) should be removed
#' prior to the calculations.
#'
#' @author Johannes Rainer
#'
#' @md
rsd <- function(x, na.rm = TRUE) {
    sd(x, na.rm = na.rm) / abs(mean(x, na.rm = na.rm))
}
rowRsd <- function(x, na.rm = TRUE)
    apply(x, MARGIN = 1, rsd, na.rm = na.rm)

#' @title Determine the proportion of missing values
#'
#' `naProp` and `rowNaProp` determine the proportion of missing values in a
#' numeric vector or in rows of a numeric matrix
#'
#' @param x `numeric` or numeric `matrix`.
#'
#' @author Johannes Rainer
#'
#' @md
naProp <- function(x) {
    sum(is.na(x)) / length(x)
}
rowNaProp <- function(x) {
    apply(x, MARGIN = 1, naProp)
}

#'@title Extract run start time stamp
#'
#' @description
#'
#' Extract the start time stamps from all mzML files of an `MSnExp` object.
#'
#' @param x `MSnExp` object
#'
#' @param format `character(1)` defining the date/time format of the time
#'     stamp. If `NULL` the time stamp will be returned as a `character`.
#' 
#' @return `character` with the start time stamps.
extractTimeStamps <- function(x, format = "%Y-%m-%dT%H:%M:%S") {
    stopifnot(inherits(x, "MSnExp"))
    ts <- sapply(fileNames(x), function(z) {
        fl <- mzR::openMSfile(z)
        run_info <- mzR::runInfo(fl)
        mzR::close(fl)
        run_info$startTimeStamp
    })
    if (length(format))
        strptime(ts, format = format)
    else ts
}

#' @title Identify features matching provided masses
#' 
#' @description
#'
#' Toy helper function to identify all potential features for a (monoisotopic)
#' mass and list of potential adducts.
#'
#' @note
#'
#' The function requires the `CompoundDb` package to be installed and loaded.
#' 
#' @param mass `numeric(1)` with the mass of compound.
#'
#' @param rt `numeric(1)` with the retention time of the compound.
#'
#' @param object `DataFrame` with feature definitions (e.g. `rowData` for a
#'     `SummarizedExperiment` of `featureDefinitions` for an `XCMSnExp` object.
#'
#' @param adduct `character` defining the adduct(s) to check (e.g. `"[M+H]+"`).
#'
#' @param ppm `numeric(1)` defining the allowed difference between expected
#'     and actual m/z (in ppm).
#' 
#' @return
#'
#' `data.frame` with all features matching the standard's m/z given
#' the provided `ppm`. Columns are:
#'
#' - `"adduct"`: the matching adduct.
#' - `"feature_idx"`: the index (row) of the matching feature in `object`.
#' - `"mz_diff"`: the difference between the feature's m/z and the adduct's
#'   m/z (feature m/z - adduct m/z).
#' - `"mz_diff_ppm"`: the (absolute) difference of the m/z values in ppm.
#' - `"rt_diff"`: the difference between the feature's rt and the provided `rt`
#'   (feature's retention time - provided `rt`, negative values thus indicate
#'   the feature to elute *earlier* than expected).
#'
#' @author Johannes Rainer
#' 
#' @noRd
match_features <- function(mass, rt = 0, object, adduct = "[M+H]+", ppm = 5) {
    ## Check we have CompoundDb loaded
    if (!require("CompoundDb", quietly = TRUE, character.only = TRUE))
        stop("Package 'CompoundDb' not available! Install from ",
             "https://github.com/EuracBiomedicalResearch/CompoundDb")
    ## Find for each standard the features matching the m/z for potential adducts. 
    mzs <- unlist(mass2mz(mass, adduct = adduct))
    mz_matches <- lapply(mzs, function(z) {
        unlist(matchWithPpm(z, y = object$mzmed, ppm = ppm))
    })
    mz_matches <- mz_matches[lengths(mz_matches) > 0]
    if (length(mz_matches)) {
        cands <- data.frame(adduct = rep(names(mz_matches), lengths(mz_matches)),
                            feature_idx = unlist(mz_matches), row.names = NULL,
                            stringsAsFactors = FALSE)
        cands$mz_diff <- object$mzmed[cands$feature_idx] -
            mzs[cands$adduct]
        cands$mz_diff_ppm <- abs(cands$mz_diff) * 1e6 / mzs[cands$adduct]
        ## Calculte difference in retention time
        cands$rt_diff <- unlist(lapply(cands$feature_idx, function(z)
            object$rtmed[z] - rt))
        cands$feature_id <- rownames(object)[cands$feature_idx]
    } else cands <- data.frame()
    cands
}
