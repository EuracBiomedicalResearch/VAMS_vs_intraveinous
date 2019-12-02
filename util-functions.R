## General utility functions for the analysis.



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
