## Utility functions to work with the output from the CAMERA package.
#' @description
#'
#' Function to group features into compounds based on adduct and/or isotope
#' information from CAMERA.
#'
#' @param x `data.frame` with results from `CAMERA` for one *pcgroup*. See
#'     example below for the expected input.
#'
#' @param prefix `character` to be appended to the newly generated
#'     *compound IDs*.
#'
#' @return `list` of `character` (length 1 or more, depending on the number of
#'     compounds a feature can possibly be assigned) with the compound IDs.
#' 
#' @author Johannes Rainer, Soren Fjelstrup
#'
#' @examples
#' 
#' df <- data.frame(pcgroup = c(1, 2, 3, 3, 3, 4, 5, 5, 3, 6, 7, 3, 7, 7, 8, 8))
#' df$isotope <- list(NA_integer_, NA_integer_, 4, 4, NA_integer_, NA_integer_,
#'                    NA_integer_, NA_integer_, NA_integer_, NA_integer_,
#'                    5, NA_integer_, 5, NA_integer_, NA_integer_, NA_integer_)
#' df$adductGroup <- list(integer(), integer(), 1:2, integer(), 1L, integer(),
#'                        integer(), integer(), 2L, integer(), integer(),
#'                       integer(), integer(), integer(), 3L, 3L)
#'
#' res <- unsplit(lapply(split.data.frame(df, df$pcgroup),
#'     FUN = group_compound), f = df$pcgroup)
group_compound <- function(x, prefix = "rt") {
    adduct_mat <- .logical_matrix(x$adductGroup)
    isotope_mat <- .logical_matrix(x$isotope)
    overlap <- .colwise_fun(adduct_mat, isotope_mat)
    ## Assign to compound groups.
    comp_grp <- lapply(seq_len(nrow(overlap)), function(z) which(overlap[z, ]))
    miss <- lengths(comp_grp) == 0
    if (any(miss))
        comp_grp[miss] <- seq((ncol(overlap) + 1), ncol(overlap) + sum(miss))
    mapply(function(y, z) paste(prefix, y, z, sep = "."), x$pcgroup, comp_grp,
           SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

#' Create a logical matrix with columns being unique elements in `x` and
#' number of rows equal to the length of `x`.
#'
#' @param x `list` of group/class assignments.
#'
#' @return `logical` `matrix`.
#'
#' @noRd
.logical_matrix <- function(x) {
    u_x <- unique(unlist(x))
    u_x <- u_x[!is.na(u_x)]
    mat <- matrix(ncol = length(u_x), nrow = length(x))
    colnames(mat) <- letters[seq_along(u_x)]
    for (i in seq_along(u_x)) {
        mat[, i] <- vapply(x, function(z) length(which(z %in% u_x[i])) > 0,
                           logical(1))
    }
    mat
}

#' Combine each column in `x` with each columns in `y` applying the function
#' `FUN` (by default `|`).
#'
#' @param `matrix` of type `logical`.
#'
#' @param `matrix` of type `logical` with the same number of rows than `x`.
#'
#' @param `FUN` function or name of function to be applied to the pairwise
#'     columns.
#' 
#' @return `matrix`, same number of rows than `x` or `y` and
#'     `ncol(x) * ncol(y)` columns.
#'
#' @noRd
#'
#' @examples
#' 
#' a <- cbind(a = c(TRUE, FALSE, FALSE, TRUE),
#'            b = c(FALSE, TRUE, TRUE, FALSE))
#'
#' b <- cbind(c = c(FALSE, TRUE, FALSE, TRUE),
#'            d = c(TRUE, FALSE, FALSE, FALSE),
#'            e = c(FALSE, TRUE, FALSE, FALSE),
#'            f = c(FALSE, FALSE, TRUE, FALSE))
#'
#' .colwise_fun(a, b)
#' .colwise_fun(b, a)
#'
#' .colwise_fun(a[, 1, drop = FALSE], b)
#'
#' .colwise_fun(a[, 1, drop = FALSE], b[, 1, drop = FALSE])
#'
#' .colwise_fun(a[, 1, drop = FALSE], matrix(c(NA, NA, NA, NA), ncol = 1))
#'
#' .colwise_fun(a[, integer(), drop = FALSE], b)
#' .colwise_fun(a, b[, integer(), drop = FALSE])
.colwise_fun <- function(x, y, FUN = "|") {
    if (!(is.matrix(x) && is.matrix(y)))
        stop("'x' and 'y' have to be a matrix")
    if (nrow(x) != nrow(y))
        stop("'x' and 'y' must have the same number of rows")
    ncx <- ncol(x)
    ncy <- ncol(y)
    if (!ncx & !ncy)
        return(matrix(NA, nrow = nrow(x)))
    if (!ncx) {
        ncx <- 1
        x <- matrix(nrow = nrow(y), ncol = 1, dimnames = list(NULL, "a"))
    }
    if (!ncy) {
        ncy <- 1
        y <- matrix(nrow = nrow(x), ncol = 1, dimnames = list(NULL, "a"))
    }
    res <- matrix(nrow = nrow(x), ncol = ncx * ncy)
    colnames(res) <- paste(rep(colnames(x), each = ncy), rep(colnames(y), ncx),
                           sep = ".")
    col <- 0
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncy)) {
            col <- col + 1
            res[, col] <- do.call(FUN, list(x[, i], y[, j]))
        }
    }
    res
}
