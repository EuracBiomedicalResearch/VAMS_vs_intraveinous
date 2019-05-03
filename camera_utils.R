## Utility functions to work with the output from the CAMERA package.

#' @description
#'
#' Function to group features into compounds based on adduct and/or isotope
#' information from CAMERA. Features are grouped into the same compound if they
#' share an adduct or isotope ID. (or share adducts or isotopes).
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
#' df <- data.frame(pcgroup = c(1, 2, 3, 3, 3, 4, 5, 5, 3, 6, 7, 3, 7, 7, 8, 8,
#'     9, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11))
#' df$isotopeGroup <- list(NA_integer_, NA_integer_, 4, 4, NA_integer_, NA_integer_,
#'                    NA_integer_, NA_integer_, NA_integer_, NA_integer_,
#'                    5, NA_integer_, 5, NA_integer_, NA_integer_, NA_integer_,
#'                    NA_integer_, NA_integer_, 6, 6, 7, 8, 8, NA_integer_,
#'                    NA_integer_, NA_integer_, 8, 8, 9)
#' df$adductGroup <- list(integer(), integer(), 1:2, integer(), 1L, integer(),
#'                        integer(), integer(), 2L, integer(), integer(),
#'                        integer(), integer(), integer(), 3L, 3L,
#'                        4:5, 4:5, integer(), integer(), integer(), 6L,
#'                        integer(), 6L, integer(), 7:8, 7L, integer(), 8L)
#'
#' res <- unsplit(lapply(split.data.frame(df, df$pcgroup),
#'     FUN = group_compound), f = df$pcgroup)
#'
#' x <- ExportTable[ExportTable$pcgroup == 1, ]
#' res <- group_compound(x, isotope_column = "isotopeGroups")
#' res[79:87]
#'
#' x <- x[79:87, ]
#' res <- group_compound(x, isotope_column = "isotopeGroups")
#' res
group_compound <- function(x, prefix = "rt", isotope_column = "isotopeGroup",
                           adduct_column = "adductGroup") {
    adduct_mat <- .logical_matrix(x[, adduct_column]) # Should we reduce this too?
    isotope_mat <- .logical_matrix(x[, isotope_column])
    ## any adduct shared with any isotope?
    ## Compare each column in adduct_mat with each column in isotope_mat and
    ## combine them with OR if they share any element.
    overlap <- .colwise_fun(adduct_mat, isotope_mat, FUN = comb_fun)
    overlap <- .drop_contained_columns(overlap, na.rm = TRUE)
    ## Assign to compound groups.
    comp_grp <- lapply(seq_len(nrow(overlap)), function(z) which(overlap[z, ]))
    miss <- lengths(comp_grp) == 0
    if (any(miss))
        comp_grp[miss] <- seq((ncol(overlap) + 1), ncol(overlap) + sum(miss))
    mapply(function(y, z) paste(prefix, y, z, sep = "."), x$pcgroup, comp_grp,
           SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

#' @description
#'
#' Reduce a `logical` matrix dropping columns which are completely within other
#' columns, i.e. if one column has values `c(TRUE, FALSE, TRUE, TRUE)` and
#' another column `c(TRUE, FALSE, FALSE, TRUE)`, only the first is reported
#' while the latter is dropped.
#'
#' @param x `matrix` of type `logical`.
#'
#' @return `matrix` of type `logical` with `rowSum` equal to 1 for all rows.
#'
#' @examples
#'
#' x <- cbind(
#'     c(TRUE, FALSE, FALSE, FALSE, TRUE),
#'     c(FALSE, TRUE, FALSE, TRUE, TRUE),
#'     c(FALSE, FALSE, TRUE, FALSE, FALSE),
#'     c(FALSE, FALSE, TRUE, FALSE, FALSE)
#'     )
#'
#' .drop_contained_colums(x)
#'
#' x <- cbind(
#'     c(TRUE, NA, FALSE, TRUE, FALSE),
#'     c(TRUE, TRUE, TRUE, TRUE, FALSE),
#'     c(TRUE, TRUE, FALSE, FALSE, FALSE),
#'     c(TRUE, FALSE, FALSE, FALSE, FALSE),
#'     c(FALSE, FALSE, FALSE, FALSE, TRUE)
#' )
#' .drop_contained_columns(x, na.rm = TRUE)
#'
#' x <- cbind(c(TRUE, FALSE, TRUE))
#' .drop_contained_columns(x)
#'
#' x <- cbind(
#'     c(FALSE, TRUE, FALSE),
#'     c(TRUE, FALSE, FALSE),
#'     c(FALSE, FALSE, TRUE)
#' )
#' .drop_contained_columns(x)
#'
#' x <- cbind(
#'     c(FALSE, FALSE, FALSE, TRUE),
#'     c(TRUE, FALSE, FALSE, TRUE),
#'     c(FALSE, TRUE, TRUE, TRUE),
#'     c(TRUE, TRUE, TRUE, TRUE)
#' )
#' .drop_contained_columns(x)
#'
#' x <- cbind(
#'     c(TRUE, FALSE, TRUE, FALSE),
#'     c(TRUE, FALSE, TRUE, FALSE),
#'     c(FALSE, FALSE, FALSE, TRUE),
#'     c(FALSE, FALSE, FALSE, TRUE)
#' )
#' .drop_contained_columns(x)
.drop_contained_columns <- function(x, na.rm = FALSE) {
    if (ncol(x) == 1)
        return(x)
    if (na.rm)
        x[is.na(x)] <- FALSE
    x <- unique(x, MARGIN = 2)
    ## Loop over columns
    res <- lapply(1:ncol(x), function(idx) {
        ovlp <- as.numeric(x[, idx] %*% x[, -idx, drop = FALSE])
        if (any(ovlp == sum(x[, idx])))
            NULL
        else x[, idx]
    })
    do.call(cbind, res)
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
    col <- 0
    res <- vector("list", ncx * ncy)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncy)) {
            col <- col + 1
            res[[col]] <- do.call(FUN, list(x[, i], y[, j]))
        }
    }
    do.call(cbind, res)
}

#' This function only combines the two input vectors with OR if they have an
#' overlap in any of their elements (i.e. `x` and `y` are `TRUE` at the same
#' position.
#'
#' @param x `logical`
#'
#' @param y `logical`
#'
#' @return `logical`, same length than `x` and `y`.
#'
#' @noRd
#'
#' @examples
#'
#' x <- c(TRUE, FALSE, TRUE, FALSE)
#' y <- c(FALSE, TRUE, TRUE, FALSE)
#'
#' comb_fun(x, y)
#'
#' y <- c(FALSE, TRUE, FALSE, FALSE)
#' comb_fun(x, y)
comb_fun <- function(x, y) {
    if (length(which(x & y)))
        x | y
    else cbind(x, y)
}
