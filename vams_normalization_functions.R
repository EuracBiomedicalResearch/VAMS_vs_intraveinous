
#' Utility functions to plot model fits for features

#' @description
#'
#' Plot model fits for a *global* model like `y ~ inj_idx * batch`
plot_feature_slopes <- function(y, yfill, data, lmod, main) {
    ##' Two plots, one per batch.
    par(mar = c(0.5, 4.5, 1, 0), mfrow = c(1, 2))
    a <- which(data$batch == data$batch[1])
    y_a <- y[a]
    yfill_a <- yfill[a]
    dta <- data[a, ]
    plot(main = paste0(main, " batch a"), xlab = "", xaxt = "n",
         ylab = expression(log[2]~abundance), ylim = range(yfill, na.rm = TRUE),
         col = paste0(col_source[dta$source], 80), x = dta$inj_idx,
         y = yfill_a)
    points(x = dta$inj_idx, y = y_a, col = paste0(col_source[dta$source], 80),
           pch = 16)
    grid()
    if (!missing(lmod)) {
        if (!is.na(lmod)) {
            if (is.finite(lmod$coefficients[1]) &
                is.finite(lmod$coefficients[2]))
                abline(lmod$coefficients[1:2], col = col_source["all"])
            else warning("Feature ", ft, " has non-finite coeffs in batch a")
        }
    }
    par(mar = c(0.5, 0, 1, 4.5))
    a <- which(data$batch != data$batch[1])
    y_a <- y[a]
    yfill_a <- yfill[a]
    dta <- data[a, ]
    plot(main = paste0(main, " batch b"), xlab = "", xaxt = "n",
         ylab = expression(log[2]~abundance), ylim = range(yfill, na.rm = TRUE),
         col = paste0(col_source[dta$source], 80), x = dta$inj_idx,
         y = yfill_a, yaxt = "n")
    points(x = dta$inj_idx, y = y_a, col = paste0(col_source[dta$source], 80),
           pch = 16)
    grid()
    if (!missing(lmod)) {
        if (!is.na(lmod)) {
            a <- sum(lmod$coefficients[c(1, 3)])
            b <- sum(lmod$coefficients[c(2, 4)])
            if (is.finite(a) & is.finite(b))
                abline(a, b, col = col_source["all"])
            else warning("Feature ", ft, " has non-finite coeffs in batch b")
        }
    }
}

#' @description
#'
#' Plot model fit for separate model fits in two batches a and b.
#'
#' @param y `numeric` vector with the intensities to plot.
#'
#' @param is_filled `logical` indicating which of the data points are filled-in.
plot_feature_slopes_batch <- function(y, is_filled, data, lmoda,
                                      lmodb, main, legend = TRUE, ...) {
    ##' Two plots, one per batch.
    par(mar = c(0.5, 4.5, 1, 0), mfrow = c(1, 2))
    a <- which(data$batch == data$batch[1])
    y_a <- y[a]
    if (!missing(is_filled))
        pch <- ifelse(is_filled[a], yes = 1, no = 16)
    else
        pch <- 16
    dta <- data[a, ]
    plot(main = paste0(main, " batch a"), xlab = "", xaxt = "n",
         ylab = expression(log[2]~abundance), ylim = range(y, na.rm = TRUE),
         col = paste0(col_source[dta$source], 80), x = dta$inj_idx,
         y = y_a, pch = pch, ...)
    grid()
    if (!missing(lmoda)) {
        if (length(lmoda) > 1) {
            if (is.finite(lmoda$coefficients[1]) &
                is.finite(lmoda$coefficients[2]))
                abline(lmoda$coefficients[1:2], col = col_source["all"])
            else warning("Feature ", ft, " has non-finite coeffs in batch a")
            ##' diff_res <- diff(quantile(residuals(lmoda), probs = c(0.1, 1)))
            diff_res <- mean(abs(residuals(lmoda)))
            if (legend)
                legend("topright",
                       legend = paste0(c("res: ", "R2: "),
                                       c(format(diff_res, digits = 3),
                                         format(summary(lmoda)$adj.r.squared,
                                                digits = 3))))
        }
    }
    par(mar = c(0.5, 0, 1, 4.5))
    a <- which(data$batch != data$batch[1])
    y_a <- y[a]
    if (!missing(is_filled))
        pch <- ifelse(is_filled[a], yes = 1, no = 16)
    else
        pch <- 16
    dta <- data[a, ]
    plot(main = paste0(main, " batch b"), xlab = "", xaxt = "n",
         ylab = expression(log[2]~abundance), ylim = range(y, na.rm = TRUE),
         col = paste0(col_source[dta$source], 80), x = dta$inj_idx,
         y = y_a, yaxt = "n", pch = pch)
    grid()
    if (!missing(lmodb)) {
        if (length(lmodb) > 1) {
            if (is.finite(lmodb$coefficients[1]) &
                is.finite(lmodb$coefficients[2]))
                abline(lmodb$coefficients[1:2], col = col_source["all"])
            else warning("Feature ", ft, " has non-finite coeffs in batch a")
            ##' diff_res <- diff(quantile(residuals(lmodb), probs = c(0.1, 1)))
            diff_res <- mean(abs(residuals(lmodb)))
            if (legend)
                legend("topright",
                       legend = paste0(c("res: ", "R2:"),
                                       c(format(diff_res, digits = 3),
                                         format(summary(lmodb)$adj.r.squared,
                                                digits = 3))))
        }
    }
}

#' @title Apply a function to sets of columns of a matrix
#'
#' @description
#'
#' Apply a function to sets of columns, i.e. to sub-matrices of selected
#' columns defined with parameter `colgroup`.
#'
#' @param x `matrix` with numeric values.
#'
#' @param colgroup `character` or `factor` defining the sets of columns.
#'
#' @param FUN `function` to be applied to the sub-matrix.
#'
#' @param simplify `logical(1)` whether `cbind` should be called on the result.
#'
#' @return `matrix` with aggregated values. Number of columns represent the
#'     number of unique groups/sets defined by `colgroup`. Number of rows is
#'     either one (for `FUN` being e.g. `mean`, `sum` etc) or equal to the
#'     number of rows of `x` (for `FUN` being e.g. `rowSums`).
#'
#' @author Johannes Rainer
apply_colgroup <- function(x, colgroup, FUN, simplify = TRUE, ...) {
    if (missing(FUN)) stop("'FUN' is missing")
    if (length(colgroup) != ncol(x))
        stop("length of 'colgroup' should match ncol of 'x'")
    grps <- unique(colgroup)
    res <- lapply(grps, function(z) {
        x_sub <- x[, colgroup == z, drop = FALSE]
        FUN(x_sub, ...)
    })
    names(res) <- grps
    if (simplify)
        do.call(cbind, res)
    else res
}

#' @description
#'
#' Function to calculate row-wise maximum ratio of abundances, i.e. the ratio of
#' the largest difference of replicated measurements.
#'
#' @param x `matrix` with values. These have to be in **natural scale**.
#'
#' @return The *MRA* or `NA` if only a single value is available.
rowMra <- function(x) {
    mra <- function(z) {
        z <- z[!is.na(z)]
        if (length(z) > 1) {
            rng <- range(z)
            rng[2]/rng[1]
        } else {
            NA
        }
    }
    apply(x, MARGIN = 1, mra)
}

#' @description
#'
#' Calculates MRA values (maximum ratio of abundances) between replicated
#' columns (defined with argument `grp` and returns a `matrix` with these.
#' 
#' @param x matrix of normalized abundances.
#'
#' @param nofill matrix, same dimensions than x, that contain only detected
#'     signal. If provided, all values in `x` that are not detected are
#'     replaced with `NA`.
#'
#' @param grp vector defining which samples are replicates.
#'
#' @param drop_samples `character` with names of samples that should be
#'     dropped.
mra_for_mat <- function(x, nofill, grp, drop_samples) {
    if (!missing(nofill))
        x[is.na(nofill)] <- NA
    mras <- apply_colgroup(x, grp, rowMra)
    mras[, !(colnames(mras) %in% drop_samples)]
}

#' @description
#'
#' Calculates column-wise summary statistic on a matrix with MRA values.
mra_summary_quant <- function(x, probs = 0.75, na.rm = TRUE, ...) {
    apply(x, 2, quantile, probs = probs, na.rm = na.rm, ...)
}
mra_summary_count <- function(x, cut = 2, na.rm = TRUE) {
    apply(x, 2, function(z) {
        if (all(is.na(z))) NA
        else sum(z > cut, na.rm = na.rm)
    })
}
mra_summary_perc <- function(x, cut = 2, na.rm = TRUE) {
    apply(x, 2, function(z) {
        if (all(is.na(z))) NA
        else sum(z > cut, na.rm = na.rm) * 100 / sum(!is.na(z))
    })
}

dobox <- function(x, only_detected = FALSE,
                  col = paste0(col_source[data_pos$source], "ff"),
                  outline = FALSE, notch = TRUE, range = 0,
                  border = paste0(col_source[data_pos$source], "60"),
                  ylab = expression(log[2]~abundance), xaxt = "n", xlab = "",
                  ...) {
    if (only_detected)
        x[is.na(fv_nofill)] <- NA
    boxplot(x, col = col, outline = outline, notch = notch, range = range,
            border = border, ylab = ylab, xaxt = xaxt, xlab = xlab, ...)
    grid(nx = NA, ny = NULL)
}

#' Calculate the mean only if we'be got enough values and report `NA`
#' otherwise
mean_if <- function(x, n = 6, na.rm = TRUE) {
    if (na.rm)
        x <- x[!is.na(x)]
    if (length(x) >= n)
        mean(x)
    else NA_real_
}
