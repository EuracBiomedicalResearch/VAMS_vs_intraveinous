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
