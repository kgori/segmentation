#' PCF - piecewise constant function fitting algorithm for 1-dimensional data.
#' 
#' @description Uses PELT (pruned exact linear time) optimisation for best case O(n), worst case O(n^2) time complexity.
#' Returns a list of segment starts, ends and widths, and within segment means.
#' Optionally, constrain segment starts to be drawn from a vector of preselected indices.
#' NB: Indices are 1-based.
#'     Segments are half-open intervals [start, end).
#' @param x Vector of data to segment
#' @param kmin Minimum segment size
#' @param gamma Penalty term added for each new segment. Higher values produce fewer segments.
#' @param bks Constrain the start indices of the result to be drawn from this vector.
#' @returns List of starts, ends and lengths of segments, and the mean value of x within each.
#' @export
pcf <- function(x, kmin, gamma, bks = NULL, integer_constraint = FALSE) {
    if (is.null(bks)) {
        bks <- integer(0)
    } else {
        bks <- bks - 1 # convert 1-based to 0-based
    }
    starts <- pelt_pcf_(x, kmin, gamma, bks, integer_constraint) + 1
    ends <- c(starts[starts > 1] - 1, length(x))
    lengths <- ends - starts + 1
    means <- sapply(seq_along(starts), function(i) {
        mean(x[starts[i]:ends[i]])
    })
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

#' MultiPCF - piecewise constant function fitting algorithm for 2-dimensional data.
#' 
#' @description Uses PELT (pruned exact linear time) optimisation for best case O(n), worst case O(n^2) time complexity.
#' Returns a list of segment starts, ends and widths, and within segment means for each column.
#' Optionally, constrain segment starts to be drawn from a vector of preselected indices.
#' NB: Input data has observations in rows, samples in columns.
#'     Indices are 1-based.
#'     Segments are half-open intervals [start, end).
#' @param x Matrix of values to segment, samples in columns
#' @param kmin Minimum segment size
#' @param gamma Penalty term added for each new segment. Higher values produce fewer segments.
#' @param bks Constrain the start indices of the result to be drawn from this vector.
#' @param w Optional vector of weights to apply to the columns of x, to up- or down-weight the contribution of each sample to the result.
#' @returns List of starts, ends and lengths of segments, and the mean value of x within each.
#' @export
multipcf <- function(x, kmin, gamma, bks = NULL, w = NULL, integer_constraint = FALSE) {
    if (!is.null(w)) {
        stopifnot(all(w > 0))
        x <- sweep(x, 2, w, "*")
    }
    if (is.null(bks)) {
        bks <- integer(0)
    } else {
        bks <- bks - 1 # convert 1-based to 0-based
    }
    starts <- pelt_multipcf_(x, kmin, gamma, bks, integer_constraint) + 1
    ends <- c(starts[starts > 1] - 1, nrow(x))
    lengths <- ends - starts + 1
    means <- t(sapply(seq_along(starts), function(i) {
        colMeans(x[starts[i]:ends[i], ])
    }))
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

#' Estimate the standard deviation of PCF input data.
#' 
#' @description Uses Median Absolute Deviation applied to the first differences
#' of the data, for robustness to any changes of location present.
#' @returns an estimate of SD for each column of the input.
#'   NB: Assumes that the deviation is constant along the input.
#' @export
estimate_sd <- function(x) {
    apply(as.matrix(x), 2, function(col) { mad(diff(col)) / sqrt(2) })
}

#' Estimate the PCF gamma parameter using BIC.
#' 
#' @description This function can provide a reasonable gamma penalty term for PCF or MultiPCF.
#' It uses BIC to adjust to the size of the data and the amount of variance in each
#' sample. It's not magic, though, and might not produce the segmentation you want.
#' 
#' NB: An alternative is to scale the input to unit variance, and then compute
#' gamma as num_samples * log(num_rows). This has the advantage of making each sample
#' contribute approximately equally to the segmentation.
#' @returns Numeric value of gamma computed as the sum of the column variances
#'  multiplied by the log of the number of data points
#' @export
estimate_gamma_bic <- function(x) {
    xmat <- as.matrix(x)
    sd_hat <- estimate_sd(xmat)
    log(nrow(xmat)) * sum(sd_hat^2)
}

#' Scale input data to unit variance.
#' @param x Numeric vector or matrix of PCF inputs
#' @returns Data of the same size and shape as the input, scaled so each column
#'   has unit variance. Vectors are treated as single column matrices.
#' @export
scale_pcf_input <- function(x) {
    orig_type_is_vector <- is.null(dim(x))
    xmat <- as.matrix(x)
    sd_hat <- estimate_sd(xmat)
    n <- sweep(xmat, 2, sd_hat, "/")
    if (orig_type_is_vector) {
        as.vector(n)
    } else {
        n
    }
}
