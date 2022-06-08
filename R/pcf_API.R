#' Segment a single sample using the PCF algorithm.
#' @param x Vector of observations
#' @param kmin Minimum number of data points permitted in a single segment
#' @param gamma Penalty parameter. Higher values produce fewer segments.
#' @export
pcf <- function(x, kmin, gamma) {
    if (any(is.na(x) | is.nan(x) | is.infinite(x))) {
        stop("Data is not permitted to have NA, NaN or infinite values")
    }
    
    if (length(x) < 24) {
        stop("Can't segment data this small")
    }
    
    if(kmin > (length(x) / 3) - 6) {
        warning("kmin is too large compared to the data size, so setting kmin=", as.integer(floor((length(x) / 3) - 6)))
        kmin <- (length(x) / 3) - 6
        #return (list(nIntervals = 0))
    }
    
    if (kmin < 1) {
        warning("Minimum allowed kmin value is 1, so setting kmin=1")
    }
    
    sd <- getMad(x, 25)
    adjusted_gamma <- gamma * sd * sd
    if (length(x) < 1000) {
        result <- runFastPcfCpp(x, kmin, adjusted_gamma, 0.15, 0.15, FALSE)
    } else if (length(x) < 15000) {
        result <- runFastPcfCpp(x, kmin, adjusted_gamma, 0.12, 0.05, FALSE)
    } else {
        result <- runPcfSubsetCpp(x, kmin, adjusted_gamma, 0.12, 0.05, FALSE)
    }
    result
}

#' Segment multiple samples simultaneously using the MultiPCF algorithm.
#' @param x Matrix of observations (samples arranged as columns)
#' @param gamma Penalty parameter. Higher values produce fewer segments.
#' @param w Optional vector of weights to apply to the data before segmenting.
#' Allows down-weighting of noisy samples to avoid them contributing excessively
#' and causing the result to be over-segmented. Must be a numeric vector the same
#' length as the number of columns of x.
#' @export
multipcf <- function(x, gamma, w = NULL) {
    if (any(is.na(x) | is.nan(x) | is.infinite(x))) {
        stop("Data is not permitted to have NA, NaN or infinite values")
    }
    
    if (nrow(x) < 15) {
        stop ("Data is too small to segment")
    }
    
    if (is.null(w)) w <- rep(1.0, ncol(x))
    stopifnot(is.numeric(w))
    if(!(length(w) == ncol(x))) {
        stop ("Weights vector wrong size - must have length equal to number of columns of x")
    }
    
    # Scale up the gamma by multiplying by the number of samples.
    # This makes gamma values provided to pcf and multipcf have similar effects
    # on the amount of segmentation.
    gamma <- gamma * ncol(x)
    
    sd <- apply(x, 2, getMad)
    x <- sweep(x, 2, sd, "/")
    x <- sweep(x, 2, w, "*")
    nrow_x <- nrow(x)
    if (nrow_x < 1000) {
        result <- runFastMultiPCFCpp(x, gamma, 15, 0.15, 0.15, FALSE)
    } else {
        if (nrow_x < 3000) {
            result <- runFastMultiPCFCpp(x, gamma, 15, 0.12, 0.1, FALSE)
        } else {
            if (nrow_x < 15000) {
                result <- runFastMultiPCFCpp(x, gamma, 15, 0.12, 0.05, FALSE)
            } else {
                result <- runMultiPcfSubsetCpp(x, gamma, 15, 0.12, 0.05, FALSE)
            }
        }
    }
    result$mean <- sweep(result$mean, 1, w, "/")
    result$mean <- sweep(result$mean, 1, sd, "*")
    result
}

#' PCF algorithm, but the algorithm can only place breakpoints at the indices
#' listed in `breakpoints`
#' @export
restricted_pcf <- function(x, breakpoints, kmin, gamma) {
    sd <- getMad(x, 25)
    adjusted_gamma <- gamma * sd * sd
    
    mark <- vector(length = length(x))
    breakpoints <- breakpoints[breakpoints > 0 & breakpoints <= length(x)]
    mark[breakpoints-1] <- TRUE
    mark[length(mark)] <- TRUE
    dense <- compactCpp(x, mark)
    result <- PottsCompactCpp(kmin, adjusted_gamma, dense$Nr, dense$Sum,
                              FALSE)
    return (result)
}

#' Reimplementation of exact PCF
#' @export
exact_pcf <- function(x, kmin, gamma) {
    sd <- getMad(x, 25)
    adjusted_gamma <- gamma * sd * sd
    starts <- exact_pcf_(x, kmin, adjusted_gamma) + 1
    ends <- c(starts[starts > 1] - 1, length(x))
    lengths <- ends - starts + 1
    means <- sapply(seq_along(starts), function(i) {
        mean(x[starts[i]:ends[i]])
    })
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

exact_multipcf <- function(x, kmin, gamma, w = NULL) {
    sds <- apply(x, 2, getMad, k = 25)
    x_ <- sweep(x, 2, sds, "/")
    if (!is.null(w)) {
        x_ <- sweep(x_, 2, w, "*")
    }
    adjusted_gamma <- ncol(x) * gamma
    starts <- exact_multipcf_(x_, kmin, adjusted_gamma) + 1
    ends <- c(starts[starts > 1] - 1, nrow(x))
    lengths <- ends - starts + 1
    means <- t(sapply(seq_along(starts), function(i) {
        colMeans(x[starts[i]:ends[i], ])
    }))
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

#' @export
mark_positions <- function(x, nmad = 1, filter_size = 4) {
    if (is.matrix(x)) {
        mark <- sort(Reduce(union,
                            apply(x, 2, mark_,
                                  nmad = nmad,
                                  filter_size = filter_size)))
    } else {
        mark <- mark_(x, nmad, filter_size)
    }
    mark + 1
}

#' Reimplementation of fast PCF
#' @export
fast_pcf <- function(x, mark, kmin, gamma) {
    mark <- mark - 1
    sd <- getMad(x, 25)
    adjusted_gamma <- gamma * sd * sd
    
    starts <- fast_pcf_(x, mark, kmin, adjusted_gamma) + 1
    ends <- c(starts[starts > 1] - 1, length(x))
    lengths <- ends - starts + 1
    means <- sapply(seq_along(starts), function(i) {
        mean(x[starts[i]:ends[i]])
    })
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

#' @export
fast_multipcf <- function(x, mark, kmin, gamma, w = NULL) {
    mark <- mark - 1
    sds <- apply(x, 2, getMad, k = 25)
    x_ <- sweep(x, 2, sds, "/")
    if (!is.null(w)) {
        x_ <- sweep(x_, 2, w, "*")
    }
    adjusted_gamma <- ncol(x) * gamma
    starts <- fast_multipcf_(x_, mark, kmin, adjusted_gamma) + 1
    ends <- c(starts[starts > 1] - 1, nrow(x))
    lengths <- ends - starts + 1
    means <- t(sapply(seq_along(starts), function(i) {
        colMeans(x[starts[i]:ends[i], ])
    }))
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

#' @export
expanding_fast_pcf <- function(x, mark, kmin, gamma) {
    mark <- mark - 1
    if (length(x) < 10000) return (fast_pcf(x, mark, kmin, gamma))
    
    sd <- getMad(x, 25)
    adjusted_gamma <- gamma * sd * sd
    
    slices <- pmin(length(x), seq(5000, length(x) + 5000 - 1, 5000))
    bks <- vector("integer")
    for (slice_end in slices) {
        bks <- fast_pcf_(x[1:slice_end],
                         c(bks, mark[mark >= slice_end - 5000 &
                                         mark < slice_end]),
                         kmin,
                         adjusted_gamma)
    }
    
    starts <- bks + 1
    ends <- c(starts[starts > 1] - 1, length(x))
    lengths <- ends - starts + 1
    means <- sapply(seq_along(starts), function(i) {
        mean(x[starts[i]:ends[i]])
    })
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}

#' @export
expanding_fast_multipcf <- function(x, mark, kmin, gamma, w = NULL) {
    mark <- mark - 1
    if (nrow(x) < 10000) return (fast_multipcf(x, mark, kmin, gamma))
    
    sds <- apply(x, 2, getMad, k = 25)
    x_ <- sweep(x, 2, sds, "/")
    if (!is.null(w)) {
        x_ <- sweep(x_, 2, w, "*")
    }
    adjusted_gamma <- ncol(x) * gamma
    
    slices <- pmin(nrow(x), seq(5000, nrow(x) + 5000 - 1, 5000))
    bks <- vector("integer")
    for (slice_end in slices) {
        bks <- fast_multipcf_(x_[1:slice_end, ],
                              c(bks, mark[mark >= slice_end - 5000 &
                                              mark < slice_end]),
                              kmin,
                              adjusted_gamma)
    }
    
    starts <- bks + 1
    ends <- c(starts[starts > 1] - 1, nrow(x))
    lengths <- ends - starts + 1
    means <- t(sapply(seq_along(starts), function(i) {
        colMeans(x[starts[i]:ends[i], ])
    }))
    list(starts = starts, ends = ends, lengths = lengths, means = means)
}
