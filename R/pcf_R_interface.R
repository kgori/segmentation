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
#' @export
multipcf <- function(x, gamma) {
    if (any(is.na(x) | is.nan(x) | is.infinite(x))) {
        stop("Data is not permitted to have NA, NaN or infinite values")
    }
    
    if (nrow(x) < 15) {
        stop ("Data is too small to segment")
    }
    
    # Scale up the gamma by multiplying by the number of samples.
    # This makes gamma values provided to pcf and multipcf have similar effects
    # on the amount of segmentation.
    gamma <- gamma * ncol(x)
    
    sd <- apply(x, 2, getMad)
    x <- sweep(x, 2, sd, "/")
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
    result$mean <- sweep(result$mean, 1, sd, "*")
    result
}
