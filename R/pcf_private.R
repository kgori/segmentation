###############
# Single PCF

runPcfSubsetCpp <- function (x, kmin, gamma, frac1 = 0.12, frac2 = 0.05, yest = FALSE)
{
    SUBSIZE <- 5000
    antGen <- length(x)
    mark <- filterMarkCpp(x, kmin, 8, 1, frac1, frac2, 0.02,
                          0.9)
    markInit <- c(mark[1:(SUBSIZE - 1)], TRUE)
    compX <- compactCpp(x[1:SUBSIZE], markInit)
    mark2 <- rep(FALSE, antGen)
    mark2[1:SUBSIZE] <- markWithPottsCpp(kmin, gamma, compX$Nr,
                                         compX$Sum, SUBSIZE)
    mark2[4 * SUBSIZE/5] <- TRUE
    start <- 4 * SUBSIZE/5 + 1
    while (start + SUBSIZE < antGen) {
        slutt <- start + SUBSIZE - 1
        markSub <- c(mark2[1:(start - 1)], mark[start:slutt])
        markSub[slutt] <- TRUE
        compX <- compactCpp(x[1:slutt], markSub)
        mark2[1:slutt] <- markWithPottsCpp(kmin, gamma, compX$Nr,
                                           compX$Sum, slutt)
        start <- start + 4 * SUBSIZE/5
        mark2[start - 1] <- TRUE
    }
    markSub <- c(mark2[1:(start - 1)], mark[start:antGen])
    compX <- compactCpp(x, markSub)
    result <- PottsCompactCpp(kmin, gamma, compX$Nr, compX$Sum,
                              yest)
    return(result)
}

runFastPcfCpp <- function (x, kmin, gamma, frac1= 0.12, frac2 = 0.05, yest = FALSE)
{
    antGen <- length(x)
    mark <- filterMarkCpp(x, kmin, 8, 1, frac1, frac2, 0.02,
                          0.9)
    mark[antGen] = TRUE
    dense <- compactCpp(x, mark)
    result <- PottsCompactCpp(kmin, gamma, dense$Nr, dense$Sum,
                              yest)
    return(result)
}

##############
# MultiPCF

runFastMultiPCFCpp <- function (x, gamma, L, frac1, frac2, yest) 
{
    mark <- rep(0, nrow(x))
    mark <- sawMarkMCpp(x, L, frac1, frac2)
    dense <- compactMultiCpp(t(x), mark)
    compPotts <- multiPCFcompactCpp(dense$Nr, dense$Sum, gamma)
    if (yest) {
        potts <- expandMultiCpp(nrow(x), ncol(x), compPotts$Lengde, 
                                compPotts$mean)
        return(list(pcf = potts, length = compPotts$Lengde, 
                    start0 = compPotts$sta, mean = compPotts$mean, nIntervals = compPotts$nIntervals))
    }
    else {
        return(list(length = compPotts$Lengde, start0 = compPotts$sta, 
                    mean = compPotts$mean, nIntervals = compPotts$nIntervals))
    }
}

runMultiPcfSubsetCpp <- function (x, gamma, L, frac1, frac2, yest)
{
    SUBSIZE <- 5000
    antGen <- nrow(x)
    
    mark <- sawMarkMCpp(x, L, frac1, frac2)
    
    markInit <- c(mark[1:(SUBSIZE - 1)], TRUE)
    compX <- compactMultiCpp(t(x[1:SUBSIZE, ]), markInit)
    
    mark2 <- rep(FALSE, antGen)
    mark2[1:SUBSIZE] <- markMultiPottsCpp(compX$Nr, compX$Sum,
                                          gamma, SUBSIZE)
    
    mark2[4 * SUBSIZE/5] <- TRUE
    start0 <- 4 * SUBSIZE/5 + 1
    times <- 0
    while (start0 + SUBSIZE < antGen) {
        times <- times + 1
        slutt <- start0 + SUBSIZE - 1
        markSub <- c(mark2[1:(start0 - 1)], mark[start0:slutt])
        markSub[slutt] <- TRUE
        compX <- compactMultiCpp(t(x[1:slutt, ]), markSub)
        
        mark2[1:slutt] <- markMultiPottsCpp(compX$Nr, compX$Sum,
                                            gamma, slutt)
        start0 <- start0 + 4 * SUBSIZE/5
        mark2[start0 - 1] <- TRUE
    }
    markSub <- c(mark2[1:(start0 - 1)], mark[start0:antGen])
    compX <- compactMultiCpp(t(x), markSub)
    compPotts <- multiPCFcompactCpp(compX$Nr, compX$Sum, gamma)
    
    if (yest) {
        potts <- expandMultiCpp(nrow(x), ncol(x), compPotts$Lengde,
                                compPotts$mean)
        return(list(pcf = potts, length = compPotts$Lengde,
                    start0 = compPotts$sta, mean = compPotts$mean, nIntervals = compPotts$nIntervals))
    }
    else {
        return(list(length = compPotts$Lengde, start0 = compPotts$sta,
                    mean = compPotts$mean, nIntervals = compPotts$nIntervals))
    }
}


###########################################################
# These functions taken verbatim from 'copynumber' package,
# license Artistic 2.0

#' @importFrom "stats" runmed
medianFilter <- function (x, k) 
{
    n <- length(x)
    filtWidth <- 2 * k + 1
    if (filtWidth > n) {
        if (n == 0) {
            filtWidth <- 1
        }
        else if (n%%2 == 0) {
            filtWidth <- n - 1
        }
        else {
            filtWidth <- n
        }
    }
    runMedian <- runmed(x, k = filtWidth, endrule = "median")
    return(runMedian)
}

#' @importFrom "stats" mad
getMad <- function (x, k = 25) {
    x <- x[x != 0]
    runMedian <- medianFilter(x, k)
    dif <- x - runMedian
    SD <- mad(dif)
    return(SD)
}
