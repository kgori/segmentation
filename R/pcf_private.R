#' Compute the residual score of a segmentation model 
score_model <- function(fit, data) {
    m <- as.matrix(data)
    calc <- meanCalculator$new(m)
    score <- list(free=0, int=0)
    for (i in seq_along(fit$starts)) {
        residual_free <- calc$ssr(fit$starts[i], fit$ends[i])
        residual_int <- calc$ssr_rounded_mean(fit$starts[i], fit$ends[i])
        score$free <- score$free + residual_free
        score$int <- score$int + residual_int
    }
    size <- length(fit$starts)
    append(lapply(score, sum), list(size=size))
}

#' Prefix sums for computing means and residuals
#' @importFrom "R6" R6Class
meanCalculator <- R6::R6Class(
    "meanCalculator",
    private=list(
        prefix_sums=NA,
        prefix_sum_of_squares=NA
    ),
    public=list(
        #' @param data Numerical vector or matrix
        initialize = function(data) {
            m <- as.matrix(data)
            self$prefix_sums <- rbind(rep(0, ncol(m)), apply(m, 2, cumsum))
            self$prefix_sum_of_squares <- rbind(rep(0, ncol(m)), apply(m^2, 2, cumsum))
        },
        
        #' @description
        #' Compute the sum of the values in the owned data
        #' for segment [start, end]
        #' @param start First index of the data segment
        #' @param end Last index of the data segment
        sum = function(start, end) {
            self$prefix_sums[end + 1,] - self$prefix_sums[start,]
        },
        
        #' @description
        #' Compute the mean of the values in the owned data
        #' for segment [start, end]
        #' @param start First index of the data segment
        #' @param end Last index of the data segment
        mean = function(start, end) {
            self$sum(start, end) / (end + 1 - start)
        },
        
        #' @description
        #' Compute the sum of squared residuals relative to the mean for
        #' segment [start, end]
        #' @param start First index of the data segment
        #' @param end Last index of the data segment
        ssr = function(start, end) {
            sum_of_squares <- self$prefix_sum_of_squares[end + 1,] - self$prefix_sum_of_squares[start,]
            square_of_sum <- self$sum(start, end)^2
            n <- end + 1 - start
            sum_of_squares - square_of_sum/n
        },
        
        #' @description
        #' Compute the sum of the squared residuals relative to an arbitrary
        #' target value for segment [start, end]
        #' @param start First index of the data segment
        #' @param end Last index of the data segment
        #' @param target Residuals are deviations from this target value
        residuals_from_target = function(start, end, target) {
            sum_of_squares <- self$prefix_sum_of_squares[end + 1,] - self$prefix_sum_of_squares[start,]
            sum <- self$sum(start, end)
            n <- end + 1 - start
            sum_of_squares - 2*target*sum + n*target*target
        },
        
        #' @description
        #' Compute the sum of the squared residuals relative to the integer
        #' rounded mean for segment [start, end]
        #' @param start First index of the data segment
        #' @param end Last index of the data segment
        ssr_rounded_mean = function(start, end) {
            target <- round(self$mean(start, end))
            self$residuals_from_target(start, end, target)
        }
    )
)

#' Collapse segments with the same integer mean state
merge_integers <- function(model) {
    rounded_means <- round(model$means)
    starts <- model$starts[1]
    ends <- model$ends[1]
    sums <- model$means[1,] * model$lengths[1]
    j <- 1
    for (i in seq_along(model$starts)) {
        if (i == 1) next
        if (all(rounded_means[i, ] == rounded_means[i-1, ])) {
            ends[j] <- model$ends[i]
            sums[j,] <- sums[j,] + model$means[i,] * model$lengths[i]
        } else {
            starts <- c(starts, model$starts[i])
            ends <- c(ends, model$ends[i])
            sums <- rbind(sums, model$means[i,] * model$lengths[i])
            j <- j + 1
        }
    }
    lengths <- ends - starts + 1
    means <- sweep(sums, 1, lengths, "/")
    dimnames(means) <- list(NULL, colnames(means))
    list(starts=starts, ends=ends, lengths=lengths, means=means)
}