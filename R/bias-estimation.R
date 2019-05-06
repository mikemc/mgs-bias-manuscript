#' Estimate bias
#'
#' The "center" method should only be used if all samples contain all taxa, and
#' in that case should give the same results as "rss"
#'
#' @param tb - data frame with Sample, Taxon, Observed, and Actual
#' @param method - "rss" (default) or "center"
#' @param bound - lower and upper bound on the log efficiencies for the "rss"
#' method
#'
#' @export
estimate_bias <- function(tb, method = "rss", bound = 10) {
    if (method == "rss") {
        # `rss` computes the residual sum of squares of the
        # log-ratio-transformed error across all samples
        rss <- function(log_bias, taxa, tb0) {
            bias <- tibble(Taxon = taxa, Bias = c(log_bias, 0))
            tb0 %>%
                left_join(bias, by = "Taxon") %>%
                mutate_by(Sample,
                    Residual = Error - Bias,
                    # Center the residuals within each sample (needed since
                    # Bias isn't centered)
                    Residual = Residual - mean(Residual)
                ) %>%
                {sum(.$Residual^2)}
        }

        taxa <- tb$Taxon %>% unique
        K <- length(taxa)

        # Rather than repeatedly take logs, we'll work with the clr'd error in
        # each sample. The clr is computed within each sample _after_
        # restricting to taxa actually in the sample.
        tb0 <- tb %>%
            filter(Actual > 0) %>%
            mutate(Error = Observed / Actual) %>%
            select(Sample, Taxon, Error) %>%
            mutate_by(Sample, Error = clr(Error))
        # initial guess for the log-bias vector relative to the last taxon
        initial <- rep(0, K-1)
        # Estimate the log bias
        log_bias_fit <- optim(initial, fn = rss, taxa = taxa, tb = tb0,
            method = "L-BFGS-B", lower = -bound, upper = bound)
        # Transform back to bias and put into a dataframe
        bias_fit <- exp(c(log_bias_fit$par, 0))
        bias_tb <- tibble(Taxon = taxa, Bias = bias_fit)
    } else if (method == "center") {
        bias_tb <- tb %>%
            group_by(Taxon) %>%
            summarize(Bias = gm_mean(Observed / Actual))
    }
    bias_tb %>% 
        mutate(Bias = center_elts(Bias)) %>%
        rename(Bias_est = Bias) %>%
        arrange(Taxon)
}

# Center (Compositional mean) -------------------------------------------------

# TODO: consider doing the matrix conversion and log in/out conversion in
# center, to reduce code replication

#' Compute the center (compositional mean) of a set of compositions
#'
#' @param data - data frame or matrix with taxa as columns
#' @param weights - sample (row) weights
#' @param method - "proj", "gm", or "rss"
#' @param in_scale - "linear" (default) or "log"
#' @param out_scale - "linear" (default) or "log"
#' @param bound - single number giving the lower and upper bound on the alr
#' efficiencies ("rss" only) 
#'
#' @export
center <- function(data, weights = rep(1, nrow(data)), method = "proj",
    in_scale = "linear", out_scale = "linear", bound = 10) {
    if (!(in_scale %in% c("linear", "log")))
        stop('`in_scale` must be "linear" or "log"')
    if (!(out_scale %in% c("linear", "log")))
        stop('`out_scale` must be "linear" or "log"')
    if (!(method %in% c("proj", "gm", "rss")))
        stop('`method` is invalid')

    if (method == "proj") {
        center.proj(data, weights, in_scale, out_scale)
    } else if (method == "gm") {
        center.gm(data, weights, in_scale, out_scale)
    } else if (method == "rss") {
        center.rss(data, weights, in_scale, out_scale, bound = bound)
    }
}

center.gm <- function(data, weights = rep(1, nrow(data)),
    in_scale = "linear", out_scale = "linear") {

    mat <- data %>% as("matrix")
    if (in_scale == "linear")
        mat <- log(mat)

    clr_bias_est <- mat %>%
        apply(2, weighted.mean, weights) %>%
        {. - mean(.)}

    # Format and return results
    if (out_scale == "log")
        return(clr_bias_est)
    else if (out_scale == "linear")
        return(exp(clr_bias_est))
    else 
        stop("Invalid `out_scale`")
}

# Method from vandenBoogaart2006; described in vandenBoogaart2013 and at
# https://core.ac.uk/download/pdf/132548286.pdf
center.proj <- function(data, weights = rep(1, nrow(data)),
    in_scale = "linear", out_scale = "linear") {

    proj_mat <- function(K, M) {
        mat <- diag(nrow = K) - 1/(K - length(M))
        mat[M,] <- 0
        mat[,M] <- 0
        mat
    }

    # Proper normalization appears to handled by the ginv matrix, so is
    # unnecessary to normalize the weights.

    mat <- data %>% as("matrix")
    K <- ncol(mat)
    if (in_scale == "linear")
        mat <- log(mat)

    P_sum <- diag(0, nrow = K)
    v_sum <- rep(0, K)
    for (i in seq(nrow(mat))) {
        M <- which(is.nan(mat[i,]))
        v <- mat[i,] %>% replace(M, 0)
        P <- proj_mat(K, M) * weights[i]
        P_sum <- P_sum + P
        v_sum <- v_sum + P %*% v
    }

    clr_bias_est <- (MASS::ginv(P_sum) %*% v_sum) %>% c
    names(clr_bias_est) <- colnames(mat)

    # Format and return results
    if (out_scale == "log")
        return(clr_bias_est)
    else if (out_scale == "linear")
        return(exp(clr_bias_est))
    else 
        stop("Invalid `out_scale`")
}

center.rss <- function(data, weights = rep(1, nrow(data)),
    in_scale = "linear", out_scale = "linear", bound = 10) {

    # This function computes the squared Aitchison norm of x on the
    # subcomposition of taxa that are observed.
    # 
    # x is a vector of the log compositional error of a sample
    row_rss <- function (x) {
        x %>%
            {. - mean(., na.rm = TRUE)} %>%
            {.^2} %>%
            sum(., na.rm = TRUE)
    }

    # mat is a matrix of the compositional errors on the log scale
    rss <- function(alr_bias, mat, weights) {
        # Substract the alr_bias from each row (sample)
        swept <- sweep(mat, 2, c(alr_bias, 0), "-")
        # Compute and sum the RSS's, weighting samples by `weights`
        apply(swept, 1, row_rss) %>%
            {. * weights} %>%
            sum
    }

    mat <- data %>% as("matrix")
    K <- ncol(mat)

    if (in_scale == "linear")
        mat <- log(mat)

    # Initial alr bias values to be passed to `par` arg of `optim()`
    initial <- rep(0, K-1)

    # Compute minimizer
    alr_bias_est <- optim(initial, fn = rss, mat = mat,
        method = "L-BFGS-B", weights = weights, lower = -bound, upper = bound)

    # Format and return results
    clr_bias_est <- c(alr_bias_est$par, 0) %>%
        {. - mean(.)}
    names(clr_bias_est) <- colnames(mat)
    if (out_scale == "log")
        return(clr_bias_est)
    else if (out_scale == "linear")
        return(exp(clr_bias_est))
    else 
        stop("Invalid `out_scale`")
}
