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
