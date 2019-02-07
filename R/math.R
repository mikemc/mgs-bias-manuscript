# TODO 
# - consider making trim and na_rm usage uniform
# - add better documentation to each function

# Element-wise functions ------------------------------------------------------

#' Convert proportions
odds <- function (x) {
    x / (1-x)
}

#' Logit (log-odds) of t
logit <- function (x) {
    log(x) - log(1-x)
}

# Summaries of relative abundance vectors--------------------------------------

#' Geometric mean of x.
g_mean <- function(x, na_rm = FALSE) {
    exp(mean(log(x), na.rm = na_rm))
}

#' Geometric standard deviation of x.
#'
#' Note, uses denominator n-1
g_sd <- function(x, na_rm = FALSE) {
    exp(sd(log(x), na.rm = na_rm))
}

#' Geometric range of x.
g_range <- function(x) {
    max(x) / min(x)
}

#' Shannon entropy of x.
entropy <- function (x) {
    x <- x[x>0]
    x <- x / sum(x)
    - sum(x * log(x))
}

# Transformations for relative abundance vectors ------------------------------

#' "Close" the elements of x to proportions
close_elts <- function (x, na_rm = FALSE) {
    x / sum(x, na.rm = na_rm)
}

#' Geometrically center the elements of x
center_elts <- function (x, na_rm = FALSE) {
    exp(log(x) - mean(log(x), na.rm = na_rm))
}

#' Centered log-ratio transform of x.
#' 
#' @param x Vector of abundances.
#' @param base Base for logarithm
clr <- function(x, base = exp(1), na_rm = FALSE) {
    log(x, base = base) - mean(log(x, base = base), na.rm = na_rm)
}

# Distance / dissimilarity between samples ------------------------------------

#' Distance or dissimilarity between relative abundance vectors x and y
#'
#' @param method : distance/dissimilarity measure
#' @param trim : should x and y be reduced to their common positive elements
#' before computing the Aitchison distance (otherwise, the distance will be
#' Inf)
#' 
#' method == "aitchison" -> Aitchison distance
#' method == "bray" -> Bray-Curtis dissimilarity between x and y. Note,
#' converts x and y to proportions before computing.
#'
# Bray method is equal to `vegan::vegdist(rbind(close_elts(x), close_elts(y)))[1]`
xydist <- function(x, y, method = "aitchison", trim = FALSE) {
    if (length(x) != length(y)) {
        stop("x and y have different lengths")
    }
    if (method == "aitchison") {
        if (trim) {
            idx <- (x > 0) & (y > 0)
            x <- x[idx]
            y <- y[idx]
        }
        sqrt( sum( (clr(x) - clr(y))^2 ) )
    } else if (method == "bray") {
        x <- close_elts(x)
        y <- close_elts(y)
        Cxy <- map2_dbl(x, y, min) %>% sum
        1 - Cxy
    }
}

#' Aitchison norm of x
#'
#' @param na_rm; remove NAs and NaNs before calculating
anorm <- function (x, na_rm = FALSE) {
    x %>%
        clr(., na_rm = na_rm) %>%
        {.^2} %>%
        sum(., na.rm = na_rm) %>%
        sqrt(.)
}

# Ratios ----------------------------------------------------------------------

#' Pairwise ratios of the elements of x
pw_ratios <- function (x) {
    tidyr::crossing(i = seq(x), j = seq(x)) %>%
        {purrr::map2_dbl(.$i, .$j, ~ x[.x] / x[.y])}
}
