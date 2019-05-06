#' Get the upper left corner of a 2-d array.
#'
#' @export
corner <- function (x, r = 5L, c = 5L) {
    if (r > nrow(x)) {
        r <- nrow(x)
    }
    if (c > ncol(x)) {
        c <- nrow(x)
    }
    x[seq(r), seq(c)]
}


# Helpers for working with "tidy" microbiome data -----------------------------

#' Mutate within groups
#'
#' @export
mutate_by <- function(.data, group_vars, ...) {
    gvs <- rlang::enquos(group_vars)
    .data %>%
        group_by_at(vars(!!!gvs)) %>%
        mutate(...) %>%
        ungroup
}

#' Create a matrix from columns of a tidy data frame
#'
#' Note: Setting `fill` will replace NAs and NaNs, along with elements
#' corresponding to missing rows, with the value of `fill`.
#'
#' @param .data A data frame with columns `rows`, `cols`, and `elts`.
#' @param rows Column that will become the rownames of the matrix.
#' @param cols Column that will become the rownames of the matrix.
#' @param elts Column that will become the matrix elements.
#' @param fill Value to use for missing elements.
#' @export
build_matrix <- function(.data, rows, cols, elts, fill = NULL) {
    rows <- rlang::enquo(rows)
    cols <- rlang::enquo(cols)
    elts <- rlang::enquo(elts)
    tb <- .data %>% select(!!rows, !!cols, !!elts)
    if (!is.null(fill)) {
        tb <- tb %>% 
            complete(!!rows, !!cols, fill = rlang::list2(!!elts := fill))
    }
    tb <- tb %>% spread(!!cols, !!elts)
    mat <- tb %>% select(-!!rows) %>% as("matrix")
    rownames(mat) <- tb %>% pull(!!rows)
    mat
}
