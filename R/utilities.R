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


# Tidy helpers ---------------------------------------------------------------

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

# mutate_by_alt <- function(.data, group_var, ...) {
#     group_var <- rlang::enquo(group_var)
#     .data %>%
#         group_by(!!group_var) %>%
#         mutate(...) %>%
#         ungroup
# }

