# Plotting --------------------------------------------------------------------     

# Taxa names will have the format `Genus_species` in data frames and matrices.
# For axes labels in plots, we instead want the format _G. species_.  For plots
# of ratios of pairs of taxa, we want the `Genus1_species1:Genus2_species2`
# label to be rendered as _G. species1 : G. species2_. We can make these
# conversions with the function
#
# See https://stackoverflow.com/questions/26781676/align-legend-text-in-ggplot
# for background on fix to legend alignment issue; use bquote without making an
# expression.

#' Formatter for taxon names in the initial form of "Genus_species"
#'
#' @export
tax_labeller <- function (name) {
    name %>%
        str_replace_all("[a-z]+_", "\\. ") %>%
        str_replace(":", " : ") %>%
        map( function (x) bquote(italic(.(x))) )
}

# To get a prettier formatting of log-scale axes, use the following labeller
# function:

#' Formatter for log axes
#'
#' @export
log_formatter <- function (x) {
    format(x, scientific = 1,
        trim = TRUE, justify  = "none", 
        drop0trailing = TRUE)
}


#' Extract the x and y values from a ggplot object
#'
#' @export
xy_values <- function(ggplot) {
    tibble(
        x = rlang::eval_tidy(ggplot$mapping$x, data = ggplot$data),
        y = rlang::eval_tidy(ggplot$mapping$y, data = ggplot$data)
        )
}

