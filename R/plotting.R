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

#' Abbreviate taxon names from "Genus_species" to "Gsp"
#'
#' @export
tax_abbrev <- function (taxa) {
    m <- str_match(taxa, "([A-Z])[a-z]+_([a-z]{2})[a-z]*")
    ifelse(is.na(m[,1]), taxa, paste0(m[,2], m[,3]))
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

# TODO: import ggplot functions

#'
#'
#' Inspired by the `boxplot()` method for the `acomp` object of the
#' `compositions` package.
#'
#' @export
plot_ratios <- function(.data) {

    samples <- rownames(.data)
    taxa <- colnames(.data)
    K <- length(taxa)

    ratios.cen <- pw_ratios(center(.data)) %>% 
        enframe("Pair", "Center") %>%
        separate(Pair, c("Taxon.x", "Taxon.y"), ":")

    ratios <- pw_ratios_matrix(.data) %>%
        filter(!is.na(Ratio)) %>%
        separate(Pair, c("Taxon.x", "Taxon.y"), ":")

    # Use facet_grid() to view the ratios for each pair of taxa
    p <- ggplot(ratios, aes(x = factor(0), y = Ratio)) +
        geom_hline(yintercept = 1, color = "grey") +
        # geom_boxplot() +
        ggbeeswarm::geom_quasirandom(groupOnX = TRUE) +
        geom_point(data = ratios.cen, aes(y = Center), shape = 3, size = 3,
            color = "darkred") +
        scale_y_log10() +
        facet_grid(Taxon.x ~ Taxon.y) +
        theme_grey() +
        theme(
            # Remove background + grid lines
            # panel.background = element_blank(),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.line = element_line(colour = "black"),
            # Remove x axis
            axis.line.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank(),
            # Remove facet strip labels
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank()
        )

    # 
    g <- ggplotGrob(p)
    # Replace the grobs corresponding to diagonal panels with grobs that have the
    # that taxon's name
    grobs_idx <- g$layout$name %in% paste0("panel-", seq(K), "-", seq(K))
    label_grobs <- map(taxa, ~ ggplotGrob(cowplot::ggdraw() +
            cowplot::draw_label(.)))
    names(label_grobs) <- paste0("panel-", seq(K), "-", seq(K))
    # Replace and draw
    g$grobs[grobs_idx] <- label_grobs
    grid::grid.draw(g)
}
