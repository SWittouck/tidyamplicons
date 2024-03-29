#' Prepare tidyamplicons object for visualization by barplot.
#' Clusters samples, adds color groups and relative abundances.
#' @param ta a tidyamplicons object
#' @param n an integer
#' @noRd
prepare_for_bp <- function(ta, n = 12, extended = TRUE) {
  # add sample_clustered if not present
  if (!"sample_clustered" %in% names(ta$samples)) {
    ta <- add_sample_clustered(ta)
  }

  # add taxon_name_color if not present
  if (!"taxon_name_color" %in% names(ta$taxa)) {
    ta <- add_taxon_name_color(ta, n = n)
  }

  # add relative abundances if not present
  if (!"rel_abundance" %in% names(ta$abundances)) {
    ta <- add_rel_abundance(ta)
  }
  # optional extension (not used by sample bp)
  if (extended) {
    ta <- ta %>% get_abundances_extended()
  }
  ta
}

#' Return a bar plot of the samples
#'
#' @export
bar_plot <- function(ta, n = 12, x = sample_clustered, geom_bar = T) {
  # convert promise to formula
  x <- enquo(x)

  # make plot and return
  plot <- prepare_for_bp(ta, n) %>%
    ggplot(aes(
      x = forcats::fct_reorder(!!x, as.integer(sample_clustered)), 
      y = rel_abundance, fill = taxon_name_color)) +
    scale_fill_brewer(palette = "Paired", name = "Taxon") +
    xlab("sample") +
    ylab("relative abundance") +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white")
    )

  # add geom_bar if requested
  if (geom_bar) {
    plot <- plot + geom_bar(stat = "identity")
  }

  plot
}

#' Return an interactive bar plot of the samples
#' @param ta a tidyamplicons object
#' @param n an integer, representing the amount of colors used to depict different taxa
#' @param x a string, representing the column name used to label and cluster the samples on.
#'
#' @export
bar_plot_ly <- function(ta, n = 12, x = sample_clustered) {
  # convert promise to formula
  x <- enquo(x)

  # wrap in eval and quosure shenannigans
  plot <- rlang::eval_tidy(rlang::quo_squash(
    quo({
      # make plot and return
      prepare_for_bp(ta, n) %>%
        plotly::plot_ly(
          x = ~forcats::fct_reorder(!!x, as.integer(sample_clustered)),
          y = ~rel_abundance,
          color = ~taxon_name_color,
          colors = palette_xgfs,
          name = ~taxon_name_color,
          hovertemplate = paste(
            "<b>%{x}</b>",
            "<br>%{y:.2%}<br>"
          ),
          type = "bar"
        ) %>%
        plotly::layout(
          barmode = "stack",
          xaxis = list(title="Sample"),
          yaxis = list(title="Relative Abundance")
        )
    })
  ))
  plot
}

#' Return an interactive pcoa plot of the samples
#' @param ta a tidyamplicons object
#' @param x a string, representing the column name used to color the sample groups on.
#' @param palette a vector of colors, used as the palette for coloring sample groups.
#'
#' @export
pcoa_plot_ly <- function(ta, x, samplenames = sample, palette = NULL, title = "PCOA plot") {
  # convert promise to formula
  x <- enquo(x)
  samplenames <- enquo(samplenames)

  # fallback to default palette
  if (is.null(palette)) {
    palette <- cols
  }

  # prepare pcoa if needed
  if (!all(c("pcoa1", "pcoa2") %in% names(ta$samples))) {
    ta <- add_pcoa(ta)
  }

  plot <- rlang::eval_tidy(rlang::quo_squash(
    quo({
      ta$samples %>%
        plotly::plot_ly(
          x = ~pcoa1,
          y = ~pcoa2,
          color = ~!!x,
          colors = palette_paired,
          text = ~ paste(!!samplenames),
          hovertemplate = paste("<i>%{text}</i>"),
          type = "scatter",
          mode = "markers"
        ) %>%
        plotly::layout(
          title = title,
          yaxis = list(zeroline = F),
          xaxis = list(zeroline = F)
        )
    })
  ))

  plot
}


#' Return a bar plot of the samples
#'
#' DEPRECATED, use \code{\link{bar_plot}}
#'
#' @export
get_bar_plot <- bar_plot

#' Return a history plot of the samples
#'
#' DEPRECATED, this function is kept for historical reasons and will probably
#' not work
#'
#' @export
history_plot <- function(ta, col = NULL) {
  # convert promise to formula
  col <- substitute(col)

  # remove lib_size if present
  ta$samples$lib_size <- NULL

  # make plot and return
  ta$lib_sizes %>%
    left_join(ta$samples, by = "sample_id") %>%
    ggplot(aes_(x = ~step, y = ~lib_size, group = ~sample_id, col = col)) +
    geom_line(size = 0.5) +
    scale_y_log10() +
    scale_color_brewer(palette = "Paired") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1),
      panel.background = element_rect(fill = "white", colour = "black")
    )
}

#' Return a history plot of the samples
#'
#' DEPRECATED, this function is kept for historical reasons and will probably
#' not work
#'
#' @export
get_history_plot <- history_plot

#' Return a visualization designed for a small number of samples
#'
#' @export
sample_plot <- function(ta, sample = sample_id, n = 15, nrow = NULL) {
  ta <- prepare_for_bp(ta, n, extended = FALSE)

  sample <- rlang::enexpr(sample)
  if (sample != rlang::expr(sample_id)) {
    ta <- change_id_samples(ta, sample_id_new = !!sample)
  }

  data <-
    ta %>%
    everything() %>%
    group_by(sample_id) %>%
    arrange(desc(rel_abundance)) %>%
    slice(1:n) %>%
    ungroup() %>%
    arrange(sample_id, rel_abundance) %>%
    mutate(row = 1:n())

  data %>%
    ggplot(aes(x = row, y = rel_abundance, fill = taxon_name_color)) +
    geom_col() +
    facet_wrap(~sample_id, scales = "free", nrow) +
    coord_flip() +
    theme_bw() +
    scale_x_continuous(
      breaks = data$row,
      labels = data$taxon_name,
      expand = c(0, 0)
    ) +
    scale_fill_brewer(palette = "Paired", name = "taxon") +
    xlab("taxon name") +
    ylab("relative abundance")
}

palette_paired <- c(
  "#e8e8e8", # light grey
  "#a6cee3", # light blue
  "#1f78b4", # dark blue
  "#b2df8a", # light green
  "#33a02c", # dark green
  "#fb9a99", # light red
  "#e31a1c", # dark red
  "#fdbf6f", # light orange
  "#ff7f00", # dark orange
  "#cab2d6", # light purple
  "#6a3d9a", # dark purple
  "#ffff99", # light brown
  "#b15928" # dark brown
)

palette_xgfs <- c(
  # source: http://tsitsul.in/blog/coloropt/
  "#bdbdbd",
  "#00a76c",
  "#878500",
  "#00c6f8",
  "#5954d6",
  "#ff9287",
  "#b24502",
  "#d163e6",
  "#00bbad",
  "#006e00",
  "#008cf9",
  "#b80058",
  "#ebac23"
)
