
get_bar_plot <- function(ta, x = sample_clustered) {

  # convert promise to formula
  x <- substitute(x)

  # add sample_clustered if not present
  if (! "sample_clustered" %in% names(ta$samples)) {
    ta <- add_sample_clustered(ta)
  }

  # add taxon_name_color if not present
  if (! "taxon_name_color" %in% names(ta$taxa)) {
    ta <- add_taxon_name_color(ta)
  }

  # make plot and return
  get_abundances_extended(ta) %>%
    ggplot(aes_(x = x, y = ~abundance, fill = ~taxon_name_color)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette = "Paired", name = "Taxon") +
    xlab("sample") + ylab("relative abundance") +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white')
    )

}

# needs lib_sizes table
get_history_plot <- function(ta, col = NULL) {

  # convert promise to formula
  col <- substitute(col)

  # remove lib_size if present
  ta$samples$lib_size <- NULL

  # make plot and return
  ta$lib_sizes %>%
    left_join(ta$samples) %>%
    ggplot(aes_(x = ~step, y = ~lib_size, group = ~sample, col = col)) +
    geom_line(size = 0.5) +
    scale_y_log10() +
    scale_color_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1),
          panel.background = element_rect(fill = "white", colour = "black"))

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
  "#b15928"  # dark brown
)
