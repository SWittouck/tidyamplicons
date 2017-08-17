
update_lib_sizes <- function(ta, step) {

  # current lib_size in tidy table
  lib_sizes_new <- ta$abundances %>%
    group_by(sample) %>%
    summarize(lib_size = sum(abundance)) %>%
    mutate(step = step)

  # make lib_sizes table if it doesn't exist
  if (is.null(ta$lib_sizes)) {
    ta$lib_sizes <- lib_sizes_new %>%
      mutate(step = factor(step))
    # update lib_sizes table if it already existed
  } else {
    levels <- levels(ta$lib_sizes$step)
    levels <- c(levels, step)
    ta$lib_sizes <- ta$lib_sizes %>%
      mutate(step = as.character(step)) %>%
      bind_rows(lib_sizes_new) %>%
      mutate(step = factor(step, levels = !! levels))
  }

  # return ta object
  ta

}
