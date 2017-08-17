
add_rel_abundance <- function(ta) {

  # add relative abundance to abundance table
  ta$abundances <- ta$abundances %>%
    group_by(sample) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup()

  # return ta object
  ta

}
