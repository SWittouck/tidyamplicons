
add_rel_abundance <- function(ta) {

  # add relative abundance to abundance table
  ta$abundances <- ta$abundances %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = abundance / sum(abundance)) %>%
    ungroup()

  # return ta object
  ta

}
