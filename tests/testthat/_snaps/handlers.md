# CLR transformation returns expected output

    Code
      urt_clr$clr_counts %>% arrange(sample_id, taxon_id)
    Output
      # A tibble: 3,272 x 3
         sample_id taxon_id count
         <chr>     <chr>    <dbl>
       1 s1        t11       3.07
       2 s1        t14       3.81
       3 s1        t146      1.25
       4 s1        t2        3.50
       5 s1        t30       3.05
       6 s1        t7        4.68
       7 s10       t1        4.20
       8 s10       t11       2.66
       9 s10       t14       1.71
      10 s10       t16       2.65
      # ... with 3,262 more rows

