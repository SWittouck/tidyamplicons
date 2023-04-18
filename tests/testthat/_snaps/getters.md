# Correct numbers reported

    Code
      ta_test %>% report_numbers()
    Message
      samples: 217
      taxa: 1957
      reads: 3873478

---

    Code
      ta_test %>% numbers()
    Output
      n_samples    n_taxa   n_reads 
            217      1957   3873478 

# Perform adonis shows stable output

    Code
      adn
    Output
      Permutation test for adonis under reduced model
      Terms added sequentially (first to last)
      Permutation: free
      Number of permutations: 999
      
      adonis2(formula = as.formula(paste("abundances_matrix", formula_RHS, sep = " ~ ")), data = metadata, permutations = permutations)
                Df SumOfSqs      R2      F Pr(>F)    
      location   1    2.800 0.03823 8.4267  0.001 ***
      method     1    0.334 0.00456 1.0056  0.409    
      Residual 211   70.118 0.95721                  
      Total    213   73.253 1.00000                  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

