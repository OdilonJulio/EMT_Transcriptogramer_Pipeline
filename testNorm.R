testNorm <-pca_result_R30[["pca_result"]][["rotation"]] [, "PC1"]

sum(testNorm*testNorm)
