source('R/import_data.R')

# region ---------
region_key <- ecoregions %>% st_set_geometry(NULL) %>% 
  dplyr::select(c(3,4,5,7)) %>% distinct() %>% 
  arrange(NA_L3NAME)
level3 <- matrix(0, 84, 84)
level2 <- matrix(0, 84, 84)
level1 <- matrix(0, 84, 84)
# rownames(corr_mat) <- region_key$NA_L3NAME
for(i in 1:84) {
  for(j in 1:84) {
    if (region_key[j, "NA_L3CODE"] == region_key[i, "NA_L3CODE"]) {
      level3[i, j] = 1 # i = j, diagonal of 1s
    } else if (region_key[j, "NA_L2CODE"] == region_key[i, "NA_L2CODE"]) {
      level1[i, j] = 1 # indicator for correlation at level 1
      level2[i, j] = 1 # indicator for correlation at level 2
    } else if (region_key[j, "NA_L1CODE"] == region_key[i, "NA_L1CODE"]) {
      level1[i, j] = 1 # indicator for correlation at level 1
    } else {
      level3[i, j] = 0
      level1[i, j] = 0
      level2[i, j] = 0
    }
  }
}
test <- 0.1*level1 + 0.5*level2 +level3
chol(test)
test[1:10, 1:10]
