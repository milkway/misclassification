library(tidyverse)

segments_class <- tibble(
  City = c("Goiana", "Goiana", "Goiana", "Santos Dumont", "Santos Dumont", "Santos Dumont"),
  Strata = c("Highly Cropland", "Cropland", "Non-Cropland", "Highly Cropland", "Cropland", "Non-Cropland"),
  Frame = c(314, 232, 376, 31, 288, 1078),
  Sample = c(42, 15, 3, 29, 16, 15)
) %>% 
  pivot_longer(-City:-Strata, names_to = "Source", values_to = "Total")


glimpse(segments_class)


