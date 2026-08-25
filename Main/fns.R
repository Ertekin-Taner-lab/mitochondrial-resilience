
# named_group_split ----
# this function split the dataframe based on your given variable and create a named list
named_group_split <- function(.tbl, ...) {
  require(tidyverse)
  
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

# meta_gen_fn ----
# this function conduct meta-analysis between two cohorts based on the effect estimates and standard error you input
meta_gen_fn = function(x, eff_idx1, eff_idx2, err_idx1, err_idx2, ...) {
  require(meta)
  
  eff = c(x[[eff_idx1]], x[[eff_idx2]]) %>% as.numeric()
  err = c(x[[err_idx1]], x[[err_idx2]]) %>% as.numeric()
  m = metagen(eff, err) %>% summary()
  return(m)
}

# get_max_contrast_color ----
# this function compute whether the text color should be black or white depending on the background color
get_max_contrast_color = function(background_color_vec, cut_off = 0.5) {
  .rgb_to_bw = function(color_chr) {
    rgb_vec = col2rgb(color_chr)
    luminance = 0
    # Some magic values from stackoverflow
    luminance = luminance + rgb_vec["red", 1, drop = TRUE] * 0.299
    luminance = luminance + rgb_vec["green", 1, drop = TRUE] * 0.587
    luminance = luminance + rgb_vec["blue", 1, drop = TRUE] * 0.114
    luminance = luminance / 255
    luminance = as.double(luminance)
    ifelse(luminance > cut_off, "black", "white")
  }
  background_color_vec %>%
    purrr::map_chr(.rgb_to_bw)
}