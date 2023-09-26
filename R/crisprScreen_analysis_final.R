

source("./crispr_screen_analysis_lib.R")

wds <- c("path_to_wd1",
         "path_to_wd2",
         "path_to_wd3")
for(i in seq_along(wds)){
  wd <- wds[i]
  res <- collect_sort_unsort_data(wd)
  sort_unsort_deseq2(res$sorted, res$unsorted)
}

