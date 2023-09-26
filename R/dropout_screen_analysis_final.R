
options(error=recover)

source("./crispr_screen_analysis_lib.R")

wds <- c("path_to_wd1",
         "path_to_wd2",
         "path_to_wd3")
hwds <- c("path_to_hwd1",
         "path_to_hwd2",
         "path_to_hwd3")

for(i in seq_along(wds)){
  wd <- wds[i]
  hwd <- hwds[i]
  scale_normalize(wd, hwd)
  qGI_analysis(wd, hwd, runGSEA = TRUE)
}
