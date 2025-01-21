
options(error=recover)
source("C:/GREENBLATT/Rscripts/CrisprScreen/R/crispr_screen_analysis_lib.R")
setup()

run_deseq2_wrapper <- function(wd, samplef, groupf, exclude){
   res <- collect_sort_unsort_gene(wd, samplef, groupf, exclude)
   sort_unsort_deseq2(res$sorted, res$unsorted, terminator=res$subject)
}

# CPEB2 ####
wd <- c("C:/GREENBLATT/Zuyao/crispr/Nov26_2019/CPEB2/Analysis_Nov2023")
setwd(wd)

exclude <- c("CPEB2_Ref_1", "CPEB2_Ref_2", "CPEB2_Ref_3", "CPEB2_Ref_4", "CPEB2_Ref_5", "CPEB2_Ref_6", "CPEB2_Ref_7", "CPEB2_Ref_8", "CPEB2_Sort_5", "CPEB2_Sort_6", "CPEB2_Sort_7", "CPEB2_Sort_8", "CPEB2_UnSort_5", "CPEB2_UnSort_6", "CPEB2_UnSort_7", "CPEB2_UnSort_8")
samplef <- c(2,4)
groupf <- c(2,3)

preprocess_mageck_terminator(wd, samplef, exclude)

run_deseq2_wrapper(wd, samplef, groupf, exclude)


# CCND2 ####
wd <- c("C:/GREENBLATT/Zuyao/crispr/Mar02_2020/CCND2/Analysis_Nov2023")
setwd(wd)

exclude <- NULL
samplef <- c(4,6)
groupf <- c(4,5)

preprocess_mageck_terminator(wd, samplef, exclude)

run_deseq2_wrapper(wd, samplef, groupf, exclude)
