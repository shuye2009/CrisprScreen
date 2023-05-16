library(pheatmap)
library(limma)
library(Rtsne)
library(psych)
library(DESeq2)
library(knitr)
library(plyr)
library(plotrix)
library(pcaMethods)

source("./R/crispr_screen_analysis_lib.R")


########### input gene name update map
color_store <- brewer.pal(n = 8, name = "Dark2")
name_map <- NULL
mapfile <- "C:/GREENBLATT/Rscripts/CrisprScreen/resources/hgnc_geneName_updates_2019.txt"
name_map <- read.table(mapfile, header=T, stringsAsFactors = FALSE)
head(name_map)
dim(name_map)


GOBP_file <- "C:/GREENBLATT/resource/GSEA_gmt/c5.go.bp.v7.5.1.symbols.gmt"
PATHWAY_file <- "C:/GREENBLATT/resource/GSEA_gmt/ReactomePathways.gmt"

corum <- read.delim("C:/GREENBLATT/Rscripts/CrisprScreen/resources/coreComplexes.txt")
idMap <- read.delim("C:/GREENBLATT/Rscripts/CrisprScreen/resources/Human_id_map.tab")

TKO3 <- read.delim("C:/GREENBLATT/Rscripts/CrisprScreen/resources/TKOv2.1-Human-Library.txt")
old <- TKO3[, "GENE"]
TKO3[, "GENE"] <- update_gene_names(name_map, TKO3[, "GENE"])
new <- TKO3[, "GENE"]

names(new) <- old

for(i in 1:length(old)){
   if(old[i] !=  new[i]){
      cat(paste(old[i], "\t", new[i], "\n"))
   }
}

TKO3_genes <- unique(TKO3[, "GENE"])
TKO3_guides <- unique(TKO3[, "SEQUENCE"])
dim(TKO3)
length(TKO3_genes)
length(TKO3_guides)

control_gene <- c("EGFP", "LacZ", "luciferase")


############################# Key parameters 
######## inputs come from: countReads_perGuide.pl $f F, run on BC cluster          
######## scale_by_sample <- FALSE:  normalized by total aligned reads per sample    
######## sizeFactor <- FA:  sizeFactore for each guide is determined with "iterate" 
######## sort_1 through sort_8: all 8 samples are included in the analysis         
######## output: sortALL_SFA_raw_upregulated_gene_table.tab is the final results       

## MAIN starts ####

wds <- c("C:/GREENBLATT/Nujhat/APAscreen/Nov28_2022/DICER1")

for (wd in wds[1]){

  setwd(wd)
  format <- "_geneRawCount.tab"
  guide_format <- "_guideRawCount.tab"
  terminator <- unlist(strsplit(wd, "/", fixed=T))[6]
  ## Load file paths and directory names
  files <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = format)
  guide_files <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = guide_format)
  length(files)
  files
  
  samples <- unlist(lapply(files, function(x)split_sample(x,1,3)))
  groups <- unlist(lapply(files, function(x)split_group(x,1,2)))
  
  data_list <- list()
  data_list_NA <- list()
  
  for(i in 1:length(files)){
   # i <- 1
    df_all <- read.table(files[i], header=T, sep="\t")
    df_all <- df_all[, c(2,4,6,8)]
    ### update gene names in row.names
    old_gene_names <- row.names(df_all)
    new_gene_names <- update_gene_names(name_map, old_gene_names)
    
    if( i == 1){
      name_df <- data.frame(old_gene_names, new_gene_names)
      name_diff <- name_df[!name_df$old_gene_names %in% name_df$new_gene_names, ]
      write.table(name_df, "updated_gene_list.tab", sep="\t")
    }
    
    row.names(df_all) <- new_gene_names
    
    df <- na.omit(df_all)
    
    df_NA <- df_all[!row.names(df_all) %in% row.names(df),]
    
    data_list[[samples[i]]] <- df
    data_list_NA[[samples[i]]] <- df_NA
  }
  
  data_df_NA <- lapply(data_list_NA, cbind)
  
  write.table(data_df_NA, "Less_than_4_guides_raw_count.tab", row.names=T, col.names=NA, sep="\t")
  #data_list[["Ref_1"]]["CPSF1",]
  #data_list[["UnSort_1"]]["CPSF1",]
  #data_list[["Sort_1"]]["CPSF1",]
  
  dim(data_list[[2]])
  
  for(i in 1:length(files)){
    print(samples[i])
    test_df <- data_list[[i]]
    test_df <- test_df[order(test_df[,3]),]
    #head(test_df)
    #print(tail(test_df))
    print(test_df[c("CPSF1","NUDT21", "FIP1", "WDR33"),])
  }
  
  guides <- paste0("guide", c(1:4))
  unsort_count_df <- NULL
  sort_count_df <- NULL 
  
  
  sort_samples <- samples[grep("_sorted_", samples)]
  unsort_samples <- samples[grep("_unsorted_", samples)]
  
  
  for(i in 1:length(sort_samples)){
    sort_sample <- sort_samples[i]
    unsort_sample <- unsort_samples[i]
    
    sort_col <- paste(sort_samples[i], guides, sep="_")
    unsort_col <- paste(unsort_samples[i], guides, sep="_")
    
    df1c <- data_list[[unsort_sample]]
    colnames(df1c) <- unsort_col
    if(is.null(unsort_count_df)){
      unsort_count_df <- df1c
    }else{
      unsort_count_df <- cbind(unsort_count_df, df1c)
    }
    
    df2c <- data_list[[sort_sample]]
    colnames(df2c) <- sort_col
    if(is.null(sort_count_df)){
      sort_count_df <- df2c
    }else{
      sort_count_df <- cbind(sort_count_df, df2c)
    }
    
  }
  
  sort_count_df <- sort_count_df[!rownames(sort_count_df) %in% c("NUDT21"),]
  unsort_count_df <- unsort_count_df[!rownames(unsort_count_df) %in% c("NUDT21"),]
  raw_count_df <- cbind(sort_count_df, unsort_count_df)
  dim(raw_count_df)
  
  plot_heatmap(raw_count_df, "raw_count")
  
  dsc <- sumstats_col(raw_count_df)
 
  write.table(dsc, "descriptive_stats_of_raw_count_perGuide.tab", row.names=T, col.names=NA, sep="\t")
 
  #boxplot.matrix(as.matrix(raw_count_df), outline=T)
  
  dsr <- sumstats_row(sort_count_df)
  dur <- sumstats_row(unsort_count_df)
  dr <- cbind(dsr, dur)
  write.table(dr, "descriptive_stats_of_sort_count_perGene.tab", row.names=T, col.names=NA, sep="\t")
  
  
  pdf(paste(terminator, "raw_count_boxplot.pdf", sep="_"), width=15, height=6)
  
  old.par <- par(mfrow=c(1,1),mar=c(10,4,2,2))
  
   boxplot(as.matrix(raw_count_df), outline=T, names=NA)
   staxlab(1,1:24,colnames(raw_count_df),srt=45)
  
  dev.off()
  par(old.par)
  
  
  sort_count_df <- as.matrix(sort_count_df)
  unsort_count_df <- as.matrix(unsort_count_df)
  head(sort_count_df)
  dim(sort_count_df)
  summary(sort_count_df)
  
  

  ################################ DEseq2 analysis ###########################
  
  if(1){
    
      
        sample_matrix <- round(cbind(sort_count_df, unsort_count_df))
        samples <- colnames(sample_matrix)
        replicates <- dim(sample_matrix)[2]/2
        conditions <- c(rep("sorted", replicates), rep("unsorted", replicates))
        s2c <- data.frame(path=samples, conditions, row.names=samples)
        #s2c$conditions <- relevel(s2c$conditions, ref="unsort")
        s2c$path <- as.character(s2c$path)
        #s2c <- s2c[order(sample_order),]
        s2c
        
        guide_stat <- sumstats_col(sample_matrix)
        guide_stat_sum <- guide_stat$Sum
        
        guide_gomean <- geometric.mean(guide_stat_sum)
        ddsMatrix <- DESeqDataSetFromMatrix(sample_matrix, 
                                            colData = s2c,
                                            design = ~ conditions)
        class(ddsMatrix)
        sizeFactors(ddsMatrix) <- guide_stat_sum/guide_gomean
        dds <- DESeq(ddsMatrix, fitType = "local")
        
        res <- results(dds, contrast=c("conditions", "sorted", "unsorted"), independentFiltering = FALSE, pAdjustMethod = "BH")
        
        contrasts <- matrix(c("conditions", "sorted", "unsorted"), nrow=1, byrow=TRUE)
        plot_DESeq2_results(dds, contrasts=contrasts, pvalue_cutoff=0.05, gsea=TRUE, verbose=TRUE)
        
        
        ## examine each guide and each sample based on zscore
        normalized_count_df <- counts(dds, normalized=TRUE)[row.names(res),]
        plot_heatmap(normalized_count_df, "normalized_count")
        pdf(paste(terminator, "normalized_count_boxplot.pdf", sep="_"), width=15, height=6)
        
        old.par <- par(mfrow=c(1,1),mar=c(10,4,2,2))
        
        boxplot(as.matrix(normalized_count_df), outline=T, names=NA)
        staxlab(1,1:24,colnames(normalized_count_df),srt=45)
        
        dev.off()
        par(old.par)
        
        sort_stats <- sumstats_row(normalized_count_df[,grep("_sorted_", colnames(normalized_count_df))])
        unsort_stats <- sumstats_row(normalized_count_df[,grep("_unsorted_", colnames(normalized_count_df))])
        
        sort_diff <- normalized_count_df[,grep("_sorted_", colnames(normalized_count_df))] - unsort_stats$Mean
        sort_zscore <- round(t(t(sort_diff)/unsort_stats$SD), 2)
        sig_guide_count <- apply(sort_zscore, 1, function(x)length(x[x>2]))  ## zscore > 2 
        
        ## for each guide, find the samples that are significant
        sig_sample_count <- list()
        sig_guide_sample_count <- list()
        for(guide in guides){
          guide_z <- sort_zscore[, grep(guide, colnames(sort_zscore))]
          print(colnames(guide_z))
          sig_sample_count[[guide]] <- apply(guide_z, 1, function(x)length(x[x>2]))
          header <- paste(guide, "samples", sep="_")
          sig_guide_sample_count[[header]] <- apply(guide_z, 1, function(x) paste(which(x>2), collapse=","))
        }
        sig_sample_count <- as.data.frame(sig_sample_count)
        sig_guide_sample_count <- as.data.frame(sig_guide_sample_count)
        
        ## for each sample, find the guides that are significant
        sig_guides <- list()
        sig_sample_guides_count <- list()
        for(samp in sort_samples){
          samp_z <- sort_zscore[, grep(samp, colnames(sort_zscore))]
          print(colnames(samp_z))
          sig_guides[[samp]] <- apply(samp_z, 1, function(x)length(x[x>2]))
          header <- paste(samp, "guides", sep="_")
          sig_sample_guides_count[[header]] <- apply(samp_z, 1, function(x) paste(which(x>2), collapse=","))
        }
        sig_guides <- as.data.frame(sig_guides)
        sig_sample_guides_count <- as.data.frame(sig_sample_guides_count)
        
        normalized_count_df_combined <- cbind(res, sig_guide_count, normalized_count_df, sort_stats, unsort_stats)
       
        allgene_zscore_combined <- cbind(res, sig_guide_count, sig_sample_count, sig_guide_sample_count, sig_guides, sig_sample_guides_count, sort_zscore)
        
        filename <- paste(terminator, "ALL_gene_normalized_data.tab", sep="_")
        write.table(normalized_count_df_combined, filename, row.names=T, col.names=NA, sep="\t")
        filename <- paste(terminator, "ALL_gene_zscore.tab", sep="_")
        write.table(allgene_zscore_combined, filename, row.names=T, col.names=NA, sep="\t")
        
  }
}


################################## guide ############################


################################## Characteristic Direction test #################################

if(0){
  library(GeoDE)
  data(AllGMTfiles)  ## load function annotations that comes with 'GeoDE'
  replacezero <- function(x) "[<-" (x, !x | is.na(x), min(x[x > 0], na.rm=TRUE)/2)
  
  sample_matrix <- replacezero(cbind(unsort_count_df, sort_count_df))
  guide_stat <- sumstats_col(sample_matrix)
  guide_stat_sum <- guide_stat$Sum
  
  guide_gomean <- geometric.mean(guide_stat_sum)
  sizeFactors <- guide_stat_sum/guide_gomean
  names(sizeFactors) <- colnames(sample_matrix)
  sizeNormalized_matrix <- scale(sample_matrix, center=FALSE, scale=sizeFactors)
  
  scaled_df <- t(scale(t(sizeNormalized_matrix), center=TRUE))
  genenames <- row.names(unsort_count_df)
  Data_ALL <- type.convert(as.data.frame(cbind(genenames, sizeNormalized_matrix)))
  head(Data_ALL)
  dim(Data_ALL)
  SC <- as.factor(c(rep(1,dim(unsort_count_df)[2]), rep(2,dim(sort_count_df)[2])))
  
  chdir_analysis <- chdirAnalysis(Data_ALL, sampleclass=SC, gammas = list(1), CalculateSig=FALSE,nnull=10)
  
  chdir_analysis$chdirprops$number_sig_genes
  
  chdir_analysis$results[[1]][1:100]
  
  multiPAEAtest <- multigmtPAEAAnalysis(chdir_analysis$chdirprops, AllGMTfiles, gammas = list(1))
  
  data(GeneOntology_BP.gmt)
  PAEAtest_GOBP <- PAEAAnalysis(chdir_analysis$chdirprops, gmt, gammas = c(1), casesensitive = FALSE, showprogress=TRUE)
  

}


################################## NOISeqBIO test #################################

if(0){
  library(NOISeq)
  
  replacezero <- function(x) "[<-" (x, !x | is.na(x), min(x[x > 0], na.rm=TRUE)/2)
  
  sample_matrix <- replacezero(cbind(unsort_count_df, sort_count_df))
  guide_stat <- sumstats_col(sample_matrix)
  guide_stat_sum <- guide_stat$Sum
  
  guide_gomean <- geometric.mean(guide_stat_sum)
  sizeFactors <- guide_stat_sum/guide_gomean
  names(sizeFactors) <- colnames(sample_matrix)
  sizeNormalized_matrix <- scale(sample_matrix, center=FALSE, scale=sizeFactors)
  
  myfactors = data.frame(treatment=c(rep("Unsort", 16), rep("Sort",16)))
  
  mydata <- readData(data = sizeNormalized_matrix, factors = myfactors)
  mydata
  
  mynoiseqbio = noiseqbio(mydata, k = 1, norm = "n", factor = "treatment",
                          lc = 0, r = 20, adj = 1.5, plot = T, a0per = 0, 
                          random.seed = 12345, filter = 1)
  res <- mynoiseqbio@results[[1]]
  res_pos <- res[res[,"theta"]>1,]
  head(res_pos[order(res_pos[,"theta"],decreasing=T),], n=50)
  
  mynoiseqbio@replicates
  mynoiseqbio@comparison
  mynoiseqbio@k
  head(mynoiseqbio@results[[1]])
  
  
  qcut = 0.8
  mynoiseqbio.deg = degenes(mynoiseqbio, q = qcut, M = NULL)
  mynoiseqbio.UP = degenes(mynoiseqbio, q = qcut, M = "up")
  mynoiseqbio.DOWN = degenes(mynoiseqbio, q = qcut, M = "down")
  
  DE.plot(mynoiseqbio, q = qcut, graphic = "expr", log.scale = TRUE)
  DE.plot(mynoiseqbio, q = qcut, graphic = "MD")
  
}


################################## EBseq test #################################

if(0){
  library(EBSeq)
  
  replacezero <- function(x) "[<-" (x, !x | is.na(x), min(x[x > 0], na.rm=TRUE)/2)
  
  sample_matrix <- replacezero(cbind(unsort_count_df, sort_count_df))
  
  Sizes=QuantileNorm(sample_matrix, 0.8)
  apply(sample_matrix, 2, quantile)
  
  guide_stat <- sumstats_col(sample_matrix)
  guide_stat_sum <- guide_stat$Sum
  
  guide_gomean <- geometric.mean(guide_stat_sum)
  sizeFactors <- guide_stat_sum/guide_gomean
  names(sizeFactors) <- colnames(sample_matrix)
  #sizeNormalized_matrix <- scale(sample_matrix, center=FALSE, scale=sizeFactors)
  
  EBOut=EBTest(Data=sample_matrix, 
               Conditions=as.factor(rep(c("Unsort", "Sort"), each=16)),
               sizeFactors=sizeFactors, 
               maxround=5,
               PoolLower = .25,
               PoolUpper = .75,
               ApproxVal = 10^-10)
  
  EBDERes=GetDEResults(EBOut, FDR=0.05, FDRMethod="hard",Threshold_FCRatio=0.2)
  str(EBDERes$DEfound)
  head(EBDERes$PPMat, n=30)
  
  GeneFC=PostFC(EBOut)
  str(GeneFC)
  as.data.frame(GeneFC)
  PlotPostVsRawFC(EBOut,GeneFC)
  
  par(mfrow=c(1,2))
  QQP(EBOut)
  DenNHist(EBOut)
  
  
  
  results_table <- cbind(as.data.frame(GeneFC), EBDERes$PPMat)
  
  results_table <- results_table[order(results_table[,"PostFC"], decreasing=T),]
  
  head(results_table, n=100)
  
  par(mfrow=c(1,1))
  sigpv <- results_table[(results_table$PPEE < 0.05),]
  colors <- rep("black", dim(results_table)[1])
  colors[row.names(results_table) %in% row.names(sigpv)] <- "red"
  plot(log(results_table$PostFC), log(-log(results_table$PPEE)), col=colors)
  
  
}





################################## t-test and wilcox test #################################

if(0){  ## NOT used, as they are inferior to DESeq2 
  ttest_results <- NULL
  wtest_results <- NULL
  for(gene in row.names(raw_count_sort_filtered_df)){
    res <- t.test(raw_count_sort_filtered_df[gene,], raw_count_unsort_filtered_df[gene,], paired=T)
    ttest_results <- rbind(ttest_results, c(res$estimate, res$statistic, res$p.value))
    
    wres <- wilcox.test(raw_count_sort_filtered_df[gene,], raw_count_unsort_filtered_df[gene,], paired=T, conf.int = TRUE)
    wtest_results <- rbind(wtest_results, c(wres$estimate, wres$statistic, wres$p.value))
  }
  
  row.names(ttest_results) <- row.names(raw_count_sort_filtered_df)
  colnames(ttest_results) <- c("diff", "t-stats", "p.value")
  ttest_results <- ttest_results[order(ttest_results[, "p.value"]),]
  ttest_results_pos <- ttest_results[ttest_results[, "diff"] > 0, ]
  
  write.table(ttest_results_pos, "filtered_samples_paired_t_test_results.tab", col.names=NA, sep="\t")
  
  
  row.names(wtest_results) <- row.names(raw_count_sort_filtered_df)
  colnames(wtest_results) <- c("diff", "w-stats", "p.value")
  wtest_results <- wtest_results[order(wtest_results[, "p.value"]),]
  wtest_results_pos <- wtest_results[wtest_results[, "diff"] > 0, ]
  
  write.table(wtest_results_pos, "filtered_samples_wilcox_test_results.tab", col.names=NA, sep="\t")
  
  ############################ test on ALL samples
  
  ttest_results <- NULL
  wtest_results <- NULL
  for(gene in row.names(sort_count_df)){
    res <- t.test(sort_count_df[gene,], unsort_count_df[gene,], paired=T)
    ttest_results <- rbind(ttest_results, c(res$estimate, res$statistic, res$p.value))
    
    wres <- wilcox.test(sort_count_df[gene,], unsort_count_df[gene,], paired=T, conf.int = TRUE)
    wtest_results <- rbind(wtest_results, c(wres$estimate, wres$statistic, wres$p.value))
  }
  
  row.names(ttest_results) <- row.names(sort_count_df)
  colnames(ttest_results) <- c("diff", "t-stats", "p.value")
  ttest_results <- ttest_results[order(ttest_results[, "p.value"]),]
  ttest_results_pos <- ttest_results[ttest_results[, "diff"] > 0, ]
  
  write.table(ttest_results_pos, "ALL_samples_paired_t_test_results.tab", col.names=NA, sep="\t")
  
  
  row.names(wtest_results) <- row.names(sort_count_df)
  colnames(wtest_results) <- c("diff", "w-stats", "p.value")
  wtest_results <- wtest_results[order(wtest_results[, "p.value"]),]
  wtest_results_pos <- wtest_results[wtest_results[, "diff"] > 0, ]
  
  write.table(wtest_results_pos, "ALL_samples_wilcox_test_results.tab", col.names=NA, sep="\t")
}


### the following sections are for exploratory analysis and plots #################

if(0){
  
  Labels <- colnames(ratio_over_ref_df)
  #pheatmap(as.matrix(ratio_over_ref_df))
  
  colors = rainbow(length(unique(Labels)))
  names(colors) = unique(Labels)
  set.seed(42)
  pdf("Tsne_analysis_of_ratio.pdf")
  tsne <- Rtsne(t(as.matrix(ratio_over_ref_df)), 
                dims = 2, 
                perplexity=10, 
                verbose=TRUE,
                num_threads = 10,
                pca=FALSE,
                normalize=FALSE, 
                max_iter = 5000)
  
  plot(tsne$Y, t='n', main="tsne: ratio over ref")
  text(tsne$Y, labels=Labels, col=colors[Labels])
  dev.off()
  
  
  Labels <- colnames(raw_count_df)
  #pheatmap(as.matrix(ratio_over_ref_df))
  
  colors = rainbow(length(unique(Labels)))
  names(colors) = unique(Labels)
  set.seed(42)
  pdf("Tsne_analysis_of_raw_count.pdf")
  tsne <- Rtsne(t(as.matrix(raw_count_df)), 
                dims = 2, 
                perplexity=10, 
                verbose=TRUE,
                num_threads = 10,
                pca=FALSE,
                normalize=FALSE, 
                max_iter = 5000)
  
  plot(tsne$Y, t='n', main="tsne: raw_count")
  text(tsne$Y, labels=Labels, col=colors[Labels])
  dev.off()
  
  
  ################################# Statistical tests ###############################
  
  ratio_over_ref_df_n <- as.matrix(type.convert(ratio_over_ref_df))
  head(ratio_over_ref_df_n)
  raw_count_df_n <- as.matrix(type.convert(raw_count_df))
  head(raw_count_df_n)
  
  med <- apply(ratio_over_ref_df, 2, median)
  batch1_sort <- c(1:16)
  batch1_unsort <- c(33:48)
  batch2_sort <- c(17:32)
  batch2_unsort <- c(49:64)
  batch_sort <- c(1:32)
  batch_unsort <- c(33:64)
  sample_list <- list("batch1"=list(batch1_sort, batch1_unsort), 
                      "batch2"=list(batch2_sort, batch2_unsort),
                      "batch"=list(batch_sort, batch_unsort))
  
  for(abatch in c("batch1", "batch2", "batch")){
    sort <- sample_list[[abatch]][[1]]
    unsort <- sample_list[[abatch]][[2]]
    sort_group <- colnames(raw_count_df_n)[sort]
    unsort_group <- colnames(raw_count_df_n)[unsort]
    sort_group
    unsort_group
    
    if(0){
      ## descriptive stats ratio
      sort_stats_ratio <- sumstats_row(ratio_over_ref_df_n[, sort_group])
      unsort_stats_ratio <- sumstats_row(ratio_over_ref_df_n[, unsort_group])
      
      stats_ratio <- cbind(sort_stats_ratio, unsort_stats_ratio)
      write.table(stats_ratio, paste(abatch, "ratio_summary_stats.tab", sep="_"), col.names=NA, sep="\t")
      
      sort_stats_ratio <- sumstats_col(ratio_over_ref_df_n[, sort_group])
      unsort_stats_ratio <- sumstats_col(ratio_over_ref_df_n[, unsort_group])
      
      stats_ratio <- rbind(sort_stats_ratio, unsort_stats_ratio)
      write.table(stats_ratio, paste(abatch, "ratio_summary_stats_byGuide.tab", sep="_"), col.names=NA, sep="\t")
      
      if(abatch == "batch"){
      
        pdf("correlation_among_guides_in_sort_ratio.pdf", width=10, height=10)
        old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
        
        sum_guide <- NULL
        for(i in c(seq(1,length(sort_group),4))){
          j <- i + 3
          pairs.panels(ratio_over_ref_df_n[, sort_group[i:j]], 
                       method = "pearson", # correlation method
                       hist.col = "#00AFBB",
                       density = TRUE,  # show density plots
                       ellipses = TRUE # show correlation ellipses
          )
          if(is.null(sum_guide)){
            sum_guide <- sumstats_row(ratio_over_ref_df_n[, sort_group[i:j]])
          }else{
            sum_guide <- cbind(sum_guide, sumstats_row(ratio_over_ref_df_n[, sort_group[i:j]]))
          }
        }
        dev.off()
        par(old.par)
        
        sort_sum_only <- sum_guide[, seq(1,length(colnames(sum_guide)),6)]
        colnames(sort_sum_only) <- sort_samples
        
        pdf("correlation_among_sumguides_in_sort_ratio.pdf", width=10, height=10)
        old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
        pairs.panels(sort_sum_only, 
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
        )
        dev.off()
        par(old.par)
        
        pdf("correlation_among_guides_in_unsort_ratio.pdf", width=10, height=10)
        old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
        sum_guide <- NULL
        for(i in seq(1,length(sort_group),4)){
          j <- i + 3
          pairs.panels(ratio_over_ref_df_n[, unsort_group[i:j]], 
                       method = "pearson", # correlation method
                       hist.col = "#00AFBB",
                       density = TRUE,  # show density plots
                       ellipses = TRUE # show correlation ellipses
          )
          if(is.null(sum_guide)){
            sum_guide <- sumstats_row(ratio_over_ref_df_n[, unsort_group[i:j]])
          }else{
            sum_guide <- cbind(sum_guide, sumstats_row(ratio_over_ref_df_n[, unsort_group[i:j]]))
          }
        }
        dev.off()
        par(old.par)
        
        unsort_sum_only <- sum_guide[, seq(1,length(colnames(sum_guide)),6)]
        colnames(unsort_sum_only) <- unsort_samples
        
        pdf("correlation_among_sumguides_in_unsort_ratio.pdf", width=10, height=10)
        old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
        pairs.panels(unsort_sum_only, 
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
        )
        dev.off()
        par(old.par)
        
        pdf("correlation_among_samples_in_sort_ratio.pdf", width=10, height=10)
        old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
        for(i in 0:3){
          s <- seq(1,length(sort_group),4)
          pairs.panels(ratio_over_ref_df_n[, sort_group[c(i+s)]], 
                       method = "pearson", # correlation method
                       hist.col = "#00AFBB",
                       density = TRUE,  # show density plots
                       ellipses = TRUE # show correlation ellipses
          )
        }
        dev.off()
        par(old.par)
        
        pdf("correlation_among_samples_in_unsort_ratio.pdf", width=10, height=10)
        old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
        for(i in 0:3){
          s <- seq(1,length(unsort_group),4)
          pairs.panels(ratio_over_ref_df_n[, unsort_group[c(i+s)]],  
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
          )
        }
        dev.off()
        par(old.par)
        
        ### all guides of a sample are summed, each sample has 1 entries
        sum_only <- cbind(sort_sum_only, unsort_sum_only)
        ttest_results <- NULL
        for(gene in row.names(sum_only)){
          res <- t.test(sum_only[gene, 1:8], sum_only[gene, 9:16], paired=T)
          ttest_results <- rbind(ttest_results, c(res$estimate, res$statistic, res$p.value))
        }
        
        row.names(ttest_results) <- row.names(sum_only)
        colnames(ttest_results) <- c("diff", "t-stats", "p.value")
        ttest_results <- ttest_results[order(ttest_results[, "p.value"]),]
        ttest_results_pos <- ttest_results[ttest_results[, "diff"] > 0, ]
        
        write.table(ttest_results_pos, paste(abatch, "ratio_paired_t_test_results_summed_guides.tab", sep="_"), col.names=NA, sep="\t")
        
      }
      
    ######################## paired t-test for ratio ###########################
    
      ### all guides together, each sample has 4 entries
      ttest_results <- NULL
      for(gene in row.names(ratio_over_ref_df_n)){
        res <- t.test(ratio_over_ref_df_n[gene, sort_group], ratio_over_ref_df_n[gene, unsort_group], paired=T)
        ttest_results <- rbind(ttest_results, c(res$estimate, res$statistic, res$p.value))
      }
      
      row.names(ttest_results) <- row.names(ratio_over_ref_df_n)
      colnames(ttest_results) <- c("diff", "t-stats", "p.value")
      ttest_results <- ttest_results[order(ttest_results[, "p.value"]),]
      ttest_results_pos <- ttest_results[ttest_results[, "diff"] > 0, ]
      
      write.table(ttest_results_pos, paste(abatch, "ratio_paired_t_test_results_all_guides.tab", sep="_"), col.names=NA, sep="\t")
      
      ### all guides are tested separately, each sample has 1 entries
      for(guide in guides){
        
        guide_sort_group <- sort_group[grep(guide, sort_group)]
        guide_unsort_group <- unsort_group[grep(guide, unsort_group)]
        
        ttest_results <- NULL
        for(gene in row.names(ratio_over_ref_df_n)){
          res <- t.test(ratio_over_ref_df_n[gene, guide_sort_group], ratio_over_ref_df_n[gene, guide_unsort_group], paired=T)
          ttest_results <- rbind(ttest_results, c(res$estimate, res$statistic, res$p.value))
        }
        
        row.names(ttest_results) <- row.names(ratio_over_ref_df_n)
        colnames(ttest_results) <- c("diff", "t-stats", "p.value")
        ttest_results <- ttest_results[order(ttest_results[, "p.value"]),]
        ttest_results_pos <- ttest_results[ttest_results[, "diff"] > 0, ]
        
        write.table(ttest_results_pos, paste(abatch, guide, "ratio_paired_t_test_results.tab", sep="_"), col.names=NA, sep="\t")
        
      }
      
    } #END if(0)
    
    ## descriptive stats raw count
    sort_stats_raw <- sumstats_row(raw_count_df_n[, sort_group])
    unsort_stats_raw <- sumstats_row(raw_count_df_n[, unsort_group])
    
    stats_raw <- cbind(sort_stats_raw, unsort_stats_raw)
    write.table(stats_raw, paste(abatch, "raw_count_summary_stats.tab", sep="_"), col.names=NA, sep="\t")
    
    sort_stats_raw <- sumstats_col(raw_count_df_n[, sort_group])
    unsort_stats_raw <- sumstats_col(raw_count_df_n[, unsort_group])
    
    stats_raw <- rbind(sort_stats_raw, unsort_stats_raw)
    write.table(stats_raw, paste(abatch, "raw_count_summary_stats_byGuide.tab", sep="_"), col.names=NA, sep="\t")
    
    if(abatch == "batch"){
      
      pdf("correlation_among_guides_in_sort_raw_count.pdf", width=10, height=10)
      old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
      
      sum_guide <- NULL
      for(i in c(seq(1,length(sort_group),4))){
        j <- i + 3
        pairs.panels(raw_count_df_n[, sort_group[i:j]], 
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
        )
        if(is.null(sum_guide)){
          sum_guide <- sumstats_row(raw_count_df_n[, sort_group[i:j]])
        }else{
          sum_guide <- cbind(sum_guide, sumstats_row(raw_count_df_n[, sort_group[i:j]]))
        }
      }
      dev.off()
      par(old.par)
      
      sort_sum_only <- sum_guide[, seq(1,length(colnames(sum_guide)),6)]
      colnames(sort_sum_only) <- sort_samples
      
      pdf("correlation_among_sumguides_in_sort_raw_count.pdf", width=10, height=10)
      old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
      pairs.panels(sort_sum_only, 
                   method = "pearson", # correlation method
                   hist.col = "#00AFBB",
                   density = TRUE,  # show density plots
                   ellipses = TRUE # show correlation ellipses
      )
      dev.off()
      par(old.par)
      
      pdf("correlation_among_guides_in_unsort_raw_count.pdf", width=10, height=10)
      old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
      sum_guide <- NULL
      for(i in seq(1,length(sort_group),4)){
        j <- i + 3
        pairs.panels(raw_count_df_n[, unsort_group[i:j]], 
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
        )
        if(is.null(sum_guide)){
          sum_guide <- sumstats_row(raw_count_df_n[, unsort_group[i:j]])
        }else{
          sum_guide <- cbind(sum_guide, sumstats_row(raw_count_df_n[, unsort_group[i:j]]))
        }
      }
      dev.off()
      par(old.par)
      
      unsort_sum_only <- sum_guide[, seq(1,length(colnames(sum_guide)),6)]
      colnames(unsort_sum_only) <- unsort_samples
      
      pdf("correlation_among_sumguides_in_unsort_raw_count.pdf", width=10, height=10)
      old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
      pairs.panels(unsort_sum_only, 
                   method = "pearson", # correlation method
                   hist.col = "#00AFBB",
                   density = TRUE,  # show density plots
                   ellipses = TRUE # show correlation ellipses
      )
      dev.off()
      par(old.par)
      
      pdf("correlation_among_samples_in_sort_raw_count.pdf", width=10, height=10)
      old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
      for(i in 0:3){
        s <- seq(1,length(sort_group),4)
        pairs.panels(raw_count_df_n[, sort_group[c(i+s)]], 
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
        )
      }
      dev.off()
      par(old.par)
      
      pdf("correlation_among_samples_in_unsort_raw_count.pdf", width=10, height=10)
      old.par <- par(mfrow=c(1,1),mar=c(2,2,2,2))
      for(i in 0:3){
        s <- seq(1,length(unsort_group),4)
        pairs.panels(raw_count_df_n[, unsort_group[c(i+s)]],  
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE # show correlation ellipses
        )
      }
      dev.off()
      par(old.par)
      
      ### all guides of a sample are summed, each sample has 1 entries
      sum_only <- as.matrix(cbind(sort_sum_only, unsort_sum_only))
      ttest_results <- NULL
      for(gene in row.names(sum_only)){
        #if(sum(sum_only[gene, 1:8]) > 10){
          res <- t.test(sum_only[gene, 1:8], sum_only[gene, 9:16], paired=T)
          ttest_results <- rbind(ttest_results, c(res$estimate, res$statistic, res$p.value))
        #}
      }
      
      row.names(ttest_results) <- row.names(sum_only)
      colnames(ttest_results) <- c("diff", "t-stats", "p.value")
      ttest_results <- ttest_results[order(ttest_results[, "p.value"]),]
      ttest_results_pos <- ttest_results[ttest_results[, "diff"] > 0, ]
      
      write.table(ttest_results_pos, paste(abatch, "rawCount_paired_t_test_results_summed_guides.tab", sep="_"), col.names=NA, sep="\t")
      
    }
    
  }
}



## END ####

