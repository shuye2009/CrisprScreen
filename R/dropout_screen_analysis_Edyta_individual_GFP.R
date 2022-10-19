### for the analysis of dropout screen data of Edyta



## setup ####
if(1){
library(ssh)
library(EnhancedVolcano)
library(corrplot)
  
source("C:/GREENBLATT/Rscripts/CrisprScreen/crispr_screen_analysis_lib.R")

merge_GFP <- FALSE ## TRUE, FALSE, NA
def_res_list <- list()  
  
name_map <- NULL
mapfile <- "C:/GREENBLATT/Rscripts/CrisprScreen/resources/hgnc_geneName_updates_2019.txt"
name_map <- read.table(mapfile, header=T, stringsAsFactors = FALSE)
head(name_map)
dim(name_map)


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
}


#### merge GFP data ####
if(!is.na(merge_GFP)){
  if(merge_GFP){
    wd <- "C:/data/raw/EDYTAdropoutScreen/all_screens" 
    
    setwd(wd)
    format1 <- ".*GFP.*_geneRawCount.tab"
    
    gfpfiles <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = format1)
    exclude <- c("Greenblatt_005_210310_GFP_T0_S5_geneRawCount.tab", 
                 "Greenblatt_006_210310_GFP_T18_A_S6_geneRawCount.tab",
                 "Greenblatt_007_210310_GFP_T18_B_S7_geneRawCount.tab",
                 "Greenblatt_008_210310_GFP_T18_C_S8_geneRawCount.tab")
    gfpfiles <- gfpfiles[!gfpfiles %in% exclude]
    length(gfpfiles)
    gfpfiles
    gfpfiles <- gfpfiles[!grepl("XXX_XXX", gfpfiles)] ## exlude existing merged files
    T0_files <- gfpfiles[grepl("T0", gfpfiles)]
    merge_geneRawCount(T0_files, "GFP_T0_merged")
    
    for(rep in c("A_", "B_", "C_")){
      rep_files <- gfpfiles[grepl(rep, gfpfiles, fixed=T) & grepl("T18_", gfpfiles, fixed=T)]
      merge_geneRawCount(rep_files, paste("GFP_",rep,"T18_merged", sep=""))
    }
    
    format2 <- ".*GFP.*_guideRawCount.tab"
    
    gfpfiles <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = format2)
    exclude <- c("Greenblatt_005_210310_GFP_T0_S5_guideRawCount.tab", 
                 "Greenblatt_006_210310_GFP_T18_A_S6_guideRawCount.tab",
                 "Greenblatt_007_210310_GFP_T18_B_S7_guideRawCount.tab",
                 "Greenblatt_008_210310_GFP_T18_C_S8_guideRawCount.tab")
    gfpfiles <- gfpfiles[!gfpfiles %in% exclude]
    length(gfpfiles)
    gfpfiles
    gfpfiles <- gfpfiles[!grepl("XXX_XXX", gfpfiles)] ## exlude existing merged files
    T0_files <- gfpfiles[grepl("T0", gfpfiles)]
    merge_guideRawCount(T0_files, "GFP_T0_merged")
    
    for(rep in c("A_", "B_", "C_")){
      rep_files <- gfpfiles[grepl(rep, gfpfiles, fixed=T) & grepl("T18_", gfpfiles, fixed=T)]
      merge_guideRawCount(rep_files, paste("GFP_",rep,"T18_merged", sep=""))
    }
    
  }
}
#### scale and normalize data ####
if(1){
  wd <- "C:/data/raw/EDYTA/dropoutScreen/all_screens" 
  
  setwd(wd)
  format2 <- "_guideRawCount.tab"
  
  guidefiles <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = format2)
  length(guidefiles)
  guidefiles
  
  if(1){
    samples <- gsub("_guideRawCount.tab", "", guidefiles)
    replicates <- unlist(lapply(guidefiles, function(x)split_group(x,4,6)))
    subjects <- unlist(lapply(guidefiles, function(x)split_group(x,4,4)))
    design <- data.frame(samples, replicates, subjects)
    rownames(design) <- guidefiles
    design <- design[order(design$subjects),]
    if(!is.na(merge_GFP)){
      if(merge_GFP){
        write.table(design, "experimental_design_merged_GFP_raw.tab", sep="\t", row.names=T, col.names=NA, quote=F)  
      }else{
        write.table(design, "experimental_design_individual_GFP_raw.tab", sep="\t", row.names=T, col.names=NA, quote=F)
      }
    }
  }
  ## modify "experimental_design_individual_GFP_raw.tab" to correct parsing errors, and save as "experimental_design_individual_GFP.tab"
  if(is.na(merge_GFP)){
    eXdesign <- data.frame(read.delim("experimental_design_only_GFP.tab", header=TRUE, sep="\t", stringsAsFactors=F))
  }else{
      if(merge_GFP){
        eXdesign <- data.frame(read.delim("experimental_design_merged_GFP.tab", header=TRUE, sep="\t", stringsAsFactors=F))
      }else if(!merge_GFP){
        eXdesign <- data.frame(read.delim("experimental_design_individual_GFP1.tab", header=TRUE, sep="\t", stringsAsFactors=F))
      }
  }
    
  rownames(eXdesign) <- eXdesign$X
  eXdesign <- eXdesign[,-1]
  eXdesign <- eXdesign[!grepl("T10|NR", eXdesign$replicates), ]  ## exclude T10 and NR samples
  eXdesign
  
  guide_data <- process_guideRawCount(eXdesign) ## output can be used for BAGEL directly
  gene_data <- process_geneRawCount(eXdesign)
  

  combined_gene_table <- NULL  ## each row is a gene
  combined_guide_table <- NULL ## each row is a guide
  combined_lfc_guide_table <- NULL ## each row is a guide
  combined_scaled_guide_table <- NULL ## each row is a guide
  combined_normalized_guide_table <- NULL ## each row is a guide
  
  #for (wd in c(wd1, wd2, wd3, wd4)){
  for(subject in unique(eXdesign$subjects)){
    
    geneRaw_table <- gene_data[[subject]]
    guideRaw_table <- guide_data[["Raw"]][[subject]]
    guideScaled_table <- guide_data[["Scaled"]][[subject]]
    guideNormalized_table <- guide_data[["Normalized"]][[subject]]
    guideLFC_table <- guide_data[["LFC"]][[subject]]
    
    head(geneRaw_table)
    dim(geneRaw_table)
    summary(geneRaw_table)
    
    head(guideRaw_table)
    dim(guideRaw_table)
    summary(guideRaw_table)
    
    head(guideScaled_table)
    dim(guideScaled_table)
    summary(guideScaled_table)
    
    geneRaw_table <- geneRaw_table[TKO3_genes, ]
    
    if(is.null(combined_gene_table)){
      combined_gene_table <- geneRaw_table
    }else{
      combined_gene_table <- cbind(combined_gene_table, geneRaw_table)  
    }
    
    
    if(is.null(combined_guide_table)){
      combined_guide_table <- guideRaw_table
    }else{
      combined_guide_table <- merge(combined_guide_table, guideRaw_table, by=c("SEQUENCE", "GENE"), all=TRUE)
    }
   
    if(is.null(combined_scaled_guide_table)){
      combined_scaled_guide_table <- guideScaled_table
    }else{
      combined_scaled_guide_table <- merge(combined_scaled_guide_table, guideScaled_table, by=c("SEQUENCE", "GENE"), all=TRUE)
    }
    
    if(is.null(combined_lfc_guide_table)){
      combined_lfc_guide_table <- guideLFC_table
    }else{
      combined_lfc_guide_table <- merge(combined_lfc_guide_table, guideLFC_table, by=c("SEQUENCE", "GENE"), all=TRUE)
    }
    
    if(is.null(combined_normalized_guide_table)){
      combined_normalized_guide_table <- guideNormalized_table
    }else{
      combined_normalized_guide_table <- merge(combined_normalized_guide_table, guideNormalized_table, by=c("SEQUENCE", "GENE"), all=TRUE)
    }
  }

    
  if(is.na(merge_GFP)){
    GFPhome_dir <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only" 
  }else{
    if(merge_GFP){
      GFPhome_dir <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged"  
    }else if(!merge_GFP){
      GFPhome_dir <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual"
    }
  }

  if(!dir.exists(GFPhome_dir)){
    dir.create(GFPhome_dir, showWarnings = F)
  }
  
  dim(combined_gene_table)
  print("writing combined table file")
  write.table(combined_gene_table, file.path(GFPhome_dir,"combined_gene_rawCount_table.tab"), row.names=T, col.names=NA, sep="\t", quote=F)
  
  dim(combined_guide_table)
  write.table(combined_guide_table, file.path(GFPhome_dir,"combined_guide_rawCount_table.tab"), row.names=F, sep="\t", quote=F)
  combined_guide_filtered <- na.omit(combined_guide_table[, 3:ncol(combined_guide_table)])
  
  dim(combined_scaled_guide_table)
  ## "combined_guide_scaledCount_table.tab" can be used as input for drugZ directly
  write.table(combined_scaled_guide_table, file.path(GFPhome_dir,"combined_guide_scaledCount_table.tab"), row.names=F, sep="\t", quote=F)
  combined_scaled_guide_filtered <- na.omit(combined_scaled_guide_table[, 3:ncol(combined_scaled_guide_table)])
  
  
  dim(combined_normalized_guide_table)
  write.table(combined_normalized_guide_table, file.path(GFPhome_dir,"combined_guide_normalizedCount_table.tab"), row.names=F, sep="\t", quote=F)
  combined_normalized_guide_filtered <- na.omit(combined_normalized_guide_table[, 3:ncol(combined_normalized_guide_table)])
  
  
  dim(combined_lfc_guide_table)
  combined_lfc_guide_table[is.na(combined_lfc_guide_table)] <- 0
  write.table(combined_lfc_guide_table, file.path(GFPhome_dir,"combined_guide_lfc_table.tab"), row.names=F, sep="\t", quote=F)
  combined_lfc_guide_filtered <- na.omit(combined_lfc_guide_table[, 3:ncol(combined_lfc_guide_table)])
  
}



#Prepare files for JACKS ####
if(1){
  wd <- "C:/data/raw/EDYTA/dropoutScreen/all_screens" 
  setwd(wd)
  if(is.na(merge_GFP)){
    eXdesign <- data.frame(read.delim("experimental_design_only_GFP.tab", header=TRUE, sep="\t", stringsAsFactors=F))
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      eXdesign <- data.frame(read.delim("experimental_design_merged_GFP.tab", header=TRUE, sep="\t", stringsAsFactors=F))
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      eXdesign <- data.frame(read.delim("experimental_design_individual_GFP1.tab", header=TRUE, sep="\t", stringsAsFactors=F))
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual"
    }
  }
  rownames(eXdesign) <- eXdesign$X
  eXdesign <- eXdesign[,-1]
  eXdesign <- eXdesign[!grepl("T10|NR", eXdesign$replicates), ]  ## exclude T10 and NR samples
  eXdesign
  
  setwd(hwd)
  ### evaluate reproducibility 
  
  scaled_count <- read.delim("combined_guide_scaledCount_table.tab", header=T, sep="\t")
  raw_count <- read.delim("combined_guide_rawCount_table.tab", header=T, sep="\t")
  normalized_count <- read.delim("combined_guide_normalizedCount_table.tab", header=T, sep="\t")
  lfc_count <- read.delim("combined_guide_lfc_table.tab", header=T, sep="\t")
  
  
  
  ## compare between-replicates correlation
  replicates <- colnames(scaled_count)[3:ncol(scaled_count)]
  samples <- unique(sort(eXdesign[eXdesign$replicates %in% replicates, "subjects"]))
  
 
  df_list <- list()
  df_list[["Scaled_to_T0"]] <- scaled_count[,c(-1,-2)]
  df_list[["Raw"]] <- raw_count[, c(-1,-2)]
  df_list[["TotalReads"]] <- normalized_count[,c(-1,-2)]
  df_list[["LFC_of_T0"]] <- lfc_count[,c(-1,-2)]
  
  if(1){
    mean_corr <- data.frame()
    for(ds in names(df_list)){
      dfs <- df_list[[ds]]
      print(ds)
      print(dim(dfs))
      dfs <- na.omit(dfs)
      
      print(dim(dfs))
      plot_heatmap(dfs, ds)
      for(screen in c(paste(samples, "T18", sep="_"))){
        subdf <- dfs[, grepl(screen, colnames(dfs))]
        cor_mean <- (sum(cor(subdf))-3)/6
        cor_mean <- round(cor_mean, digits=2)
        screen <- gsub("_T18", "", screen)
        mean_corr <- rbind(mean_corr, c(ds, screen, cor_mean))
      }
      
    }
    
    colnames(mean_corr) <- c("Normalization", "Replicate", "Correlation")
    p <- ggbarplot(data=mean_corr, x="Replicate", y="Correlation", fill="Normalization", 
                   xlab=FALSE, position = position_dodge(0.9)) %>% 
      ggpar(font.tickslab = 10)
    ggexport(p, filename="comparison_of_normalization_methods.pdf")
  
  }
  
  if(1){
  ## create tables for JACKS
  input_for_JACKS <- scaled_count[,-2]
  input_for_JACKS <- input_for_JACKS[, !grepl("T0_", colnames(input_for_JACKS))]
  dim(input_for_JACKS)
  summary(input_for_JACKS)
  input_for_JACKS[is.na(input_for_JACKS)] <- 0
  summary(input_for_JACKS)
  head(input_for_JACKS)
  
  
  replicates <- colnames(input_for_JACKS)[2:ncol(input_for_JACKS)]
  samples <- unlist(lapply(replicates, function(x) gsub("_A|_B|_C", "", x)))
  if(merge_GFP){
    controls <- rep("GFPm_T18", length(samples))
    controls[grepl("_T10", samples)] <- "GFPm_T10"
  }else{
    controls <- rep("GFP1_T18", length(samples))
    controls[grepl("_T10", samples)] <- "GFP1_T10"
  }
  #samples <- unlist(lapply(replicates, function(x) unlist(strsplit(x, split="_"))[1]))
  #samples <- eXdesign[eXdesign$replicates %in% replicates, "subjects"]
  #controls <- eXdesign[eXdesign$replicates %in% replicates, "controls"]

  
  design <- data.frame(replicates, samples, controls)
  design
  
  guide_gene <- scaled_count[,1:2]
  head(guide_gene)
  
  
  control_gene <- c("EGFP", "LacZ", "luciferase")
  control_guide <- guide_gene[guide_gene$GENE %in% control_gene, 1]
  control_count <- input_for_JACKS[input_for_JACKS$SEQUENCE %in% control_guide,]
  dim(control_count)
  
  jacks_dir <- file.path(hwd,"JACKS")
  
  if(!dir.exists(jacks_dir)){
    dir.create(jacks_dir, showWarnings = F)
  }
  setwd(jacks_dir)
  
  write.table(input_for_JACKS, "countfile.tab", row.names=F, sep="\t", quote=F)
  write.table(design, "replicatemapfile.tab", row.names=F, sep="\t", quote=F)
  write.table(guide_gene, "guidegenemapfile.tab", row.names=F, sep="\t", quote=F)
  writeLines(control_gene, "controlGenefile.txt")
  
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  if(is.na(merge_GFP)){
    system("scp countfile.tab replicatemapfile.tab guidegenemapfile.tab controlGenefile.txt shuyepu@dc01.ccbr.utoronto.ca:./JACKS/only_GFP")
  }else{
    if(merge_GFP){
      scp_upload(sshSession, c("countfile.tab", "replicatemapfile.tab", "guidegenemapfile.tab", "controlGenefile.txt"), "./JACKS/merged_GFP")
      
    }else{
      scp_upload(sshSession, c("countfile.tab", "replicatemapfile.tab", "guidegenemapfile.tab", "controlGenefile.txt"), "./JACKS/individual_GFP")
    }
  }
  ## Run JACKS on dc cluster here
  ssh_disconnect(sshSession)
  }
}

## prepare for drugZ analysis scaled ####
if(1){
  ## scp "combined_guide_scaledCount_table.tab" to dc01.ccbr.utoronto.ca
  dataset <- "drugZ"
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
    setwd(hwd)
    system("scp combined_guide_scaledCount_table.tab shuyepu@dc01.ccbr.utoronto.ca:./drugz/only_GFP")
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
      setwd(hwd)
      scp_upload(sshSession, "combined_guide_scaledCount_table.tab", "./drugz/merged_GFP")
      
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual"
      setwd(hwd)
      scp_upload(sshSession, "combined_guide_scaledCount_table.tab", "./drugz/individual_GFP") 
    }
  }
  
  ssh_disconnect(sshSession)
  
  wd <- paste(hwd, dataset, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  ## run drugz on dc1 cluster
}


## JACKS analysis ####
if(1){
  
  dataset <- "JACKS"
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
    setwd(file.path(hwd, dataset))
    system("scp shuyepu@dc01.ccbr.utoronto.ca:./JACKS/only_GFP/*_JACKS_results.txt .")
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged"
      setwd(file.path(hwd, dataset))
      scp_download(sshSession, "./JACKS/merged_GFP", ".")
      system("mv ./merged_GFP/*_JACKS_results.txt .")
      #system("scp shuyepu@dc01.ccbr.utoronto.ca:./JACKS/merged_GFP/*_JACKS_results.txt .")
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual"
      setwd(file.path(hwd, dataset))
      scp_download(sshSession, "./JACKS/individual_GFP", ".")
      system("mv ./individual_GFP/*_JACKS_results.txt .")
      #system("scp shuyepu@dc01.ccbr.utoronto.ca:./JACKS/individual_GFP/*_JACKS_results.txt .")
    }
  }
  ssh_disconnect(sshSession)
  
  jacks_list <- list()
  jacks_results_neg <- read.delim("negative_screen_gene_JACKS_results.txt", header=T) 
  jacks_results_pos <- read.delim("positive_screen_gene_JACKS_results.txt", header=T)
  jacks_pvalue_neg <- read.delim("negative_screen_gene_pval_JACKS_results.txt", header=T)
  jacks_pvalue_pos <- read.delim("positive_screen_gene_pval_JACKS_results.txt", header=T)
  
  #plot_heatmap(jacks_results_neg[,-1], "JACKS_negative")
  #plot_heatmap(jacks_results_pos[,-1], "JACKS_positive")
  
  jacks_pvalue_neg[is.na(jacks_pvalue_neg)] <- 1
  jacks_padj_neg <- jacks_pvalue_neg
  for(i in 2:ncol(jacks_pvalue_neg)){
    jacks_padj_neg[,i] <- p.adjust(jacks_pvalue_neg[,i])
  }
  
  jacks_pvalue_pos[is.na(jacks_pvalue_pos)] <- 1
  jacks_padj_pos <- jacks_pvalue_pos
  for(i in 2:ncol(jacks_pvalue_pos)){
    jacks_padj_pos[,i] <- p.adjust(jacks_pvalue_pos[,i])
  }
  
  #head(jacks_pvalue_neg[order(jacks_pvalue_neg$N),])
  #head(jacks_padj_neg[order(jacks_padj_neg$N),])
  
  gene_universe <- TKO3_genes
  
  screens <- colnames(jacks_pvalue_neg)[2:ncol(jacks_pvalue_neg)]
  for(i in 1:length(screens)){
    ascreen <- screens[i]
    res <- NULL
    
    for(df in list(jacks_results_neg, jacks_pvalue_neg, jacks_padj_neg, jacks_results_pos, jacks_pvalue_pos, jacks_padj_pos)){
      rownames(df) <- df[,1]
      res <- cbind(res, df[gene_universe, ascreen])
    }
    rownames(res) <- gene_universe;
    colnames(res) <- c("neg", "neg_p", "neg_padj", "pos", "pos_p", "pos_padj");
    
    jacks_list[[ascreen]] <- res
  }
  
  
  
  def_res_list[[dataset]] <- jacks_list
  
} ## end if(1)



## drugZ analysis scaled ####
if(1){
  ## scp "combined_guide_scaledCount_table.tab" to dc01.ccbr.utoronto.ca
  
  library(car)
  dataset <- "drugZ"
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
    setwd(file.path(hwd, dataset))
    system("scp shuyepu@dc01.ccbr.utoronto.ca:./drugz/only_GFP/T0_scaled_*_drugz-output.txt .")
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
      setwd(file.path(hwd, dataset))
      scp_download(sshSession, files="./drugz/merged_GFP", to=".")
      system("mv ./merged_GFP/T0_scaled_*_T1*_drugz-output.txt .")
      #system("scp -i C:/cygwin64/home/greenblatt/.ssh/id_rsa.pub shuyepu@dc01.ccbr.utoronto.ca:./drugz/merged_GFP/T0_scaled_*_T1*_drugz-output.txt .")
      
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
      setwd(file.path(hwd, dataset))
      scp_download(sshSession, "./drugz/individual_GFP/T0_scaled_*_T1*_drugz-output.txt", ".")
      system("mv ./individual_GFP/T0_scaled_*_T1*_drugz-output.txt .")
      #system("scp shuyepu@dc01.ccbr.utoronto.ca:./drugz/individual_GFP/T0_scaled_*_drugz-output.txt .")
    }
  }
  ssh_disconnect(sshSession)
  
  drugZ_analysis=list.files(pattern="^T0_scaled_.+_T1.+_drugz-output.txt")
  
  screens <- unlist(lapply(drugZ_analysis, function(x) paste(unlist(strsplit(x, split="_"))[3:4], collapse="_")))
  screens
  data_list <- list()
  
  for(i in 1:length(screens)){
    ascreen <- screens[i]
    afile <- drugZ_analysis[i]
    data_list[[ascreen]] <- read.delim(file.path(hwd, dataset,afile), header=T)
  }
  
  def_res_list[[dataset]] <- data_list
  
  pvalue_cutoff <- 0.05
  
  sig_list <- list()
  
  
  for(treat in names(data_list)){
    treat <- "N_T18"
    df <- data_list[[treat]]
    downlist <- df[df$pval_synth < pvalue_cutoff, "GENE"]
    uplist <- df[df$pval_supp < pvalue_cutoff, "GENE"]
    
    sig_list[[paste(treat, "up", sep="_")]] <- uplist
    sig_list[[paste(treat, "down", sep="_")]] <- downlist
    
    df <- na.omit(df)
    rownames(df) <- df$GENE
    #df <- mutate(df, pvalue=min(pval_synth, pval_supp))
    df$pvalue <- apply(df[, c("pval_synth", "pval_supp")], 1, min)
    df$fdr <- apply(df[, c("fdr_synth", "fdr_supp")], 1, min)
    
    ggscatter(df, x = "rank_synth", y = "normZ",
              size = abs(df$normZ), color = "normZ",
              label = df$GENE, label.select = df$GENE[1:10], repel = TRUE) + gradient_color(c("blue", "white", "red"))
    
    
    pdf(paste(treat, "vocano_plot_allgenes.pdf", sep="_"), width=8, height=8)
    p <- EnhancedVolcano(df,
                         lab = df$GENE,
                         x = 'normZ',
                         y = 'pvalue',
                         #selectLab = df$GENE[2:11],
                         title = paste("drugz analysis", treat),
                         subtitle = "",
                         caption = "NormZ score cutoff = +/-2; pvalue cutoff = 0.01",
                         legendLabels = c("NS", "NormZ", "P-value", "P-value and NormZ"),
                         xlab = "NormZ score",
                         pCutoff = 0.01,
                         FCcutoff = 2,
                         labSize = 4,
                         pointSize = 3.0,
                         drawConnectors = FALSE,
                         xlim = c(-7.5, 7.5),
                         ylim = c(0, 6))
    print(p)
    dev.off()
    
    normZ <- df$normZ
    names(normZ) <- df$GENE
    normZ <- normZ[TKO3_genes]
    normZ <- na.omit(normZ)
    
    print(length(normZ))
    
    
    while(!is.null(dev.list())){
      dev.off()
    }
    if(0){
      run_gseGO_simpleList(normZ, treat)
    }
    if(0){
      run_enrichGO_simpleList(uplist, "BP", paste(treat,"upGenes", sep="_"))
      run_enrichGO_simpleList(downlist, "BP", paste(treat,"downGenes", sep="_"))
      #run_enrichDO_simpleList(uplist, paste(treatSF,"upGenes", sep="_"))
      #run_enrichDO_simpleList(downlist, paste(treatSF,"downGenes", sep="_"))
      #run_gseGO_simpleList(normZ, treat)
      
      complexes_up <- overlap_with_CORUM(uplist, idMap, corum)
      complexes_down <- overlap_with_CORUM(downlist, idMap, corum)
      write.table(complexes_up, paste(treat, "up_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
      write.table(complexes_down, paste(treat, "down_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
    } ## end if
    
  } ## end for
  if(0){
    gene_lists_ENTREZ <- lapply(sig_list, function(x)y<- bitr(x, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID)
    gene_lists_ENTREZ <- lapply(gene_lists_ENTREZ, noquote)
    
    plot_clusterProfile_GO(gene_lists_ENTREZ, "drugZ_modulated_genes")
  } ## end if
  
  
  
  
} ## end if(1)


## Limma analysis on LFC ####

if(1){
  
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  dataset <- "Limma_on_lfc"
  wd <- paste(hwd, dataset, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  
  
  lfc_combined <- read.delim(file.path(hwd, "combined_guide_lfc_table.tab"), header=T, sep="\t")
  
  
  dim(lfc_combined)
  head(lfc_combined)
  summary(lfc_combined)
  
  
  NAs <- apply(lfc_combined, 1, function(x) sum(is.na(x)))
  NAss <- which(NAs > 0)
  length(NAss)
  
  ## expand rows for each gene using guide id
  guides <- c("g1", "g2", "g3", "g4")
  
  allGenes <- sort(unique(lfc_combined$GENE))
  allGenes[1:10]
  lfc_expanded <- NULL
  c <- 0
  valid_gene <- NULL
  for(gene in allGenes){
    sub <- as.matrix(lfc_combined[lfc_combined$GENE %in% gene, 3:ncol(lfc_combined)])
    if(nrow(sub) < 4){
      
      if(nrow(sub) == 3){
        arow <- apply(sub, 2, mean) ## replace missing guide with mean of existing guides
        sub <- rbind(sub, arow)
        valid_gene <- c(valid_gene, gene)
        subvector <- as.vector(sub)
        lfc_expanded <- rbind(lfc_expanded, subvector)
      }else{
        c <- c + 1
        print(paste(c, gene, nrow(sub)))
      }
      
      
    }else if(nrow(sub) == 4){
      valid_gene <- c(valid_gene, gene)
      subvector <- as.vector(sub)
      lfc_expanded <- rbind(lfc_expanded, subvector)
    }
    
  }
  
  print(paste("genes with less than 3 guides ", c))
  
  cn <- unlist(lapply(colnames(lfc_combined)[3:ncol(lfc_combined)], function(x) paste(x, guides, sep="_")))
  colnames(lfc_expanded) <- cn
  rownames(lfc_expanded) <- valid_gene
  
  dim(lfc_expanded)
  head(lfc_expanded)
  summary(lfc_expanded)
  
  
  lfc_table <- lfc_expanded
  lfc_table[is.na(lfc_table)] <- 0
  lfc_mat <- as.matrix(lfc_table)
  dim(lfc_mat)
  
  varVec <- sort(apply(lfc_mat, 1, var), decreasing=T)
  
  top <- length(varVec)
  lfcMatTopVar <- lfc_mat[names(varVec[1:top]), ]
  
  screens <- unlist(lapply(colnames(lfcMatTopVar), function(x) paste(unlist(strsplit(x, split="_"))[1:2], collapse="_")))
  labels <- factor(
    screens,
    levels = unique(screens)
  )
  
  design <- model.matrix(~ -1 + labels)
  colnames(design) <- levels(labels)
  rownames(design) <- colnames(lfc_mat)
  
  par(mfrow = c(2,2))
  for(i in 1:4){
    plot(lfcMatTopVar[i, ] ~ labels, main=rownames(lfcMatTopVar)[i], ylim = c(-10, 10))
  }
  par(mfrow = c(1,1))
  
  plotMDS(lfcMatTopVar)
  
  eset <- ExpressionSet(assayData = lfcMatTopVar)
  
  if(is.na(merge_GFP)){
    contrast.matrix <- makeContrasts(
      GFP1_GFPm = GFP1 - GFPm,
      GFP2_GFPm = GFP2 - GFPm,
      GFP3_GFPm = GFP3 - GFPm,
      GFP4_GFPm = GFP4 - GFPm,
      levels = design
    )
  }else{
    if(!merge_GFP){
      contrast.matrix <- makeContrasts(
        N_T18_GFP = N_T18 - GFP1_T18,
        NR_T18_GFP = NR_T18 - GFP1_T18,
        NSP12_T18_GFP = NSP12_T18 - GFP1_T18,
        N_T10_GFP = N_T10 - GFP1_T10,
        NR_T10_GFP = NR_T10 - GFP1_T10,
        NSP12_T10_GFP = NSP12_T10 - GFP1_T10,
        S_T18_GFP = S_T18 - GFP1_T18,
        NSP9_T18_GFP = NSP9_T18 - GFP1_T18,
        NSP10_T18_GFP = NSP10_T18 - GFP1_T18,
        ORF9B_T18_GFP = ORF9B_T18 - GFP1_T18,
        ORF3A_T18_GFP = ORF3A_T18 - GFP1_T18,
        ORF8_T18_GFP = ORF8_T18 - GFP1_T18,
        NSP2_T18_GFP = NSP2_T18 - GFP1_T18,
        NSP3_T18_GFP = NSP3_T18 - GFP1_T18,
        NSP2_T10_GFP = NSP2_T10 - GFP1_T10,
        NSP3_T10_GFP = NSP3_T10 - GFP1_T10,
        M_T18_GFP = M_T18 - GFP1_T18,
        levels = design
      )
    }else{
      contrast.matrix <- makeContrasts(
        N_T18_GFP = N_T18 - GFPm_T18,
        NR_T18_GFP = NR_T18 - GFPm_T18,
        NSP12_T18_GFP = NSP12_T18 - GFPm_T18,
        N_T10_GFP = N_T10 - GFPm_T10,
        NR_T10_GFP = NR_T10 - GFPm_T10,
        NSP12_T10_GFP = NSP12_T10 - GFPm_T10,
        S_T18_GFP = S_T18 - GFPm_T18,
        NSP9_T18_GFP = NSP9_T18 - GFPm_T18,
        NSP10_T18_GFP = NSP10_T18 - GFPm_T18,
        ORF9B_T18_GFP = ORF9B_T18 - GFPm_T18,
        ORF3A_T18_GFP = ORF3A_T18 - GFPm_T18,
        ORF8_T18_GFP = ORF8_T18 - GFPm_T18,
        NSP2_T18_GFP = NSP2_T18 - GFPm_T18,
        NSP3_T18_GFP = NSP3_T18 - GFPm_T18,
        NSP2_T10_GFP = NSP2_T10 - GFPm_T10,
        NSP3_T10_GFP = NSP3_T10 - GFPm_T10,
        M_T18_GFP = M_T18 - GFPm_T18,
        levels = design
      )
    }
  }
  
  fit <- lmFit(eset, design)
  fit.cont <- contrasts.fit(fit, contrast.matrix)
  fit.eb   <- eBayes(fit.cont)
  
  genas(fit.eb, coef=c(1,2), subset="all", plot=TRUE, alpha=0.4)
  
  res_list <- list()
  modulated_gene_list <- list()
  pvalue_cutoff <- 0.01
  for(i in 1:length(colnames(fit.eb$contrasts))){
    #i <- 1
    contrast <- colnames(fit.eb$contrasts)[i]
    
    term1 <- rownames(fit.eb$contrasts)[fit.eb$contrasts[,contrast] == 1]
    term2 <- rownames(fit.eb$contrasts)[fit.eb$contrasts[,contrast] == -1]
    
    tpall <- topTable(fit.eb, number=top, coef=i, p.value=1)
    res_list[[contrast]] <- tpall
    
    tn <- propTrueNull(tpall$P.Value, method="lfdr", nbins=20)
    
    #write.table(tpall, paste(contrast, "limma_table.tab", sep="_"), sep="\t", quote=F, col.names=NA)
    tp <- tpall[tpall$P.Val < pvalue_cutoff,]
    tp_up <- rownames(tp[tp$logFC>0,])
    tp_down <- rownames(tp[tp$logFC<0,])
    
    rankedlist <- tpall$t
    names(rankedlist) <- rownames(tpall)
    
    if(0){
    run_enrichGO_simpleList(tp_up, ont="BP", paste(contrast,"_limma_up", sep=""))
    run_enrichGO_simpleList(tp_down, ont="BP", paste(contrast,"_limma_down", sep=""))
    run_gseGO_simpleList(rankedlist, contrast)
    
    complexes_up <- overlap_with_CORUM(tp_up, idMap, corum)
    complexes_down <- overlap_with_CORUM(tp_down, idMap, corum)
    write.table(complexes_up, paste(contrast, "limma_up_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
    write.table(complexes_down, paste(contrast, "limma_down_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
    }
    
    modulated_gene_list[[paste(contrast, "up", sep="_")]] <- tp_up
    modulated_gene_list[[paste(contrast, "down", sep="_")]] <- tp_down
    
    pdf(paste(contrast, "_limma_volvanoplot.pdf", sep=""), height=8, width=8)
    volcanoplot(fit.eb,coef=i, highlight=30, xlab="Log fold change", main=contrast, names=rownames(fit.eb$coefficients))
    dev.off()
    
    pdf(paste(contrast, "_limma_scatterplot.pdf", sep=""), height=8, width=8)
    plot(x=fit$coefficients[,term2], y=fit$coefficients[,term1], col="white", xlab=paste(term2, "lfc"), ylab=paste(term1, "lfc"))
    colors <- rep("black", dim(fit$coefficients)[1])
    names(colors) <- rownames(fit$coefficients)
    colors[tp_up] <- "red"
    colors[tp_down] <- "blue"
    #plot(x=broom_fit[,3], y=broom_fit[,2], col=colors, xlab="GFP_lfc", ylab=TREAT)
    text(rownames(fit$coefficients), x=fit$coefficients[,term2], y=fit$coefficients[,term1], col=colors, cex=0.5)
    
    dev.off()
  }
  
  def_res_list[[dataset]] <- res_list
  
  
  if(0){
  gene_lists_ENTREZ <- lapply(modulated_gene_list, function(x)y<- bitr(x, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID)
  gene_lists_ENTREZ <- lapply(gene_lists_ENTREZ, noquote)
  
  plot_clusterProfile_GO(gene_lists_ENTREZ, "limma_modulated_genes")
  }
  #GOA <- goana(gene_lists_ENTREZ)
  #topup <- topGO(GOA, sort = "N_GFP_up")
  #topdown <- topGO(GOA, sort = "N_GFP_down")
  
  # Mean-difference plot
  plotMD(fit.eb,column=1)
  
  # Q-Q plot of moderated t-statistics
  qqt(fit.eb$t[,1],df=fit.eb$df.residual+fit.eb$df.prior)
  abline(0,1)
  
}

#DESeq2 on scaled data ####

if(1){
  
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  dataset <- "deseq2"
  wd <- paste(hwd, dataset, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  scaledcombined <- read.delim(file.path(hwd, "combined_guide_scaledCount_table.tab"), header=T, sep="\t")
  
  
  dim(scaledcombined)
  head(scaledcombined)
  summary(scaledcombined)
  
  NAs <- apply(scaledcombined, 1, function(x) sum(is.na(x)))
  NAss <- which(NAs > 0)
  length(NAss)
  
  ## expand rows for each gene using guide id
  guides <- c("g1", "g2", "g3", "g4")
  
  allGenes <- sort(unique(scaledcombined$GENE))
  allGenes[1:10]
  scaledexpanded <- NULL
  c <- 0
  valid_gene <- NULL
  for(gene in allGenes){
    sub <- as.matrix(scaledcombined[scaledcombined$GENE %in% gene, 3:ncol(scaledcombined)])
    if(nrow(sub) < 4){
      
      if(nrow(sub) == 3){
        arow <- apply(sub, 2, mean) ## replace missing guide with mean of existing guides
        sub <- rbind(sub, arow)
        valid_gene <- c(valid_gene, gene)
        subvector <- as.vector(sub)
        scaledexpanded <- rbind(scaledexpanded, subvector)
      }else{
        c <- c + 1
        print(paste(c, gene, nrow(sub)))
      }
      
      
    }else if(nrow(sub) == 4){
      valid_gene <- c(valid_gene, gene)
      subvector <- as.vector(sub)
      scaledexpanded <- rbind(scaledexpanded, subvector)
    }
    
  }
  
  print(paste("genes with less than 3 guides ", c))
  
  cn <- unlist(lapply(colnames(scaledcombined)[3:ncol(scaledcombined)], function(x) paste(x, guides, sep="_")))
  colnames(scaledexpanded) <- cn
  rownames(scaledexpanded) <- valid_gene
  
  dim(scaledexpanded)
  head(scaledexpanded)
  summary(scaledexpanded)
  

  data_table <- scaledexpanded
  data_table[is.na(data_table)] <- 0
  data_mat <- as.matrix(data_table)
  dim(data_mat)
  
  varVec <- sort(apply(data_mat, 1, var), decreasing=T)

  top <- length(varVec)
  data_table <- data_mat[names(varVec[1:top]), ]
  
  samples <- colnames(data_table)
  ncol <- length(samples)
  conditions <- unlist(lapply(colnames(data_table), function(x) paste(unlist(strsplit(x, split="_"))[1:2], collapse="_")))
  
  s2c <- data.frame(path=samples, conditions, row.names=samples)
  #s2c$conditions <- relevel(s2c$conditions, ref="T0")
  s2c$path <- as.character(s2c$path)
  s2c$conditions <- as.factor(s2c$conditions)
  #s2c <- s2c[order(sample_order),]
  s2c
  
  guide_stat <- sumstats_col(data_table)
  guide_stat_sum <- guide_stat$Sum
  guide_gomean <- geometric.mean(guide_stat_sum)
  
  
  
  ddsMatrix <- DESeqDataSetFromMatrix(round(data_table), 
                                      colData = s2c,
                                      design = ~ conditions)
  class(ddsMatrix)
  
  
  
  #sizeFactors(ddsMatrix) <- guide_stat_sum/guide_gomean
  dds <- DESeq(ddsMatrix, fitType = "local") #  , sfType = "iterate"
  rnames <- resultsNames(dds)
  rnames
  
  if(is.na(merge_GFP)){
    contrasts <- cbind(c(rep("conditions", 4)), c('GFP1', 'GFP2', "GFP3", 'GFP4'), c(rep("GFPm", 4)))
  }else{
    if(merge_GFP){
      contrasts <- cbind(c(rep("conditions", 17)), 
                         c("N_T18", "NSP12_T18", "S_T18", "NSP9_T18", "NSP10_T18", "ORF9B_T18", "NR_T18", "ORF3A_T18", "ORF8_T18", "NSP2_T18", "NSP3_T18", "M_T18", "NSP2_T10", "NSP3_T10", "NR_T10", "N_T10", "NSP12_T10"), 
                         c(rep(c("GFPm_T18", "GFPm_T10"), times=c(12, 5))))
    }else{
      contrasts <- cbind(c(rep("conditions", 17)), 
                         c("N_T18", "NSP12_T18", "S_T18", "NSP9_T18", "NSP10_T18", "ORF9B_T18", "NR_T18", "ORF3A_T18", "ORF8_T18", "NSP2_T18", "NSP3_T18", "M_T18", "NSP2_T10", "NSP3_T10", "NR_T10", "N_T10", "NSP12_T10"), 
                         c(rep(c("GFP1_T18", "GFP1_T10"), times=c(12, 5))))
    }
  }
  
  res_list <- list()
  for(i in 1:nrow(contrasts)){
    contrast <- contrasts[i,]
    print(paste("Processing ", paste(contrast, collapse=" ")))
    res <- results(dds, contrast=contrast, independentFiltering = FALSE, pAdjustMethod = "BH")
    
    summary(res)
    res <- na.omit(res)
    mcols(res)
    
    allup <- res[(res$log2FoldChange > 0),]
    alldown <- res[(res$log2FoldChange < 0),]
    
    print(dim(res)[1])
    print(dim(allup)[1])
    print(dim(alldown)[1])
    
    sigpv <- res[(res$pvalue < pvalue_cutoff),]
    sig <- res[(res$pvalue < pvalue_cutoff),]
    nc <- res[(res$pvalue >= pvalue_cutoff),]
    sigup <- sig[(sig$log2FoldChange > 0),]
    sigdown <- sig[(sig$log2FoldChange < 0),]
    
    uplist <- rownames(sigup)
    downlist <- rownames(sigdown)
    nclist <- rownames(nc)
    
    
    treatSF <- paste(contrast[2], "vs", contrast[3], sep="_")
    
    
    res <- res[order(-res$stat),]
    res_list[[treatSF]] <- res
    ress <- cbind(row.names(res), res$stat)
    
    filename <- paste(treatSF, "DESeq2_results_table.tab", sep="_")
    write.table(res, filename, row.names=T, col.names=NA, quote=F, sep="\t")
    
    write.table(as.matrix(sizeFactors(dds)), paste(treatSF, "Size_factors.txt", sep="_"), row.names = T)
    
    rld <- vst(dds)
    png(paste(treatSF, "PCA_plot.png", sep="_"))
    print(plotPCA(rld, intgroup=c("conditions")))
    dev.off()
    
    if(0){
    run_enrichGO_simpleList(uplist, "BP", paste(treatSF,"upGenes", sep="_"))
    run_enrichGO_simpleList(downlist, "BP", paste(treatSF,"downGenes", sep="_"))
    #run_enrichDO_simpleList(uplist, paste(treatSF,"upGenes", sep="_"))
    #run_enrichDO_simpleList(downlist, paste(treatSF,"downGenes", sep="_"))
    run_gseGO_simple(res, treatSF)
    
    complexes_up <- overlap_with_CORUM(uplist, idMap, corum)
    complexes_down <- overlap_with_CORUM(downlist, idMap, corum)
    write.table(complexes_up, paste(treatSF, "deseq2_up_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
    write.table(complexes_down, paste(treatSF, "deseq2_down_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
    }
    
    filename <- paste(treatSF, "upregulated_gene_list.txt", sep="_")
    writeLines(uplist, filename)
    filename <- paste(treatSF, "downregulated_gene_list.txt", sep="_")
    writeLines(downlist, filename)
    filename <- paste(treatSF, "noChange_gene_list.txt", sep="_")
    writeLines(nclist, filename)
    
    filename <- paste(treatSF, "upregulated_gene_table.tab", sep="_")
    sigup <- sigup[order(sigup$padj),]
    write.table(sigup, filename, row.names=T, col.names=NA, sep="\t")
    filename <- paste(treatSF, "downregulated_gene_table.tab", sep="_")
    sigdown <- sigdown[order(sigdown$padj),]
    write.table(sigdown, filename, row.names=T, col.names=NA, sep="\t")
    filename <- paste(treatSF, "noChange_gene_table.tab", sep="_")
    nc <- nc[order(nc$padj),]
    write.table(nc, filename, row.names=T, col.names=NA, sep="\t")
    
    pdf(paste(treatSF, "MA_plot.pdf", sep="_"))
    plotMA(res, ylim=c(-2,2))
    dev.off()
    
    pdf(paste(treatSF, "raw_count_plot_upgenes.pdf", sep="_"), width=8, height=8)
    
    old.par <- par(mfrow=c(2,2),mar=c(2,4,2,2))
    for (gene in row.names(sigup)){
      plotCounts(dds, gene, intgroup = "conditions")
      p <- sigup[gene, "padj"]
      legend("top", legend=paste("p",round(p, 4), sep="="))
    }
    par(old.par)
    dev.off()
    
    pdf(paste(treatSF, "raw_count_plot_downgenes.pdf", sep="_"), width=8, height=8)
    
    old.par <- par(mfrow=c(2,2),mar=c(2,4,2,2))
    for (gene in row.names(sigdown)){
      plotCounts(dds, gene, intgroup = "conditions")
      p <- sigdown[gene, "padj"]
      legend("top", legend=paste("p",round(p, 4), sep="="))
    }
    par(old.par)
    dev.off()
    
    pdf(paste(treatSF, "vocano_plot_allgenes.pdf", sep="_"), width=8, height=8)
    
    colors <- rep("black", dim(res)[1])
    colors[row.names(res) %in% row.names(sigpv)] <- "cyan"
    colors[row.names(res) %in% row.names(sig)] <- "red"
    plot(res$log2FoldChange, log10(-log10(res$pvalue)), col=colors, ylim=c(-4, 4), 
         xlab="Log2(FoldChange)",
         ylab="Log10(-Log10(Pvalue))")
    par(old.par)
    dev.off()
  }
  
  def_res_list[[dataset]] <- res_list
  
}


## save results as Rdata #####
if(1){
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  names(def_res_list)
  save(def_res_list, file=file.path(hwd, "def_res_list.Rdata"))
  summary(def_res_list)
}



 
## integrative GO analysis ####

if(1){
  library(org.Hs.eg.db)
  library(GO.db)
  library(RColorBrewer)
  
  functionAnalysis <- TRUE
  GOpvalue_cutoff <- 0.05
  pvalue_cutoffs <- c(0.1, 0.05, 0.02, 0.01)
  
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  wd <- paste(hwd, "integrated_results", sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  load(file.path(hwd, "def_res_list.Rdata"))
  
  ana_names <- names(def_res_list)
  print(ana_names)
  
  Limma_on_lfc <- def_res_list[["Limma_on_lfc"]] # list of dataframe
  drugZ <- def_res_list[["drugZ"]] # list of dataframe
  deseq2_on_reads <- def_res_list[["deseq2"]] # list of dataframe
  jacks <- def_res_list[["JACKS"]] # list of dataframe
  
  names(jacks)
  names(Limma_on_lfc)
  names(deseq2_on_reads)
  names(drugZ)
  
  names(jacks) <- paste(names(jacks), "_", sep="")
  names(Limma_on_lfc) <- gsub("_GFP.*", "_", names(Limma_on_lfc))
  names(deseq2_on_reads) <- unlist(lapply(names(deseq2_on_reads), function(x)strsplit(x,split="vs")[[1]][1]))
  names(drugZ) <- paste(names(drugZ), "_", sep="")
  
  tests <- c("drugZ", "Limma", "Deseq2", "jacks")
  all_screens_ <- names(jacks)
  all_screens <- gsub("_$", "", all_screens_)
  
  synth_sig_GO_term <- NULL
  supp_sig_GO_term <- NULL
  
  synth_sig_count_pval <- list()
  supp_sig_count_pval <- list()
  
  synth_sig_count_FDR <- list()
  supp_sig_count_FDR <- list()
  
  
  for(screen in all_screens_){
    
    Limma_on_lfc_df <- Limma_on_lfc[[names(Limma_on_lfc)[grep(screen, names(Limma_on_lfc))]]]
    jacks_df <- as.data.frame(jacks[[names(jacks)[grep(screen, names(jacks))]]])
    drugZ_df <- drugZ[[names(drugZ)[grep(screen, names(drugZ))]]]
    rownames(drugZ_df) <- drugZ_df$GENE
    
    deseq2_on_reads_df <- deseq2_on_reads[[names(deseq2_on_reads)[grep(screen, names(deseq2_on_reads))]]]
    
    print(screen)
    print(dim(Limma_on_lfc_df))
    print(dim(drugZ_df))
    print(dim(deseq2_on_reads_df))
    print(dim(jacks_df))
    head(Limma_on_lfc_df)
    head(drugZ_df)
    head(deseq2_on_reads_df)
    head(jacks_df)
    

    ## pvalue integration for GO enrichment analysis ##
    
    gene_at_pval_synth <- list()
    gene_at_pval_supp <- list()
    
    gene_at_FDR_synth <- list()
    gene_at_FDR_supp <- list()
    
    for(pvalue_cutoff in pvalue_cutoffs){
      print(pvalue_cutoff)
      drugZPval_synth <- rownames(drugZ_df[drugZ_df$pval_synth < pvalue_cutoff,])
      limmaPval_synth <- rownames(Limma_on_lfc_df[Limma_on_lfc_df$P.Value < pvalue_cutoff & Limma_on_lfc_df$t < 0,])
      deseq2Pval_synth <- rownames(deseq2_on_reads_df[deseq2_on_reads_df$pvalue < pvalue_cutoff & deseq2_on_reads_df$stat < 0,])
      jacksPval_synth <- rownames(jacks_df[jacks_df$neg_p < pvalue_cutoff,])
      
      drugZFDR_synth <- rownames(drugZ_df[drugZ_df$fdr_synth < pvalue_cutoff,])
      limmaFDR_synth <- rownames(Limma_on_lfc_df[Limma_on_lfc_df$adj.P.Val < pvalue_cutoff & Limma_on_lfc_df$t < 0,])
      deseq2FDR_synth <- rownames(deseq2_on_reads_df[deseq2_on_reads_df$padj < pvalue_cutoff & deseq2_on_reads_df$stat < 0,])
      jacksFDR_synth <- rownames(jacks_df[jacks_df$neg_padj < pvalue_cutoff,])
      
      
      alist <- list("drugZ"=drugZPval_synth, "Limma"=limmaPval_synth, "Deseq2"=deseq2Pval_synth, "jacks"=jacksPval_synth)
      atitle <- paste("synthetic lethal genes at pvalue < ", pvalue_cutoff, sep="")
      plotVenn(alist, paste(screen, pvalue_cutoff, "synth_pval_venn_plot.png",sep=""), atitle)
      
      gene_at_pval_synth[["drugZ"]][[as.character(pvalue_cutoff)]] <- drugZPval_synth
      gene_at_pval_synth[["Limma"]][[as.character(pvalue_cutoff)]] <- limmaPval_synth
      gene_at_pval_synth[["Deseq2"]][[as.character(pvalue_cutoff)]] <- deseq2Pval_synth
      gene_at_pval_synth[["jacks"]][[as.character(pvalue_cutoff)]] <- jacksPval_synth 
      
      gene_at_FDR_synth[["drugZ"]][[as.character(pvalue_cutoff)]] <- drugZFDR_synth
      gene_at_FDR_synth[["Limma"]][[as.character(pvalue_cutoff)]] <- limmaFDR_synth
      gene_at_FDR_synth[["Deseq2"]][[as.character(pvalue_cutoff)]] <- deseq2FDR_synth
      gene_at_FDR_synth[["jacks"]][[as.character(pvalue_cutoff)]] <- jacksFDR_synth
      
      
      if(functionAnalysis){
        run_enrichGO_simpleList(drugZPval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_drugZ_synth", sep=""))
        run_enrichGO_simpleList(limmaPval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_limma_synth", sep=""))
        run_enrichGO_simpleList(deseq2Pval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_deseq2_synth", sep=""))
        run_enrichGO_simpleList(jacksPval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_jacks_synth", sep=""))
      }
      
     
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_jacks_synth_genes.tab",sep=""))){
        GO_jacks <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_jacks_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_jacks[GO_jacks$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_drugZ_synth_genes.tab",sep=""))){
        GO_drugZ <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_drugZ_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_drugZ[GO_drugZ$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_deseq2_synth_genes.tab",sep=""))){
        GO_deseq2 <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_deseq2_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_deseq2[GO_deseq2$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_limma_synth_genes.tab",sep=""))){
        GO_limma <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_limma_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_limma[GO_limma$p.adjust < GOpvalue_cutoff, "ID"])
      }
      
      
      
      
      ## supp genes ##
      drugZPval_supp <- rownames(drugZ_df[drugZ_df$pval_supp < pvalue_cutoff,])
      limmaPval_supp <- rownames(Limma_on_lfc_df[Limma_on_lfc_df$P.Value < pvalue_cutoff & Limma_on_lfc_df$t > 0,])
      deseq2Pval_supp <- rownames(deseq2_on_reads_df[deseq2_on_reads_df$pvalue < pvalue_cutoff & deseq2_on_reads_df$stat > 0,])
      jacksPval_supp <- rownames(jacks_df[jacks_df$pos_p < pvalue_cutoff,])
      
      drugZFDR_supp <- rownames(drugZ_df[drugZ_df$fdr_supp < pvalue_cutoff,])
      limmaFDR_supp <- rownames(Limma_on_lfc_df[Limma_on_lfc_df$adj.P.Value < pvalue_cutoff & Limma_on_lfc_df$t > 0,])
      deseq2FDR_supp <- rownames(deseq2_on_reads_df[deseq2_on_reads_df$padh < pvalue_cutoff & deseq2_on_reads_df$stat > 0,])
      jacksFDR_supp <- rownames(jacks_df[jacks_df$pos_padj < pvalue_cutoff,])
      
      alist <- list("drugZ"=drugZPval_supp, "Limma"=limmaPval_supp, "Deseq2"=deseq2Pval_supp, "jacks"=jacksPval_supp)
      atitle <- paste("suppression genes at pvalue < ",pvalue_cutoff,sep="")
      plotVenn(alist, paste(screen, pvalue_cutoff, "supp_pval_venn_plot.png",sep=""), atitle)
      
      gene_at_pval_supp[["drugZ"]][[as.character(pvalue_cutoff)]] <- drugZPval_supp
      gene_at_pval_supp[["Limma"]][[as.character(pvalue_cutoff)]] <- limmaPval_supp
      gene_at_pval_supp[["Deseq2"]][[as.character(pvalue_cutoff)]] <- deseq2Pval_supp
      gene_at_pval_supp[["jacks"]][[as.character(pvalue_cutoff)]] <- jacksPval_supp
      
      gene_at_FDR_supp[["drugZ"]][[as.character(pvalue_cutoff)]] <- drugZFDR_supp
      gene_at_FDR_supp[["Limma"]][[as.character(pvalue_cutoff)]] <- limmaFDR_supp
      gene_at_FDR_supp[["Deseq2"]][[as.character(pvalue_cutoff)]] <- deseq2FDR_supp
      gene_at_FDR_supp[["jacks"]][[as.character(pvalue_cutoff)]] <- jacksFDR_supp
      
      
      if(functionAnalysis){
        run_enrichGO_simpleList(drugZPval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_drugZ_supp", sep=""))
        run_enrichGO_simpleList(limmaPval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_limma_supp", sep=""))
        run_enrichGO_simpleList(deseq2Pval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_deseq2_supp", sep=""))
        run_enrichGO_simpleList(jacksPval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_jacks_supp", sep=""))
      }
      
      
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_jacks_supp_genes.tab",sep=""))){
        GO_jacks <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_jacks_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_jacks[GO_jacks$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_drugZ_supp_genes.tab",sep=""))){
        GO_drugZ <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_drugZ_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_drugZ[GO_drugZ$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_deseq2_supp_genes.tab",sep=""))){
        GO_deseq2 <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_deseq2_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_deseq2[GO_deseq2$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_limma_supp_genes.tab",sep=""))){
        GO_limma <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_limma_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_limma[GO_limma$p.adjust < GOpvalue_cutoff, "ID"])
      }
      
    } # end for pvalue_cutoff
    
    
    for(test in tests){
      synth_pval_list <- gene_at_pval_synth[[test]]
      supp_pval_list <- gene_at_pval_supp[[test]]
      
      synth_sig_count_pval[[test]][[screen]] <- lapply(synth_pval_list, length)
      supp_sig_count_pval[[test]][[screen]] <- lapply(supp_pval_list, length)
      
      synth_FDR_list <- gene_at_FDR_synth[[test]]
      supp_FDR_list <- gene_at_FDR_supp[[test]]
      
      synth_sig_count_FDR[[test]][[screen]] <- lapply(synth_FDR_list, length)
      supp_sig_count_FDR[[test]][[screen]] <- lapply(supp_FDR_list, length)
      
      if(0){
      gene_lists_ENTREZ <- lapply(synth_pval_list, function(x)y<- bitr(x, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID)
      gene_lists_ENTREZ <- lapply(gene_lists_ENTREZ, noquote)
      
      try(
      plot_clusterProfile_GO(gene_lists_ENTREZ, paste(screen, test, "_pval_synth_genes", sep=""))
      )
      
      gene_lists_ENTREZ <- lapply(supp_pval_list, function(x)y<- bitr(x, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID)
      gene_lists_ENTREZ <- lapply(gene_lists_ENTREZ, noquote)
      try(
      plot_clusterProfile_GO(gene_lists_ENTREZ, paste(screen, test, "_pval_supp_genes", sep=""))
      )
      } # end if(functionAnalysis)
    } # end for test
    
    
    system('rm *.log')
  } ## for screen
  synth_sig_GO_terms <- na.omit(synth_sig_GO_term)
  supp_sig_GO_terms <- na.omit(supp_sig_GO_term)
  sig_GO_terms <- list("synth"=synth_sig_GO_terms, "supp"=supp_sig_GO_terms)
  sig_count_pval <- list("synth"=synth_sig_count_pval, "supp"=supp_sig_count_pval)
  sig_count_FDR <- list("synth"=synth_sig_count_FDR, "supp"=supp_sig_count_FDR)
  
  
  save(sig_GO_terms, file=file.path(hwd, "sig_GO_terms.Rdata"))
  save(sig_count_pval, file=file.path(hwd, "sig_count_pval.Rdata"))
  save(sig_count_FDR, file=file.path(hwd, "sig_count_FDR.Rdata"))
  
}
  
######## plotting GO results ####
if(1){
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  load(file.path(hwd, "sig_GO_terms.Rdata"))
  load(file.path(hwd, "sig_count_pval.Rdata"))
  load(file.path(hwd, "sig_count_FDR.Rdata"))
  
  
  nscreen <- length(all_screens)
  npval <- length(pvalue_cutoffs)
  ntest <- length(tests)
  
  pdf(paste("sig_gene_counts_by_test.pdf", sep="_"), width=12, height=10)
  par(mfrow=c(2, ntest))
  cols <- brewer.pal(nscreen,"Paired")
  for(ex in c("synth", "supp")){
    for(test in tests){
      sig_counts_pval <- matrix(unlist(sig_count_pval[[ex]][[test]]), nrow=npval, ncol=nscreen, byrow=F, dimnames=list(pvalue_cutoffs, all_screens))
      
      matplot(x=pvalue_cutoffs, y=sig_counts_pval, lwd=2, xlab="p value", ylab="Number of genes", lty=1, col=cols, ylim=c(0,2000), xlim=c(max(pvalue_cutoffs), min(pvalue_cutoffs)), type="b", pch=2, main=paste(ex,test), log="x")
      legend("top", legend=all_screens, lwd=2, col=cols, bty="n")
    }
    for(test in tests){
      sig_counts_FDR <- matrix(unlist(sig_count_FDR[[ex]][[test]]), nrow=npval, ncol=nscreen, byrow=F, dimnames=list(pvalue_cutoffs, all_screens))
      
      matplot(x=pvalue_cutoffs, y=sig_counts_FDR, lwd=2, xlab="FDR", ylab="Number of genes", lty=1, col=cols, ylim=c(0,300), xlim=c(max(pvalue_cutoffs), min(pvalue_cutoffs)), type="b", pch=2, main=paste(ex,test), log="x")
      legend("top", legend=all_screens, lwd=2, col=cols, bty="n")
    }
  }
  dev.off()
  
  pdf(paste("sig_gene_counts_by_screen.pdf", sep="_"), width=20, height=10)
  par(mfrow=c(2, nscreen))
  cols <- brewer.pal(ntest,"Paired")
  for(ex in c("synth", "supp")){
    for(screen in all_screens){
      count_pval <- NULL
      for(test in tests){
        sig_counts_pval <- matrix(unlist(sig_count_pval[[ex]][[test]]), nrow=npval, ncol=nscreen, byrow=F, dimnames=list(pvalue_cutoffs, all_screens))
        count_pval <- cbind(count_pval, sig_counts_pval[, screen])      
      }
      
      matplot(x=pvalue_cutoffs, y=count_pval, lwd=2, xlab="p value", ylab="Number of genes", lty=1, col=cols, ylim=c(0,2000), xlim=c(max(pvalue_cutoffs), min(pvalue_cutoffs)), type="b", pch=2, main=paste(ex,screen), log="x")
      legend("top", legend=tests, lwd=2, col=cols, bty="n")
    }
    for(screen in all_screens){
      count_FDR <- NULL
      for(test in tests){
        sig_counts_FDR <- matrix(unlist(sig_count_FDR[[ex]][[test]]), nrow=npval, ncol=nscreen, byrow=F, dimnames=list(pvalue_cutoffs, all_screens))
        count_FDR <- cbind(count_FDR, sig_counts_FDR[, screen]) 
      }
      
      matplot(x=pvalue_cutoffs, y=count_FDR, lwd=2, xlab="FDR", ylab="Number of genes", lty=1, col=cols, ylim=c(1,300), xlim=c(max(pvalue_cutoffs), min(pvalue_cutoffs)), type="b", pch=2, main=paste(ex,screen), log="x")
      legend("top", legend=tests, lwd=2, col=cols, bty="n")
    }
  }
  dev.off()
  
  lapply(sig_GO_terms, length)
  GO_term_padj <- list()
  
  for(ex in c("synth", "supp")){
    terms <- sig_GO_terms[[ex]]
    
    
    for(screen in all_screens_){
      
      for(term in terms){
        #go_id = GOID( GOTERM[ Term(GOTERM) == term])
        #allegs = get(go_id, org.Hs.egGO2ALLEGS)
        #icr_genes = unlist(mget(allegs,org.Hs.egSYMBOL))
        description <- Term(GOTERM[GOID(GOTERM) == term])
        print(paste("processing ", ex, screen, description))
        
        padj_GO <- NULL
        
        for(pvalue_cutoff in pvalue_cutoffs){
          
          GO_drugZ_termp <- GO_limma_termp <- GO_deseq2_termp <- GO_jacks_termp <- 0
          
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_jacks_",ex,"_genes.tab",sep=""))){
            GO_jacks <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_jacks_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_jacks$ID){
              GO_jacks_termp <- -log(GO_jacks[GO_jacks$ID == term, "p.adjust"])
            }
          }
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_drugZ_",ex,"_genes.tab",sep=""))){
            GO_drugZ <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_drugZ_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_drugZ$ID){
              GO_drugZ_termp <- -log(GO_drugZ[GO_drugZ$ID == term, "p.adjust"])
            }
          }
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_limma_",ex,"_genes.tab",sep=""))){
            GO_limma <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_limma_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_limma$ID){
              GO_limma_termp <- -log(GO_limma[GO_limma$ID == term, "p.adjust"])
            }
          }
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_deseq2_",ex,"_genes.tab",sep=""))){
            GO_deseq2 <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_deseq2_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_deseq2$ID){
              GO_deseq2_termp <- -log(GO_deseq2[GO_deseq2$ID == term, "p.adjust"])
            }
          }
          
          padj_GO <- rbind(padj_GO, c(GO_drugZ_termp, GO_limma_termp, GO_deseq2_termp, GO_jacks_termp))
          
        } ## for pvalue_cutoff 
        
        GO_term_padj[[ex]][[screen]][[term]] <- padj_GO
        
      } ## for term
      
    } ## for screen
  } ## for ex
  
  GO_term_padj[["synth"]][["NSP9_"]]
  
  for(ex in c("synth", "supp")){
    terms <- sig_GO_terms[[ex]]
    
    pdf(paste(ex,"GO_padj_profile.pdf", sep="_"), width=15, height=10)
    par(mfrow=c(3, nscreen))
    cols <- brewer.pal(ntest,"Paired")
    for(term in terms){
      description <- Term(GOTERM[GOID(GOTERM) == term])
      for(screen in all_screens_){
        pval_df <- GO_term_padj[[ex]][[screen]][[term]]
        maintext <- gsub("_", "", screen)
        maintext <- paste( " ", maintext, sep="\n")
        if(screen == all_screens_[1]){
          maintext <- paste(description, screen, sep="\n")
        }
        matplot(x=pvalue_cutoffs, y=pval_df[, 1:ntest], lwd=2, lty=1, col=cols,  ylab="-log(pvalue)", ylim=c(0,20), xlim=c(max(pvalue_cutoffs), min(pvalue_cutoffs)), type="b", pch=2, main=maintext, log="x")
        abline(h=-log(GOpvalue_cutoff), lty=2) # -log(0.01)
        legend("top", legend=c(tests), lwd=2, col=cols, bty="n")
      }
    }
    
    dev.off()
    
    for(test in tests){
      pdf(paste(ex,"GO_padj_profile", test, ".pdf", sep="_"), width=12, height=10)
      par(mfrow=c(2, 2))
      cols <- brewer.pal(nscreen,"Paired")
      terms_padjmean <- NULL
      descriptions <- NULL
      for(term in terms){
        description <- Term(GOTERM[GOID(GOTERM) == term])
        descriptions <- c(descriptions, description)
        test_df <- NULL
        for(screen in all_screens_){
          pval_df <- GO_term_padj[[ex]][[screen]][[term]]
          colnames(pval_df) <- tests
          test_df <- cbind(test_df, pval_df[,test])
          
        }
        test_df_colmean <- apply(test_df, 2, mean)
        test_df_colmedian <- apply(test_df, 2, median)
        for(ind in 1:length(test_df_colmean)){
          test_df_colmean[ind] <- max(test_df_colmean[ind], test_df_colmedian[ind])
        }
        terms_padjmean <- rbind(terms_padjmean, test_df_colmean)
        test_df[is.na(test_df)] <- 0
        if(max(test_df) > -log(GOpvalue_cutoff)){
          matplot(x=pvalue_cutoffs, y=test_df[, 1:nscreen], lwd=2, lty=1, col=cols, ylab="-log(pvalue)", ylim=c(0,20), xlim=c(max(pvalue_cutoffs), min(pvalue_cutoffs)), type="b", pch=2, main=description, log="x")
          abline(h=-log(GOpvalue_cutoff), lty=2) # -log(0.01)
          legend("top", legend=all_screens, lwd=2, col=cols, bty="n")
        }
      }
      
      dev.off()
      
      rownames(terms_padjmean) <- terms
      colnames(terms_padjmean) <- all_screens
      
      terms_padjmean <- na.omit(terms_padjmean)
      
      terms_padjmean
      terms_padjmean <- terms_padjmean[apply(terms_padjmean, 1, max) > -log(GOpvalue_cutoff), ]
      
      if(!is.null(dim(terms_padjmean)) && nrow(terms_padjmean) > 3){
      
        rowOrder <- arrangeGOterms(rownames(terms_padjmean))
        
        
        terms_padjmean <- terms_padjmean[names(rowOrder),]
        rownames(terms_padjmean) <- rowOrder
        
        while(!is.null(dev.list())){
          dev.off()
        }
        pdf(paste(ex,"GO_padj_heatmap", test, ".pdf", sep="_"), width=10, height=12)
        hm <- pheatmap(terms_padjmean, cluster_rows = F, silent=T)
        grid.draw(hm)
        
        dev.off()
        
        png(paste(ex,"GO_padj_heatmap", test, ".png", sep="_"), width=10, height=12, res=300, units="in")
        hm <- pheatmap(terms_padjmean, cluster_rows = F, silent=T)
        grid.draw(hm)
        
        dev.off()
      }
    } ## for test
    
  } ## for ex
  
} ## for (1)

#merge_GFP <- TRUE
## integrative overlap analysis #####
for(merge_GFP in c(FALSE, TRUE)){
  
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  wd <- paste(hwd, "Comparative_analysis", sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  
  load(file.path(hwd, "def_res_list.RData"))
  
  ana_names <- names(def_res_list)
  print(ana_names)
  
  Limma_on_lfc <- def_res_list[["Limma_on_lfc"]] # list of dataframe
  drugZ <- def_res_list[["drugZ"]] # list of dataframe
  deseq2_on_reads <- def_res_list[["deseq2"]] # list of dataframe
  jacks <- def_res_list[["JACKS"]] # list of dataframe
  
  names(jacks)
  names(Limma_on_lfc)
  names(deseq2_on_reads)
  names(drugZ)
  
  names(jacks) <- paste(names(jacks), "_", sep="")
  names(Limma_on_lfc) <- gsub("_GFP.*", "_", names(Limma_on_lfc))
  names(deseq2_on_reads) <- unlist(lapply(names(deseq2_on_reads), function(x) unlist(strsplit(x, split="vs"))[1]))
  names(drugZ) <- paste(names(drugZ), "_", sep="")
  
  Limma_on_lfc_t_table <- NULL
  deseq2_on_reads_stat_table <- NULL
  jacks_neg_table <- NULL
  drugz_normz_table <- NULL
  
  tests <- c("drugZ", "Limma", "Deseq2", "jacks")
  all_screens_ <- names(jacks)
  all_screens <- gsub("_$", "", all_screens_)
  top <- 10
  
  top_gene_list <- list() ## by screen
  top_gene_list_by_test <- list() ## by test
  top_gene_list_by_screentest <- list() ## by screen and test
  
  pval_df_list <- list()
  for(screen in all_screens_){
    #screen <- all_screens_[1]  
    Limma_on_lfc_df <- Limma_on_lfc[[names(Limma_on_lfc)[grep(screen, names(Limma_on_lfc))]]]
    jacks_df <- as.data.frame(jacks[[names(jacks)[grep(screen, names(jacks))]]])
    drugZ_df <- drugZ[[names(drugZ)[grep(screen, names(drugZ))]]]
    rownames(drugZ_df) <- drugZ_df$GENE
    
    deseq2_on_reads_df <- deseq2_on_reads[[names(deseq2_on_reads)[grep(screen, names(deseq2_on_reads))]]]
    
    print(screen)
    print(dim(Limma_on_lfc_df))
    print(dim(drugZ_df))
    print(dim(deseq2_on_reads_df))
    print(dim(jacks_df))
    head(Limma_on_lfc_df)
    head(drugZ_df)
    head(deseq2_on_reads_df)
    head(jacks_df)
    
    Limma_on_lfc_df <- Limma_on_lfc_df %>%
        mutate(P.value_synth=P.Value, P.value_supp=P.Value)
    Limma_on_lfc_df$P.value_synth[Limma_on_lfc_df$t>0] <- 1 - Limma_on_lfc_df$P.Value[Limma_on_lfc_df$t>0]
    Limma_on_lfc_df$P.value_supp[Limma_on_lfc_df$t<=0] <- 1 - Limma_on_lfc_df$P.Value[Limma_on_lfc_df$t<=0]
    
    deseq2_on_reads_df <- as.data.frame(deseq2_on_reads_df) %>%
      mutate(P.value_synth=pvalue, P.value_supp=pvalue)
    deseq2_on_reads_df$P.value_synth[deseq2_on_reads_df$stat>0] <- 1 - deseq2_on_reads_df$pvalue[deseq2_on_reads_df$stat>0]
    deseq2_on_reads_df$P.value_supp[deseq2_on_reads_df$stat<=0] <- 1 - deseq2_on_reads_df$pvalue[deseq2_on_reads_df$stat<=0]
    
    synth_pval_combined <- cbind(drugZ_df[TKO3_genes, "pval_synth"],
                                 Limma_on_lfc_df[TKO3_genes, "P.value_synth"],
                                 deseq2_on_reads_df[TKO3_genes, "P.value_synth"],
                                 jacks_df[TKO3_genes, "neg_p"])
    supp_pval_combined <- cbind(drugZ_df[TKO3_genes, "pval_supp"],
                                Limma_on_lfc_df[TKO3_genes, "P.value_supp"],
                                deseq2_on_reads_df[TKO3_genes, "P.value_supp"],
                                jacks_df[TKO3_genes, "pos_p"])
    
    colnames(synth_pval_combined) <- colnames(supp_pval_combined) <- tests
    rownames(synth_pval_combined) <- rownames(supp_pval_combined) <- TKO3_genes
    
    synth_pval_combined_filtered <- na.omit(synth_pval_combined)
    supp_pval_combined_filtered <- na.omit(supp_pval_combined)
    
    #pheatmap(synth_pval_combined_filtered)
    #pheatmap(supp_pval_combined_filtered)
    
    pval_df_list[["supp"]][[screen]] <- supp_pval_combined_filtered
    pval_df_list[["synth"]][[screen]] <- synth_pval_combined_filtered
    
    for(test in tests){
      #test <- "drugZ"
      screentest <- paste(screen, test, sep="")
      top_gene_list[["supp"]][[screen]][[test]] <- rownames(supp_pval_combined_filtered[order(supp_pval_combined_filtered[,test]), ])[1:top]
      top_gene_list[["synth"]][[screen]][[test]] <- rownames(synth_pval_combined_filtered[order(synth_pval_combined_filtered[,test]), ])[1:top]
      top_gene_list_by_test[["supp"]][[test]][[screen]] <- rownames(supp_pval_combined_filtered[order(supp_pval_combined_filtered[,test]), ])[1:top]
      top_gene_list_by_test[["synth"]][[test]][[screen]] <- rownames(synth_pval_combined_filtered[order(synth_pval_combined_filtered[,test]), ])[1:top]
      top_gene_list_by_screentest[["supp"]][[screentest]] <- rownames(supp_pval_combined_filtered[order(supp_pval_combined_filtered[,test]), ])[1:top]
      top_gene_list_by_screentest[["synth"]][[screentest]] <- rownames(synth_pval_combined_filtered[order(synth_pval_combined_filtered[,test]), ])[1:top]
    }
    
    ## collect stats for each test
    Limma_on_lfc_t_table <- cbind(Limma_on_lfc_t_table, Limma_on_lfc_df[TKO3_genes, "t"])
    deseq2_on_reads_stat_table <- cbind(deseq2_on_reads_stat_table, deseq2_on_reads_df[TKO3_genes, "stat"])
    jacks_neg_table <- cbind(jacks_neg_table, jacks_df[TKO3_genes, "neg"])
    drugz_normz_table <- cbind(drugz_normz_table, drugZ_df[TKO3_genes, "normZ"])
  }
  
  stat_table_list <- list()
  stat_table_list[["deseq2_stat"]] <- deseq2_on_reads_stat_table
  stat_table_list[["Limma_t"]] <- Limma_on_lfc_t_table
  stat_table_list[["jacks_neg"]] <- jacks_neg_table
  stat_table_list[["drugz_norm"]] <- drugz_normz_table
  
  for(aname in names(stat_table_list)){
    atable <- stat_table_list[[aname]]
    dimnames(atable) <- list(TKO3_genes, all_screens)
    plot_heatmap(na.omit(atable), aname)
  }
  
  
  for(ex in c("synth", "supp")){
    top_list <- top_gene_list_by_screentest[[ex]]
    top_mat <- union2mat(top_list)
    pheatmap(top_mat[,!grepl("Sum", colnames(top_mat))])
    write.table(top_mat, paste(ex, "_allScreenTest_top_gene_matrix.tab",sep=""), sep="\t", col.names=NA, quote=F)
    for(screen in all_screens_){
      #ex <- "synth"
      #screen <- "N_T10_"
      pval_df <- pval_df_list[[ex]][[screen]]
      pval_df["ZNF629",]
      top_list <- top_gene_list[[ex]][[screen]]
      top_mat <- union2mat(top_list)
      pval_mat <- pval_df[rownames(top_mat), colnames(top_mat)[1:length(tests)]]
      geo_mean <- apply(pval_mat, 1, geometric.mean)
      pval_mat <- as.data.frame(cbind(pval_mat, geo_mean))
      
      top_mat <- cbind(top_mat, pval_mat[rownames(top_mat),])
      write.table(top_mat, paste(screen, ex, "_top_gene_matrix.tab",sep=""), sep="\t", col.names=NA, quote=F)
      plotVenn(top_list, paste(screen, ex, "_top_gene_venn_plot.png",sep=""), paste(screen, ex, "_top_genes",sep=""), "text")
    }
    system('rm *.log')
  }
    
  
  for(ex in c("synth", "supp")){
    for(test in tests){
      #ex <- "synth"
      #test <- "drugZ"
      if(!is.na(merge_GFP)){
        screen_pairs_list <- list("N_NR_T18"=c("N_T18_","NR_T18_"), 
                                  "N_NR_T10"=c("N_T10_","NR_T10_"),
                                  "NSP12"=c("NSP12_T18_", "NSP12_T10_"), 
                                  "NSP2"=c("NSP2_T18_", "NSP2_T10_"),
                                  "NSP3"=c("NSP3_T18_", "NSP3_T10_"),
                                  "N"=c("N_T18_", "N_T10_"),
                                  "NR"=c("NR_T18_", "NR_T10_"),
                                  "N_NR"=c("N_T18_", "N_T10_", "NR_T18_", "NR_T10_"))
      }else{
        screen_pairs_list <- list("GFP1_2"=c("GFP1_","GFP2_"), "GFP1_3"=c("GFP1_","GFP3_"), "GFP1_4"=c("GFP1_","GFP4_"))
      }
      
      for(apair in names(screen_pairs_list)){
        #apair <- "N"
        common_list <- list()
        screen_pair <- screen_pairs_list[[apair]]
        for(screen in screen_pair){
          
          top_100 <- top_gene_list[[ex]][[screen]][[test]]
          common_list[[screen]] <- top_100
          
        }
        top_mat <- union2mat(common_list)
        write.table(top_mat, paste(apair, ex, test, "top_gene_matrix.tab",sep="_"), sep="\t", col.names=NA, quote=F)
        plotVenn(common_list, paste(apair, ex, test, "top_gene_venn_plot.png",sep="_"), paste(apair, ex, "top_genes",sep="_"), "text")
        
      }  
      system('rm *.log')
    }
  }
  
 
  
  for(ex in c("synth", "supp")){
    for(test in tests){ 
      top_genes <- top_gene_list_by_test[[ex]][[test]]
      top_mat <- union2mat(top_genes)
      cn <- colnames(top_mat)
      for(screen in all_screens_){
        #ex <- "synth"
        #screen <- "N_"
        pval_df <- pval_df_list[[ex]][[screen]]
        top_mat <- cbind(top_mat, pval_df[rownames(top_mat), test])
        cn <- c(cn, paste(screen, "pval", sep=""))
      }
        
      colnames(top_mat) <- cn
      write.table(top_mat, paste(ex, test, "top_gene_in_allScreens_matrix.tab",sep="_"), sep="\t", col.names=NA, quote=F)
      if(length(top_genes) < 6){
        plotVenn(top_genes, paste(ex, test, "top_gene_in_allScreens_venn_plot.png",sep="_"), paste(ex, test, "top_genes",sep="_"), "text")
      }else{
        five_screens <- top_genes[1:5]
        plotVenn(five_screens, paste(ex, test, "top_gene_in_5Screens_venn_plot.png",sep="_"), paste(ex, test, "top_genes",sep="_"), "text")
      }
      system('rm *.log')
      
    }
    
  }
  if(0){
    
    pdf(paste(screenn, "sig_gene_data_plot_supp.png", sep="_"), width=11, height=8)
    par(mfrow=c(3,3))
    for(agene in top_gene_supp){
      inspect_data_plot(combined_table, agene, control_prefix, treat_prefix, "Raw read count")
      inspect_data_plot(T18_expanded, agene, control_prefix, treat_prefix, "scaled read count")
      inspect_data_plot(lfc_expanded, agene, control_prefix, treat_prefix, "Log2(folderChange)")
    }
    par(mfrow=c(1,1))
    dev.off()
    
    pdf(paste(screenn, "sig_gene_data_plot_synth.pdf", sep="_"), width=11, height=8)
    par(mfrow=c(3,3))
    for(agene in top_gene_synth){
      inspect_data_plot(combined_table, agene, control_prefix, treat_prefix, "Raw read count")
      inspect_data_plot(T18_expanded, agene, control_prefix, treat_prefix, "scaled read count")
      inspect_data_plot(lfc_expanded, agene, control_prefix, treat_prefix, "Log2(folderChange)")
    }
    par(mfrow=c(1,1))
    dev.off()
  }
    
  
}


## DRACO ####
if(1){
  merge_GFP <- TRUE
  if(is.na(merge_GFP)){
    hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_only"
  }else{
    if(merge_GFP){
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged" 
    }else{
      hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual" 
    }
  }
  
  dataset <- "DRACO"
  wd <- paste(hwd, dataset, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  
  raw_combined <- read.delim(file.path(hwd, "combined_guide_rawCount_table.tab"), header=T, sep="\t")
  lfc_combined <- read.delim(file.path(hwd, "combined_guide_lfc_table.tab"), header=T, sep="\t")
  boxplot(raw_combined[, c(-1,-2)])
  boxplot(lfc_combined[, c(-1,-2)])
  dim(lfc_combined)
  head(lfc_combined)
  summary(lfc_combined)
  
  replicates <- colnames(lfc_combined)[3:ncol(lfc_combined)]
  subjects <- unique(unlist(lapply(replicates, function(x)unlist(strsplit(x, split="_"))[1])))
  control <- "GFPm"
  subjects <- subjects[!subjects %in% control]
  
  lfcOverControl <- lfc_combined[, colnames(lfc_combined)[!grepl(control, colnames(lfc_combined))]]
  for(screen in subjects){
    for(rep in c("T18_A", "T18_B", "T18_C")){
      control_col <- replicates[grepl(paste(control, rep, sep="_"), replicates)]
      screen_col <- replicates[grepl(paste(screen, rep, sep="_"), replicates)]
      lfcOverControl[, screen_col] <- lfc_combined[, screen_col] - lfc_combined[, control_col]
    }
  }
}

#N_col <- replicates[grepl("N_", replicates)]
#chart.Correlation(log2fc[,N_col])
#chart.Correlation(norm_combined[,N_col])
summary(lfcOverControl)
boxplot(lfcOverControl[, c(-1,-2)])
l## expand rows for each gene using guide id
guides <- c("g1", "g2", "g3", "g4")

allGenes <- sort(unique(log2fc$GENE))
allGenes[1:10]
lfc_expanded <- NULL
c <- 0
valid_gene <- NULL
for(gene in allGenes){
  gene <- allGenes[10]
  sub <- as.matrix(log2fc[log2fc$GENE %in% gene, 3:ncol(log2fc)])
  if(nrow(sub) < 4){
    
    if(nrow(sub) == 3){
      arow <- apply(sub, 2, mean) ## replace missing guide with mean of existing guides
      sub <- rbind(sub, arow)
      valid_gene <- c(valid_gene, gene)
      subvector <- as.vector(sub)
      lfc_expanded <- rbind(lfc_expanded, subvector)
    }else{
      c <- c + 1
      print(paste(c, gene, nrow(sub)))
    }
    
    
  }else if(nrow(sub) == 4){
    valid_gene <- c(valid_gene, gene)
    subvector <- as.vector(sub)
    lfc_expanded <- rbind(lfc_expanded, subvector)
  }
  
}

print(paste("genes with less than 3 guides ", c))

cn <- unlist(lapply(colnames(log2fc)[3:ncol(log2fc)], function(x) paste(x, guides, sep="_")))
colnames(lfc_expanded) <- cn
rownames(lfc_expanded) <- valid_gene

dim(lfc_expanded)
head(lfc_expanded)
summary(lfc_expanded)

## Merged_vs_individual ####

hwd_merged <- "C:/data/raw/EDYTA/dropoutScreen/GFP_merged"
hwd_individual <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual"

load(file.path(hwd_merged, "def_res_list.Rdata"))
res_merged <- def_res_list
load(file.path(hwd_individual, "def_res_list.Rdata"))
res_individual <- def_res_list
top100_list <- list()
top100_list[["merged"]] <- head(res_merged$drugZ$S_T18$GENE, 100)
top100_list[["individual"]] <- head(res_individual$drugZ$S_T18$GENE, 100)
venn(top100_list)

combined_table <- merge(res_merged$drugZ$S_T18, res_individual$drugZ$N_T18, by=c("GENE"))

ggscatterhist(combined_table, x="normZ.x", y="normZ.y", add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x = -5, label.sep = "\n"),
              margin.plot = "histogram"
)

## BAGEL analysis ####

if(1){
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_individual"

  dataset <- "bagel"
  wd <- file.path(hwd, dataset)
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  gsP <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/training_essentials_CRISPR.txt")
  gsN <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/training_nonessential.txt")
  
  ## this table is a supp material for the paper
  ## Evaluation and Design of Genome-Wide CRISPR/SpCas9 Knockout Screens, Hart et al.
  ## published in G3 Genes|Genomes|Genetics, August 2017
  print("Processing G3 table")
  
  G3_table3 <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/2719TableS3.txt", header=T)
  rownames(G3_table3) <- G3_table3$GENE
  G3_table3 <- G3_table3[, -1]
  G3_table3_NA <- G3_table3
  G3_table3[is.na(G3_table3)] <- 0 ## set NA to 0
  brks <- c(-3500, seq(-50, 50, by=5), 1000)
  for(cn in colnames(G3_table3)){
    hist(as.numeric(G3_table3[,cn]), breaks=brks, xlim=c(-50, 50), main=cn)
  }
  
  
  gsP <- update_gene_names(name_map, gsP$Gene)
  gsN <- update_gene_names(name_map, gsN$Gene)
  
  
  gsp <- intersect(gsP, rownames(G3_table3))
  gsn <- intersect(gsN, rownames(G3_table3))
  
  ## G3_table3 not related for now
  if(0){
  fdr_stats <- calc_FDR_cutoff(gsn, gsp, G3_table3)
  
  ## use ROC to find the dataset having low quality, and exclude them from further analysis
  pdf("G3_ROC_plot.pdf", height=8, width=8)
  par(mfrow = c(3,3))
  colors <- terrain.colors(dim(G3_table3)[2])
  for(i in 1:dim(G3_table3)[2]){
    plot(fdr_stats$FDR[[i]]$Precision, fdr_stats$FDR[[i]]$Recall, type="l", xlim=c(0.4,1), ylim=c(0,1), col=colors[1], lwd=2, main=names(fdr_stats$FDR)[i], xlab="Precision", ylab="Recall")
    abline(v=0.95, lty=1, lwd=1)
    abline(h=0.8, lty=1, lwd=1)
  }
  dev.off()
  par(mfrow = c(1,1))
  
  exclude <- c("BF_HELA", "BF_RPE1")
  G3_essential_list <- within(fdr_stats$essentialGene, rm("BF_HELA"))
  G3_essential_list <- within(G3_essential_list, rm("BF_RPE1"))
  
  unlist(lapply(G3_essential_list, length))
  unlist(lapply(fdr_stats$essentialGene, length))
  
  G3_essential_union <- G3_essential_list[[1]]
  G3_essential_intersect <- G3_essential_list[[1]]
  for(i in 2:length(G3_essential_list)){
    G3_essential_union <- union(G3_essential_union, G3_essential_list[[i]])
    G3_essential_intersect <- intersect(G3_essential_intersect, G3_essential_list[[i]])
  }
  
  length(G3_essential_union)
  length(G3_essential_intersect)
  
  
  
  Frequency <- unlist(lapply(G3_essential_union, function(x)countOccurance(x, G3_essential_list)))
  
  G3_table3_15 <- G3_table3_NA[, names(G3_essential_list)]
  tail(G3_table3_15)
  
  valid_assays <- apply(G3_table3_15, 1, function(x) sum(!is.na(x)))
  Total_assays <- valid_assays[names(Frequency)]
  Ratio <- Frequency/Total_assays
  union_stats <- cbind(Frequency, Total_assays, Ratio)
  union_stats <- union_stats[order(union_stats[,"Ratio"], decreasing=T),]
  union_stats_core <- union_stats[union_stats[,"Ratio"] > 0.5 & union_stats[,"Total_assays"] > 5,]
  dim(union_stats_core)
  tail(union_stats_core)
  write.table(union_stats, "essential_genes_from_15assays.tab", sep="\t", quote=F, col.names=NA)
  }
  
  
  print("Processing BAGEL results")
  
  bagel_results_list <- list()
  
  
  bagel_files <- list.files(pattern="_bagel.bf")
  
  for (bagel_file in bagel_files){
    bf <- read.delim(bagel_file)
    subject <- unlist(strsplit(bagel_file, split="_"))[1]
    bagel_results_list[[subject]] <- bf
  }
  
  bagel_combined <- bagel_results_list[[1]]
  for(i in 2:length(bagel_results_list)){
    print(dim(bagel_results_list[[i]]))
    bagel_combined <- merge(bagel_combined, bagel_results_list[[i]], by=c("GENE"), all=T)
  }
  
  dim(bagel_combined)
  head(bagel_combined)
  
  header <-c("GENE", unlist(lapply(names(bagel_results_list), function(x)paste(x, c("BF", "STD", "NumObs"), sep="_"))))
  colnames(bagel_combined) <- header
  
  bf_table <- bagel_combined[, grepl("_BF", colnames(bagel_combined))]
  rownames(bf_table) <- bagel_combined$GENE
  screens <- unlist(lapply(colnames(bf_table), function(x) unlist(strsplit(x, split="_"))[1]))
  colnames(bf_table) <- screens
  dim(bf_table)
  NAs <- apply(bf_table, 1, function(x) sum(is.na(x))>0)
  NA_table <- bf_table[NAs,]
  bf_table[is.na(bf_table)] <- 0
  
  
  pdf("BF_all_screen_correlation.pdf")
  col0 <- colorRampPalette(c("cyan4", "red3"))
  corrplot(cor(bf_table), method="number", is.corr = F, col=col0(15), type = "upper", order = "hclust", cl.lim=c(0.7, 1),
           tl.col = "black", tl.srt = 45)
  dev.off()
  
  #run_gseGO_SimpleList(GFP_BF, "GFP_bf")
  #run_gseGO_SimpleList(N_BF, "N_bf")
  #run_gseGO_SimpleList(NSP12_BF, "NSP12_bf")
  
  
  #cor.plot(bf_table)
  head(bf_table)
  
  
  gsp <- intersect(gsP, rownames(bf_table))
  gsn <- intersect(gsN, rownames(bf_table))
  
  fdr_stats_bf <- calc_FDR_cutoff(gsn, gsp, bf_table)
  
  essential_list <- fdr_stats_bf$essentialGene
  names(essential_list)
  #essential_list[["G3_Union"]] <- G3_essential_union
  essential_mat <- union2mat(essential_list)
  dim(essential_mat)
  tail(essential_mat)
  
  
  gghistogram(essential_mat, x="Sum")
  screen_count_essential <- as.data.frame(table(essential_mat$Sum))
  colnames(screen_count_essential) <- c("Number of screens", "Number of essential genes")
  pdf("Essential_gene_overlap_plot.pdf", height=8, width=8)
  ggbarplot(screen_count_essential, 
            x="Number of screens", 
            y="Number of essential genes",
            fill = "Number of screens",
            palette = "grey"
            ) + rremove("x.axis") + rremove("x.ticks")  + rremove("legend") +
            font("xy.title", size=16) + font("xy.text", size=16)
  dev.off()
  
  pdf("ROC_plot.pdf", height=8, width=8)
  screens <- unlist(lapply(names(fdr_stats_bf$FDR), function(x) unlist(strsplit(x, split="_"))[1]))
  colors <- terrain.colors(length(screens))
  plot(fdr_stats_bf$FDR$GFP1_BF$Precision, fdr_stats_bf$FDR$GFP1_BF$Recall, type="l", col=colors[1], lwd=2, main="BF score ROC", xlab="Precision", ylab="Recall")
  for(i in 2:length(screens)){
    lines(fdr_stats_bf$FDR[[i]]$Precision, fdr_stats_bf$FDR[[i]]$Recall, col=colors[i], lwd=2)
  }
  abline(v=0.95, lty=2, lwd=1)
  abline(h=0.90, lty=2, lwd=1)
  legend("bottomleft", legend=screens, col=colors, lwd=2)
  dev.off()
  
  if(0){
  pdf("number_of_essential_genes.pdf", height=8, width=8)
  n_essential_gens <- list("G3"=unlist(lapply(G3_essential_list, length)), "Screen"=unlist(lapply(fdr_stats_bf$essentialGene, length)))
  boxplot2(n_essential_gens, ylim=c(0,2000), ylab="Number of genes")
  dev.off()
  
  n_essential_gens_table <- as.data.frame(unlist(n_essential_gens))
  write.table(as.data.frame(unlist(n_essential_gens)), "number_of_essential_genes.tab", col.names=NA, sep="\t", quote=F)
  
  alist <- fdr_stats_bf$essentialGene[c("N_BF", "NSP12_BF", "S_BF", "GFP_BF", "HEK293T_BF")]
  filename <- "venn_all5.png"
  atitle <- "Essential genes in all five screens"
  plotVenn(alist, filename, atitle)
  
  alist <- list()
  alist <- fdr_stats_bf$essentialGene[c("GFP2_BF", "GFP_BF", "HEK293T_BF")]
  
  filename <- "venn_GFP2_GFP_HEK293T.png"
  atitle <- "Essential genes in GFP and HEK293T"
  plotVenn(alist, filename, atitle)
  
  alist <- list()
  alist <- fdr_stats_bf$essentialGene[c("GFP_BF", "HEK293T_BF")]
  alist[["Essential_union"]] <- G3_essential_union
  alist[["Essential_core"]] <- rownames(union_stats_core)
  
  filename <- "venn_HEK293T_GFP_G3Union.png"
  atitle <- "Essential genes in GFP and union of Hart"
  plotVenn(alist, filename, atitle)
  
  
  nlist <- alist
  nlist[["N_BF"]] <- fdr_stats_bf$essentialGene$N_BF
  filename <- "venn_N_GFP_G3Union.png"
  atitle <- "Essential genes in N, GFP and union of Hart"
  plotVenn(nlist, filename, atitle)
  
  nsp12list <- alist
  nsp12list[["NSP12_BF"]] <- fdr_stats_bf$essentialGene$NSP12_BF
  filename <- "venn_NSP12_GFP_G3Union.png"
  atitle <- "Essential genes in NSP12, GFP and union of Hart"
  plotVenn(nsp12list, filename, atitle)
  
  slist <- alist
  slist[["S_BF"]] <- fdr_stats_bf$essentialGene$S_BF
  filename <- "venn_S_GFP_G3Union.png"
  atitle <- "Essential genes in S, GFP and union of Hart"
  plotVenn(slist, filename, atitle)
  
  alist <- fdr_stats_bf$essentialGene[c("N_BF", "NSP12_BF", "S_BF")]
  filename <- "venn_N_NSP12_S.png"
  atitle <- "Essential genes in N, NSP12, S"
  plotVenn(alist, filename, atitle)
  
  system("rm *.log")
  
  #grid.draw(venn.plot);
  grid.newpage();
  control_union <- union(G3_essential_union, fdr_stats_bf$essentialGene$GFP_BF)
  control_union <- union(control_union, fdr_stats_bf$essentialGene$GFP2_BF)
  
  N_wo_control <- fdr_stats_bf$essentialGene$N_BF[!fdr_stats_bf$essentialGene$N_BF %in% control_union]
  length(N_wo_control)
  
  
  bf_table$N_GFP <- bf_table$N_BF - bf_table$GFP_BF
  bf_table <- bf_table[order(bf_table$N_GFP, decreasing=T),]
  write.table(bf_table, paste("bf_all_table.tab", sep="_"), col.names=NA, row.names=T, sep="\t", quote=F)
  
  run_enrichGO_simpleList(N_wo_control, ont="BP", "N_wo_control")
  N_wo_control_complex <- overlap_with_CORUM(N_wo_control, idMap, corum)
  write.table(N_wo_control_complex, paste("N_wo_control_overlap_with_Corum.tab", sep="_"), row.names=F, sep="\t", quote=F)
  }
  
}

### END ####
