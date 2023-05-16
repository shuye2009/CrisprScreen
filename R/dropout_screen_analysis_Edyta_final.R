### for the analysis of dropout screen data of Edyta



## setup ####
if(1){
library(RColorBrewer)
library(ssh)
library(EnhancedVolcano)
library(corrplot)
library(patchwork)
library(UpSetR)
library(ggsci)
library(AnnotationDbi)
  
source("C:/GREENBLATT/Rscripts/CrisprScreen/R/crispr_screen_analysis_lib.R")
color_store <- brewer.pal(n = 8, name = "Dark2")
name_map <- NULL
mapfile <- "C:/GREENBLATT/Rscripts/CrisprScreen/resources/hgnc_geneName_updates_2019.txt"
name_map <- read.table(mapfile, header=T, stringsAsFactors = FALSE)
head(name_map)
dim(name_map)

#color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

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
}


## scale and normalize data ####
if(1){
   #wd <- "C:/data/raw/EDYTA/dropoutScreen/all_screens" 
   wd <- "C:/data/raw/EDYTA/ebvScreen/all_screens" 
   
   setwd(wd)
   format2 <- "_guideRawCount.tab"
   
   guidefiles <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = format2)
   length(guidefiles)
   guidefiles
   
   if(1){
      samples <- gsub("_guideRawCount.tab", "", guidefiles)
      replicates <- unlist(lapply(guidefiles, function(x)split_group(x,3,7))) #4,6 for dropoutScreen
      subjects <- unlist(lapply(guidefiles, function(x)split_group(x,3,4))) #4,4 for dropoutScreen
      design <- data.frame(samples, replicates, subjects)
      rownames(design) <- guidefiles
      design <- design[order(design$subjects),]
      design <- design[!grepl("T10|NR", design$replicates), ]  ## exclude T10 and NR samples
      if(!file.exists("experimental_design_final.tab")){
         write.table(design, "experimental_design_final.tab", sep="\t", row.names=T, col.names=NA, quote=F)
         message("EDIT experimental_design_final.tab NOW!")
      }else{
         eXdesign <- data.frame(read.delim("experimental_design_final.tab", header=TRUE, sep="\t", stringsAsFactors=F))
      }
   }
   
   rownames(eXdesign) <- eXdesign$X
   eXdesign <- eXdesign[,-1]
   
   eXdesign
   
   sink("process_guideRawCount_log.txt")
   guide_data <- process_guideRawCount(eXdesign) ## output can be used for BAGEL directly
   sink()
   gene_data <- process_geneRawCount(eXdesign)
   
   
   combined_gene_table <- NULL  ## each row is a gene
   combined_guide_table <- NULL ## each row is a guide
   combined_lfc_guide_table <- NULL ## each row is a guide
   combined_scaled_guide_table <- NULL ## each row is a guide
   combined_normalized_guide_table <- NULL ## each row is a guide
  
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
      
      head(guideNormalized_table)
      dim(guideNormalized_table)
      summary(guideNormalized_table)
      
      head(guideLFC_table)
      dim(guideLFC_table)
      summary(guideLFC_table)
      
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
   
   
   #GFPhome_dir <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
   GFPhome_dir <- "C:/data/raw/EDYTA/ebvScreen/GFP_final" 
   
   
   if(!dir.exists(GFPhome_dir)){
      dir.create(GFPhome_dir, showWarnings = F)
   }
   
   setwd(GFPhome_dir)
   
   dim(combined_gene_table)
   
   print("writing combined table file")
   write.table(combined_gene_table, file.path(GFPhome_dir,"combined_gene_rawCount_table.tab"), row.names=T, col.names=NA, sep="\t", quote=F)
   
   dim(combined_guide_table)
   NAc <- apply(combined_guide_table, 1, function(x) sum(is.na(x)))
   NAs <- NAc == ncol(combined_guide_table)-2
   combined_guide_table <- combined_guide_table[!NAs,]
   write.table(combined_guide_table, file.path(GFPhome_dir,"combined_guide_rawCount_table.tab"), row.names=F, sep="\t", quote=F)
   combined_guide_filtered <- na.omit(combined_guide_table[, 3:ncol(combined_guide_table)])
   plot_heatmap(combined_guide_filtered, "guide_rawCount")
   
   dim(combined_scaled_guide_table)
   NAc <- apply(combined_scaled_guide_table, 1, function(x) sum(is.na(x)))
   NAs <- NAc == ncol(combined_scaled_guide_table)-2
   combined_scaled_guide_table <- combined_scaled_guide_table[!NAs,]
   ## "combined_guide_scaledCount_table.tab" can be used as input for DrugZ directly
   write.table(combined_scaled_guide_table, file.path(GFPhome_dir,"combined_guide_scaledCount_table.tab"), row.names=F, sep="\t", quote=F)
   combined_scaled_guide_filtered <- na.omit(combined_scaled_guide_table[, 3:ncol(combined_scaled_guide_table)])
   plot_heatmap(combined_scaled_guide_filtered, "guide_scaledCount")
   
   
   dim(combined_normalized_guide_table)
   NAc <- apply(combined_normalized_guide_table, 1, function(x) sum(is.na(x)))
   NAs <- NAc == ncol(combined_normalized_guide_table)-2
   combined_normalized_guide_table <- combined_normalized_guide_table[!NAs,]
   write.table(combined_normalized_guide_table, file.path(GFPhome_dir,"combined_guide_normalizedCount_table.tab"), row.names=F, sep="\t", quote=F)
   combined_normalized_guide_filtered <- na.omit(combined_normalized_guide_table[, 3:ncol(combined_normalized_guide_table)])
   plot_heatmap(combined_normalized_guide_filtered, "guide_normalizedCount")
   
   
   dim(combined_lfc_guide_table)
   NAc <- apply(combined_lfc_guide_table, 1, function(x) sum(is.na(x)))
   NAs <- NAc == ncol(combined_lfc_guide_table)-2
   combined_lfc_guide_table <- combined_lfc_guide_table[!NAs,]
   write.table(combined_lfc_guide_table, file.path(GFPhome_dir,"combined_guide_lfc_table.tab"), row.names=F, sep="\t", quote=F)
   combined_lfc_guide_filtered <- na.omit(combined_lfc_guide_table[, 3:ncol(combined_lfc_guide_table)])
   plot_heatmap(combined_lfc_guide_filtered, "guide_lfcCount")
   
}



## qGI analysis ####
if(1){
   #wd <- "C:/data/raw/EDYTA/dropoutScreen/all_screens" 
   wd <- "C:/data/raw/EDYTA/ebvScreen/all_screens"
   setwd(wd)
   
   eXdesign <- data.frame(read.delim("experimental_design_final.tab", header=TRUE, sep="\t", stringsAsFactors=F))
   
   rownames(eXdesign) <- eXdesign$X
   eXdesign <- eXdesign[,-1]
   
   eXdesign
   
   #hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final" 
   hwd <- "C:/data/raw/EDYTA/ebvScreen/GFP_final"
  
   setwd(hwd)
   combined_lfc_guide_table <- read.delim("combined_guide_lfc_table.tab", header=T, sep="\t")
   
   qGI_dir <- file.path(hwd, "qGI") # "qGI_diff" when loess is set to FALSE in 'compute_lfc_residuals'
   
   if(!dir.exists(qGI_dir)){
      dir.create(qGI_dir, showWarnings = F)
   }
   setwd(qGI_dir)
   
   dim(combined_lfc_guide_table)
   summary(combined_lfc_guide_table) 
   
   lfc_residuals <- compute_lfc_residuals(eXdesign, combined_lfc_guide_table, loess=TRUE, plot=TRUE)
   
   head(lfc_residuals)
   
   fiteb_results <- modt_lfc_residuals(eXdesign, lfc_residuals)
   
   qGI_results <- plot_qGI_results(test_results = fiteb_results, pvalue_cutoff=0.05, gsea=FALSE)
   
   save(qGI_results, file=file.path(qGI_dir, "results_list.Rdata"))
}


#Prepare files for JACKS ####
if(1){
  wd <- "C:/data/raw/EDYTA/dropoutScreen/all_screens" 
  setwd(wd)
  eXdesign <- data.frame(read.delim("experimental_design_final.tab", header=TRUE, sep="\t", stringsAsFactors=F))
  rownames(eXdesign) <- eXdesign$X
  eXdesign <- eXdesign[,-1]
  dim(eXdesign)
  eXdesign
  
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
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
  
  if(0){
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
  ## create tables for JACKS using scaled_count
  input_for_JACKS <- scaled_count[,-2]
  dim(input_for_JACKS)
  summary(input_for_JACKS)
  input_for_JACKS[is.na(input_for_JACKS)] <- 0
  summary(input_for_JACKS)
  head(input_for_JACKS)
  
  
  replicates <- colnames(input_for_JACKS)[2:ncol(input_for_JACKS)]
  samples <- unlist(lapply(replicates, function(x) unlist(strsplit(x, split="_"))[1]))
  controls <- rep("GFP1", length(samples))
  

  #samples <- unlist(lapply(replicates, function(x) gsub("_A|_B|_C", "", x)))
  #samples <- eXdesign[eXdesign$replicates %in% replicates, "subjects"]
  #controls <- eXdesign[eXdesign$replicates %in% replicates, "controls"]

  
  design <- data.frame(replicates, samples, controls)
  design
  
  guide_gene <- scaled_count[,1:2]
  head(guide_gene)
  
  
  control_guide <- guide_gene[guide_gene$GENE %in% control_gene, 1]
  control_count <- input_for_JACKS[input_for_JACKS$SEQUENCE %in% control_guide,]
  dim(control_count)
  
  JACKS_dir <- file.path(hwd,"JACKS")
  
  if(!dir.exists(JACKS_dir)){
    dir.create(JACKS_dir, showWarnings = F)
  }
  setwd(JACKS_dir)
  
  write.table(input_for_JACKS, "countfile.tab", row.names=F, sep="\t", quote=F)
  write.table(design, "replicatemapfile.tab", row.names=F, sep="\t", quote=F)
  write.table(guide_gene, "guidegenemapfile.tab", row.names=F, sep="\t", quote=F)
  writeLines(control_gene, "controlGenefile.txt")
  
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  
  scp_upload(sshSession, c("countfile.tab", "replicatemapfile.tab", "guidegenemapfile.tab", "controlGenefile.txt"), "./JACKS/final_GFP")
   
  ## Run JACKS on dc cluster here
  ssh_disconnect(sshSession)
  }
}

## prepare for DrugZ analysis scaled ####
if(1){
  ## scp "combined_guide_scaledCount_table.tab" to dc01.ccbr.utoronto.ca
  dataset <- "DrugZ"
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  setwd(hwd)
  scp_upload(sshSession, "combined_guide_scaledCount_table.tab", "./DrugZ/final_GFP") 
    
  ssh_disconnect(sshSession)
  
  wd <- paste(hwd, dataset, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  ## run DrugZ on dc1 cluster
}


## JACKS analysis ####
if(1){
  
  dataset <- "JACKS"
 
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  setwd(file.path(hwd, dataset))
  
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  
  scp_download(sshSession, "./JACKS/final_GFP", ".")
  system("mv ./final_GFP/*_JACKS_results.txt .")
  
  ssh_disconnect(sshSession)
  
  JACKS_results <- list()
  JACKS_results_neg <- read.delim("negative_screen_gene_JACKS_results.txt", header=T) 
  JACKS_results_pos <- read.delim("positive_screen_gene_JACKS_results.txt", header=T)
  JACKS_pvalue_neg <- read.delim("negative_screen_gene_pval_JACKS_results.txt", header=T)
  JACKS_pvalue_pos <- read.delim("positive_screen_gene_pval_JACKS_results.txt", header=T)
  
  #plot_heatmap(JACKS_results_neg[,-1], "JACKS_negative")
  #plot_heatmap(JACKS_results_pos[,-1], "JACKS_positive")
  
  JACKS_pvalue_neg[is.na(JACKS_pvalue_neg)] <- 1
  JACKS_padj_neg <- JACKS_pvalue_neg
  for(i in 2:ncol(JACKS_pvalue_neg)){
    JACKS_padj_neg[,i] <- p.adjust(JACKS_pvalue_neg[,i])
  }
  
  JACKS_pvalue_pos[is.na(JACKS_pvalue_pos)] <- 1
  JACKS_padj_pos <- JACKS_pvalue_pos
  for(i in 2:ncol(JACKS_pvalue_pos)){
    JACKS_padj_pos[,i] <- p.adjust(JACKS_pvalue_pos[,i])
  }
  
  #head(JACKS_pvalue_neg[order(JACKS_pvalue_neg$N),])
  #head(JACKS_padj_neg[order(JACKS_padj_neg$N),])
  
  
  screens <- colnames(JACKS_pvalue_neg)[2:ncol(JACKS_pvalue_neg)]
  for(i in 1:length(screens)){
    ascreen <- screens[i]
    res <- NULL
    
    for(df in list(JACKS_results_neg, JACKS_pvalue_neg, JACKS_padj_neg, JACKS_results_pos, JACKS_pvalue_pos, JACKS_padj_pos)){
      rownames(df) <- df[,1]
      res <- cbind(res, df[TKO3_genes, ascreen])
    }
    rownames(res) <- TKO3_genes;
    colnames(res) <- c("neg", "neg_p", "neg_padj", "pos", "pos_p", "pos_padj");
    
    JACKS_results[[ascreen]] <- res
  }
  
  save(JACKS_results, file="results_list.Rdata")
  
  plot_JACKS_results(test_results=JACKS_results, pvalue_cutoff=0.05, gsea=F)
  #def_res_list[[dataset]] <- JACKS_list
  
} ## end if(1)



## DrugZ analysis scaled ####
if(1){
  ## scp "combined_guide_scaledCount_table.tab" to dc01.ccbr.utoronto.ca
  
  library(car)
  library(GSEABase)
  
  dataset <- "DrugZ"
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final" 
  setwd(file.path(hwd, dataset))
  
  if(1){
  sshSession <- ssh_connect(host="shuyepu@dc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  
  scp_download(sshSession, files="./DrugZ/final_GFP", to=".")
  system("mv ./final_GFP/T0_scaled_*_T18_DrugZ-output.txt .")

  ssh_disconnect(sshSession)
  }
  
  DrugZ_analysis=list.files(pattern="^T0_scaled_.+_T18_DrugZ-output.txt")
  
  screens <- unlist(lapply(DrugZ_analysis, function(x) paste(unlist(strsplit(x, split="_"))[3], collapse="_")))
  screens
  DrugZ_results <- list()
  
  for(i in 1:length(screens)){
    ascreen <- screens[i]
    afile <- DrugZ_analysis[i]
    DrugZ_results[[ascreen]] <- read.delim(file.path(hwd, dataset, afile), header=T)
  }
  
  #def_res_list[[dataset]] <- DrugZ_results

  save(DrugZ_results, file="results_list.Rdata")
  
  plot_DrugZ_results(test_results=DrugZ_results, pvalue_cutoff=0.05, gsea=F)
  
  ## following code should be simplified
    
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
    
 
  
  for(pv in c(0.01, 0.05)){
    screen_gsea_results <- list()
    for(i in 1:length(screens)){
      ascreen <- screens[i]
      afile <- paste(ascreen, "GO_BP_gsea_analysis.tab", sep="_")
      if(file.exists(afile)){
        gsea <- read.delim(file.path(hwd, dataset,afile), header=T)
        gsea <- gsea[gsea$p.adjust < pv,]
        if(!is.null(gsea)){
          screen_gsea_results[[ascreen]] <- gsea  
        }
          
      }
      
    }
    names(screen_gsea_results)
    all_sig_terms <- screen_gsea_results[[1]]$Description
    names(all_sig_terms) <- screen_gsea_results[[1]]$ID
    for(i in 2:length(screen_gsea_results)){
      descr <- screen_gsea_results[[i]]$Description
      names(descr) <- screen_gsea_results[[i]]$ID
      all_sig_terms <- c(all_sig_terms, descr)
    }
    
    unique_sig_terms <- unique(names(all_sig_terms))
    all_sig_terms <- all_sig_terms[unique_sig_terms]
    
    
    # try to map to GOslim, did not work properly
    if(0){
    hsGO <- godata('org.Hs.eg.db', ont="BP")
    goIC <- hsGO@IC
    
    slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_metagenomics.obo")
   
    slim_mapping <- list()
    for(i in 1:length(unique_sig_terms)){
      myCollection <- GOCollection(unique_sig_terms[i])
      gs <- goSlim(myCollection, slim, "BP")
      aterm <- all_sig_terms[unique_sig_terms[i]]
      ids <- rownames(gs[gs[,"Count"] > 0,])
      maxIC <- which.max(goIC[ids])
      description <- Term(GOTERM[GOID(GOTERM) == names(maxIC)])
      
      names(aterm) <- NULL
      
      print(aterm)
      print(paste("##     ", description))
      
      slim_mapping[[aterm]] <- description
    }
    sorted_slim_mapping <- slim_mapping[order(unlist(slim_mapping))]
    
    term_df <- as.data.frame(do.call(rbind, sorted_slim_mapping))
    colnames(term_df) <- "GOslim"
    term_df$GOslim <- as.factor(term_df$GOslim)
    term_count <- table(term_df[,"GOslim"])
    myColors = colorRampPalette(c("red", "orange", "blue"), space = "rgb")(length(term_count))
    annotation_colors <- rep(myColors, term_count)
    names(annotation_colors) <- rownames(term_df)
    }
    
    sig_term_matrix <- matrix(0, ncol=length(screens), nrow=length(all_sig_terms), dimnames=list(all_sig_terms, screens))
    for(ascreen in names(screen_gsea_results)){
      descriptions <- screen_gsea_results[[ascreen]]$Description
      enrichment <- screen_gsea_results[[ascreen]]$enrichmentScore
      sig_term_matrix[descriptions, ascreen] <- enrichment
    }
    write.table(sig_term_matrix, paste("GSEA_GOBP_terms_table_", pv, ".tsv", sep=""), col.names = NA, sep="\t", quote=F)
    
    pdf(paste("GSEA_GOBP_terms_heatmap_", pv, ".pdf", sep=""), width=10, height=16)
    pheatmap(sig_term_matrix, cluster_rows=TRUE, fontsize_row=10, border_color="grey60") 
             #annnotation_row=term_df, annotation_colors = annotation_colors)
    dev.off()
  }
  
  ## examine DrugZ results for certain GO terms
  for( aterm in c("anaphase-promoting complex-dependent catabolic process", "proteasome complex")){
  
    go_id = GOID(GOTERM[Term(GOTERM) == aterm]) 
    allegs = get(go_id, org.Hs.egGO2ALLEGS)
    
    agenes = unlist(mget(allegs,org.Hs.egSYMBOL))
    

    #agenes <- unlist(strsplit(screen_gsea_results[["S"]]$core_enrichment[match(aterm, screen_gsea_results[["S"]]$Description)], split="/", fixed=T))
    S_data <- data_list[["S"]]
    S_normZ <- S_data[S_data$GENE %in% agenes, c("GENE", "pval_synth")]
    NSP9_data <- data_list[["NSP9"]]
    NSP9_normZ <- NSP9_data[NSP9_data$GENE %in% agenes, c("GENE", "pval_synth")]
    S_NSP9_df <- merge(S_normZ, NSP9_normZ, by="GENE")
    rownames(S_NSP9_df) <- S_NSP9_df$GENE
    colnames(S_NSP9_df) <- c("GENE", "S", "NSP9")
    fn <- paste(aterm, "DrugZ_pval_synth_heatmap.pdf", sep="_")
    pheatmap(S_NSP9_df[, c(2,3)], cluster_cols=FALSE, display_numbers = T, fontsize_row=7, filename = fn)
    
  }
  
  
  if(0){
    gene_lists_ENTREZ <- lapply(sig_list, function(x)y<- bitr(x, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID)
    gene_lists_ENTREZ <- lapply(gene_lists_ENTREZ, noquote)
    
    plot_clusterProfile_GO(gene_lists_ENTREZ, "DrugZ_modulated_genes")
  } ## end if
  
  
  
  
} ## end if(1)


## Limma analysis on LFC ####

if(1){
  
  
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
 
  dataset <- "Limma"
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
        #arow <- apply(sub, 2, mean, na.rm=TRUE) ## replace missing guide with mean of existing guides
        arow <- rep(NA, ncol(sub)) ## replace missing guide with NA
        sub <- rbind(sub, arow)
        valid_gene <- c(valid_gene, gene)
        subvector <- as.vector(sub)
        lfc_expanded <- rbind(lfc_expanded, subvector)
      }else if(nrow(sub) == 2){
         #arow <- apply(sub, 2, mean, na.rm=TRUE) ## replace missing guide with mean of existing guides
         arow <- rep(NA, ncol(sub)) ## replace missing guide with NA
         sub <- rbind(sub, arow)
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
  
  print(paste("genes with less than 2 guides ", c))
  
  cn <- unlist(lapply(colnames(lfc_combined)[3:ncol(lfc_combined)], function(x) paste(x, guides, sep="_")))
  colnames(lfc_expanded) <- cn
  rownames(lfc_expanded) <- valid_gene
  
  dim(lfc_expanded)
  head(lfc_expanded)
  summary(lfc_expanded[,1:15])
  
  NAs <- apply(lfc_expanded, 1, function(x) sum(is.na(x)))
  NAs[which(NAs == max(NAs))]
  NAss <- which(NAs > 0)
  length(NAss)
  
  lfc_mat <- as.matrix(lfc_expanded)
  dim(lfc_mat)
  
  varVec <- sort(apply(lfc_mat, 1, var, na.rm=TRUE), decreasing=T)
  
  top <- length(varVec)
  lfcMatTopVar <- lfc_mat[names(varVec[1:top]), ]
  
  #lfcMatTopVar <- normalizeCyclicLoess(lfcMatTopVar, span=0.4, iterations = 3, method="fast")
  #lfcMatTopVar <- normalizeQuantiles(lfcMatTopVar)
  #lfcMatTopVar <- normalizeVSN(lfcMatTopVar)
 # meanSdPlot(lfcMatTopVar)
  
  screens <- unlist(lapply(colnames(lfcMatTopVar), function(x) paste(unlist(strsplit(x, split="_"))[1], collapse="_")))
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
  
  #plotMDS(lfcMatTopVar)
  
  eset <- ExpressionSet(assayData = lfcMatTopVar)
  
  
   contrast.matrix <- makeContrasts(
     N_GFP = N - GFP1,
     NSP12_GFP = NSP12 - GFP1,
     S_GFP = S - GFP1,
     NSP9_GFP = NSP9 - GFP1,
     NSP10_GFP = NSP10 - GFP1,
     ORF9B_GFP = ORF9B - GFP1,
     ORF3A_GFP = ORF3A - GFP1,
     ORF8_GFP = ORF8 - GFP1,
     NSP2_GFP = NSP2 - GFP1,
     NSP3_GFP = NSP3 - GFP1,
     M_GFP = M - GFP1,
     levels = design
   )
    
  
  fit <- lmFit(eset, design)
  fit.cont <- contrasts.fit(fit, contrast.matrix)
  fit.eb   <- eBayes(fit.cont)
  
  
  Limma_results <- plot_Limma_results(fit.eb, fit, pvalue_cutoff=0.05, gsea=F)
 
  save(Limma_results, file="results_list.Rdata")
  
  if(0){
     genas(fit.eb, coef=c(1,2), subset="all", plot=TRUE, alpha=0.4)
     modulated_gene_list <- list()
     gene_lists_ENTREZ <- lapply(modulated_gene_list, function(x)y<- bitr(x, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID)
     gene_lists_ENTREZ <- lapply(gene_lists_ENTREZ, noquote)
     
     plot_clusterProfile_GO(gene_lists_ENTREZ, "Limma_modulated_genes")
     
     #GOA <- goana(gene_lists_ENTREZ)
     #topup <- topGO(GOA, sort = "N_GFP_up")
     #topdown <- topGO(GOA, sort = "N_GFP_down")
     
     # Mean-difference plot
     plotMD(fit.eb,column=3)
     
     # Q-Q plot of moderated t-statistics
     qqt(fit.eb$t[,1],df=fit.eb$df.residual+fit.eb$df.prior)
     abline(0,1)
  
  }
  
}

#DESeq2 on scaled data ####

if(1){
  
  
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  
  
  dataset <- "DESeq2"
  wd <- paste(hwd, dataset, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  scaledcombined <- read.delim(file.path(hwd, "combined_guide_scaledCount_table.tab"), header=T, sep="\t")
  
  NAs <- apply(scaledcombined, 1, function(x) sum(is.na(x)))
  dim(scaledcombined[which(NAs == (ncol(scaledcombined)-2)),])
  
  dim(scaledcombined)
  head(scaledcombined)
  summary(scaledcombined)
  
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
          #arow <- apply(sub, 2, mean, na.rm=TRUE) ## replace missing guide with mean of existing guides
          arow <- rep(NA, ncol(sub)) ## replace missing guide with NA
          sub <- rbind(sub, arow)
          valid_gene <- c(valid_gene, gene)
          subvector <- as.vector(sub)
          scaledexpanded <- rbind(scaledexpanded, subvector)
       }else if(nrow(sub) == 2){
          #arow <- apply(sub, 2, mean, na.rm=TRUE) ## replace missing guide with mean of existing guides
          arow <- rep(NA, ncol(sub)) ## replace missing guide with NA
          sub <- rbind(sub, arow)
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
  
  print(paste("genes with less than 2 guides ", c))
  
  cn <- unlist(lapply(colnames(scaledcombined)[3:ncol(scaledcombined)], function(x) paste(x, guides, sep="_")))
  colnames(scaledexpanded) <- cn
  rownames(scaledexpanded) <- valid_gene
  
  dim(scaledexpanded)
  head(scaledexpanded)
  summary(scaledexpanded[, 1:10])
  
  
  data_table <- scaledexpanded
  data_table[is.na(data_table)] <- 0
  data_mat <- as.matrix(data_table)
  dim(data_mat)
  
  varVec <- sort(apply(data_mat, 1, var), decreasing=T)

  top <- length(varVec)
  data_table <- data_mat[names(varVec[1:top]), ]
  
  samples <- colnames(data_table)
  ncol <- length(samples)
  conditions <- unlist(lapply(colnames(data_table), function(x) paste(unlist(strsplit(x, split="_"))[1], collapse="_")))
  
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
  
  
  contrasts <- cbind(c(rep("conditions", 11)), 
                         c("N", "NSP12", "S", "NSP9", "NSP10", "ORF9B", "ORF3A", "ORF8", "NSP2", "NSP3", "M"), 
                         c(rep(c("GFP1"), times=c(11))))
  
  DESeq2_results <- plot_DESeq2_results(test_results=dds, contrasts, pvalue_cutoff=0.05, gsea=F, verbose=FALSE)
  
  save(DESeq2_results, file="results_list.Rdata")
  
}


## load results as Rdata #####
if(1){
  
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  datasets <- c("qGI", "DrugZ", "JACKS", "DESeq2", "Limma")
  
  for(x in datasets){
     #x <- "DrugZ"
     load(file.path(hwd, x, "results_list.Rdata"), verbose=F)
  }
  
  Limma <- Limma_results # list of dataframe
  DrugZ <- DrugZ_results # list of dataframe
  DESeq2 <- DESeq2_results # list of dataframe
  JACKS <- JACKS_results # list of dataframe
  qGI <- qGI_results 
  
  names(JACKS)
  names(Limma)
  names(DESeq2)
  names(DrugZ)
  names(qGI)
  
  names(JACKS) <- paste(names(JACKS), "_", sep="")
  names(Limma) <- gsub("_GFP.*", "_", names(Limma))
  names(DESeq2) <- unlist(lapply(names(DESeq2), function(x) unlist(strsplit(x, split="vs"))[1]))
  names(DrugZ) <- paste(names(DrugZ), "_", sep="")
  names(qGI) <- paste(names(qGI), "_", sep="")
  
  lapply(list(Limma, DrugZ, DESeq2, JACKS, qGI), function(x) lapply(x, dim))
  
  DrugZ <- lapply(DrugZ, function(x){
     rownames(x) <- x$GENE
     x <- select(x, -GENE)
  })
  
  tests <- c("DrugZ", "Limma", "qGI", "DESeq2", "JACKS")
  all_screens_ <- names(JACKS)
  all_screens <- gsub("_$", "", all_screens_)
  
  Limma_t_table <- list()
  DESeq2_stat_table <- list()
  JACKS_neg_table <- list()
  DrugZ_normz_table <- list()
  qGI_t_table <- list()
  
  pval_df_list <- list()
  stat_table_list <- list()
  
  for(screen in all_screens_){
     #screen <- all_screens_[1]  
     qGI_df <- qGI[[names(qGI)[grep(screen, names(qGI))]]]
     Limma_df <- Limma[[names(Limma)[grep(screen, names(Limma))]]]
     JACKS_df <- as.data.frame(JACKS[[names(JACKS)[grep(screen, names(JACKS))]]])
     DrugZ_df <- DrugZ[[names(DrugZ)[grep(screen, names(DrugZ))]]]
     DESeq2_df <- DESeq2[[names(DESeq2)[grep(screen, names(DESeq2))]]]
     
     print(screen)
     head(Limma_df)
     head(DrugZ_df)
     head(DESeq2_df)
     head(JACKS_df)
     head(qGI_df)
     
     Limma_df[is.na(Limma_df)] <- 0
     qGI_df[is.na(qGI_df)] <- 0
     
     Limma_df <- Limma_df %>%
        mutate(P.value_synth=P.Value, P.value_supp=P.Value)
     Limma_df$P.value_synth[Limma_df$t>0] <- 1 - Limma_df$P.Value[Limma_df$t>0]
     Limma_df$P.value_supp[Limma_df$t<=0] <- 1 - Limma_df$P.Value[Limma_df$t<=0]
     
     qGI_df <- qGI_df %>%
        mutate(P.value_synth=P.Value, P.value_supp=P.Value)
     qGI_df$P.value_synth[qGI_df$t>0] <- 1 - qGI_df$P.Value[qGI_df$t>0]
     qGI_df$P.value_supp[qGI_df$t<=0] <- 1 - qGI_df$P.Value[qGI_df$t<=0]
     
     DESeq2_df <- as.data.frame(DESeq2_df) %>%
        mutate(P.value_synth=pvalue, P.value_supp=pvalue)
     DESeq2_df$P.value_synth[DESeq2_df$stat>0] <- 1 - DESeq2_df$pvalue[DESeq2_df$stat>0]
     DESeq2_df$P.value_supp[DESeq2_df$stat<=0] <- 1 - DESeq2_df$pvalue[DESeq2_df$stat<=0]
     
     synth_pval_combined <- cbind(DrugZ_df[TKO3_genes, "pval_synth"],
                                  Limma_df[TKO3_genes, "P.value_synth"],
                                  qGI_df[TKO3_genes, "P.value_synth"],
                                  DESeq2_df[TKO3_genes, "P.value_synth"],
                                  JACKS_df[TKO3_genes, "neg_p"])
     supp_pval_combined <- cbind(DrugZ_df[TKO3_genes, "pval_supp"],
                                 Limma_df[TKO3_genes, "P.value_supp"],
                                 qGI_df[TKO3_genes, "P.value_supp"],
                                 DESeq2_df[TKO3_genes, "P.value_supp"],
                                 JACKS_df[TKO3_genes, "pos_p"])
     
     colnames(synth_pval_combined) <- colnames(supp_pval_combined) <- tests
     rownames(synth_pval_combined) <- rownames(supp_pval_combined) <- TKO3_genes
     
     synth_pval_combined_filtered <- na.omit(synth_pval_combined)
     supp_pval_combined_filtered <- na.omit(supp_pval_combined)
     
     #pheatmap(synth_pval_combined_filtered)
     #pheatmap(supp_pval_combined_filtered)
     
     pval_df_list[["supp"]][[screen]] <- supp_pval_combined_filtered
     pval_df_list[["synth"]][[screen]] <- synth_pval_combined_filtered
     
     
     ## collect stats for each test
     Limma_t_table[[screen]] <- Limma_df[TKO3_genes, "t"]
     qGI_t_table[[screen]] <- qGI_df[TKO3_genes, "t"]
     DESeq2_stat_table[[screen]] <- DESeq2_df[TKO3_genes, "stat"]
     JACKS_neg_table[[screen]] <- JACKS_df[TKO3_genes, "neg"]
     DrugZ_normz_table[[screen]] <- DrugZ_df[TKO3_genes, "normZ"]
  }
 
  stat_table_list[["DESeq2"]] <- as.data.frame(DESeq2_stat_table, row.names=TKO3_genes)
  stat_table_list[["Limma"]] <- as.data.frame(Limma_t_table, row.names=TKO3_genes)
  stat_table_list[["JACKS"]] <- as.data.frame(JACKS_neg_table, row.names=TKO3_genes)
  stat_table_list[["DrugZ"]] <- as.data.frame(DrugZ_normz_table, row.names=TKO3_genes)
  stat_table_list[["qGI"]] <- as.data.frame(qGI_t_table, row.names=TKO3_genes)
  
  lapply(stat_table_list, dim)
}
  
 
## integrative GO analysis ####
if(1){
hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"

wd <- paste(hwd, "integrated_GO_results", sep="/")
if(!dir.exists(wd)){
   dir.create(wd, showWarnings = F)
}
print(wd)
setwd(wd)


tests <- c("qGI", "DrugZ", "DESeq2", "Limma", "JACKS")

all_screens <- c("S", "N", "M", "NSP2", "NSP3", "NSP9", "NSP10", "NSP12", "ORF3A", "ORF8", "ORF9B")
all_screens_ <- paste0(all_screens, "_")

### GO term gene heatmap ####

#terms <- c("GOCC_PROTEASOME_COMPLEX", "GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS")
terms <- c("proteasome complex")
data_list <- pval_df_list[["synth"]]


pdf("qGI_synthPval_GOTERM_heatmap.pdf", width=8, height=8)
plot_GO_term_heatmap_pval(terms, genes=all_screens_, arrangeBy="S", tests="qGI", data_list)
dev.off()

pdf("qGI_moderatedT_GOTERM_heatmap.pdf", height=8, width=8)
plot_GO_term_heatmap(terms, "S", stat_table_list[["qGI"]], measurement="moderated t")
dev.off()


## GSEA analysis all data ####
if(0){
   ## this part takes very long time, probably a couple of days
   pval <- 0.05 
   gseaRe <- NULL
   for(test in tests){
      for(screen in all_screens_){
         gene_list <- stat_table_list[[test]][[screen]]
         names(gene_list) <- rownames(stat_table_list[[test]])
         dataName <- paste0(screen, test)
         if(test == "qGI" && screen == "S_"){
            for(ont in c("MF", "CC")){
               run_gseGO_simpleList(gene_list, dataName=dataName, adjp_cutoff=pval, simplify=F, ont=ont, GO_file=NULL)
            }
            run_gseGO_simpleList(gene_list, dataName=dataName, adjp_cutoff=pval, ont="BP", GO_file=NULL)
            gseaRe <- run_gseGO_simpleList(gene_list, dataName=dataName, adjp_cutoff=pval, ont="REACTOME_Pathway", GO_file=PATHWAY_file)
         }else{
            run_gseGO_simpleList(gene_list, dataName=dataName, adjp_cutoff=pval, ont="BP", GO_file=NULL)
            run_gseGO_simpleList(gene_list, dataName=dataName, adjp_cutoff=pval, ont="REACTOME_Pathway", GO_file=PATHWAY_file)
         }
      }
   }
   
   ## inspect the 'S_Reactome_pathway_gsea_analysis.tab', get indices of viral and immune related terms
   viral_subset <- c(21, 27, 29, 33, 35, 47, 89, 94, 112, 132, 1147, 159, 161, 180, 195, 196, 209, 222, 230, 
                     235, 236, 239, 249, 281,282,338,253,264,270,272,284,290,291,295,297,299,327,348,369,374,417) - 1
   
   pdf("S_qGI_REACTOME_Pathway_gsea_analysis_viral.pdf", height=8, width=10)
   clusterProfiler_GSEA_plot(gseaRe, geneList=gene_list, labFormat=75, pval=0.05, subsetDesc=viral_subset) 
   dev.off()
}

## Number of terms heatmap for all screens ####
padj_filter <- 0.001
gsea_GOBP <- list()
gsea_Reactome <- list()
GO_term_count <- matrix(rep(0, length(all_screens)*length(tests)), nrow=length(tests), dimnames=list(tests, all_screens))
Reactome_term_count <- matrix(rep(0, length(all_screens)*length(tests)), nrow=length(tests), dimnames=list(tests, all_screens))


for(test in tests){
   #setwd(file.path(hwd, test))
   BPfiles <- list.files(path=wd, pattern=paste0(test, "_BP_gsea_analysis.tab"))
   print(BPfiles)
   Reactomefiles <- list.files(path=wd, pattern=paste0(test,"_REACTOME_Pathway_gsea_analysis.tab"))
   print(Reactomefiles)
   for(s in all_screens){
      for(afile in BPfiles){
         if(grepl(paste0(s,"_"), afile)){
            gsea_GOBP[[test]][[s]] <- read.delim(file.path(wd, afile), header=T, sep="\t") %>%
               filter(p.adjust < padj_filter)
            print(dim(gsea_GOBP[[test]][[s]]))
            GO_term_count[test, s] <- nrow(gsea_GOBP[[test]][[s]])
         }
      }
      for(afile in Reactomefiles){
         if(grepl(paste0(s,"_"), afile)){
            gsea_Reactome[[test]][[s]] <- read.delim(file.path(wd, afile), header=T, sep="\t") %>%
               filter(p.adjust < padj_filter)
            Reactome_term_count[test, s] <- nrow(gsea_Reactome[[test]][[s]])
         }
      }
   }
}

pdf("Number_of_GSEA_terms_heatmap.pdf", height=4, width=8)
pheatmap(GO_term_count, display_numbers = T, number_format = "%.0f", fontsize_number = 12, main="Number of GOBP terms")
pheatmap(Reactome_term_count, display_numbers = T, number_format = "%.0f", fontsize_number = 12, main="Number of Reactome pathway terms")
dev.off()

## shared and unique terms in qGI ####
names(gsea_GOBP)
names(gsea_GOBP[["qGI"]])
sapply(gsea_GOBP[["qGI"]], dim)
colnames(gsea_GOBP[["qGI"]][["S"]])

qGI_gsea <- list()
qGI_gsea[["GOBP"]] <- gsea_GOBP[["qGI"]]
qGI_gsea[["Reactome"]] <- gsea_Reactome[["qGI"]]


for(geneset in c("GOBP", "Reactome")){
   enrichment <- qGI_gsea[[geneset]]
   top_supp_terms <- NULL
   top_synth_terms <- NULL
   top_supp_by_qGI <- list()
   top_synth_by_qGI <- list()
   
   all_supp_terms <- NULL
   all_synth_terms <- NULL
   all_supp_by_qGI <- list()
   all_synth_by_qGI <- list()
   top <- 20
   for(screen in names(enrichment)){
      df <- enrichment[[screen]]
      df_pos <- df %>%
         filter(NES > 0) %>%
         arrange(p.adjust)
      top_supp_by_qGI[[screen]] <- df_pos[1:min(top, nrow(df_pos)), c("Description", "p.adjust")]
      top_supp_terms <- c(top_supp_terms, df_pos[1:min(top, nrow(df_pos)), "Description"])
      all_supp_by_qGI[[screen]] <- df_pos[, c("Description", "p.adjust")]
      all_supp_terms <- c(all_supp_terms, df_pos[, "Description"])
      df_neg <- df %>%
         filter(NES < 0) %>%
         arrange(p.adjust)
      top_synth_by_qGI[[screen]] <- df_neg[1:min(top, nrow(df_neg)), c("Description", "p.adjust")]
      top_synth_terms <- c(top_synth_terms, df_neg[1:min(top, nrow(df_neg)), "Description"])
      all_synth_by_qGI[[screen]] <- df_neg[, c("Description", "p.adjust")]
      all_synth_terms <- c(all_synth_terms, df_neg[, "Description"])
   }
   
   top_supp_terms <- top_supp_terms[!is.na(top_supp_terms)]
   top_synth_terms <- top_synth_terms[!is.na(top_synth_terms)]
   length(unique(top_supp_terms))
   length(unique(top_synth_terms))
   
   supp_term_pmat <- matrix(rep(0, length(all_screens)*length(unique(top_supp_terms))), nrow=length(unique(top_supp_terms)), byrow=T, dimnames=list(unique(top_supp_terms), all_screens))
   synth_term_pmat <- matrix(rep(0, length(all_screens)*length(unique(top_synth_terms))), nrow=length(unique(top_synth_terms)), byrow=T, dimnames=list(unique(top_synth_terms), all_screens))
   
   length(unique(all_supp_terms))
   length(unique(all_synth_terms))
   
   all_supp_term_pmat <- matrix(rep(0, length(all_screens)*length(unique(all_supp_terms))), nrow=length(unique(all_supp_terms)), byrow=T, dimnames=list(unique(all_supp_terms), all_screens))
   all_synth_term_pmat <- matrix(rep(0, length(all_screens)*length(unique(all_synth_terms))), nrow=length(unique(all_synth_terms)), byrow=T, dimnames=list(unique(all_synth_terms), all_screens))
   
   for(screen in all_screens){
      dfp <- top_supp_by_qGI[[screen]]
      for(term in unique(top_supp_terms)){
         if(term %in% dfp$Description) supp_term_pmat[term, screen] <- -log10(dfp[dfp$Description==term, "p.adjust"])
      }
      dfn <- top_synth_by_qGI[[screen]]
      for(term in unique(top_synth_terms)){
         if(term %in% dfn$Description) synth_term_pmat[term, screen] <- -log10(dfn[dfn$Description==term, "p.adjust"])
      }
      
      all_dfp <- all_supp_by_qGI[[screen]]
      for(term in unique(all_supp_terms)){
         if(term %in% all_dfp$Description) all_supp_term_pmat[term, screen] <- -log10(all_dfp[all_dfp$Description==term, "p.adjust"])
      }
      all_dfn <- all_synth_by_qGI[[screen]]
      for(term in unique(all_synth_terms)){
         if(term %in% all_dfn$Description) all_synth_term_pmat[term, screen] <- -log10(all_dfn[all_dfn$Description==term, "p.adjust"])
      }
   }
   
   pdf(paste0("qGI_", geneset, "_terms_heatmap.pdf"), height=10, width=10)
   
   if(geneset == "GOBP"){
      suppTermsOrdered <- reorderGOTerms(unique(top_supp_terms), k=5)
      synthTermsOrdered <- reorderGOTerms(unique(top_synth_terms), k=8)
      
      h1 <- Heatmap(supp_term_pmat[names(suppTermsOrdered),], 
              name = "-log(adj.p)", 
              width = unit(4, "inch"),
              row_order = seq_along(suppTermsOrdered),
              column_title = paste("Top", top, geneset, "terms enrich in suppression genes"),
              show_heatmap_legend = T,
              row_names_max_width = unit(5, "inch"),
              row_names_gp = gpar(col = color_store[suppTermsOrdered], fontsize = 10))
      
      h2 <- Heatmap(synth_term_pmat[names(synthTermsOrdered),], 
              name = "-log(adj.p)", 
              width = unit(4, "inch"),
              row_order = seq_along(synthTermsOrdered),
              column_title = paste("Top", top, geneset, "terms enrich in synthetic lethal genes"),
              show_heatmap_legend = T,
              row_names_max_width = unit(5, "inch"),
              row_names_gp = gpar(col = color_store[synthTermsOrdered], fontsize = 10))
      
      suppTermsOrdered <- reorderGOTerms(unique(all_supp_terms), k=5)
      synthTermsOrdered <- reorderGOTerms(unique(all_synth_terms), k=8)
      
      h3 <- Heatmap(all_supp_term_pmat[names(suppTermsOrdered),], 
              name = "-log(adj.p)", 
              width = unit(4, "inch"),
              row_order = seq_along(suppTermsOrdered),
              column_title = paste("All", geneset, "terms enrich in suppression genes"),
              show_heatmap_legend = T,
              row_names_max_width = unit(5, "inch"),
              row_names_gp = gpar(col = color_store[suppTermsOrdered], fontsize = 10))
      
      
      h4 <- Heatmap(all_synth_term_pmat[names(synthTermsOrdered),], 
              name = "-log(adj.p)", 
              width = unit(4, "inch"),
              row_order = seq_along(synthTermsOrdered),
              column_title = paste("All", geneset, "terms enrich in synthetic lethal genes"),
              show_heatmap_legend = T,
              row_names_max_width = unit(5, "inch"),
              row_names_gp = gpar(col = color_store[synthTermsOrdered], fontsize = 10))
    }else{
       
       h1 <- Heatmap(supp_term_pmat, 
                     name = "-log(adj.p)", 
                     width = unit(4, "inch"),
                     column_title = paste("Top", top, geneset, "terms enrich in suppression genes"),
                     show_heatmap_legend = T,
                     row_names_max_width = unit(5, "inch"),
                     row_names_gp = gpar(fontsize = 10))
       
       h2 <- Heatmap(synth_term_pmat, 
                     name = "-log(adj.p)", 
                     width = unit(4, "inch"),
                     column_title = paste("Top", top, geneset, "terms enrich in synthetic lethal genes"),
                     show_heatmap_legend = T,
                     row_names_max_width = unit(5, "inch"),
                     row_names_gp = gpar(fontsize = 10))
       
       h3 <- Heatmap(all_supp_term_pmat, 
                     name = "-log(adj.p)", 
                     width = unit(4, "inch"),
                     column_title = paste("All", geneset, "terms enrich in suppression genes"),
                     show_heatmap_legend = T,
                     row_names_max_width = unit(5, "inch"),
                     row_names_gp = gpar(fontsize = 10))
       
       h4 <- Heatmap(all_synth_term_pmat, 
                     name = "-log(adj.p)", 
                     width = unit(4, "inch"),
                     column_title = paste("All", geneset, "terms enrich in synthetic lethal genes"),
                     show_heatmap_legend = T,
                     row_names_max_width = unit(5, "inch"),
                     row_names_gp = gpar(fontsize = 10))
                     
    }
   
   draw(h1, heatmap_legend_side="left")
   draw(h2, heatmap_legend_side="left")
   draw(h3, heatmap_legend_side="left")
   draw(h4, heatmap_legend_side="left")
      
   dev.off()
   
   screen_total <- apply(all_synth_term_pmat, 1, function(x){length(x[x>0])})
   all_synth_term_pmat <- cbind(all_synth_term_pmat, screen_total)
   all_synth_term_pmat <- all_synth_term_pmat[order(all_synth_term_pmat[,"screen_total"], decreasing=T),]
   terms_total <- apply(all_synth_term_pmat, 2, function(x){length(x[x>0])})
   all_synth_term_pmat <- rbind(all_synth_term_pmat, terms_total)
   
   write.table(all_synth_term_pmat, paste0("qGI_all_", geneset, "_terms_synth.tab"), sep="\t", col.names=NA, quote=F)
   
   screen_total <- apply(all_supp_term_pmat, 1, function(x){length(x[x>0])})
   all_supp_term_pmat <- cbind(all_supp_term_pmat, screen_total)
   all_supp_term_pmat <- all_supp_term_pmat[order(all_supp_term_pmat[,"screen_total"], decreasing=T),]
   terms_total <- apply(all_supp_term_pmat, 2, function(x){length(x[x>0])})
   all_supp_term_pmat <- rbind(all_supp_term_pmat, terms_total)
   
   write.table(all_supp_term_pmat, paste0("qGI_all_", geneset, "_terms_supp.tab"), sep="\t", col.names=NA, quote=F)
}
}

if(0){ ## obsolete
  library(org.Hs.eg.db)
  library(GO.db)
  library(RColorBrewer)
  
  functionAnalysis <- TRUE
  GOpvalue_cutoff <- 0.05
  pvalue_cutoffs <- c(0.1, 0.05, 0.02, 0.01)

  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  
  
  wd <- paste(hwd, "integrated_results", sep="/")
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  
  load(file.path(hwd, "def_res_list.Rdata"))
  
  ana_names <- names(def_res_list)
  print(ana_names)
  
  Limma <- def_res_list[["Limma"]] # list of dataframe
  DrugZ <- def_res_list[["DrugZ"]] # list of dataframe
  DESeq2_on_reads <- def_res_list[["DESeq2"]] # list of dataframe
  JACKS <- def_res_list[["JACKS"]] # list of dataframe
  
  names(JACKS)
  names(Limma)
  names(DESeq2_on_reads)
  names(DrugZ)
  
  names(JACKS) <- paste(names(JACKS), "_", sep="")
  names(Limma) <- gsub("_GFP.*", "_", names(Limma))
  names(DESeq2_on_reads) <- unlist(lapply(names(DESeq2_on_reads), function(x)strsplit(x,split="vs")[[1]][1]))
  names(DrugZ) <- paste(names(DrugZ), "_", sep="")
  
  tests <- c("DrugZ", "Limma", "DESeq2", "JACKS")
  all_screens_ <- names(JACKS)
  all_screens <- gsub("_$", "", all_screens_)
  
  synth_sig_GO_term <- NULL
  supp_sig_GO_term <- NULL
  
  synth_sig_count_pval <- list()
  supp_sig_count_pval <- list()
  
  synth_sig_count_FDR <- list()
  supp_sig_count_FDR <- list()
  
  
  for(screen in all_screens_){
    
    Limma_df <- Limma[[names(Limma)[grep(screen, names(Limma))]]]
    JACKS_df <- as.data.frame(JACKS[[names(JACKS)[grep(screen, names(JACKS))]]])
    DrugZ_df <- DrugZ[[names(DrugZ)[grep(screen, names(DrugZ))]]]
    rownames(DrugZ_df) <- DrugZ_df$GENE
    
    DESeq2_on_reads_df <- DESeq2_on_reads[[names(DESeq2_on_reads)[grep(screen, names(DESeq2_on_reads))]]]
    
    print(screen)
    print(dim(Limma_df))
    print(dim(DrugZ_df))
    print(dim(DESeq2_on_reads_df))
    print(dim(JACKS_df))
    head(Limma_df)
    head(DrugZ_df)
    head(DESeq2_on_reads_df)
    head(JACKS_df)
    

    ## pvalue integration for GO enrichment analysis ##
    
    gene_at_pval_synth <- list()
    gene_at_pval_supp <- list()
    
    gene_at_FDR_synth <- list()
    gene_at_FDR_supp <- list()
    
    for(pvalue_cutoff in pvalue_cutoffs){
      print(pvalue_cutoff)
      DrugZPval_synth <- rownames(DrugZ_df[DrugZ_df$pval_synth < pvalue_cutoff,])
      LimmaPval_synth <- rownames(Limma_df[Limma_df$P.Value < pvalue_cutoff & Limma_df$t < 0,])
      DESeq2Pval_synth <- rownames(DESeq2_on_reads_df[DESeq2_on_reads_df$pvalue < pvalue_cutoff & DESeq2_on_reads_df$stat < 0,])
      JACKSPval_synth <- rownames(JACKS_df[JACKS_df$neg_p < pvalue_cutoff,])
      
      DrugZFDR_synth <- rownames(DrugZ_df[DrugZ_df$fdr_synth < pvalue_cutoff,])
      LimmaFDR_synth <- rownames(Limma_df[Limma_df$adj.P.Val < pvalue_cutoff & Limma_df$t < 0,])
      DESeq2FDR_synth <- rownames(DESeq2_on_reads_df[DESeq2_on_reads_df$padj < pvalue_cutoff & DESeq2_on_reads_df$stat < 0,])
      JACKSFDR_synth <- rownames(JACKS_df[JACKS_df$neg_padj < pvalue_cutoff,])
      
      
      alist <- list("DrugZ"=DrugZPval_synth, "Limma"=LimmaPval_synth, "DESeq2"=DESeq2Pval_synth, "JACKS"=JACKSPval_synth)
      atitle <- paste("synthetic lethal genes at pvalue < ", pvalue_cutoff, sep="")
      plotVenn(alist, paste(screen, pvalue_cutoff, "synth_pval_venn_plot.png",sep=""), atitle)
      
      gene_at_pval_synth[["DrugZ"]][[as.character(pvalue_cutoff)]] <- DrugZPval_synth
      gene_at_pval_synth[["Limma"]][[as.character(pvalue_cutoff)]] <- LimmaPval_synth
      gene_at_pval_synth[["DESeq2"]][[as.character(pvalue_cutoff)]] <- DESeq2Pval_synth
      gene_at_pval_synth[["JACKS"]][[as.character(pvalue_cutoff)]] <- JACKSPval_synth 
      
      gene_at_FDR_synth[["DrugZ"]][[as.character(pvalue_cutoff)]] <- DrugZFDR_synth
      gene_at_FDR_synth[["Limma"]][[as.character(pvalue_cutoff)]] <- LimmaFDR_synth
      gene_at_FDR_synth[["DESeq2"]][[as.character(pvalue_cutoff)]] <- DESeq2FDR_synth
      gene_at_FDR_synth[["JACKS"]][[as.character(pvalue_cutoff)]] <- JACKSFDR_synth
      
      
      if(functionAnalysis){
        run_enrichGO_simpleList(DrugZPval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_DrugZ_synth", sep=""))
        run_enrichGO_simpleList(LimmaPval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_Limma_synth", sep=""))
        run_enrichGO_simpleList(DESeq2Pval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_DESeq2_synth", sep=""))
        run_enrichGO_simpleList(JACKSPval_synth, ont="BP", paste(screen, pvalue_cutoff, "_pval_JACKS_synth", sep=""))
      }
      
     
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_JACKS_synth_genes.tab",sep=""))){
        GO_JACKS <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_JACKS_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_JACKS[GO_JACKS$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DrugZ_synth_genes.tab",sep=""))){
        GO_DrugZ <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DrugZ_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_DrugZ[GO_DrugZ$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DESeq2_synth_genes.tab",sep=""))){
        GO_DESeq2 <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DESeq2_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_DESeq2[GO_DESeq2$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_Limma_synth_genes.tab",sep=""))){
        GO_Limma <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_Limma_synth_genes.tab",sep=""), header=T, sep="\t")
        synth_sig_GO_term <- union(synth_sig_GO_term, GO_Limma[GO_Limma$p.adjust < GOpvalue_cutoff, "ID"])
      }
      
      
      
      
      ## supp genes ##
      DrugZPval_supp <- rownames(DrugZ_df[DrugZ_df$pval_supp < pvalue_cutoff,])
      LimmaPval_supp <- rownames(Limma_df[Limma_df$P.Value < pvalue_cutoff & Limma_df$t > 0,])
      DESeq2Pval_supp <- rownames(DESeq2_on_reads_df[DESeq2_on_reads_df$pvalue < pvalue_cutoff & DESeq2_on_reads_df$stat > 0,])
      JACKSPval_supp <- rownames(JACKS_df[JACKS_df$pos_p < pvalue_cutoff,])
      
      DrugZFDR_supp <- rownames(DrugZ_df[DrugZ_df$fdr_supp < pvalue_cutoff,])
      LimmaFDR_supp <- rownames(Limma_df[Limma_df$adj.P.Value < pvalue_cutoff & Limma_df$t > 0,])
      DESeq2FDR_supp <- rownames(DESeq2_on_reads_df[DESeq2_on_reads_df$padh < pvalue_cutoff & DESeq2_on_reads_df$stat > 0,])
      JACKSFDR_supp <- rownames(JACKS_df[JACKS_df$pos_padj < pvalue_cutoff,])
      
      alist <- list("DrugZ"=DrugZPval_supp, "Limma"=LimmaPval_supp, "DESeq2"=DESeq2Pval_supp, "JACKS"=JACKSPval_supp)
      atitle <- paste("suppression genes at pvalue < ",pvalue_cutoff,sep="")
      plotVenn(alist, paste(screen, pvalue_cutoff, "supp_pval_venn_plot.png",sep=""), atitle)
      
      gene_at_pval_supp[["DrugZ"]][[as.character(pvalue_cutoff)]] <- DrugZPval_supp
      gene_at_pval_supp[["Limma"]][[as.character(pvalue_cutoff)]] <- LimmaPval_supp
      gene_at_pval_supp[["DESeq2"]][[as.character(pvalue_cutoff)]] <- DESeq2Pval_supp
      gene_at_pval_supp[["JACKS"]][[as.character(pvalue_cutoff)]] <- JACKSPval_supp
      
      gene_at_FDR_supp[["DrugZ"]][[as.character(pvalue_cutoff)]] <- DrugZFDR_supp
      gene_at_FDR_supp[["Limma"]][[as.character(pvalue_cutoff)]] <- LimmaFDR_supp
      gene_at_FDR_supp[["DESeq2"]][[as.character(pvalue_cutoff)]] <- DESeq2FDR_supp
      gene_at_FDR_supp[["JACKS"]][[as.character(pvalue_cutoff)]] <- JACKSFDR_supp
      
      
      if(functionAnalysis){
        run_enrichGO_simpleList(DrugZPval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_DrugZ_supp", sep=""))
        run_enrichGO_simpleList(LimmaPval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_Limma_supp", sep=""))
        run_enrichGO_simpleList(DESeq2Pval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_DESeq2_supp", sep=""))
        run_enrichGO_simpleList(JACKSPval_supp, ont="BP", paste(screen, pvalue_cutoff, "_pval_JACKS_supp", sep=""))
      }
      
      
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_JACKS_supp_genes.tab",sep=""))){
        GO_JACKS <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_JACKS_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_JACKS[GO_JACKS$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DrugZ_supp_genes.tab",sep=""))){
        GO_DrugZ <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DrugZ_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_DrugZ[GO_DrugZ$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DESeq2_supp_genes.tab",sep=""))){
        GO_DESeq2 <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DESeq2_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_DESeq2[GO_DESeq2$p.adjust < GOpvalue_cutoff, "ID"])
      }
      if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_Limma_supp_genes.tab",sep=""))){
        GO_Limma <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_Limma_supp_genes.tab",sep=""), header=T, sep="\t")
        supp_sig_GO_term <- union(supp_sig_GO_term, GO_Limma[GO_Limma$p.adjust < GOpvalue_cutoff, "ID"])
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
 
  hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  
  
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
          
          GO_DrugZ_termp <- GO_Limma_termp <- GO_DESeq2_termp <- GO_JACKS_termp <- 0
          
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_JACKS_",ex,"_genes.tab",sep=""))){
            GO_JACKS <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_JACKS_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_JACKS$ID){
              GO_JACKS_termp <- -log(GO_JACKS[GO_JACKS$ID == term, "p.adjust"])
            }
          }
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DrugZ_",ex,"_genes.tab",sep=""))){
            GO_DrugZ <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DrugZ_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_DrugZ$ID){
              GO_DrugZ_termp <- -log(GO_DrugZ[GO_DrugZ$ID == term, "p.adjust"])
            }
          }
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_Limma_",ex,"_genes.tab",sep=""))){
            GO_Limma <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_Limma_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_Limma$ID){
              GO_Limma_termp <- -log(GO_Limma[GO_Limma$ID == term, "p.adjust"])
            }
          }
          if(file.exists(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DESeq2_",ex,"_genes.tab",sep=""))){
            GO_DESeq2 <- read.table(paste("GO_enrichment_",screen, pvalue_cutoff, "_pval_DESeq2_",ex,"_genes.tab",sep=""), header=T, sep="\t")
            if(term %in% GO_DESeq2$ID){
              GO_DESeq2_termp <- -log(GO_DESeq2[GO_DESeq2$ID == term, "p.adjust"])
            }
          }
          
          padj_GO <- rbind(padj_GO, c(GO_DrugZ_termp, GO_Limma_termp, GO_DESeq2_termp, GO_JACKS_termp))
          
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
if(1){
 
   hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
   
   wd <- paste(hwd, "Comparative_analysis", sep="/")
   if(!dir.exists(wd)){
      dir.create(wd, showWarnings = F)
   }
   print(wd)
   setwd(wd)
 
  
  for(aname in names(stat_table_list)){
    atable <- stat_table_list[[aname]]
    plot_heatmap(na.omit(atable), aname)
  }
  
  if(0){
     top <- 1000
     ex <- "synth"
     screen <- "S_"
     test <- "qGI"
  }
   
  top_count_list <- list()
  tops <- c(1000, 500, 100, 50, 20, 10)
  for(top in tops){
   
    top_gene_list <- list() ## by screen
    top_gene_list_by_test <- list() ## by test
    top_gene_list_by_screentest <- list() ## by screen and test
    print(top)
    
    for(ex in c("supp", "synth")){
      for(screen in all_screens_){
         
        pval_combined_filtered <- pval_df_list[[ex]][[screen]]
        for(test in tests){
          screentest <- paste(screen, test, sep="")
          top_gene_list[[ex]][[screen]][[test]] <- rownames(pval_combined_filtered[order(pval_combined_filtered[,test]), ])[1:top]
          top_gene_list_by_test[[ex]][[test]][[screen]] <- rownames(pval_combined_filtered[order(pval_combined_filtered[,test]), ])[1:top]
          top_gene_list_by_screentest[[ex]][[screentest]] <- rownames(pval_combined_filtered[order(pval_combined_filtered[,test]), ])[1:top]
        }
      }
    }
  
    ## overlap among tests
    for(ex in c("synth", "supp")){
      #ex <- "synth"
      top_list <- top_gene_list_by_screentest[[ex]]
      top_mat <- union2mat(top_list)
      #upset(top_mat)
      #pheatmap(top_mat[,!grepl("Sum", colnames(top_mat))])
      write.table(top_mat, paste(ex, "allScreenTest_top", top, "gene_matrix.tab",sep="_"), sep="\t", col.names=NA, quote=F)
      for(screen in all_screens_){
        #ex <- "synth"
        #screen <- "N_"
        pval_df <- pval_df_list[[ex]][[screen]]
        pval_df["ZNF629",]
        top_list <- top_gene_list[[ex]][[screen]]
        top_mat <- union2mat(top_list)
        top_count_list[[ex]][[screen]][[as.character(top)]] <- table(top_mat[, "Sum"])
        
        top_genes <- rownames(top_mat[top_mat$Sum>0, ])
        #run_enrichGO_simpleList(top_genes, ont="BP", paste(screen, ex, "top", top, sep="_"))
        
        pdf(paste(ex, screen, "top", top, "gene_upset_plot_by_test.pdf", sep="_"), height=8, width=10)
        p <- upset(top_mat, nsets = length(tests), nintersects=NA)
        print(p)
        dev.off()
        
        pval_mat <- pval_df[rownames(top_mat), colnames(top_mat)[1:length(tests)]]
        geo_mean <- apply(pval_mat, 1, geometric.mean)
        pval_mat <- as.data.frame(cbind(pval_mat, geo_mean))
        
        top_mat <- cbind(top_mat, pval_mat[rownames(top_mat),])
        write.table(top_mat, paste(screen, ex, "top", top, "gene_matrix.tab",sep="_"), sep="\t", col.names=NA, quote=F)
        plotVenn(top_list, paste(screen, ex, "top", top, "gene_venn_plot.png",sep="_"), paste(screen, ex, "_top_genes",sep=""), "text")
      }
      system('rm *.log')
    }
      
   
    ## overlap among screens
    for(ex in c("synth", "supp")){
      for(test in tests){
        #ex <- "synth"
        #test <- "DrugZ"
        top_genes <- top_gene_list_by_test[[ex]][[test]]
        top_mat <- union2mat(top_genes)
        colnames(top_mat) <- gsub("_", "", colnames(top_mat))
        pdf(paste(ex, test, "top", top, "gene_upset_plot.pdf", sep="_"), height=8, width=20)
        p <- upset(top_mat, nsets = length(all_screens_), nintersects = NA, set_size.show=F, point.size = 1)
        print(p)
        rowlim <- min(50, nrow(filter(top_mat, Sum>1)))
        collim <- ncol(top_mat)-1
        h_mat <- top_mat[1:rowlim, 1:collim]
        h <- pheatmap(t(h_mat), cluster_cols = F)
        dev.off()
        cn <- colnames(top_mat)
        for(screen in all_screens_){
          #ex <- "synth"
          #screen <- "N_"
          pval_df <- pval_df_list[[ex]][[screen]]
          top_mat <- cbind(top_mat, pval_df[rownames(top_mat), test])
          cn <- c(cn, paste(screen, "pval", sep="_"))
        }
          
        colnames(top_mat) <- cn
        write.table(top_mat, paste(ex, test, "top", top, "gene_in_allScreens_matrix.tab",sep="_"), sep="\t", col.names=NA, quote=F)
        if(length(top_genes) < 6){
          plotVenn(top_genes, paste(ex, test, "top", top, "gene_in_allScreens_venn_plot.png",sep="_"), paste(ex, test, "top_genes",sep="_"), "text")
        }else{
          five_screens <- top_genes[1:5]
          plotVenn(five_screens, paste(ex, test, "top", top, "gene_in_5Screens_venn_plot.png",sep="_"), paste(ex, test, "top_genes",sep="_"), "text")
        }
        system('rm *.log')
        
      }
      
    }
    
    if(1){
       screen_pairs_list <- combn(all_screens_, 2, simplify = F)
       test_pairs_list <- combn(tests, 2, simplify = F)
       for(ex in c("synth", "supp")){
         
         pdf(paste0("pairwise_screen_percentOverlap_top_", top, "_", ex, ".pdf"), height=8, width=10)
         par(mfrow=c(2,3))
        for(test in tests){
          #ex <- "synth"
          #test <- "qGI"
          
          pair_overlap_mat <- matrix(0, nrow=length(all_screens_), ncol=length(all_screens_), dimnames = list(all_screens_, all_screens_))
          for( i in seq_along(screen_pairs_list)){
             apair <- screen_pairs_list[[i]]
             common <- intersect(top_gene_list[[ex]][[apair[1]]][[test]], top_gene_list[[ex]][[apair[2]]][[test]])
             pair_overlap_mat[apair[1], apair[2]] <- 100*length(common)/length(top_gene_list[[ex]][[apair[1]]][[test]])
             pair_overlap_mat[apair[2], apair[1]] <- 100*length(common)/length(top_gene_list[[ex]][[apair[2]]][[test]])
          }
          colnames(pair_overlap_mat) <- gsub("_", "", colnames(pair_overlap_mat))
          rownames(pair_overlap_mat) <- gsub("_", "", rownames(pair_overlap_mat))
          corrplot(pair_overlap_mat, is.corr=F, method = 'square', type = 'lower', order = 'alphabet', diag = FALSE, 
                   tl.col = 'black', tl.srt = 0, col.lim = c(0, 100), addCoef.col = 'grey10')
          
          text(x=7, y=9, labels=test, cex=1.5, pos=1)
        
          
        }
         
         dev.off()
         
         
         pdf(paste0("pairwise_test_percentOverlap_top_", top, "_", ex, ".pdf"), height=8, width=10)
         par(mfrow=c(3,4))
         for(screen in all_screens_){
            
            pair_overlap_mat <- matrix(0, nrow=length(tests), ncol=length(tests), dimnames = list(tests, tests))
            for( i in seq_along(test_pairs_list)){
               apair <- test_pairs_list[[i]]
               common <- intersect(top_gene_list[[ex]][[screen]][[apair[1]]], top_gene_list[[ex]][[screen]][[apair[2]]])
               pair_overlap_mat[apair[1], apair[2]] <- 100*length(common)/length(top_gene_list[[ex]][[screen]][[apair[1]]])
               pair_overlap_mat[apair[2], apair[1]] <- 100*length(common)/length(top_gene_list[[ex]][[screen]][[apair[2]]])
            }
            colnames(pair_overlap_mat) <- gsub("_", "", colnames(pair_overlap_mat))
            rownames(pair_overlap_mat) <- gsub("_", "", rownames(pair_overlap_mat))
            corrplot(pair_overlap_mat, is.corr=F, method = 'square', type = 'lower', order = 'alphabet', diag = FALSE, 
                     tl.col = 'black', tl.srt = 0, col.lim = c(0, 100), addCoef.col = 'grey10')
            
            text(x=3.5, y=4, labels=gsub("_", "",screen), cex=1.5, pos=1, offset=1.2)
         }
         
         dev.off()
         par(mfrow=c(1,1))
      }
    }
  }
  
  df <- data.frame()
  for(ex in c("supp", "synth")){
    count_list <- top_count_list[[ex]]
    overlap5 <- matrix(0, ncol=length(all_screens_), nrow=length(tops), dimnames=list(as.character(tops),all_screens_))
    for(screen in all_screens_){
      for(top in as.character(tops)){
        t <- count_list[[screen]][[top]]
        if(is.element("5", names(t))){
          overlap5[top, screen] <- t["5"]
          count <- t["5"]
        }else{
          count <- 0
        }
        percent <- count*100/as.numeric(top)
        ascreen <- gsub("_", "", screen)
        df <- rbind(df, c(ex, ascreen, top, count, percent))
      }
    }
    overlap5[is.na(overlap5)] <- 0
    colnames(overlap5) <- gsub("_", "", colnames(overlap5))
    overlap5_percent <- overlap5*100/tops
    barplot(overlap5_percent, )
    #print(overlap4_percent)
  }
  colnames(df) <- c("Ex", "Screen", "Top", "Count", "Percent")
  df <- df %>% mutate_at(c("Screen"), factor) %>%
    mutate_at(c("Count", "Percent"), as.numeric) %>%
    mutate(Top=factor(Top, levels=c(10,20,50,100,500, 1000))) %>%
    mutate(Ex=factor(Ex, levels=c("synth", "supp")))   
  
  
  p <- ggplot(df, aes(fill=Top, y=Percent, x=Top)) +
    scale_fill_manual(values=color_store[1:length(tops)]) +
    geom_bar(position="dodge", stat="identity") +
     theme_bw() +
    facet_grid(Ex~Screen) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p
  pdf("percent_common_to_5methods_in_top_genes.pdf", width=12, height=8)
  print(p)
  dev.off()
  
  if(0){
    
    pdf(paste(screenn, "sig_gene_data_plot_supp.pdf", sep="_"), width=11, height=8)
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
if(0){
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
top100_list[["merged"]] <- head(res_merged$DrugZ$S_T18$GENE, 100)
top100_list[["individual"]] <- head(res_individual$DrugZ$S_T18$GENE, 100)
venn(top100_list)

combined_table <- merge(res_merged$DrugZ$S_T18, res_individual$DrugZ$N_T18, by=c("GENE"))

ggscatterhist(combined_table, x="normZ.x", y="normZ.y", add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x = -5, label.sep = "\n"),
              margin.plot = "histogram"
)

## BAGEL analysis ####

if(1){
  gsP <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/training_essentials_CRISPR.txt")
  gsN <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/training_nonessential.txt")
  
  gsP <- update_gene_names(name_map, gsP$Gene)
  gsN <- update_gene_names(name_map, gsN$Gene)
  
  #hwd <- "C:/data/raw/EDYTA/dropoutScreen/GFP_final"
  hwd <- "C:/data/raw/EDYTA/ebvScreen/GFP_final"
  setwd(hwd)
  lfc_table <- read.delim("combined_guide_lfc_table.tab")
  mat <- lfc_table[,3:ncol(lfc_table)]
  mat[is.na(mat)] <- 0
  na_replaced_table <- cbind(lfc_table[, 1:2], mat)
  write.table(na_replaced_table, "combined_guide_lfc_table_for_bagel.tab", row.names=FALSE, quote=FALSE, sep="\t")
  ## copy file to rc cluster for bagel scoring
  sshSession <- ssh_connect(host="shuyepu@rc01.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  #scp_upload(sshSession, "combined_guide_lfc_table.tab", "./bagel-for-knockout-screens-code/data") 
  scp_upload(sshSession, "combined_guide_lfc_table_for_bagel.tab", "./bagel-for-knockout-screens-code/ebvScreen")
  
  ssh_disconnect(sshSession)
  ###
  ### run run_BAGEL_oneInput.sh on rc cluster now
  ###

  dataset <- "bagel"
  wd <- file.path(hwd, dataset)
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = F)
  }
  print(wd)
  setwd(wd)
  

  print("Processing BAGEL results")
  
  sshSession <- ssh_connect(host="shuyepu@bc2.ccbr.utoronto.ca", keyfile="C:/cygwin64/home/greenblatt/.ssh/id_rsa")
  
  scp_download(sshSession, files="./bagel-for-knockout-screens-code/ebvScreen", to=".")
  system("mv ./ebvScreen/*_bagel.bf .")
  
  ssh_disconnect(sshSession)
  
  
  ## start processing bf data
  bagel_results_list <- list()
  
  
  bagel_files <- list.files(pattern="_bagel.bf")
  
  for (bagel_file in bagel_files){
    bf <- read.delim(bagel_file)
    subject <- gsub("_bagel.bf", "", bagel_file, fixed=TRUE)
    bagel_results_list[[subject]] <- bf
  }
  
  HEK293T <- update_gene_names(name_map, bagel_results_list$HEK293T$GENE)
  bagel_results_list$HEK293T$GENE <- HEK293T
  
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
  screens <- unlist(lapply(colnames(bf_table), function(x) gsub("_BF", "", x, fixed=TRUE)))
  colnames(bf_table) <- screens
  dim(bf_table)
  NAs <- apply(bf_table, 1, function(x) sum(is.na(x))>0)
  NA_table <- bf_table[NAs,]
  bf_table[is.na(bf_table)] <- 0
  
  
  pdf("BF_all_screen_correlation.pdf")
  col0 <- colorRampPalette(c("cyan4", "red3"))
  corrplot(cor(bf_table), method="number", is.corr = F, col=col0(20), type = "upper", order = "hclust", cl.lim=c(0.65, 1),
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
  #screens <- unlist(lapply(names(fdr_stats_bf$FDR), function(x) gsub("_BF", "", x, fixed=TRUE)))
  colors <- terrain.colors(length(screens))
  colors[which(screens == "HEK293T")] <-  "black"
  plot(fdr_stats_bf$FDR$HEK293T$Precision, fdr_stats_bf$FDR$HEK293T$Recall, type="l", col=colors[1], lwd=2, main="BF score ROC", xlab="Precision", ylab="Recall")
  for(i in 2:length(screens)){
    lines(fdr_stats_bf$FDR[[i]]$Precision, fdr_stats_bf$FDR[[i]]$Recall, col=colors[i], lwd=2)
  }
  abline(v=0.95, lty=2, lwd=1)
  abline(h=0.90, lty=2, lwd=1)
  legend("bottomleft", legend=screens, col=colors, lwd=2)
  dev.off()
  
  setwd("C:/data/raw/EDYTA/dropoutScreen")
  cfile <- "common_essentials_Depmap_public_2022Q2.tab"
  hfile <- "ListofEssentialGenesHEK293.tab"
  op <- "Overlap_of_essential"
  plot_overlap_genes(c(cfile, hfile), c(1,1),  outPrefix = op)
  
  
  ## G3_table3 not relevant for now
  if(0){
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
    
    
    gsp <- intersect(gsP, rownames(G3_table3))
    gsn <- intersect(gsN, rownames(G3_table3))
    
    
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
  
  
}

### END ####
