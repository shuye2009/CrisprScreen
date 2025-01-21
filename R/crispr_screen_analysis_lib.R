

setup <- function(){
  library(pheatmap)
  library(gplots)
  library(corrplot)
  library(ggpubr)
  library(patchwork)
  library(Rtsne)
  library(psych)
  library(DESeq2)
  library(knitr)
  library(dplyr)
  library(plotrix)
  library(pcaMethods)
  library(PerformanceAnalytics)
  library(VennDiagram)
  library(limma)
  library(GO.db)
  library(GOSemSim)
  library(ggpointdensity)
  library(viridis)
  library(EnhancedVolcano)
  library(GenomicPlot)
  library(RColorBrewer)
  
  setwd("C:/GREENBLATT/Rscripts/CrisprScreen")
  source("../RNAseq/Function_analysis_for_RNAseq_lib.R")
  color_store <<- brewer.pal(n = 8, name = "Dark2")
  mapfile <- "./resources/hgnc_geneName_updates_2019.txt"
  name_map <<- read.table(mapfile, header=T, stringsAsFactors = FALSE)
  head(name_map)
  dim(name_map)
  
  #color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")
  
  GOBP_file <<- "C:/GREENBLATT/resource/GSEA_gmt/c5.go.bp.v7.5.1.symbols.gmt"
  PATHWAY_file <<- "C:/GREENBLATT/resource/GSEA_gmt/ReactomePathways.gmt"
  
  corum <<- read.delim("./resources/coreComplexes.txt")
  idMap <<- read.delim("./resources/Human_id_map.tab")
  
  TKO3 <<- read.delim("./resources/TKOv2.1-Human-Library.txt")
  oldName <- TKO3[, "GENE"]
  TKO3[, "GENE"] <- update_gene_names(name_map, TKO3[, "GENE"])
  newName <<- TKO3[, "GENE"]
  
  names(newName) <<- oldName
  
  for(i in 1:length(oldName)){
    if(oldName[i] !=  newName[i]){
      cat(paste(oldName[i], "\t", newName[i], "\n"))
    }
  }
  
  TKO3_genes <<- unique(TKO3[, "GENE"])
  TKO3_guides <<- unique(TKO3[, "SEQUENCE"])
  dim(TKO3)
  length(TKO3_genes)
  length(TKO3_guides)
  
  control_gene <<- c("EGFP", "LacZ", "luciferase")
  
}

arrangeGOterms <- function(terms, cat="BP"){
  #terms <- rownames(terms_padjmean)
  #cat <- "BP"
  if(length(terms) < 3){
    return(NULL)
  }
  
    d <- godata('org.Hs.eg.db', ont=cat, computeIC=FALSE)
    
    descriptions <- NULL
    simMat <- matrix(0, nrow=length(terms), ncol=length(terms))
    
    for(i in 1:length(terms)){
      term1 <- terms[i]
      description <- Term(GOTERM[GOID(GOTERM) == term1])
      descriptions <- c(descriptions, description)
      for(j in 1:length(terms)){
        term2 <- terms[j]
        simMat[i,j] <- goSim(term1, term2, semData=d, measure="Wang")
        if(i == j){
          simMat[i,j] <- 0
        }
      }
      print(paste(i, description))
    }
    simMat
    rownames(simMat) <- descriptions
    colnames(simMat) <- descriptions
    
    maxNeighbore <- apply(simMat, 1, function(x) colnames(simMat)[which.max(x)])
    
    usedTerms <- NULL
    unUsedTerms <- descriptions
    names(unUsedTerms) <- NULL
    n_usedTerms <- 0
    
    startTerm <- unUsedTerms[1]
    
    restart <-0 
    
    while(n_usedTerms < length(descriptions)){
      
      usedTerms <- c(usedTerms, startTerm)
      n_usedTerms <- n_usedTerms + 1
      unUsedTerms <- setdiff(unUsedTerms, startTerm)
      aNeighbore <- maxNeighbore[startTerm]
      names(aNeighbore) <- NULL
      if(!aNeighbore %in% usedTerms){
        startTerm <- aNeighbore
      }else{
        for(iter in 1:3){
          for(aTerm in unUsedTerms){
            xNeighbore <- maxNeighbore[aTerm]
            if(xNeighbore %in% usedTerms){
              usedTerms <- c(usedTerms, aTerm)
              n_usedTerms <- n_usedTerms + 1
              unUsedTerms <- setdiff(unUsedTerms, aTerm)
            }
          }
        }
        if(length(unUsedTerms) > 0){
          startTerm <- unUsedTerms[1]
          
          restart <- restart + 1
          print(paste(restart, "restarting ...", startTerm))
          
        }
      }
    }
    
      
    reordered <- usedTerms
      
    
    seq <- NULL
    for(i in 1:length(reordered)){
      for(j in 1:length(descriptions)){
        if(reordered[i] == descriptions[j]){
          seq[i] <- names(descriptions[j])
        }
      }
    }
    
    names(reordered) <- seq
    
    #for(mea in c("Wang", "Jiang", "Rel", "Lin", "Resnik")){
    #  id1 <- GOID(GOTERM[Term(GOTERM) == "telomere maintenance"])
    #  id2 <- GOID(GOTERM[Term(GOTERM) == "telomere organization"])
    #  score <- goSim(id1, id2, semData=d, measure=mea) 
    #  print(paste(mea, score))
    #}
    
    return(reordered)
  
  
}

sumstats_row <- function(z){
  Sum <- apply(z, 1, sum)
  Mean <- apply(z, 1, mean)
  Median <- apply(z, 1, median)
  SD <- apply(z, 1, sd)
  SE <- apply(z, 1, function(x) sd(x)/sqrt(length(x)))
  CV <- apply(z, 1, function(x) sd(x)/mean(x))
  nonZero <- apply(z, 1, function(x) length(x[x != 0]))
  N <- apply(z, 1, length)
  result <- data.frame(Sum, Mean, Median, SD, SE, CV, nonZero, N)
  return(result)
}

sumstats_col <- function(z){
  Sum <- apply(z, 2, sum)
  Mean <- apply(z, 2, mean)
  Median <- apply(z, 2, median)
  SD <- apply(z, 2, sd)
  SE <- apply(z, 2, function(x) sd(x)/sqrt(length(x)))
  CV <- apply(z, 2, function(x) sd(x)/mean(x))
  nonZero <- apply(z, 2, function(x) length(x[x != 0]))
  N <- apply(z, 2, length)
  result <- data.frame(Sum, Mean, Median, SD, SE, CV, nonZero, N)
  return(result)
}

split_sample <- function(x, a, b){
  y<-unlist(strsplit(x, split="_")) 
  z <- y[a:b] ## 4:6 for Mar02_2020, 2:4 for Nov26_2019
  z <- paste(z, collapse="_")
  return(z)
}

split_group <- function(x, a, b){
  y<-unlist(strsplit(x, split="_")) 
  z <- y[a:b]  ## 4:5 for Mar02_2020, 2:3 for Nov26_2019
  z <- paste(z, collapse="_")
  return(z)
}


## produce an updated list of gene names, the output list preserves the order. To ensure one-to-one mapping,
## if the multiple old names map to a new name, only the first old name is updated, the rest remain the same
## if one old name maps to multiple new names, only the first new name is used
update_gene_names <- function(nameMap, nameList){
  #nameMap <- name_map
  #nameList <- old
  
  
  nameList_uptodate <- unique(nameList[nameList %in% nameMap$New])
  nameList_outdate <- unique(nameList[!nameList %in% nameMap$New])
  
  nameList_mapped <- NULL  ## to be output
  
  n_uptodate <- 0
  n_updated <- 0
  n_notupdated <- 0
  
  new_old_mapping <- list() # if more than 1 old map to 1 new, only allow 1 old to map, other old will remain
  
  for(aname in nameList){
    if(is.element(aname, nameList_uptodate)){
      nameList_mapped <- c(nameList_mapped, aname)
      n_uptodate <- n_uptodate + 1
    }else if(is.element(aname, nameMap$Old)){
      new_names <- nameMap[nameMap$Old %in% aname, "New"]
      new_name <- new_names[1] 
      
      if(is.null(new_old_mapping[[new_name]])){
        new_old_mapping[[new_name]] <- aname 
        ## make sure no other old map to this new, only allow 1 old name map to a new name
      }  
      
      ## mappable is the allowed old name, this also allows multiple instances of an old name to map to the same new name
      mappable <- new_old_mapping[[new_name]] 
      
      ## if the input list is unique, to maintain uniqueness, a new name is not allowed to collide with uptodate names
      if(aname == mappable && !is.element(new_name, nameList_uptodate)){ 
        nameList_mapped <- c(nameList_mapped, new_name)
        n_updated <- n_updated + 1
        #print(paste(aname, "is updated to", new_name))
      }else{
        nameList_mapped <- c(nameList_mapped, aname) 
        n_notupdated <- n_notupdated + 1
        #print(paste("      ", aname, "cannot be updated, due to name colliding"))
      }
      
      
    }else{  ## aname is not in the nameMap, cannot be updated, use the original name
      nameList_mapped <- c(nameList_mapped, aname) 
      n_notupdated <- n_notupdated + 1
      #print(paste("      ", aname, "cannot be updated, not in the name map"))
    }
  }

  print(paste("n_uptodate", n_uptodate))
  print(paste("n_updated", n_updated))
  print(paste("n_notupdated", n_notupdated))
  
  return(nameList_mapped)
}

format_for_BAGEL <- function(df, cn){
  #df <- df_all
  #cn <- "TEST"
  out_table <- NULL
  
  df <- t(df)
  for(i in 1:dim(df)[2]){
    #i <- 1
    gene <- rep(colnames(df)[i], length(rownames(df)))
    tmp <- data.frame(gene, df[,i])
    rownames(tmp) <- NULL
    
    out_table <- rbind(out_table, tmp, stringsAsFactors=F)
  }
  
  colnames(out_table) <- c("GENE", cn)
  return(out_table)
}

false_discovery_rate <- function(gsn, gsp, namedList, desc = TRUE){
  #namedList <- GFP_BF[GFP_BF > 0]
  namedList <- namedList[order(namedList, decreasing=desc)]
  
  fdr_list <- list()
  TP <- 0
  FP <- 0
  Precision <- 0
  Recall <- 0
  for(i in 1:length(namedList)){
    sublist <- namedList[i]
    gene <- names(sublist)
    names(sublist) <- NULL
    TP <- length(intersect(gsp, gene)) + TP
    FP <- length(intersect(gsn, gene)) + FP
    
    if(TP > 0){
      fdr_sub <- FP/(FP + TP)
    }else{
      fdr_sub <- 0
    }
    
    Precision <- 1 - fdr_sub
    Recall <- TP/length(gsp)
    arow <- c(i, sublist, gene, TP, FP, fdr_sub, Precision, Recall)
    names(arow) <- c("Index", "Score", "Gene", "TP", "FP", "FDR", "Precision", "Recall")
    fdr_list[[i]] <- arow
    
  }
  fdr <- dplyr::bind_rows(fdr_list)
  
  return(fdr)
}

## based on doi: 10.1534/g3.117.041277, cutoff_for_essential=6 
calc_FDR_cutoff <- function(gsn, gsp, bf_mat, cutoff_for_essential=6, ge = TRUE){
  fdr_list <- list()
  cutoff_list <- list()
  essential_gene_list <- list()
  for(i in 1:dim(bf_mat)[2]){
    subject <- colnames(bf_mat)[i]
    bf <- bf_mat[,i]
    names(bf) <- rownames(bf_mat)
    
    if(ge){
      fdr <- false_discovery_rate(gsn, gsp, bf)
      fdr_lt_p <- fdr[as.numeric(fdr$Score) >= cutoff_for_essential, ] 
    }else{
      fdr <- false_discovery_rate(gsn, gsp, bf, desc = FALSE)
      fdr_lt_p <- fdr[as.numeric(fdr$Score) <= cutoff_for_essential, ] 
    }
     
    cutoff_bf <- as.numeric(tail(fdr_lt_p, 1)[, "Score"])
    gene_bf <- fdr_lt_p$Gene
    
    fdr_list[[subject]] <- fdr
    cutoff_list[[subject]] <- cutoff_bf
    essential_gene_list[[subject]] <- gene_bf
  }
  return(list("FDR"=fdr_list, "cutoff"=cutoff_list, "essentialGene"=essential_gene_list))
}


essential_ROC <- function(rank_df, cutoff = 0.1, ge = FALSE){
  gsP <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/training_essentials_CRISPR.txt")
  gsN <- read.delim("C:/data/raw/EDYTA/bagel-for-knockout-screens-code/training_nonessential.txt")
  
  gsP <- update_gene_names(name_map, gsP$Gene)
  gsN <- update_gene_names(name_map, gsN$Gene)
  
  gsp <- intersect(gsP, rownames(rank_df))
  gsn <- intersect(gsN, rownames(rank_df))
  
  fdr_stats_bf <- calc_FDR_cutoff(gsn, gsp, rank_df, cutoff_for_essential=cutoff, ge = ge)
  
  essential_list <- fdr_stats_bf$essentialGene
  names(essential_list)
  
  essential_mat <- union2mat(essential_list)
  dim(essential_mat)
  tail(essential_mat)
  
  colors <- c("#1B85ED", "#A7D400", "#00BFFF", "#00CC00", "#FFD700",  "#636363", "#FFA500")
  pdf("Essential_gene_overlap_plot.pdf", height=8, width=8)
  p <- UpSetR::upset(essential_mat, main.bar.color = colors[3], sets.bar.color = colors[2])
  print(p)
  dev.off()
  
  pdf("ROC_plot.pdf", height=8, width=8)

  
  plot(1-as.numeric(fdr_stats_bf$FDR[[1]]$Precision), fdr_stats_bf$FDR[[1]]$Recall, type="l", col=colors[1], lwd=2, main="Essential gene ROC", xlab="1-Precision", ylab="Recall")
  for(i in 2:length(essential_list)){
    lines(1-as.numeric(fdr_stats_bf$FDR[[i]]$Precision), fdr_stats_bf$FDR[[i]]$Recall, col=colors[i], lwd=2)
  }
  abline(v=0.05, lty=2, lwd=1)
  abline(h=0.95, lty=2, lwd=1)
  legend("bottomright", legend=names(essential_list), col=colors[1:length(essential_list)], lwd=2)
  dev.off()
  
  return(essential_list)
}

overlap_with_CORUM <- function(inputList, idMap, corum){
  
  dim(corum)
  corum <- corum[corum$Organism == "Human", ]
  complexList <- as.list(c(corum$subunits.UniProt.IDs.))
  complexListGenes <- lapply(complexList, function(x){y<-unlist(strsplit(x, split=";")); TranslateWithSymbol(y, idMap=idMap)})
  names(complexListGenes) <- corum$ComplexName
  length(complexList)
  
  complexListGenes
  
  
  overlap_table <- NULL
  
    for(complex in unique(names(complexListGenes))){
      #complex <- "CCC complex"
      subunits <- complexListGenes[[complex]]
      
      common <- intersect(inputList, subunits)
      geneRatio <- length(common)/length(subunits)
      if(length(common) > 0){
        overlap_table <- rbind(overlap_table, c(geneRatio, paste(inputList, collapse=";"), complex, paste(subunits,collapse=";"), paste(common,collapse=";")))
      }
    }
  if(!is.null(overlap_table)){
    overlap_table <- data.frame(overlap_table)
    colnames(overlap_table) <- c("GeneRatio", "Genes", "corumComplex", "Subunits", "Common")
    head(overlap_table)
    dim(overlap_table)
    overlap_table <- overlap_table[order(overlap_table$GeneRatio, decreasing = T),]
  }
  return(overlap_table)
}


## translate uniprot sp|id|name to gene name
TranslateWithSymbol <- function(idList, idMap){
  #idList <- top
  #idMap <- geneName
  
  luni <- lapply(idList, function(prot) {
    uniprot <- prot
    if(grepl("sp\\|", prot)) {
      uniprot <- unlist(strsplit(prot, "|", fixed=TRUE))[2]
    }
    c(prot, uniprot)
  })
  
  ids <- as.data.frame(do.call(rbind, luni))
  names(ids) <- c("protein", "uniprot")
  head(ids)
  
  
  annotations.id <- merge(ids, idMap, by.x="uniprot", by.y="uniprot_acc", all.x=T)
  annotations.id <- unique(annotations.id)
  row.names(annotations.id) <- annotations.id$protein
  annotations.id <- annotations.id[idList,]
  annotations.id[is.na(annotations.id$uniprot_name), "uniprot_name"] <- annotations.id[is.na(annotations.id$uniprot_name), "uniprot"]
  
  return(annotations.id$uniprot_name)
}

countOccurance <- function(x, alist){ 
  c<-0; 
  for(i in 1:length(alist)){
    if(is.element(x, alist[[i]])){
      c<- c+1
    } 
  }
  names(c) <- x
  return(c)
}

scale_by_T0 <- function(df, operation="normalize", mincount = 40, maxcount = 10000){
  #df <- input_table_final
  #operation <- "scale"
  #head(df)
 
  T0_replicate <- names(df)[grepl("T0", names(df))]
  df <- df[!is.na(df[, T0_replicate]), ] ## remove guides that counts for T0 is missing
  df[is.na(df)] <- 0 ## missing guides are considered as 0 for non-T0 samples
  nc <- ncol(df)
  cs <- colSums(df[,3:nc])  ## column sum of samples
  normalized <- t(t(df[,3:nc])*10000000/cs) ## normalize to 10 million reads in each sample
  print(colSums(normalized, na.rm=TRUE)) ## check total normalized reads
  
  filters <- normalized[, T0_replicate] >= mincount & normalized[, T0_replicate] < maxcount
  normalized <- normalized[filters,]

  info_col <- df[filters, 1:2]
  
  if(operation == "scale"){ ## find out if T0 are started at equal values, what the T18 values will be
    T0mean <- mean(normalized[,T0_replicate], na.rm=TRUE) ## mean counts for T0
    scales <- normalized[,T0_replicate]/T0mean ## derive scale factor for each guide
    scaled <- normalized[,2:ncol(normalized)]/scales ## scale each sample
    normalized <-scaled 
  }else if(operation == "normalize"){
    message(paste0("Normalized to reads per 10 million, min T0 count: ", mincount,
                   "max T0 count: ", maxcount))
    print(colSums(normalized, na.rm=TRUE)) ## check total normalized reads
  }else if(operation == "lfc"){
    lfc <- log2(normalized[,2:ncol(normalized)]+1) - log2(normalized[,1]+1)  ## add peudocount 1, divide each T18 by T0
    normalized <- lfc
  }else{
    stop("operation not supported")
  }
  
  out <- cbind(info_col, normalized)
  print(operation)
  print(summary(out))
  
  #chart.Correlation(normalized)
  
  return(out)
}

## convert a list of vectors to a binary matrix of presence
union2mat <- function(alist){
  #alist <- synth
  u <- alist[[1]]
  for(i in 2:length(alist)){
    u <- union(u, alist[[i]])
  }
  
  mat <- NULL
  for(i in 1:length(alist)){
    v <- as.numeric(u %in% alist[[i]])
    mat <- cbind(mat, v)
  }
  mat <- as.data.frame(mat)
  rownames(mat) <- u
  colnames(mat) <- names(alist)
  mat$Sum <- apply(mat, 1, sum)
  mat <- mat[order(mat$Sum, decreasing=T),]
  
  return(mat)
}

## provide a visual inspection of data difference for significant genes
inspect_data_plot <- function(data_mat, gene, ctl, exp, ylabel){
  #data_mat <- combined_table
  #gene <- "OGA"
  #ctl <- "GFP_"
  #exp <- "N_";
  #ylabel <- "Log2(fold change)"
  colors <- c("#1B85ED", "#A7D400", "#00BFFF", "#00CC00", "#FFD700",  "#636363", "#FFA500")
  data_mat <- as.data.frame(data_mat)
  
  submat <- data_mat[data_mat$GENE==gene, c(ctl, exp)]
  submat[is.na(submat)] <- 0
  rownames(submat) <- data_mat[data_mat$GENE==gene, "SEQUENCE"]
  submat <- type.convert(t(submat), as.is = TRUE)
  
  cols <- colors[1:nrow(submat)]
  
  matplot(submat, type="b", pch=18, lwd=2, lty=c(1,1), xaxt="n", col=cols, main=gene, ylab=ylabel, ylim = c(min(submat), 1.1*max(submat)))
  staxlab(1, 1:length(rownames(submat)), rownames(submat), srt=45)
  legend("topright", col=cols, lwd = 2, legend=colnames(submat), title = "Guide", bty="n", bg="transparent")
 
  
}

## guide raw count for each guide are combined, the output "guideRawCount_", subject, ".tab"
## can be used as input for BAGEL directly
process_guideRawCount <- function(design, mincount = 40, maxcount = 10000){
  #wd <- "C:\\GREENBLATT\\Edyta\\dropoutScreen\\all_screens"
  #setwd(wd)
  #design <- eXdesign

  
  guidefiles <- rownames(design)
  length(guidefiles)
  guidefiles
  
  
  ## collected guide raw count data from guidefiles
  input_table_list <- list()
  
  print("Per guide input table size")
  for(i in 1:length(guidefiles)){
    #i <- 1
    input_table <- read.delim(guidefiles[i], header=F, sep="\t")
    prefix <- gsub("_guideRawCount.tab", "", guidefiles[i])
    replicate <- design[prefix, "replicates"]
    subject <- design[prefix, "subjects"]
    colnames(input_table) <- c("SEQUENCE", "GENE", replicate)
    
    print(paste(replicate, nrow(input_table)))
    
    old_gene_names_guide <- input_table[,2]
    new_gene_names_guide <- newName[old_gene_names_guide]
   
    input_table[,2] <- new_gene_names_guide
    
    input_table_list[[subject]][[replicate]] <- input_table
    
  }
  
  final_raw_table <- list()
  final_normalized_table <- list()
  final_lfc_table <- list()
  final_scaled_table <- list()
  
  print("Per subject input table size")
  for(subject in unique(design$subjects)){
    #subject <- "N"
    subject_input <- input_table_list[[subject]]
    T0_replicate <- names(subject_input)[grepl("T0", names(subject_input))]
    T18_replicates <- names(subject_input)[!grepl("T0", names(subject_input))]
    input_table_final <- data.frame(subject_input[[T0_replicate]])
    for(i in T18_replicates){
      input_table_final <- merge(input_table_final, subject_input[[i]], by=c("SEQUENCE", "GENE"), all=TRUE)
    }
    print(subject)
    
    print(dim(input_table_final))
    #print(head(input_table_final))
    print(summary(input_table_final))
    if(0){
      #input_table_final[is.na(input_table_final)] <- 0 ## set missing value to 0
      input_table_final <- input_table_final[!is.na(input_table_final[,T0_replicate]),] ## filter read count
      T0_cutoff <- median(input_table_final[,T0_replicate], na.rm=TRUE) * 0.05 ### use 5% of median as noise level cutoff
      print("T0 cutoff")
      print(paste(subject, T0_cutoff))
      print("Number of guides with T0 count less than T0 cutoff: ")
      print(dim(input_table_final[input_table_final[,T0_replicate] < T0_cutoff,]))
      print("Number of guides with T0 count greater than 100000: ")
      print(dim(input_table_final[input_table_final[,T0_replicate] > 100000,]))
      
      input_table_final <- input_table_final[input_table_final[,T0_replicate] >= T0_cutoff,] ## filter read count
      input_table_final <- input_table_final[input_table_final[,T0_replicate] <= 100000,] ## filter read count
    }
    
    final_raw_table[[subject]] <- input_table_final
    
    ## normalized data output, each samples is normalized to 10 million reads (include T0 and T18)
    normalized_table <- scale_by_T0(input_table_final, operation="normalize", mincount, maxcount)
    final_normalized_table[[subject]] <- normalized_table
    
    ## Log2 fold change output, each samples is normalized to 10 million reads (include T0 and T18), then log2(T18/T0)
    lfc_table <- scale_by_T0(input_table_final, operation="lfc", mincount, maxcount)
    final_lfc_table[[subject]] <- lfc_table
    
    ## scaled data output, each T18 samples is scaled to T0 sample, then normalized to 10 million ( output T18 only)
    scaled_table <- scale_by_T0(input_table_final, operation="scale", mincount, maxcount)
    final_scaled_table[[subject]] <- scaled_table
  }
  
  return(list("Raw"=final_raw_table, "Normalized"=final_normalized_table, "Scaled"=final_scaled_table, "LFC"=final_lfc_table))
  
}

process_geneRawCount <- function(design){
  
  #nameStart <- 4
  #nameEnd <- 6
  #design <- eXdesign
  
  genefiles <- rownames(design)
  genefiles <- gsub("guide", "gene", genefiles)
  length(genefiles)
  genefiles
  
  
  allGenes <- NULL
  input_table_list <- list()
  for(i in 1:length(genefiles)){
    input_table <- read.delim(genefiles[i], header=T, sep="\t")
    input_table <- input_table[, grepl("Count", colnames(input_table))]
    
    prefix <- gsub("_geneRawCount.tab", "", genefiles[i])
    replicate <- design[prefix, "replicates"]
    subject <- design[prefix, "subjects"]
    print(paste(replicate, nrow(input_table)))
    
    old_gene_names_gene <- rownames(input_table)
    new_gene_names_gene <- newName[old_gene_names_gene]
    
    length(unique(old_gene_names_gene))
    length(unique(new_gene_names_gene))
    
    rownames(input_table) <- new_gene_names_gene
    colnames(input_table) <- paste(replicate, colnames(input_table), sep="_")
    input_table_list[[subject]][[replicate]] <- input_table
    allGenes <- union(allGenes, rownames(input_table))
  }
  
  length(unique(allGenes))
  
  final_raw_table <- list()
  
  for(subject in unique(design$subjects)){
    #subject <- "N"
    subject_input <- input_table_list[[subject]]
    T0_replicate <- names(subject_input)[grepl("T0", names(subject_input))]
    T18_replicates <- names(subject_input)[!grepl("T0", names(subject_input))]
    input_table_final <- data.frame(subject_input[[T0_replicate]][allGenes, ])
    for(i in T18_replicates){
      input_table_final <- cbind(input_table_final, subject_input[[i]][allGenes, ])
    }
    head(input_table_final)
    dim(input_table_final)
    summary(input_table_final)
  
    final_raw_table[[subject]] <- input_table_final
  
  }
  
  return(final_raw_table)
}


merge_guideRawCount <- function(guidefiles, groupName){
  #wd <- "C:/RSYNC/Edyta/dropoutScreen/Mar18_2021/ORF9B"
  #guidefiles <- T0_files
  
  
  input_table_list <- list()
  for(i in 1:length(guidefiles)){
    input_table <- read.delim(guidefiles[i], header=F, sep="\t")
    colnames(input_table) <- c("SEQUENCE", "GENE", paste("rep", i, sep="_"))
    
    input_table_list[[i]] <- input_table
  }

  input_table_final <- data.frame(input_table_list[[1]])
  for(i in 2:length(input_table_list)){
    input_table_final <- merge(input_table_final, input_table_list[[i]], by=c("SEQUENCE", "GENE"), all=TRUE)
  }
  head(input_table_final)
  dim(input_table_final)
  summary(input_table_final)
  
  input_table_final[is.na(input_table_final)] <- 0
  mean_table <- cbind(input_table_final[,1:2], apply(input_table_final[,3:ncol(input_table_final)], 1, median))
  
  write.table(mean_table, paste("Greenblatt_XXX_XXX", groupName, "guideRawCount.tab", sep="_"), sep="\t", row.names=F, col.names=F, quote=F)
}

merge_geneRawCount <- function(genefiles, groupName){
  #wd <- "C:/RSYNC/Edyta/dropoutScreen/Mar18_2021/ORF9B"
  #genefiles <- T0_files
  
  allGenes <- NULL
  input_table_list <- list()
  for(i in 1:length(genefiles)){
    input_table <- read.delim(genefiles[i], header=T, sep="\t")
    #input_table <- input_table[, grepl("Count", colnames(input_table))]
    #colnames(input_table) <- c(paste("rep", i, colnames(input_table), sep="_"))
    
    input_table_list[[i]] <- input_table
    allGenes <- union(allGenes, rownames(input_table))
  }
  
  input_table_final <- data.frame(input_table_list[[1]][allGenes,])
  input_table_final[is.na(input_table_final)] <- 0
  for(i in 2:length(input_table_list)){
    tmp_table <- input_table_list[[i]][allGenes, ]
    tmp_table[is.na(tmp_table)] <- 0
    for(coln in colnames(input_table_final)[grepl("Count", colnames(input_table))]){
      input_table_final[, coln] <- input_table_final[, coln] + tmp_table[, coln]
    }
  }
  
  for(coln in colnames(input_table_final)[grepl("Count", colnames(input_table))]){
    input_table_final[, coln] <- input_table_final[, coln]/length(input_table_list)
  }
  head(input_table_final)
  dim(input_table_final)
  summary(input_table_final)
  
  input_table_final[is.na(input_table_final)] <- 0
  
  
  write.table(input_table_final, paste("Greenblatt_XXX_XXX", groupName, "geneRawCount.tab", sep="_"), sep="\t", row.names=T, col.names=T, quote=F)
}

plot_GO_term_heatmap_pval <- function(terms, genes, arrangeBy, tests, data_list){
   ## examine DrugZ results for certain GO terms
   if(0){
      terms <- c("anaphase-promoting complex-dependent catabolic process", "proteasome complex")
      genes <- c("S_", "NSP9_")
      tests <- c("DrugZ", "Limma", "DESeq2", "jacks")
      data_list <- pval_df_list[["synth"]]
      aterm <- terms[1]
      test <- tests[1]
   }
   
   gene_data_list <- data_list[genes]
   names(gene_data_list) <- genes
   
   for(test in tests){
      df <- sapply(gene_data_list, function(x)x[,test])
      plot_GO_term_heatmap(terms, arrangeBy, df)
   }
   
}

plot_GO_term_heatmap <- function(terms, arrangeBy, df, measurement=NULL){
   ## examine DrugZ results for certain GO terms
   if(0){
      terms <- c("anaphase-promoting complex-dependent catabolic process", "proteasome complex")
      genes <- c("S_", "NSP9_")
      tests <- c("DrugZ", "Limma", "DESeq2", "jacks")
      data_list <- pval_df_list[["synth"]]
      aterm <- terms[1]
      test <- tests[1]
   }
   #pGO <- build_term2gene_table("c5")
   for( aterm in terms){
      if(1){
         go_id = GOID(GOTERM[Term(GOTERM) == aterm]) 
         allegs = get(go_id, org.Hs.egGO2ALLEGS)
         
         agenes = unique(unlist(mget(allegs,org.Hs.egSYMBOL)))
      }
      
      if(0) agenes <- as.character(pGO[pGO$term==aterm, "gene"])
         
      rn <- rownames(df)[rownames(df) %in% agenes]
      dfrn <- as.data.frame(df[rn, ])
      colnames(dfrn) <- gsub("_$", "", colnames(df))
      dfrn <- dfrn %>%
         arrange(.data[[arrangeBy]]) %>%
         as.matrix()
     
      if(min(dfrn) < 0){
         med <- 0
      }else{
         med <- (max(dfrn)-min(dfrn))/2
      }
      
      p <- ComplexHeatmap::Heatmap(dfrn,
              name=measurement,
              col = circlize::colorRamp2(c(min(dfrn), med, max(dfrn)), c("darkblue", "white", "darkred")),
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows=FALSE,
              column_title = aterm,
              column_title_side = "top",
              row_names_side="left",
              row_names_gp = gpar(fontsize = 7))
      print(p)
      #print(pheatmap(dfrn, color=circlize::colorRamp2(c(min(dfrn), med, max(dfrn)), c("#3e26a7", "#13beb7", "#f9fb15")),
       #              cluster_cols=T, cluster_rows=F, display_numbers = F, fontsize_row=7, main=aterm))
      
   }
}

reorderGOTerms <- function(terms, k=6, ont="BP"){
   
   go_id = Term(GOTERM)[Term(GOTERM) %in% terms]
   
   Sim <- termSim(names(go_id), names(go_id), method="Wang", semData=SemList[[ont]])
   colnames(Sim) <- rownames(Sim) <- go_id
   
   d <- dist(Sim, method = "euclidean") # distance matrix
   fit <- hclust(d, method="average") 
   
   h <- Heatmap(Sim, show_column_names = F, clustering_distance_rows = "euclidean",
                clustering_distance_columns = "euclidean",
                clustering_method_rows = "average", clustering_method_columns = "average")
   #h <- draw(h)
   #row_order(h)
  
   # trim the tree, k for number of groups, h for height
   fit.k <- cutree(fit, k=k)
   fit.k <- sort(fit.k[row_order(h)]) 
   
   return(fit.k)
}


plot_GSEA <- function(screens, tests, data_list, ont="BP"){
   
   ## examine DrugZ results for certain GO terms
   if(0){
      goALLfile <- "C:/GREENBLATT/resource/GSEA_gmt/all.go.v7.5.1.symbols.gmt"
      goIMMUNEfile <- "C:/GREENBLATT/resource/GSEA_gmt/immune.all.v7.5.1.symbols.gmt"
      goBPfile <- "C:/GREENBLATT/resource/GSEA_gmt/c5.go.bp.v7.5.1.symbols.gmt"
      goCCfile <- "C:/GREENBLATT/resource/GSEA_gmt/c5.go.cc.v7.5.1.symbols.gmt"
      goMFfile <- "C:/GREENBLATT/resource/GSEA_gmt/c5.go.mf.v7.5.1.symbols.gmt"
      
      screens <- c("S_")
      tests <- c("DrugZ_norm")
      data_list <- stat_table_list
      screen <- screens[1]
      test <- tests[1]
      ont <- "BP"
      pval <- 0.05
   }
   
   for(test in tests){
      df <- stat_table_list[[test]]
      for( screen in screens){
         acol <- df[, screen]
         names(acol) <- rownames(df)
         acol[is.na(acol)] <- 0
         plot(density(acol))
         run_gseGO_simpleList(acol, paste(screen, test, sep=""), ont=ont)
         
         res <- GSEA(acol, goBPfile, pval)
         res$Results
         res$Plot
      }
   }
}

#design <- eXdesign
#lfc_table <- combined_lfc_guide_table
compute_lfc_residuals <- function(design, lfc_table, loess=TRUE, plot=FALSE){
   
   subjects <- unique(design$subjects)
   controls <- unique(design$controls)
   rownames(design) <- NULL
   subject_control_map <- design %>%
      select(c(subjects, controls)) %>%
      filter(subjects != controls) %>%
      unique()
   
   mean_lfc <- sapply(subjects, function(x){
      df <- lfc_table[, grepl(paste0(x, "_"), colnames(lfc_table))]
      m <- rowMeans(df, na.rm=TRUE)
   })
   mean_lfc <- as.data.frame(mean_lfc)
   noNArows <- !apply(mean_lfc,1, function(x)any(is.na(x)))
   mean_lfc <- mean_lfc[noNArows,] 
   ## remove row with any  NA
   
   diff_lfc <- sapply(seq.int(nrow(subject_control_map)), function(x){
      treat <- subject_control_map[x, "subjects"]
      control <- subject_control_map[x, "controls"]
      residual <- mean_lfc[,treat] - mean_lfc[, control]
   })
   colnames(diff_lfc) <- subject_control_map$subjects
   
   diff_lfc <- cbind(lfc_table[noNArows, 1:2], diff_lfc)
   
   loess_fitted_lfc <- sapply(seq.int(nrow(subject_control_map)), function(x){
      treat <- subject_control_map[x, "subjects"]
      control <- subject_control_map[x, "controls"]
      predicted <- loessFit(y=diff_lfc[,treat], x=mean_lfc[, control], span=0.4)
      predicted$fitted
   })
   loess_fitted_lfc <- as.data.frame(loess_fitted_lfc)
   colnames(loess_fitted_lfc) <- subject_control_map$subjects
   
   loess_residuals_lfc <- sapply(seq.int(nrow(subject_control_map)), function(x){
      treat <- subject_control_map[x, "subjects"]
      control <- subject_control_map[x, "controls"]
      predicted <- loessFit(y=diff_lfc[,treat], x=mean_lfc[, control], span=0.4)
      predicted$residuals
   })
   loess_residuals_lfc <- as.data.frame(loess_residuals_lfc)
   colnames(loess_residuals_lfc) <- subject_control_map$subjects
   
   loess_residuals_lfc <- cbind(lfc_table[noNArows, 1:2], loess_residuals_lfc)
   
   ## plot LOESS regression with respect to GFP1 control
   if(plot){
      pdf("qGI_diagnostic_plots.pdf", height = 8, width = 8)
      if(ncol(diff_lfc) > 3){
        chart.Correlation(diff_lfc[, 3:ncol(diff_lfc)], main="diff of lfc", 
                          line.main=1.5, oma=c(2,2,3,2))
        chart.Correlation(loess_residuals_lfc[, 3:ncol(loess_residuals_lfc)], 
                          main="loess residual of lfc", line.main=1.5, 
                          oma=c(2,2,3,2))
      }
      for(i in 1:nrow(subject_control_map)){
         subject <- subject_control_map[i, "subjects"]
         control <- subject_control_map[i, "controls"]
         if(subject != control){
            p1 <- ggplot(data = mean_lfc, 
                         mapping = aes(x = .data[[control]], y = .data[[subject]])) +
               geom_pointdensity(show.legend = F) +
               scale_color_viridis() + 
               geom_smooth(method = "loess", span=0.1, size = 1.5, se=FALSE) +
               theme_bw() + 
               labs(x=paste0("mean lfc(",control,")"), 
                    y=paste0("mean lfc(", subject, ")")) +
               ggtitle(paste0("Mean lfc(", subject, ") ~ Mean lfc(", control, ")"))
            
            
            #plot(mean_lfc[,"GFP1"], mean_lfc[,"S"])
            #points(loess.smooth(y=mean_lfc[,"S"], x=mean_lfc[,"GFP1"], span=0.4), col="red")
            
            p2 <- ggplot(data = data.frame(acontrol=mean_lfc[,control], asubject=diff_lfc[,subject]), mapping = aes(x = acontrol, y = asubject)) +
               geom_pointdensity(show.legend = F) +
               scale_color_viridis() + 
               geom_smooth(method = "loess", span=0.1, size = 1.5, se=FALSE) +
               theme_bw() + 
               labs(x=paste0("mean lfc(",control,")"), y=paste0("diff lfc(", subject, ")")) +
               ggtitle(paste0("Diff lfc(", subject, ") ~ Mean lfc(", control, ")")) +
               geom_hline(yintercept=0, color="red", linetype="dashed", size=1)
            
            
            #plot(mean_lfc[,"GFP1"], diff_lfc[,"S"])
            #points(mean_lfc[,"GFP1"], loess_fitted_lfc[,"S"], col="red3")
            
            #plot(mean_lfc[,"GFP1"], loess_residuals_lfc[,"S"])
            #points(mean_lfc[,"GFP1"], loess_fitted_lfc[,"S"], col="red3")
            
            p3 <- ggplot(data = data.frame(acontrol=mean_lfc[,control], asubject=loess_residuals_lfc[,subject]), mapping = aes(x = acontrol, y = asubject)) +
               geom_pointdensity(show.legend = F) +
               scale_color_viridis() + 
               geom_smooth(method = "loess", span=0.1, size = 1.5, se=FALSE) +
               theme_bw() + 
               labs(x=paste0("mean lfc(",control,")"), y=paste0("residual lfc(", subject, ")")) +
               ggtitle(paste0("Residual lfc(", subject, ") ~ Mean lfc(", control, ")")) +
               geom_hline(yintercept=0, color="red", linetype="dashed", size=1)
            
            
            #plot(mean_lfc[,"GFP1"], loess_residuals_lfc[,"S"])
            #points(loess.smooth(y=loess_residuals_lfc[,"S"], x=mean_lfc[,"GFP1"], span=0.4), col="red")
            
            p4 <- ggplot(data = data.frame(diff=diff_lfc[,subject], asubject=loess_residuals_lfc[,subject]), mapping = aes(x = diff, y = asubject)) +
               geom_pointdensity(show.legend = F) +
               scale_color_viridis() + 
               geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1) +
               theme_bw() + 
               labs(x=paste0("diff lfc(",subject,")"), y=paste0("residual lfc(", subject, ")")) +
               ggtitle(paste0("Residual lfc(", subject, ") ~ Diff lfc(", subject, ")"))
            
            outp <- plot_grid(p1, p2, p3, p4, ncol = 2, align = 'v', axis="b")
            print(outp)
            #plot(loess_fitted_lfc[,"S"], loess_residuals_lfc[,"S"])
            #points(loess.smooth(y=loess_residuals_lfc[,"S"], x=loess_fitted_lfc[,"S"], span=0.4), col="red")
            
            ## plot LOESS regression with respect to subject S
            #plot(mean_lfc[,"S"], diff_lfc[,"S"])
            #points(mean_lfc[,"S"], loess_fitted_lfc[,"S"], col="red3")
            
            #plot(mean_lfc[,"S"], loess_residuals_lfc[,"S"])
            #points(loess.smooth(y=loess_residuals_lfc[,"S"], x=mean_lfc[,"S"], span=0.4), col="red")
            
         }
      }
      dev.off()
   }
      
   if(loess){
      return(loess_residuals_lfc)
   }else{
      return(diff_lfc)
   }
}

modt_lfc_residuals <- function(design, lfc_residuals){
   #design <- eXdesign
   #lfc_residuals
   
   genes <- unique(lfc_residuals$GENE)
   subjects <- unique(design$subjects)
   controls <- unique(design$controls)
   subjects <- subjects[!subjects %in% controls]
   
   guides <- c("g1", "g2", "g3", "g4")
   
   allGenes <- sort(unique(lfc_residuals$GENE))
   allGenes[1:10]
   lfc_expanded <- NULL
   c <- 0 #genes with less than 2 guides
   valid_gene <- NULL
   for(gene in allGenes){
      sub <- as.matrix(lfc_residuals[lfc_residuals$GENE %in% gene, 3:ncol(lfc_residuals)])
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
   
   cn <- unlist(lapply(colnames(lfc_residuals)[3:ncol(lfc_residuals)], function(x) paste(x, guides, sep="_")))
   colnames(lfc_expanded) <- cn
   rownames(lfc_expanded) <- valid_gene
   
   dim(lfc_expanded)
   
   test_results <- lapply(subjects, function(s){
      #s <- "S"
      X <- lfc_expanded[, grepl(paste0(s,"_"), colnames(lfc_expanded))]
      designX <- matrix(1, nrow = ncol(X), ncol = 1)
      colnames(designX) <- s
      fit1 <- lmFit(X, designX)
      fit2 <- eBayes(fit1)
      
      return(fit2)
   })
   names(test_results) <- subjects
   return(test_results)
}

plot_qGI_results <- function(test_results, pvalue_cutoff=0.05, gsea=FALSE,
                             topN = 100, bottomN = 100){
   
   type <- "qGI"
   
   qGI_results <- lapply(names(test_results), function(s){
      print(s)
      fit.eb <- test_results[[s]]
      
      #pdf(paste(s, "_qGI_volvanoplot.pdf", sep=""), height=8, width=8)
      #volcanoplot(fit.eb,coef=1, highlight=30, xlab="Log fold change", main=s, names=rownames(fit.eb$coefficients))
      #dev.off()
      
      Y <- topTable(fit.eb, coef = 1, number = Inf, confint = TRUE, sort.by = "p") %>%
         arrange(t)
      write.table(Y, paste(s, type, "moderated_t_test_results.tab", sep="_"), sep="\t", col.names=NA, quote=F)
      
      df <- cbind(GENE=rownames(Y), Y, rankT=rank(Y$t)) %>%
         na.omit %>%
         arrange(t)
      
      pdf(paste(s, type, "rank_plot.pdf", sep="_"), height=8, width=8)
      
      p1 <- ggscatter(df, x = "rankT", y = "t", xlab="Rank", ylab="Moderated t", title=s,
                      size = abs(df$t), color = "t") + 
         gradient_color(c("cyan4", "white", "red3")) +
         rremove("legend")
      
      p2 <- as_ggplot(text_grob(paste(rev(head(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      p3 <- as_ggplot(text_grob(paste(rev(tail(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      
      p <- p1 + 
         inset_element(p2, left = 0.1, bottom = 0.1, right = 0.3, top = 0.3) +
         inset_element(p3, left = 0.7, bottom = 0.7, right = 0.95, top = 0.95)
      print(p)
      
      p <- EnhancedVolcano(df,
                           lab = df$GENE,
                           x = 'logFC',
                           y = 'P.Value',
                           #selectLab = df$GENE[2:11],
                           title = paste(type, "analysis of", s),
                           subtitle = "",
                           caption = paste("LFC = +/-1; pvalue cutoff =", pvalue_cutoff),
                           legendLabels = c("NS", "LFC", "P-value", "P-value and LFC"),
                           xlab = "LFC",
                           pCutoff = pvalue_cutoff,
                           FCcutoff = 1,
                           labSize = 4,
                           pointSize = 3.0,
                           drawConnectors = FALSE)
      print(p)
      dev.off()
      
      if(gsea){
         Yt <- Y$t
         names(Yt) <- rownames(Y)
         
         gse1 <- run_gseGO_simpleList(Yt, s, ont="BP", GO_file=NULL, simplify = TRUE)
         gse2 <- run_gseGO_simpleList(Yt, s, ont="Reactome_pathway", GO_file=PATHWAY_file)
         
         downlist <- names(head(Yt, n = topN))
         uplist <- names(tail(Yt, n = bottomN))
         run_enrichGO_simpleList(uplist, "BP", paste(s, topN, "upGenes", sep="_"))
         run_enrichGO_simpleList(downlist, "BP", paste(s, bottomN,"downGenes", sep="_"))
      }
      
      return(Y)
   })
   
   
   names(qGI_results) <- names(test_results)
   
   ## plot t value correlation
   if(length(qGI_results)> 1) {
     t_table <- sapply(names(qGI_results), function(x){
        res <- qGI_results[[x]]
        y <- res$t; 
        names(y) <- row.names(res); 
        y <- y[sort(names(y))]
   
        return(y)
     })
     
     pdf(paste(type, "correlation_of_tvalue.pdf", sep="_"), height=8, width=8)
     chart.Correlation(t_table)
     dev.off()
     
     ## plot top gene overlap
     
     named_vector_list <- list()
     for (x in names(qGI_results)){
        res <- qGI_results[[x]]
        y <- res$t; 
        names(y) <- row.names(res); 
        
        named_vector_list[[x]] <- y
     }
     
     names(named_vector_list) <- names(qGI_results)
     tops <- c(10, 50, 100, 500, 1000)
     pdf(paste(type, "positive_top_genes_overlap.pdf", sep="_"), height=8, width=8)
     plot_overlap_among_top_genes(geneLists=named_vector_list, tops, decreasing=TRUE)
     dev.off()
     
     pdf(paste(type, "negative_top_genes_overlap.pdf", sep="_"), height=8, width=8)
     plot_overlap_among_top_genes(geneLists=named_vector_list, tops, decreasing=FALSE)
     dev.off()
   }
   
   return(qGI_results)
}

plot_DESeq2_results <- function(dds, contrasts, pvalue_cutoff=0.05, gsea=FALSE, 
                                label=FALSE, verbose=FALSE){

   DESeq2_results <- list()
   
   for(i in 1:nrow(contrasts)){
      contrast <- contrasts[i,]
      treatSF <- paste(contrast[2], "vs", contrast[3], sep="_")
      print(paste("Processing ", paste(contrast, collapse=" ")))
      
      res <- results(dds, contrast=contrast, independentFiltering = FALSE, pAdjustMethod = "BH")
      
      res <- res[order(-res$stat),]
      DESeq2_results[[treatSF]] <- res
      
      sigpv <- res[(res$padj < pvalue_cutoff+0.2),] ## low confidence cutoff
      sig <- res[(res$padj < pvalue_cutoff),] ## high confidence cutoff
      nc <- res[(res$padj >= pvalue_cutoff),]
      sigup <- sig[(sig$log2FoldChange > 0),]
      sigdown <- sig[(sig$log2FoldChange < 0),]
      
      uplist <- rownames(sigup)
      downlist <- rownames(sigdown)
      nclist <- rownames(nc)
         
      filename <- paste(treatSF, "DESeq2_results_table.tab", sep="_")
      write.table(res, filename, row.names=T, col.names=NA, quote=F, sep="\t")
     
         
      if(verbose){
         allup <- res[(res$log2FoldChange > 0),]
         alldown <- res[(res$log2FoldChange < 0),]
         
         print(dim(res)[1])
         print(dim(allup)[1])
         print(dim(alldown)[1])
         
         rld <- vst(dds)
         pdf(paste(treatSF, "PCA_plot.pdf", sep="_"))
         print(plotPCA(rld, intgroup=c("conditions")))
         dev.off()
         
         pdf(paste(treatSF, "MA_plot.pdf", sep="_"))
         DESeq2::plotMA(res, ylim=c(-2,2))
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
      }
         
      df <- as.data.frame(cbind(GENE=rownames(res), res, rankT=rank(res$stat))) %>%
         arrange(stat)
      pdf(paste(treatSF, "vocano_plot_allgenes.pdf", sep="_"), width=8, height=8)
      
      colors <- rep("black", dim(res)[1])
      colors[row.names(res) %in% row.names(sigpv)] <- "cyan"
      colors[row.names(res) %in% row.names(sig)] <- "red"
      texts <- rep("", dim(res)[1])
      texts[row.names(res) %in% row.names(sig)] <- row.names(res)[row.names(res) %in% row.names(sig)]
      plot(res$log2FoldChange, -log10(res$pvalue), col=colors, ylim=c(0, 30), 
           xlab="Log2(FoldChange)",
           ylab="-Log10(Pvalue)")
      if(label) text(res$log2FoldChange, -log10(res$pvalue), labels=texts, cex=0.5)
      
      plot(res$log2FoldChange, log(-log10(res$pvalue)), col=colors, ylim=c(-4, 4), 
           xlab="Log2(FoldChange)",
           ylab="Log(-Log10(Pvalue))")
      if(label) text(res$log2FoldChange, log(-log10(res$pvalue)), labels=texts, cex=0.5)
      
      p <- EnhancedVolcano(df,
                           lab = df$GENE,
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           #selectLab = df$GENE[2:11],
                           title = paste("DESeq2 analysis of ", treatSF),
                           subtitle = "",
                           caption = paste("LFC = +/-1; pvalue cutoff =", pvalue_cutoff),
                           legendLabels = c("NS", "LFC", "P-value", "P-value and LFC"),
                           xlab = "LFC",
                           pCutoff = pvalue_cutoff,
                           FCcutoff = 1,
                           labSize = 4,
                           pointSize = 3.0,
                           drawConnectors = FALSE)
      print(p)
      
      dev.off()
      
      
      pdf(paste(treatSF, "_DESeq2_rank_plot.pdf", sep=""), height=8, width=8)
      
      p1 <- ggscatter(df, x = "rankT", y = "stat", xlab="Rank", ylab="Wald test statistics", title=treatSF,
                      size = abs(df$stat), color = "stat") + 
         gradient_color(c("cyan4", "white", "red3")) +
         rremove("legend")
      
      p2 <- as_ggplot(text_grob(paste(rev(head(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      p3 <- as_ggplot(text_grob(paste(rev(tail(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      
      p <- p1 + 
         inset_element(p2, left = 0.1, bottom = 0.1, right = 0.3, top = 0.3) +
         inset_element(p3, left = 0.7, bottom = 0.7, right = 0.95, top = 0.95)
      print(p)
      
      dev.off()
      
      if(gsea){
         ress <- res$stat
         names(ress) <- rownames(res)
         
         run_enrichGO_simpleList(uplist, "BP", paste(treatSF,"upGenes", sep="_"))
         run_enrichGO_simpleList(downlist, "BP", paste(treatSF,"downGenes", sep="_"))
         
         run_gseGO_simpleList(ress, treatSF, ont="BP", GO_file=GOBP_file)
         run_gseGO_simpleList(ress, treatSF, ont="Reactome_pathway", GO_file=PATHWAY_file)
      }
   }

   return(DESeq2_results)
}


plot_DrugZ_results <- function(test_results, pvalue_cutoff=0.05, gsea=FALSE){
   #test_results <- qGI_results
   
   for(treat in names(test_results)){
     #treat <- "S"
      df <- test_results[[treat]] %>%
         filter(!GENE %in% control_gene) %>%
         na.omit() %>%
         arrange(normZ)
      
      rownames(df) <- df$GENE

      downlist <- df[df$pval_synth < pvalue_cutoff, "GENE"]
      uplist <- df[df$pval_supp < pvalue_cutoff, "GENE"]
      
      #df <- mutate(df, pvalue=min(pval_synth, pval_supp))
      df$pvalue <- apply(df[, c("pval_synth", "pval_supp")], 1, min)
      df$fdr <- apply(df[, c("fdr_synth", "fdr_supp")], 1, min)
      
      pdf(paste(treat, "rank_plot_allgenes.pdf", sep="_"), width=4, height=8)
      
      p1 <- ggscatter(df, x = "rank_synth", y = "normZ", xlab="Rank", ylab="NormZ", title=treat,
                      size = abs(df$normZ), color = "normZ") + 
         gradient_color(c("cyan4", "white", "red3")) +
         rremove("legend")
      
      p2 <- as_ggplot(text_grob(paste(rev(head(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      p3 <- as_ggplot(text_grob(paste(rev(tail(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      
      p <- p1 + 
         inset_element(p2, left = 0.1, bottom = 0.1, right = 0.3, top = 0.3) +
         inset_element(p3, left = 0.7, bottom = 0.7, right = 0.95, top = 0.95)
      print(p)
      dev.off()
      
      pdf(paste(treat, "vocano_plot_allgenes.pdf", sep="_"), width=8, height=8)
      p <- EnhancedVolcano(df,
                           lab = df$GENE,
                           x = 'normZ',
                           y = 'pvalue',
                           #selectLab = df$GENE[2:11],
                           title = paste("DrugZ analysis", treat),
                           subtitle = "",
                           caption = paste("NormZ score cutoff = +/-3; pvalue cutoff =", pvalue_cutoff),
                           legendLabels = c("NS", "NormZ", "P-value", "P-value and NormZ"),
                           xlab = "NormZ score",
                           pCutoff = pvalue_cutoff,
                           FCcutoff = 3,
                           labSize = 4,
                           pointSize = 3.0,
                           drawConnectors = FALSE)
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
      
      if(gsea){
         #run_gseGO_simpleList(normZ, treat, ont="MSIGDB_GOBP", GO_file=GOBP_file)
         run_gseGO_simpleList(normZ, treat, ont="BP", GO_file=NULL)
         run_gseGO_simpleList(normZ, treat, ont="Reactome_pathway", GO_file=PATHWAY_file)
      }
   }
}

plot_JACKS_results <- function(test_results, pvalue_cutoff=0.05, gsea=FALSE){
   #test_results <- jacks_results
   
   for(treat in names(test_results)){
      #treat <- "S"
      res <- as.data.frame(test_results[[treat]]) %>%
         na.omit()
      
      df <- cbind(GENE=rownames(res), res, rankT=rank(res$neg)) %>%
         filter(!GENE %in% control_gene) %>%
         arrange(rankT)
      
      #df <- mutate(df, pvalue=min(pval_synth, pval_supp))
      df$pvalue <- apply(df[, c("neg_p", "pos_p")], 1, min)
      df$fdr <- apply(df[, c("neg_padj", "pos_padj")], 1, min)
      
      pdf(paste(treat, "rank_plot_allgenes.pdf", sep="_"), width=4, height=8)
      
      p1 <- ggscatter(df, x = "rankT", y = "neg", xlab="Rank", ylab="JACKS score", title=treat,
                      size = abs(df$neg), color = "neg") + 
         gradient_color(c("cyan4", "white", "red3")) +
         rremove("legend")
      
      p2 <- as_ggplot(text_grob(paste(rev(head(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      p3 <- as_ggplot(text_grob(paste(rev(tail(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      
      p <- p1 + 
         inset_element(p2, left = 0.1, bottom = 0.1, right = 0.3, top = 0.3) +
         inset_element(p3, left = 0.7, bottom = 0.7, right = 0.95, top = 0.95)
      print(p)
      dev.off()
      
      pdf(paste(treat, "vocano_plot_allgenes.pdf", sep="_"), width=8, height=8)
      p <- EnhancedVolcano(df,
                           lab = df$GENE,
                           x = 'neg',
                           y = 'pvalue',
                           #selectLab = df$GENE[2:11],
                           title = paste("JACKS analysis", treat),
                           subtitle = "",
                           caption = paste("JACKS score cutoff = +/-1; pvalue cutoff =", pvalue_cutoff),
                           legendLabels = c("NS", "Score", "P-value", "P-value and Score"),
                           xlab = "JACKS score",
                           pCutoff = pvalue_cutoff,
                           FCcutoff = 1,
                           labSize = 4,
                           pointSize = 3.0,
                           drawConnectors = FALSE)
      print(p)
      while(!is.null(dev.list())){
         dev.off()
      }
      
      neg <- df$neg
      names(neg) <- df$GENE
      neg <- neg[TKO3_genes]
      neg <- na.omit(neg)
      
      print(length(neg))
      
      if(gsea){
         run_gseGO_simpleList(neg, treat, ont="BP", GO_file=NULL)
         run_gseGO_simpleList(neg, treat, ont="Reactome_pathway", GO_file=PATHWAY_file)
      }
   }
}


plot_Limma_results <- function(fit.eb, fit, pvalue_cutoff=0.05, gsea=FALSE, type="Limma"){
   
   limma_results <- list()
   
   for(i in 1:length(colnames(fit.eb$contrasts))){
      #i <- 3
      contrast <- colnames(fit.eb$contrasts)[i]
      
      term1 <- rownames(fit.eb$contrasts)[fit.eb$contrasts[,contrast] == 1]
      term2 <- rownames(fit.eb$contrasts)[fit.eb$contrasts[,contrast] == -1]
      
      tpall <- topTable(fit.eb, number = Inf, coef=i, p.value=1)
      limma_results[[contrast]] <- tpall
      
      tp_up <- rownames(tpall[tpall$P.Value<pvalue_cutoff & tpall$logFC > 1, ])
      tp_down <- rownames(tpall[tpall$P.Value<pvalue_cutoff & tpall$logFC < -1, ])
     
      pdf(paste(type, contrast, "volcano_plot.pdf", sep="_"), height=8, width=8)
      
      limma::volcanoplot(fit.eb,coef=i, highlight=30, xlab="Log fold change", main=contrast, names=rownames(fit.eb$coefficients))
      limma::plotMA(fit.eb,coef=i, xlab = "Average expression")
      abline(h=0, col="red2")
      plot(x=fit$coefficients[,term2], y=fit$coefficients[,term1], col="white", xlab=paste(term2, "lfc"), ylab=paste(term1, "lfc"))
      colors <- rep("black", dim(fit$coefficients)[1])
      names(colors) <- rownames(fit$coefficients)
      colors[tp_up] <- "red"
      colors[tp_down] <- "blue"
      #plot(x=broom_fit[,3], y=broom_fit[,2], col=colors, xlab="GFP_lfc", ylab=TREAT)
      text(rownames(fit$coefficients), x=fit$coefficients[,term2], y=fit$coefficients[,term1], col=colors, cex=0.5)
      lines(loess.smooth(x=fit$coefficients[,term2], y=fit$coefficients[,term1], span=0.4), col="cyan", lwd=1)
      
      dev.off()
   
      
      write.table(tpall, paste(type, contrast, "moderated_t_test_results.tab", sep="_"), sep="\t", col.names=NA, quote=F)
      
      df <- cbind(GENE=rownames(tpall), tpall, rankT=rank(tpall$t)) %>%
         na.omit %>%
         arrange(t)
      
      pdf(paste(type, contrast, "rank_plot.pdf", sep="_"), height=8, width=8)
      
      p1 <- ggscatter(df, x = "rankT", y = "t", xlab="Rank", ylab="Moderated t", title=contrast,
                      size = abs(df$t), color = "t") + 
         gradient_color(c("cyan4", "white", "red3")) +
         rremove("legend")
      
      p2 <- as_ggplot(text_grob(paste(rev(head(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      p3 <- as_ggplot(text_grob(paste(rev(tail(df$GENE, 10)),collapse="\n"), face = "italic", color = "steelblue"))
      
      p <- p1 + 
         inset_element(p2, left = 0.1, bottom = 0.1, right = 0.3, top = 0.3) +
         inset_element(p3, left = 0.7, bottom = 0.7, right = 0.95, top = 0.95)
      print(p)
      
      p <- EnhancedVolcano(df,
                           lab = df$GENE,
                           x = 'logFC',
                           y = 'P.Value',
                           #selectLab = df$GENE[2:11],
                           title = paste(type, "analysis of", term1),
                           subtitle = "",
                           caption = paste("LFC = +/-1; pvalue cutoff =", pvalue_cutoff),
                           legendLabels = c("NS", "LFC", "P-value", "P-value and LFC"),
                           xlab = "LFC",
                           pCutoff = pvalue_cutoff,
                           FCcutoff = 1,
                           labSize = 4,
                           pointSize = 3.0,
                           drawConnectors = FALSE)
      print(p)
      dev.off()
      
      if(gsea){
         rankedlist <- tpall$t
         names(rankedlist) <- rownames(tpall)
         
         run_gseGO_simpleList(rankedlist, contrast, ont="BP", GO_file=NULL)
         run_gseGO_simpleList(rankedlist, contrast, ont="Reactome_pathway", GO_file=PATHWAY_file)
      }
      
   }
 
   return(limma_results)
}

geneLists <- list()
for(i in 1:4){
   s <- sample(20)
   names(s) <- LETTERS[1:20]
   geneLists[[i]] <- s
}

tops <- c(5, 10, 15)

plot_overlap_among_top_genes <- function(geneLists, tops, decreasing=TRUE){
   overlap_count <- vapply(tops, FUN=function(x){
      top_list <- lapply(geneLists, function(y){
         y <- sort(y, decreasing=decreasing)
         return(names(y[1:x]))
      })
      common <- Reduce(intersect, top_list)
      return(length(common))
   }, FUN.VALUE = integer(1), USE.NAMES = TRUE)
   
   df <- data.frame(TOP=as.factor(tops), COUNT=overlap_count, PERCENT=round(overlap_count*100/tops,2))
   
   p <- ggplot2::ggplot(df, aes(x=TOP, y=PERCENT)) +
      geom_bar(stat="identity", position="nudge") +
      theme_classic()
   
   print(p)
   
   grid.newpage()
   lapply(tops, FUN=function(x){
      top_list <- lapply(geneLists, function(y){
         y <- sort(y, decreasing=decreasing)
         return(names(y[1:x]))
      })
      pairs <- combn(top_list, 2, simplify = FALSE)
      
      lapply(pairs, GenomicPlot::overlap_pair, intersect)
      
      if(length(top_list)>2){
         triples <- combn(top_list, 3, simplify = FALSE)
         lapply(triples, overlap_triple, intersect)
         if(length(top_list)>3){
            quads <- combn(top_list, 4, simplify = FALSE)
            lapply(quads, overlap_quad, intersect)
         }
      }
      
   })
}


scale_normalize <- function(wd, hwd, samplef = c(5,8), mincount = 40, maxcount = 10000){
  
  setwd(wd)
  format2 <- "_guideRawCount.tab"
  
  guidefiles <- list.files(path = getwd(), recursive = TRUE, 
                           full.names = FALSE, pattern = format2)
  length(guidefiles)
  guidefiles
  
  if(1){
    samples <- gsub("_guideRawCount.tab", "", guidefiles)
    replicates <- unlist(lapply(guidefiles, function(x)split_group(x,samplef[1],samplef[2]))) #4,6 for dropoutScreen
    subjects <- unlist(lapply(guidefiles, function(x)split_group(x,samplef[1],samplef[2]))) #4,4 for dropoutScreen
    design <- data.frame(samples, replicates, subjects)
    rownames(design) <- guidefiles
    design <- design[order(design$subjects),]
    design <- design[!grepl("T10|NR", design$replicates), ]  ## exclude T10 and NR samples
    if(!file.exists("experimental_design_final.tab")){
      write.table(design, "experimental_design_final.tab", sep="\t", row.names=T, col.names=NA, quote=F)
      stop("EDIT experimental_design_final.tab NOW!")
    }else{
      eXdesign <- data.frame(read.delim("experimental_design_final.tab", 
                                        header=TRUE, 
                                        sep="\t", 
                                        stringsAsFactors=F))
    }
  }
  
  rownames(eXdesign) <- eXdesign$X
  eXdesign <- eXdesign[,-1]
  
  eXdesign
  
  #sink("process_guideRawCount_log.txt")
  guide_data <- process_guideRawCount(eXdesign, mincount, maxcount) ## output can be used for BAGEL directly
  #sink()
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
  
  
  if(!dir.exists(hwd)){
    dir.create(hwd, showWarnings = F)
  }
  
  setwd(hwd)
  
  dim(combined_gene_table)
  
  print("writing combined table file")
  write.table(combined_gene_table, file.path(hwd,"combined_gene_rawCount_table.tab"), row.names=T, col.names=NA, sep="\t", quote=F)
  
  dim(combined_guide_table)
  NAc <- apply(combined_guide_table, 1, function(x) sum(is.na(x)))
  NAs <- NAc == ncol(combined_guide_table)-2
  combined_guide_table <- combined_guide_table[!NAs,] ## remove all_NA rows
  combined_guide_table[is.na(combined_guide_table)] <- 0 ## convert NA in single sample to 0
  write.table(combined_guide_table, file.path(hwd,"combined_guide_rawCount_table.tab"), row.names=F, sep="\t", quote=F)
  combined_guide_filtered <- na.omit(combined_guide_table[, 3:ncol(combined_guide_table)])
  plot_heatmap(combined_guide_filtered, "guide_rawCount")
  
  ## scale by mutiple T0, to make nonT0 samples are comparable across samples ##
  T0_scaled_guide_table <- combined_guide_table
  T0_scaled_guide_table[is.na(T0_scaled_guide_table)] <- 0 
  T0_scaled_guide_table[,3:ncol(T0_scaled_guide_table)] <- T0_scaled_guide_table[,3:ncol(T0_scaled_guide_table)] + 1 ## to avoid devide by 0
  for(control in unique(eXdesign$controls)){
    sub_df <- dplyr:: filter(eXdesign, controls == control)
    control_T0 <- sub_df %>%
      dplyr::filter(grepl("T0", replicates) & subjects == control) %>%
      dplyr::pull(replicates)
    treat_T0s <- sub_df %>%
      dplyr::filter(grepl("T0", replicates) & subjects != control) %>%
      dplyr::pull(replicates)
    for(treat_T0 in treat_T0s){
      subect <- sub_df %>% dplyr::filter(replicates == treat_T0) %>%
        dplyr::pull(subjects)
      treat_columns <- sub_df %>% dplyr::filter(subjects == subject & replicates != treat_T0) %>%
        dplyr::pull(replicates)
      T0_scaled_guide_table[, treat_T0] <- T0_scaled_guide_table[, treat_T0]/T0_scaled_guide_table[, control_T0]
      for(treat_col in treat_columns){
        T0_scaled_guide_table[, treat_col] <- T0_scaled_guide_table[, treat_col]/T0_scaled_guide_table[, treat_T0]
      }
    }
  }
  T0_scaled_guide_table[is.na(T0_scaled_guide_table)] <- 0
  write.table(T0_scaled_guide_table, file.path(hwd,"T0_scaled_guide_rawCount_table_for_MAGeCK.tab"), row.names=F, sep="\t", quote=F)
  
  
  dim(combined_scaled_guide_table)
  NAc <- apply(combined_scaled_guide_table, 1, function(x) sum(is.na(x)))
  NAs <- NAc == ncol(combined_scaled_guide_table)-2
  combined_scaled_guide_table <- combined_scaled_guide_table[!NAs,]
  ## "combined_guide_scaledCount_table.tab" can be used as input for DrugZ directly
  write.table(combined_scaled_guide_table, file.path(hwd,"combined_guide_scaledCount_table.tab"), row.names=F, sep="\t", quote=F)
  combined_scaled_guide_filtered <- na.omit(combined_scaled_guide_table[, 3:ncol(combined_scaled_guide_table)])
  plot_heatmap(combined_scaled_guide_filtered, "guide_scaledCount")
  
  
  dim(combined_normalized_guide_table)
  NAc <- apply(combined_normalized_guide_table, 1, function(x) sum(is.na(x)))
  NAs <- NAc == ncol(combined_normalized_guide_table)-2
  combined_normalized_guide_table <- combined_normalized_guide_table[!NAs,]
  write.table(combined_normalized_guide_table, file.path(hwd,"combined_guide_normalizedCount_table.tab"), row.names=F, sep="\t", quote=F)
  combined_normalized_guide_filtered <- na.omit(combined_normalized_guide_table[, 3:ncol(combined_normalized_guide_table)])
  plot_heatmap(combined_normalized_guide_filtered, "guide_normalizedCount")
  
  
  dim(combined_lfc_guide_table)  ## used for Bagel analysis
  NAc <- apply(combined_lfc_guide_table, 1, function(x) sum(is.na(x)))
  NAs <- NAc == ncol(combined_lfc_guide_table)-2
  combined_lfc_guide_table <- combined_lfc_guide_table[!NAs,] ## remove all_NA rows
  write.table(combined_lfc_guide_table, file.path(hwd,"combined_guide_lfc_table.tab"), row.names=F, sep="\t", quote=F)
  combined_lfc_guide_filtered <- na.omit(combined_lfc_guide_table[, 3:ncol(combined_lfc_guide_table)])
  plot_heatmap(combined_lfc_guide_filtered, "guide_lfcCount")
  
}

preprocess_bagel_dropout <- function(hwd) {
  setwd(hwd)
  lfc_df <- read.delim2("combined_guide_lfc_table.tab")
  lfc_df[is.na(lfc_df)] <- 0
  write.table(lfc_df, file.path(hwd,"BAGEL", "combined_guide_lfc_for_bagel.tab"), row.names=F, sep="\t", quote=F)
}

collect_BAGEL_results<- function(wd){
  setwd(wd)
  
  ## start processing bf data
  bagel_results_list <- list()
  
  
  bagel_files <- list.files(pattern="_bagel.bf")
  
  for (bagel_file in bagel_files){
    bf <- read.delim(bagel_file)
    subject <- gsub("_bagel.bf", "", bagel_file, fixed=TRUE)
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
  samples <- unlist(lapply(colnames(bf_table), function(x) gsub("_BF", "", x, fixed=TRUE)))
  colnames(bf_table) <- samples
  dim(bf_table)
  NAs <- apply(bf_table, 1, function(x) sum(is.na(x))>0)
  NA_table <- bf_table[NAs,]
  bf_table[is.na(bf_table)] <- 0
  
  
  pdf("BF_all_sample_correlation.pdf")
  col0 <- colorRampPalette(c("cyan4", "red3"))
  corrplot::corrplot(cor(bf_table), method="number", is.corr = F, col=col0(20), type = "upper", order = "hclust", cl.lim=c(0.65, 1), tl.col = "black", tl.srt = 45)
  dev.off()
  
  return(bf_table)
}

collect_MAGeCK_results <- function(wd){
  setwd(wd)
  MAGeCK_essentials <- list.files(pattern = "_mageck_essential.gene_summary.txt")
  samples <- unlist(lapply(MAGeCK_essentials, function(x) 
    gsub("_mageck_essential.gene_summary.txt", "", x, fixed=TRUE)))
  df_list <- lapply(MAGeCK_essentials, function(x){
    df <- read.delim(x, header = TRUE) 
    df <- df[, c(1, 5)]
    colnames(df) <- c("gene", "fdr")
    df
  })
  
  fdr_df <- merge(df_list[[1]], df_list[[2]], by = c("gene"))
  rn <- fdr_df[,1]
  fdr_df <- fdr_df[,2:ncol(fdr_df)]
  colnames(fdr_df) <- samples
  rownames(fdr_df) <- rn
  
  pdf("MAGeCK_all_sample_correlation.pdf")
  col0 <- colorRampPalette(c("cyan4", "red3"))
  corrplot::corrplot(cor(fdr_df), method="number", is.corr = F, col=col0(20), type = "upper", order = "hclust", cl.lim=c(0.65, 1), tl.col = "black", tl.srt = 45)
  dev.off()
  
  return(fdr_df)
}


qGI_analysis <- function(wd, hwd, runGSEA = FALSE){
  setwd(wd)
  
  eXdesign <- data.frame(read.delim("experimental_design_final.tab", header=TRUE, sep="\t", stringsAsFactors=F))
  
  rownames(eXdesign) <- eXdesign$X
  eXdesign <- eXdesign[,-1]
  
  eXdesign
  
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
  
  qGI_results <- plot_qGI_results(test_results = fiteb_results, pvalue_cutoff=0.05, gsea=runGSEA)
  
  save(qGI_results, file=file.path(qGI_dir, "results_list.Rdata"))
}

preprocess_mageck_terminator <- function(wd, samplef=c(2,4), exclude=NULL){
  setwd(wd)
  dir.create( "MAGeCK")
  guide_format <- "_guideRawCount.tab"
  terminator <- unlist(strsplit(wd, "/", fixed=T))[6]
  ## Load file paths and directory names
 
  files <- list.files(path = wd, recursive = TRUE, full.names = FALSE, pattern = guide_format)
  files
  
  samples <- unlist(lapply(files, function(x)split_sample(x,samplef[1], samplef[2])))
  
  if(!is.null(exclude)){
    files <- files[!samples %in% exclude]
    samples <- samples[!samples %in% exclude]
  }
  
  samples <- toupper(samples)
  data_list <- list()
  
  
  for(i in 1:length(files)){
    # i <- 1
    df_all <- read.table(files[i], header=F, sep="\t")
    colnames(df_all) <- c("Guide", "Sequence", samples[i])
    data_list[[samples[i]]] <- df_all
    
  }
  sapply(data_list, dim)
  count_df <- Reduce(function(x, y) merge(x, y, all=TRUE), data_list)
  
  ### update gene names in row.names
  old_gene_names <- count_df[,2]
  new_gene_names <- update_gene_names(name_map, old_gene_names)
  
  
  name_df <- data.frame(old_gene_names, new_gene_names)
  name_diff <- name_df[!name_df$old_gene_names %in% name_df$new_gene_names, ]
  write.table(name_diff, "updated_gene_list.tab", sep="\t")
 
  count_df[,2] <- new_gene_names
  colnames(count_df)[1:2] <- c("sgRNA", "gene")
  
  count_df[is.na(count_df)] <- 0
  
  
  write.table(count_df, file.path("MAGeCK", paste0(terminator, "_sgRNA_read_count.tab")), 
              row.names=F, col.names=T, sep="\t", quote = FALSE)
  
  
  for(i in 1:length(files)){
    print(samples[i])
    test_df <- data_list[[samples[i]]]
    test_df <- test_df[order(test_df[,3]),]
    #head(test_df)
    #print(tail(test_df))
    print(test_df[test_df[,2] %in% c("CPSF1","NUDT21", "FIP1", "WDR33"),])
  }
  
  raw_count_df <- count_df[,3:ncol(count_df)]
  plot_heatmap(raw_count_df, paste0(terminator, "_guide"))
  
  dsc <- sumstats_col(raw_count_df)
  
  write.table(dsc, "descriptive_stats_of_raw_count_guide.tab", row.names=T, col.names=NA, sep="\t")
  
  
  pdf(paste(terminator, "raw_count_Guide_boxplot.pdf", sep="_"), width=15, height=6)
  
  old.par <- par(mfrow=c(1,1),mar=c(10,4,2,2))
  
  boxplot(as.matrix(raw_count_df), outline=T, names=NA)
  staxlab(1,1:ncol(raw_count_df),colnames(raw_count_df),srt=45)
  
  dev.off()
  par(old.par)
  
  design_mat <- data.frame(Samples = colnames(raw_count_df),
                           baseline = 1,
                           SORTED = ifelse(grepl("_SORT", colnames(raw_count_df)), 1, 0))
  
  write.table(design_mat, file.path("MAGeCK", paste0(terminator, "_design_matrix.tab")), 
                                    row.names = FALSE, col.names = TRUE,
                                    sep = "\t", quote = FALSE)
  
  sgRNA <- TKO3[, c("GUIDE_ID", "SEQUENCE", "GENE")]
  write.table(sgRNA, file.path("MAGeCK", "sgRNA_library.tab"), col.names = FALSE, row.names = FALSE, 
              quote = FALSE, sep = "\t")
  negative_sgRNA <- TKO3[TKO3$GENE %in% control_gene, "GUIDE_ID"]
  writeLines(negative_sgRNA, file.path("MAGeCK", "negative_control_sgRNA.txt"))
  writeLines(control_gene, file.path("MAGeCK", "negative_control_gene.txt"))
  
  file.copy(GOBP_file, file.path("MAGeCK", "GOBP.gmt"))
  
}

preprocess_mageck_dropout <- function(wd){
  setwd(wd)
  dir.create("MAGeCK")
  
 
  sgRNA <- TKO3[, c("GUIDE_ID", "SEQUENCE", "GENE")]
  write.table(sgRNA, file.path("MAGeCK", "sgRNA_library.tab"), col.names = FALSE, row.names = FALSE, 
              quote = FALSE, sep = "\t")
  negative_sgRNA <- TKO3[TKO3$GENE %in% control_gene, "GUIDE_ID"]
  writeLines(negative_sgRNA, file.path(wd, "MAGeCK", "negative_control_sgRNA.txt"))
  writeLines(control_gene, file.path(wd, "MAGeCK", "negative_control_gene.txt"))
  
  file.copy(GOBP_file, file.path(wd, "MAGeCK", "GOBP.gmt"), overwrite = TRUE)
  
  Guide_raw_count_file <- file.path(wd, "final_analysis/combined_guide_rawCount_table.tab")
  file.copy(Guide_raw_count_file, file.path(wd, "MAGeCK", "results_sgRNA_read_count.tab"),, overwrite = TRUE)
  
  Guide_T0_scaled_count_file <- file.path(wd, "final_analysis/T0_scaled_guide_rawCount_table_for_MAGeCK.tab")
  file.copy(Guide_T0_scaled_count_file, file.path(wd, "MAGeCK", "T0_scaled_sgRNA_read_count.tab"), overwrite = TRUE)
  
}

preprocess_chronos_dropout <- function(hwd, treat = "AGS_EBV", control = "AGS", day = 21){
  setwd(hwd)
  dir.create("chronos")
  
  
  sgRNA <- TKO3[, c("GUIDE_ID", "SEQUENCE", "GENE")]
  colnames(sgRNA) <- c("sgrna", "sequence", "gene")
  write.csv(sgRNA, file.path("chronos", "sgRNA_library.csv"), 
              row.names = FALSE, quote = FALSE)
  negative_sgRNA <- TKO3[TKO3$GENE %in% control_gene, "GUIDE_ID"]
  writeLines(negative_sgRNA, file.path("chronos", "negative_control_sgRNA.txt"))
  writeLines(control_gene, file.path("chronos", "negative_control_gene.txt"))
  
  Guide_raw_count_df <- read.delim2("combined_guide_rawCount_table.tab")
  Guide_raw_count_mat <- Guide_raw_count_df[, 3:ncol(Guide_raw_count_df)]
  rownames(Guide_raw_count_mat) <- Guide_raw_count_df[,1]
  Guide_raw_count_mat <- t(Guide_raw_count_mat)
  write.table(Guide_raw_count_mat, file.path("chronos", "sgRNA_raw_read_count.tab"),
              sep = "\t", col.names = NA, quote = FALSE)
  
  sequence_ID <- rownames(Guide_raw_count_mat)
  cell_line_name <- rep(control, length(sequence_ID))
  cell_line_name[grepl(treat, sequence_ID)] <- treat
  cell_line_name[grepl("T0", sequence_ID)] <- "pDNA"
  pDNA_batch <- rep("batch1", length(sequence_ID))
  pDNA_batch[grepl(treat, sequence_ID)] <- "batch2"
  
  sequence_map <- data.frame(sequence_ID, cell_line_name, days = day, pDNA_batch)
  write.csv(sequence_map, file.path("chronos/sample_map.csv"), row.names = FALSE)
  
  
}


collect_sort_unsort_gene <- function(wd, samplef=c(2,4), groupf=c(2,3), exlude=NULL){
  
  format <- "_geneRawCount.tab"
  terminator <- unlist(strsplit(wd, "/", fixed=T))[6]
  
  ## Load file paths and directory names
  files <- list.files(path = getwd(), recursive = TRUE, full.names = FALSE, pattern = format)
  
  length(files)
  files
  
  samples <- unlist(lapply(files, function(x)split_sample(x,samplef[1],samplef[2])))
  groups <- unlist(lapply(files, function(x)split_group(x,groupf[1],groupf[2])))
  
  if(!is.null(exclude)){
    groups <- groups[!samples %in% exclude]
    files <- files[!samples %in% exclude]
    samples <- samples[!samples %in% exclude]
  }
  
  samples <- toupper(samples)
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
      write.table(name_diff, "updated_gene_list.tab", sep="\t")
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
    test_df <- data_list[[samples[i]]]
    test_df <- test_df[order(test_df[,3]),]
    #head(test_df)
    #print(tail(test_df))
    print(test_df[c("CPSF1","NUDT21", "FIP1", "WDR33"),])
  }
  
  guides <- paste0("guide", c(1:4))
  unsort_count_df <- NULL
  sort_count_df <- NULL 
  
  sort_samples <- samples[grep("_SORT", samples)]
  unsort_samples <- samples[grep("_UNSORT", samples)]
  
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
  
  #sort_count_df <- sort_count_df[!rownames(sort_count_df) %in% c("NUDT21"),]
  #unsort_count_df <- unsort_count_df[!rownames(unsort_count_df) %in% c("NUDT21"),]
  raw_count_df <- cbind(sort_count_df, unsort_count_df)
  dim(raw_count_df)
  
  plot_heatmap(raw_count_df, terminator)
  
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
  staxlab(1,1:ncol(raw_count_df),colnames(raw_count_df),srt=45)
  
  dev.off()
  par(old.par)
  
  
  sort_count_df <- as.matrix(sort_count_df)
  unsort_count_df <- as.matrix(unsort_count_df)
  head(sort_count_df)
  dim(sort_count_df)
  summary(sort_count_df)
  
  return(list(sorted = sort_count_df, unsorted = unsort_count_df, subject = terminator))
}

# DEseq2 analysis 

sort_unsort_deseq2 <- function(sort_count_df, unsort_count_df, terminator) {
  
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
  staxlab(1,1:ncol(normalized_count_df),colnames(normalized_count_df),srt=45)
  
  dev.off()
  par(old.par)
  
  sort_stats <- sumstats_row(normalized_count_df[,grep("_SORT", colnames(normalized_count_df))])
  unsort_stats <- sumstats_row(normalized_count_df[,grep("_UNSORT", colnames(normalized_count_df))])
  
  sort_diff <- normalized_count_df[,grep("_SORT", colnames(normalized_count_df))] - unsort_stats$Mean
  sort_zscore <- round(t(t(sort_diff)/unsort_stats$SD), 2)
  sig_guide_count <- apply(sort_zscore, 1, function(x)length(x[x>2]))  ## zscore > 2 
  
  ## for each guide, find the samples that are significant
  sig_sample_count <- list()
  sig_guide_sample_count <- list()
  Guides <- paste0("guide", seq_len(4))
  for(guide in Guides){
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
  sort_samples <- samples[grep("_SORT", samples)]
  sort_samples <- unique(gsub("_guide[1-4]", "", sort_samples))
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
