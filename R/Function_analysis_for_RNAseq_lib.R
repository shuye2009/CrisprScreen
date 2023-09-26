## collection of functions for plotting functional enrichment of DESeq2 results
## or simple gene lists

library(dplyr)
library(pheatmap)
library(PCAtools)
library(cowplot)
library(ggplotify)
library(DESeq2)
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(VennDiagram)
library(limma)
library(DOSE)
library(GOSemSim)
library(viridis)

ap_cutoff <- 0.05
sim_cutoff <- 0.7
SemList <- list()
SemList[["BP"]] <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
SemList[["CC"]] <- godata('org.Hs.eg.db', ont="CC", computeIC=TRUE)
SemList[["MF"]] <- godata('org.Hs.eg.db', ont="MF", computeIC=TRUE)

plot_heatmap <- function(countData, dataName=NULL, h=10, w=10){

  while(!is.null(dev.list())){dev.off()}

  ## plot correlation matrix
  corMat <- cor(countData)
  dn <- TRUE
  if(ncol(corMat) > 25){
     dn <- FALSE
  }
  pdf(paste(dataName, "_correlation_heatmap.pdf", sep=""), height=h, width=w)
  plotMDS(corMat)
  grid.newpage()
  grid.draw(pheatmap::pheatmap(corMat, silent=T, display_numbers=dn))

  dev.off()

  pdf(paste(dataName, "_rawValue_heatmap.pdf", sep=""), height=h, width=w)

  topN <- max(round(nrow(countData)/2), 1000)

  plotMDS(countData, top=topN)
  if(ncol(countData) < 7){
    chart.Correlation(countData)
  }
  if(nrow(countData) > 10000){
    grid.newpage()
    grid.draw(pheatmap(countData[sample.int(nrow(countData), 10000, replace=FALSE),], silent=T, display_numbers=FALSE))
    message("the data is sampled to a limit of 10000 rows")
  }else{
     grid.newpage()
     grid.draw(pheatmap::pheatmap(countData, silent=T, display_numbers=FALSE))
  }

  dev.off()

}


plot_PCA <- function(countData, metadata, featureForColor, colorKey, featureForShape, shapeKey){
  prc <- pca(countData, metadata=metadata)

  pbiplot <- biplot(prc, x='PC1', y="PC2",
                    colby = featureForColor, colkey = colorKey,
                    shape = featureForShape, shapekey = shapeKey,
                    hline = 0, vline = c(0),
                    legendPosition = 'none',
                    gridlines.major = FALSE, gridlines.minor = FALSE,
                    pointSize = 5, axisLabSize = 12,
                    drawConnectors = TRUE,
                    title = 'PCA bi-plot',
                    subtitle = 'PC1 versus PC2',
                    returnPlot = FALSE)

  pscree <- screeplot(prc,
                      gridlines.major = FALSE, gridlines.minor = FALSE,
                      axisLabSize = 10,
                      returnPlot = FALSE)
  ppairs <- pairsplot(prc,
                      components = getComponents(prc, c(1:5)),
                      triangle = TRUE, trianglelabSize = 12,
                      hline = 0, vline = 0,
                      pointSize = 2,
                      gridlines.major = FALSE, gridlines.minor = FALSE,
                      colby = featureForColor, colkey = colorKey,
                      shape = featureForShape, shapekey = shapeKey,
                      title = '', titleLabSize=16,
                      plotaxes = FALSE,
                      margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
                      returnPlot = FALSE)

  ploadings <- plotloadings(prc,
                            rangeRetain = 0.002,
                            labSize = 2.0,
                            title = 'Loadings plot',
                            caption = 'Top 0.2% variables',
                            shape = 18,
                            shapeSizeRange=c(2,2),
                            col = c('limegreen', 'black', 'red3'),
                            drawConnectors = TRUE,
                            legendPosition = 'none',
                            axisLabSize = 10,
                            returnPlot = FALSE)

  # only works for continuous variables
  if(0){
    peigencor <- eigencorplot(prc,
                              components = getComponents(prc, 1:5),
                              metavars = c('conditions','covar'),
                              #col = c('royalblue', '', 'gold', 'forestgreen', 'darkgreen'),
                              cexCorval = 0.6,
                              fontCorval = 2,
                              posLab = 'all',
                              rotLabX = 45,
                              scale = TRUE,
                              main = "PC group correlates",
                              cexMain = 1.5,
                              plotRsquared = FALSE,
                              corFUN = 'pearson',
                              corUSE = 'pairwise.complete.obs',
                              signifSymbols = c('****', '***', '**', '*', ''),
                              signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                              returnPlot = FALSE)
  }

  ## plot multipanels

  top_row <- plot_grid(pscree, pbiplot,
                       ncol = 2,
                       labels = c('A', 'B'),
                       label_fontfamily = 'serif',
                       label_fontface = 'bold',
                       label_size = 22,
                       align = 'h',
                       rel_widths = c(1.5, 1.5))

  bottom_row <- plot_grid(ploadings,
                          ppairs,
                          ncol = 2,
                          labels = c('C', 'D Pairs plot'),
                          label_fontfamily = 'serif',
                          label_fontface = 'bold',
                          label_size = 22,
                          align = 'h',
                          rel_widths = c(1.5, 1.5))

  pdf("PCA results plot.pdf", height=10, width = 10)
  print(plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1.0, 1.0)))
  dev.off()
}



runDESeq2 <- function(e, metadata, contrast) {
  #e <- subMats[["N"]]
  #metadata <- subMetas[["N"]]
  #contrast <- contrasts[["N"]]

  dds <- DESeqDataSetFromMatrix(round(e), metadata, ~ CONDITION)
  dds <- DESeq(dds,quiet=TRUE)
  print(resultsNames(dds))
  res <- results(dds, contrast=contrast)
  write.table(res, paste(paste(contrast,collapse="_"), "deseq2_results.tab", sep="_"), col.names=NA, sep="\t", quote=F)
  return(res)
}

processDESeq2 <- function(res, contrast){
  beta_cutoff <- 0
  beta <- res$stat
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  padj[is.na(padj)] <- 1
  colors <- rep("black", length(beta))
  colors[beta < -beta_cutoff & padj < 0.05] <- "blue"
  colors[beta > beta_cutoff & padj < 0.05] <- "red"
  comparison <- paste(contrast[2:3], collapse="_vs_")
  print(plot(beta, -log(padj), pch=8, col="white", main=comparison))
  print(text(beta, -log(padj), row.names(res), col=colors))
  ups <- row.names(res)[beta > beta_cutoff & padj < 0.05]
  downs <- row.names(res)[beta < -beta_cutoff & padj < 0.05]
  writeLines(ups, paste(paste(contrast,collapse="_"), "_deseq2_up_genes.txt", sep=""))
  writeLines(downs, paste(paste(contrast,collapse="_"), "_deseq2_down_genes.txt", sep=""))

  return(list(comparison=comparison, ups=ups, downs=downs, result=res))
}


plot_DESEQ2 <- function(subMatList, subMetaList, contrastList){
  pdf("DESeq2_volcano_plot.pdf", height=10, width=10)
  #oldpar <- par(mfrow=c(2,2))
  treats <- names(subMatList)
  res_collect <- list()

  for(treat in treats){
    subMat <- subMatList[[treat]]
    subMeta <- subMetaList[[treat]]
    contrast <- contrastList[[treat]]
    res <- runDESeq2(subMat, subMeta, contrast)

    res_collect[[treat]] <- processDESeq2(res, contrast)  ## plot volcano plot
  }
  #par(oldpar)
  dev.off()

  return(res_collect)
}

run_enrichGO <- function(dge, adjp_cutoff=ap_cutoff){
  #treat <- "N"
  #dge <- res_collect[[treat]]
  print("run_enrichGO")

  upgenes <- dge$ups
  downgenes <-dge$downs
  dge_table <- dge$result

  comparison <- dge$comparison

  print(comparison)

  if(length(upgenes) > 5){
    y_up1 <- enrichGO(upgenes, 'org.Hs.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff=0.05)
    if(!is.null(y_up1)){
      print(dim(y_up1@result))
      y_up <- simplify(y_up1, cutoff=sim_cutoff, by="p.adjust", select_fun=min, measure="Wang", semData=SemList[[ont]])
      print(dim(y_up@result))

      result <- na.omit(y_up@result)
      result <- result[result$p.adjust < adjp_cutoff, ]
      print(dim(result))

      if(dim(result)[1] > 0){
        write.table(result, paste("GO_enrichment_", comparison, "_up_genes.tab", sep=""), sep="\t", quote=F)
      }else{
        print(paste(comparison, "No enrichment for up genes"))
      }
      if(dim(result)[1] > 2){
        oldpar <- par(mfrow=c(2,3))
        pdf(paste("GO_enrichment_", comparison, "_up_genes.pdf", sep=""), height=10, width=10)

        #plotGOgraph(y_up)
        lfc <- dge_table[upgenes, "stat"]
        names(lfc) <- upgenes
        print(cnetplot(y_up, colorEdge=T, foldChange=lfc))
        print(goplot(y_up))
        print(emapplot(y_up))
        print(dotplot(y_up))
        print(upsetplot(y_up))
        print(barplot(y_up))
        print(heatplot(y_up))
        #ggtable(y_up)

        dev.off()
        par(oldpar)
      }
    }
  }

  if(length(downgenes) > 5){
    y_down1 <- enrichGO(downgenes, 'org.Hs.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff=0.05)
    if(!is.null(y_down1)){
      print(dim(y_down1@result))
      y_down <- simplify(y_down1, cutoff=sim_cutoff, by="p.adjust", select_fun=min, measure="Wang", semData=SemList[[ont]])
      print(dim(y_down@result))

      result <- na.omit(y_down@result)
      result <- result[result$p.adjust < adjp_cutoff, ]
      print(dim(result))


      if(dim(result)[1] > 0){
        write.table(result, paste("GO_enrichment_", comparison, "_down_genes.tab", sep=""), sep="\t", quote=F)
      }else{
        print(paste(comparison, "No enrichment for down genes"))
      }
      if(dim(result)[1] > 2){
        oldpar <- par(mfrow=c(2,3))
        pdf(paste("GO_enrichment_", comparison, "_down_genes.pdf", sep=""), height=10, width=10)

        #plotGOgraph(y)
        lfc <- dge_table[downgenes, "stat"]
        names(lfc) <- downgenes
        print(cnetplot(y_down, colorEdge=T, foldChange=lfc))
        print(goplot(y_down))
        print(emapplot(y_down))
        print(dotplot(y_down))
        print(upsetplot(y_down))
        print(barplot(y_down))
        print(heatplot(y_down))
        #ggtable(y)

        dev.off()
        par(oldpar)
      }
    }
  }
}

run_enrichGO_simpleList <- function(geneList, ont, dataName, adjp_cutoff=ap_cutoff){

  print(paste(dataName, "run_enrichGO_simpleList", ont))

  if(length(geneList) > 2){
    y_res1 <- enrichGO(geneList, 'org.Hs.eg.db', keyType='SYMBOL', ont=ont, pvalueCutoff=0.05)
    if(!is.null(y_res1)){
      cat("before ")
      print(dim(y_res1@result))
      y_res <- simplify(y_res1, cutoff=sim_cutoff, by="p.adjust", select_fun=min, measure="Wang", semData=SemList[[ont]])
      cat("after ")
      print(dim(y_res@result))

      result <- na.omit(y_res@result)

      coln <- colnames(result)
      coln[coln=="ID"] <- "termID" ## fix the SYLK error when open with Excel, which caused by 'ID" as the first entry
      colnames(result) <- coln
      write.table(result, paste("GO_enrichment_", dataName, "_genes.tab", sep=""), row.name=F, sep="\t", quote=F)
      result <- result[result$p.adjust < adjp_cutoff, ]
      cat(paste("padj <", adjp_cutoff, " "))
      print(dim(result))

      if(dim(result)[1] > 0){
       # write.table(result, paste("GO_enrichment_", dataName, "_genes_adj.tab", sep=""), row.name=F, sep="\t", quote=F)
      }else{
        print(paste(dataName, "No GO enrichment found after p value adjust!"))
      }
      if(dim(result)[1] > 5){

          oldpar <- par(mfrow=c(2,3))
          pdf(paste("GO_enrichment_", dataName, "_genes.pdf", sep=""), height=10, width=15)


          (print(cnetplot(y_res, colorEdge=T)))
          (print(goplot(y_res)))
          #(print(emapplot(y_res)))
          (print(dotplot(y_res)))
          (print(upsetplot(y_res)))
          (print(barplot(y_res)))
          (print(heatplot(y_res)))


          dev.off()
          par(oldpar)

      }
    }
  }
}

#pathways <- gmtPathways("C:/RSYNC/GSEA/BaderHuman_feb2020/Human_GO_bp_no_GO_iea_symbol.gmt")


run_gseGO <- function(dge, adjp_cutoff=ap_cutoff){
  #treat <- "N"
  #dge <- res_collect[[treat]]
  print("run_gseGO")
  dge_table <- dge$result
  ranked_table <- dge_table[order(dge_table$stat, decreasing=T),]
  ranked_genes <- ranked_table$stat
  names(ranked_genes) <- row.names(ranked_table)

  comparison <- dge$comparison
  print(comparison)

  y1 <- gseGO(geneList=ranked_genes,
             OrgDb=org.Hs.eg.db,
             keyType='SYMBOL',
             ont="BP",
             pvalueCutoff=adjp_cutoff,
             by='fgsea')


  if(!is.null(y1)){
    print(dim(y1@result))
    y <- simplify(y1, cutoff=sim_cutoff, by="p.adjust", select_fun=min, measure="Wang", semData=SemList[[ont]])
    print(dim(y@result))

    head(y@result)
    if(dim(y@result)[1] > 0){
      pdf(paste(comparison, "GO_BP_gsea_analysis_plot.pdf", sep="_"), height=10, width=15)
      for(i in 1:dim(y@result)[1]){
        print(gseaplot2(y, geneSetID=i, title=y@result$Description[i]))
      }
      dev.off()
      write.table(y@result, paste("GO_gsea_analysis_", comparison, ".tab", sep=""), sep="\t", quote=F)
    }else{
      print(paste(comparison, "No gse enrichment for genes"))
    }
    if(dim(y@result)[1] > 1){
      pdf(paste("GO_gsea_analysis_", comparison, ".pdf", sep=""), height=10, width=15)
      print(ridgeplot(y))
      print(cnetplot(y, colorEdge=T, foldChange=ranked_genes))
      print(emapplot(y))
      print(dotplot(y))
      print(upsetplot(y))
      print(heatplot(y))
      #print(barplot(y))
      dev.off()
    }
  }
}

run_gseGO_simple <- function(dge_table, dataName=NULL, adjp_cutoff=ap_cutoff){
  #dge_table <- res
  #adjp_cutoff <- 0.05
  #dataName <- "S_vs_GFP"
  print(paste(dataName, "run_gseGO_Simple"))

  ranked_table <- dge_table[order(dge_table$stat, decreasing=T),]
  ranked_genes <- ranked_table$stat
  names(ranked_genes) <- row.names(ranked_table)

  y1 <- gseGO(geneList=ranked_genes,
             OrgDb=org.Hs.eg.db,
             keyType='SYMBOL',
             ont="BP",
             pvalueCutoff=adjp_cutoff,
             by='fgsea')

  if(!is.null(y1)){
    print(dim(y1@result))
    y <- simplify(y1, cutoff=sim_cutoff, by="p.adjust", select_fun=min, measure="Wang", semData=SemList[[ont]])
    print(dim(y@result))

    if(dim(y@result)[1] > 0){
      pdf(paste(dataName, "GO_BP_gsea_analysis_plot.pdf", sep="_"), height=10, width=15)
      for(i in 1:dim(y@result)[1]){
        print(gseaplot2(y, geneSetID=i, title=y@result$Description[i]))
      }
      dev.off()
      write.table(y@result, paste(dataName, "GO_BP_gsea_analysis.tab", sep="_"), row.names=F, sep="\t", quote=F)
    }else{
      print(paste(dataName, "No gse enrichment"))
    }
    if(dim(y@result)[1] > 1){
      pdf(paste(dataName, "GO_BP_gsea_analysis.pdf", sep="_"), height=10, width=15)
      print(ridgeplot(y))
      print(cnetplot(y, colorEdge=T, foldChange=ranked_genes))
      print(emapplot(y))
      print(dotplot(y))
      print(upsetplot(y))
      print(heatplot(y))
      #print(barplot(y))
      dev.off()
    }
  }
}

gsea_NES_plot <- function(gseaRes, pval=0.05, topn=30, label_format=75) {

   fgRes <- gseaRes@result %>%
      as.data.frame() %>%
      dplyr::filter(p.adjust < !!pval) %>%
      top_n(-topn, p.adjust) %>%
      arrange(desc(NES)) %>%
      mutate(Description = stringr::str_wrap(as.character(Description), label_format, 0.0) )

   #upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))(nrow(fgRes[fgRes$NES>0,]))
   #downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))(nrow(fgRes[fgRes$NES<0,]))
   #colos = c(upcols, downcols)


   g = ggplot(fgRes, aes(reorder(Description, NES), NES)) +
      geom_point(aes(color=NES, size=-log(p.adjust))) +
      scale_color_gradient2(
         low = "blue4",
         mid = "white",
         high = "red4",
         midpoint = 0
      )+
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="GSEA") +
      theme_minimal() +
      theme(axis.title = element_text(size = 16)) +
      theme(axis.text = element_text(size = 12))

   output = list("Results" = fgRes, "Plot" = g)
   return(output)
}

build_term2gene_table <- function(collection=NULL, subcollection=NULL){
   library(msigdb)
   library(ExperimentHub)
   library(GSEABase)

   #collection <- NULL
   #subcollection <- "GO:BP"

   if(is.null(collection) & is.null(subcollection)) stop("collection or subcollection has to be provided!")

   msigdb.hs <- getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
   msigdb.hs <- appendKEGG(msigdb.hs)
   listCollections(msigdb.hs)
   SC <- subsetCollection(msigdb.hs, collection=collection, subcollection=subcollection)
   listSubCollections(SC)
   msigdb_ids <- geneIds(SC)
   unique(sapply(names(msigdb_ids), function(x)unlist(strsplit(x, split="_"))[1]))

   listOfdf <- lapply(names(msigdb_ids), function(aterm){
      y <- msigdb_ids[[aterm]]
      adf <- as.data.frame(cbind(term=rep(aterm, length(y)), gene=y))
   })

   df <- bind_rows(listOfdf)

   return(df)
}

clusterProfiler_GSEA_plot <- function(gseaRes, geneList, labFormat=75, pval=0.05, subsetDesc=NULL){
   if(!is.null(subsetDesc)){
      gseaRes@result <- gseaRes@result[subsetDesc, ]
   }

   try(print(enrichplot::ridgeplot(gseaRes, showCategory=30, label_format=labFormat)))
   try(print(gsea_NES_plot(gseaRes, pval=pval, topn=30, label_format=labFormat)$Plot))
   try(print(cnetplot(gseaRes, colorEdge=T, foldChange=geneList)))
   #try(print(emapplot(gseaRes)))
   try(print(dotplot(gseaRes, showCategory=30, label_format=labFormat)))
   try(print(upsetplot(gseaRes)))
   try(print(heatplot(gseaRes, showCategory=30, foldChange=geneList)))
}

run_gseGO_simpleList <- function(gene_list, dataName=NULL, adjp_cutoff=0.05, simplify=TRUE, ont="BP", GO_file=NULL, collection=NULL, subcollection=NULL){
  set.seed(54321)
  print(paste(dataName, "run_gseGO_simpleList", ont))

   if(any( duplicated(names(gene_list)))  ) {
      warning("Duplicates in gene names")
      gene_list = gene_list[!duplicated(names(gene_list))]
   }
   if (!all(order(gene_list, decreasing = TRUE) == 1:length(gene_list))){
      warning("Gene list not sorted")
      gene_list = sort(gene_list, decreasing = TRUE)
   }

   y <- NULL
  if(ont %in% c("BP", "MF", "CC")){
     #d <- godata('org.Hs.eg.db', ont=ont, computeIC=FALSE)
     y <- clusterProfiler::gseGO(geneList=gene_list,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont=ont,
                exponent = 1,
                minGSSize = 10,
                maxGSSize = 500,
                eps = 0, # 1e-50,
                pvalueCutoff=adjp_cutoff,
                seed=TRUE,
                by='fgsea')
     y@readable <- FALSE
  }else if(!is.null(GO_file)){
     pGO <- read.gmt(GO_file)

     y <- clusterProfiler::GSEA(
        geneList=gene_list,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 0, #1e-50,
        pvalueCutoff = adjp_cutoff,
        pAdjustMethod = "BH",
        TERM2GENE = pGO,
        by = "fgsea",
        seed = TRUE
     )

  }else if(!(is.null(collection) && is.null(subcollection))){
     pGO <- build_term2gene_table(collection, subcollection)

     y <- clusterProfiler::GSEA(
        geneList=gene_list,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 0, #1e-50,
        pvalueCutoff = adjp_cutoff,
        pAdjustMethod = "BH",
        TERM2GENE = pGO,
        by = "fgsea",
        seed = TRUE
     )

  }else{
     stop("Neither valid ontology nor GO file is provided!")
  }


  if(!is.null(y)){
    print(nrow(y@result))

     if(ont %in% c("BP", "CC", "MF") && simplify){
        write.table(df, paste(dataName, ont, "gsea_analysis_unsimplified.tsv", sep="_"), row.names=F, sep="\t", quote=F)
        y <- clusterProfiler::simplify(y, cutoff=sim_cutoff, by="setSize", select_fun=max, measure="Wang", semData=SemList[[ont]])
        print(paste("simplified", nrow(y@result)))
     }

     df <- y@result ## reset df
     coln <- colnames(df)
     coln[coln=="ID"] <- "termID" ## fix the SYLK error when open with Excel, which caused by 'ID" as the first entry
     colnames(df) <- coln

    if(nrow(y@result) > 0){
      pdf(paste(dataName, ont, "gsea_analysis_plot.pdf", sep="_"), height=10, width=15)
      for(i in 1:nrow(y@result)){
        print(gseaplot2(y, geneSetID=i, title=y@result$Description[i]))
      }
      dev.off()
      write.table(df, paste(dataName, ont, "gsea_analysis.tab", sep="_"), row.names=F, sep="\t", quote=F)
    }else{
      print(paste(dataName, "No gsea enrichment"))
    }

    if(nrow(y@result) > 1){
      pdf(paste(dataName, ont, "gsea_analysis.pdf", sep="_"), height=10, width=15)
      clusterProfiler_GSEA_plot(y, gene_list, labFormat=100, pval=adjp_cutoff, subsetDesc=NULL)
      #print(barplot(y))
      dev.off()
    }
  }
  return(y)
}

run_enrichDO <- function(dge, adjp_cutoff=ap_cutoff){
  #treat <- "N"
  #dge <- res_collect[[treat]]
  print("run_enrichDO")

  upgenes <- bitr(dge$ups, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID
  downgenes <-bitr(dge$downs, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID
  dge_table <- dge$result

  comparison <- dge$comparison
  print(comparison)



  if(length(upgenes) > 5){

    y_up <- enrichDO(upgenes, ont="DO", pvalueCutoff=0.05, readable=TRUE)

    if(!is.null(y_up)){
      print(dim(y_up@result))

      result <- na.omit(y_up@result)
      result <- result[result$p.adjust < adjp_cutoff, ]
      print(dim(result))


      if(dim(result)[1] > 0){
        write.table(result, paste("Disease_gsea_analysis_", comparison, "_up_gene.tab", sep=""), sep="\t", quote=F)
      }else{
        print("No disease enrichment for up genes")
      }
      if(dim(result)[1] > 1){
        oldpar <- par(mfrow=c(2,3))
        pdf(paste("Disease_enrichment_", comparison, "_up_genes.pdf", sep=""), height=10, width=10)

        #plotGOgraph(y_up)
        lfc <- dge_table[dge$ups, "stat"]
        names(lfc) <- dge$ups
        print(cnetplot(y_up, colorEdge=T, foldChange=lfc))
        print(heatplot(y_up))
        print(emapplot(y_up))
        print(dotplot(y_up))
        print(upsetplot(y_up))
        print(barplot(y_up))
        #ggtable(y_up)

        dev.off()
        par(oldpar)
      }
    }
  }

  if(length(downgenes) > 5){
    y_down <- enrichDO(downgenes, ont="DO", pvalueCutoff=0.05, readable=TRUE)

    if(!is.null(y_down)){
      print(dim(y_down@result))
      result <- na.omit(y_down@result)
      result <- result[result$p.adjust < adjp_cutoff, ]
      print(dim(result))

      if(dim(result)[1] > 0){
        write.table(result, paste("Disease_gsea_analysis_", comparison, "_down_gene.tab", sep=""), sep="\t", quote=F)
      }else{
        print("No disease enrichment for down genes")
      }
      if(dim(result)[1] > 1){
        oldpar <- par(mfrow=c(2,3))
        pdf(paste("Disease_enrichment_", comparison, "_down_genes.pdf", sep=""), height=10, width=10)

        #plotGOgraph(y)
        lfc <- dge_table[dge$downs, "stat"]
        names(lfc) <- dge$downs
        print(cnetplot(y_down, colorEdge=T, foldChange=lfc))
        print(heatplot(y_down))
        print(emapplot(y_down))
        print(dotplot(y_down))
        print(upsetplot(y_down))
        print(barplot(y_down))
        #ggtable(y)

        dev.off()
        par(oldpar)
      }
    }
  }
}

run_enrichDO_simpleList <- function(geneList, dataName, adjp_cutoff=ap_cutoff){

  print(paste(dataName, "run_enrichDO_simpleList"))
  if(length(geneList) > 0){
    idList <- bitr(geneList, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID

    y_res <- enrichDO(idList, ont="DO", pvalueCutoff=0.05, readable=TRUE)

    if(!is.null(y_res)){
      print(dim(y_res@result))


      result <- na.omit(y_res@result)
      result <- result[result$p.adjust < adjp_cutoff, ]
      print(dim(result))


      if(dim(result)[1] > 0){
        write.table(result, paste("Disease_gsea_analysis_", dataName, ".tab", sep=""), sep="\t", quote=F)
      }else{
        print(paste(dataName, "No disease enrichment found!"))
      }
      if(dim(result)[1] > 1){
        oldpar <- par(mfrow=c(2,3))
        pdf(paste("Disease_enrichment_", dataName, "_genes.pdf", sep=""), height=10, width=10)

        print(cnetplot(y_res, colorEdge=T))
        print(emapplot(y_res))
        print(dotplot(y_res))
        print(upsetplot(y_res))
        print(barplot(y_res))
        print(heatplot(y_res))
        #ggtable(y_res)

        dev.off()
        par(oldpar)
      }
    }
  }
}


run_gseDO <- function(dge, adjp_cutoff=ap_cutoff){
  #treat <- "N"
  #dge <- res_collect[[treat]]
  print("run_gseDO")
  dge_table <- dge$result
  ranked_table <- dge_table[order(dge_table$stat, decreasing=T),]
  ranked_genes <- ranked_table$stat
  names(ranked_genes) <- bitr(row.names(ranked_table), "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID
  #names(ranked_genes) <- row.names(ranked_table)
  comparison <- dge$comparison
  print(comparison)

  gdo <- gseDO(ranked_genes,
               minGSSize     = 10,
               pvalueCutoff  = adjp_cutoff,
               pAdjustMethod = "BH",
               verbose       = TRUE)


  if(!is.null(gdo)){
    gdo <- setReadable(gdo, 'org.Hs.eg.db')
    print(dim(gdo@result))
    if(dim(gdo@result)[1] > 0){
      write.table(gdo@result, paste("Disease_gsea_analysis_", comparison, ".tab", sep=""), sep="\t", quote=F)
    }else{
      print(paste(comparison, "No disease enrichment for genes"))
    }

    if(dim(gdo@result)[1] > 1){
      pdf(paste("Disease_gsea_analysis_", comparison, ".pdf", sep=""), height=10, width=10)
      print(ridgeplot(gdo))
      print(cnetplot(gdo, colorEdge=T, foldChange=ranked_genes))
      print(emapplot(gdo))
      print(dotplot(gdo))
      print(upsetplot(gdo))
      print(heatplot(gdo))
      dev.off()
    }
  }
}

plot_enrichment <- function(dge_results){
  regulated_list <- list()

  for(treat in names(dge_results)){

    dge <- dge_results[[treat]]


    upgenes <- bitr(dge$ups, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID
    downgenes <-bitr(dge$downs, "SYMBOL", "ENTREZID", 'org.Hs.eg.db')$ENTREZID
    comparison <- dge$comparison

    up <- paste(comparison, "up", sep="_")
    down <- paste(comparison, "down", sep="_")

    regulated_list[[up]] <- upgenes
    regulated_list[[down]] <- downgenes

    run_enrichGO(dge)
    run_gseGO(dge)
    run_enrichDO(dge)
    run_gseDO(dge)
  }

  return(regulated_list)
}


## input: a named list of gene names
plot_clusterProfile_DO <- function(regulated_list, dataName){

  ck_DO <- compareCluster(geneCluster = regulated_list, fun = "enrichDO")

  if(dim(as.data.frame(ck_DO))[1] > 0){
    write.table(as.data.frame(ck_DO), paste("Disease_enrichment_profile_",  dataName, ".tab", sep=""), sep="\t", quote=F, row.names=F)

    pdf(paste("Disease_enrichment_profile_", dataName, ".pdf", sep=""), height=10, width=15)
    print(dotplot(ck_DO))
    dev.off()
  }else{
    print("No disease enrichment")
  }
}

plot_clusterProfile_GO <- function(regulated_list, dataName){

  ck_GO1 <- compareCluster(geneCluster = regulated_list, OrgDb=org.Hs.eg.db, fun = "enrichGO", ont="BP", pvalueCutoff=0.05)
  if(dim(as.data.frame(ck_GO1))[1] > 0){
    ck_GO <- simplify(ck_GO1, cutoff=sim_cutoff, by="p.adjust", select_fun=min, measure="Wang", semData=SemList[[ont]])
  }

  if(dim(as.data.frame(ck_GO))[1] > 0){

    write.table(as.data.frame(ck_GO), paste("GO_enrichment_profile_",  dataName, ".tab", sep=""), sep="\t", quote=F, row.names=F)

    pdf(paste("GO_enrichment_profile_", dataName, ".pdf", sep=""), height=15, width=20)

    p <- dotplot(ck_GO, showCategory = 20)
    grid.draw(p);
    grid.newpage();
    dev.off()
  }else{
    print("No GO BP enrichment")
  }

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

plotVenn <- function(alist, fileName, title, pos="outer"){
  colors <- terrain.colors(5)
  venn.plot <- venn.diagram(alist,
               lty = "blank",
               fill = colors[1:length(alist)],
               fontface = "bold",
               cat.default.pos = pos,
               cat.col = rep("darkred",length(alist)),
               cat.cex = 1.2,
               #cat.pos = 0,
               cat.dist = 0.1,
               cat.fontface = "bold",
               main = title,
               main.pos = c(0.5, 1.0), main.fontface = "bold",
               main.fontfamily = "serif", main.col = "black",
               main.cex = 1.5, main.just = c(0, 0),
               alpha = 0.5,
               euler.d = TRUE,
               filename=NULL,
               height = 500,
               width = 500,
               units = "px",
               resolution = 500,
               margin = 0.1
               )

  #grid.newpage()
  if(grepl(".svg", fileName, fixed=T)){
    svg(fileName);
    grid.draw(venn.plot);
    dev.off()
  }else if(grepl(".png", fileName, fixed=T)){
    png(fileName);
    grid.draw(venn.plot);
    dev.off()
  }

}



## obtained from https://bioinformaticsbreakdown.com/how-to-gsea/
GSEA = function(gene_list, GO_file, pval) {
   set.seed(54321)
   library(dplyr)
   library(gage)
   library(fgsea)

   if ( any( duplicated(names(gene_list)) )  ) {
      warning("Duplicates in gene names")
      gene_list = gene_list[!duplicated(names(gene_list))]
   }
   if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
      warning("Gene list not sorted")
      gene_list = sort(gene_list, decreasing = TRUE)
   }
   myGO = fgsea::gmtPathways(GO_file)

   fgRes <- fgsea::fgsea(pathways = myGO,
                         stats = gene_list,
                         minSize=15,
                         maxSize=600,
                         nperm=10000) %>%
      as.data.frame() %>%
      dplyr::filter(padj < !!pval)
   print(dim(fgRes))

   ## Filter FGSEA by using gage results. Must be significant and in same direction to keep
   gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))

   ups = as.data.frame(gaRes$greater) %>%
      tibble::rownames_to_column("Pathway") %>%
      dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
      dplyr::select("Pathway")

   downs = as.data.frame(gaRes$less) %>%
      tibble::rownames_to_column("Pathway") %>%
      dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
      dplyr::select("Pathway")

   print(dim(rbind(ups,downs)))
   ## Define up / down pathways which are significant in both tests
   keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
   keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]

   fgRes = fgRes[ !is.na(match(fgRes$pathway,
                               c( keepups$pathway, keepdowns$pathway))), ] %>%
      arrange(desc(NES))
   fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")

   fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
   filtRes = rbind(head(fgRes, n = 10),
                   tail(fgRes, n = 10 ))


   upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
   downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
   colos = c(upcols, downcols)
   names(colos) = 1:length(colos)
   filtRes$Index = as.factor(1:nrow(filtRes))

   g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
      geom_col( aes(fill = Index )) +
      scale_fill_manual(values = colos ) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="GSEA - Biological Processes") +
      theme_minimal()


   output = list("Results" = fgRes, "Plot" = g)
   return(output)
}


xfun <- function (x, showCategory = 30, fill = "p.adjust", core_enrichment = TRUE,
                  label_format = 30, orderBy = "NES", decreasing = FALSE)
{
   if (!is(x, "gseaResult"))
      stop("currently only support gseaResult")
   if (fill == "qvalue") {
      fill <- "qvalues"
   }
   if (!fill %in% colnames(x@result)) {
      stop("'fill' variable not available ...")
   }
   if (orderBy != "NES" && !orderBy %in% colnames(x@result)) {
      message("wrong orderBy parameter; set to default `orderBy = \"NES\"`")
      orderBy <- "NES"
   }
   n <- showCategory
   if (core_enrichment) {
      gs2id <- geneInCategory(x)[seq_len(n)]
   }else {
      gs2id <- x@geneSets[x$ID[seq_len(n)]]
   }
   if (x@readable) {
      id <- match(names(x@geneList), names(x@gene2Symbol))
      names(x@geneList) <- x@gene2Symbol[id]
   }
   gs2val <- lapply(gs2id, function(id) {
      res <- x@geneList[id]
      res <- res[!is.na(res)]
   })
   nn <- names(gs2val)
   i <- match(nn, x$ID)
   nn <- x$Description[i]
   j <- order(x@result[[orderBy]][i], decreasing = decreasing)
   len <- sapply(gs2val, length)
   gs2val.df <- data.frame(category = rep(nn, times = len),
                           color = rep(x[i, fill], times = len), value = unlist(gs2val))
   colnames(gs2val.df)[2] <- fill
   gs2val.df$category <- factor(gs2val.df$category, levels = nn[j])
   label_func <- default_labeller(label_format)
   if (is.function(label_format)) {
      label_func <- label_format
   }
   ggplot(gs2val.df, aes_string(x = "value", y = "category", fill = fill)) +
      ggridges::geom_density_ridges() +
      scale_fill_continuous(low = "red",  high = "blue", name = fill, guide = guide_colorbar(reverse = TRUE)) +
      scale_y_discrete(labels = label_func) + xlab(NULL) +
      ylab(NULL) + theme_dose()
}
