
options(error=recover)
project_dir <- "C:/GREENBLATT/Rscripts/CrisprScreen"
source(file.path(project_dir, "R/crispr_screen_analysis_lib.R"))
setup()

wd <- c("C:/GREENBLATT/Ahmed/ebv/Nov24_2023/results")
hwd <- c("C:/GREENBLATT/Ahmed/ebv/Nov24_2023/results/final_analysis")

scale_normalize(wd, hwd, mincount = 10, maxcount = 10000)
qGI_analysis(wd, hwd, runGSEA = TRUE)

#############################################################
#####           process chronos results                #####
############################################################# 


#############################################################
#####           process MAGeCK results                  #####
#############################################################      

setwd(wd)
preprocess_mageck_dropout(wd) ## move data to rc cluster and run MAGeCK test
fdr_df <- collect_MAGeCK_results(file.path(wd, "MAGeCK"))
mageck_essential_list <- essential_ROC(fdr_df, cutoff = 0.05, ge = FALSE)


#############################################################
#####           process BAGEL results                   #####
#############################################################
preprocess_bagel_dropout(hwd) 
bf_df <- collect_BAGEL_results(file.path(hwd, "BAGEL"))
bagel_essential_list <-  essential_ROC(bf_df, cutoff = 6, ge = TRUE)

#############################################################
#####           overlap with DepMap EG                  #####
#############################################################
agsfile <- "C:/GREENBLATT/Ahmed/depmap/CRISPRGeneDependency.csv"
cfile <- "C:/GREENBLATT/Ahmed/depmap/CRISPR_common_essentials_22Q2.csv"
hfile <- "C:/data/raw/EDYTA/dropoutScreen/ListofEssentialGenesHEK293.tab"

ags_mat <- read.csv(agsfile)
ags_all<- ags_mat %>%
  filter(X == "ACH-000880")
ags_all <- ags_all[, -1]
ags_essential <-  colnames(ags_all)[ags_all[1,] > 0.9]
ags_essential <- unlist(lapply(strsplit(ags_essential, split="..", fixed = TRUE), 
                               function(x)x[1]))

public_essential_list <-lapply(c(cfile, hfile), function(x){
  genes_df <- read.delim(x)
  genes_df[,1]
})
names(public_essential_list) <- c("depmap_common", "Hek293")

tmp <- unlist(lapply(strsplit(public_essential_list[["depmap_common"]], split=" "), 
                     function(x)x[1]))
public_essential_list[["depmap_common"]] <- tmp

public_essential_list[["depmap_AGS"]] <- ags_essential[!is.na(ags_essential)]


combined_list <- c(mageck_essential_list, bagel_essential_list, public_essential_list)
names(combined_list) <- c(paste0("mageck_", names(mageck_essential_list)),
                          paste0("bagel_", names(bagel_essential_list)),
                          names(public_essential_list))

pdf(file.path(wd, "Essential_gene_Venndiagram.pdf"), height=8, width=8)
overlap_pair(combined_list[1:2], intersect, "MAGeCK_essential_gene_overlap")
overlap_pair(combined_list[3:4], intersect, "BAGEL_essential_gene_overlap")
overlap_pair(combined_list[c("depmap_common", "depmap_AGS")], intersect, "DepMap_essential_gene")
overlap_quad(combined_list[1:4], intersect, "MAGeCK_BAGEL_essential_gene_overlap")
overlap_quad(combined_list[c("bagel_AGS", "mageck_AGS", "depmap_common", "depmap_AGS")], intersect, "MAGeCK_BAGEL_DepMap_essential_gene_overlap")
overlap_quad(combined_list[c("bagel_AGSEBV", "mageck_AGSEBV", "depmap_common", "depmap_AGS")], intersect, "MAGeCK_BAGEL_DepMap_essential_gene_overlap_EBV")
dev.off()

combined_df <- union2mat(combined_list)
EBV_df <- combined_df %>%
  dplyr::filter(mageck_AGSEBV==1 &	mageck_AGS==0 &	bagel_AGS==0 & bagel_AGSEBV ==1)
write.table(EBV_df, file.path(wd, "EBV_essential_gene_table.tab"),
            sep= "\t", col.names = NA)

pdf(file.path(wd, "overlap_with_common_essential_genes.pdf"), height=8, width=10)
p <- UpSetR::upset(combined_df, nsets=7, nintersects = 100)
print(p)
dev.off()

#############################################################
#####             Integrate qGI with EG                 #####
#############################################################

qGI_df <- read.delim2(file.path(hwd, "qGI", "AGSEBV_qGI_moderated_t_test_results.tab"))
combined_df[["GENE"]] <- rownames(combined_df)
merged_df <- merge(qGI_df, combined_df, by.x="X", by.y="GENE", all.x=TRUE, sort = FALSE)
write.table(merged_df, file.path(hwd, "qGI", "AGSEBV_qGI_moderated_t_test_results_with_essentiality.tab"),
            col.names = NA, sep="\t")

#############################################################
#####           Inspect individual guides               #####
#############################################################

setwd(hwd)

top_gene_list <- read.delim2(file.path(hwd, "qGI", "AGSEBV_qGI_moderated_t_test_results.tab")) %>%
  dplyr::filter(adj.P.Val < 0.2, sign(as.numeric(logFC)) == -1) %>%
  dplyr::pull(X)
normalized_df <- read.delim2("combined_guide_normalizedCount_table.tab")
rawcount_df <- read.delim2("combined_guide_rawCount_table.tab")
lfc_df <- read.delim2("combined_guide_lfc_table.tab")
df_list <- list("Raw count"=rawcount_df, "Normalize count"=normalized_df, "Log fold change"=lfc_df)

pdf("guide_details_plot.pdf", height=10, width=5)
oldpar <- par(mfcol = c(3,1))
lapply(top_gene_list, function(agene){
  lapply(names(df_list), function(aname){
    df <- df_list[[aname]]
    control_col <- colnames(df)[grepl("AGS_", colnames(df))]
    treat_col <- colnames(df)[grepl("AGSEBV_", colnames(df))]
    
    inspect_data_plot(data_mat=df, gene=agene, ctl=control_col, 
                    exp=treat_col, ylabel=aname)
  })
})
par(oldpar)
dev.off()
