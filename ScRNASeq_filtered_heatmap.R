
library(pheatmap)

cytokines <- c("ALL") #set treatments

red2blue <- colorRampPalette(c("#2066b4", "#ffffff", "#b4182a"))(100) #set colour scale

  for (cytokine in cytokines) {
    res_means <-read.delim(paste0("significant_data/",cytokine,"_deg_table_1L2FC_nobaseMean_cutoff_004.csv"), 
                           header = TRUE, row.names = 1, sep = ",")%>%
      mutate_all(~replace(., .=="", "no")) # read in DEGs and remove blanks
    
    potential_targets <- read.delim(paste0("reference_data/cocktail_drug_target_screen_with_percentiles.csv"), header = TRUE, sep = ",") %>%
      as.data.frame() #read in filter list
    potential_targets <- potential_targets %>%
      filter(GSE13592 == 1)
    
    merged_dat <- merge(res_means, potential_targets, by.x = "symbol", by.y = "gene") #transect DEG and filter list
    
    print(paste0(cytokine," filtered list generated")) #print message to show it worked
    
    #write.csv(pheno_list[[cytokine]], file = paste0(out_path,cytokine,"_",phenotype,"_filtered_genes_1L2FC.csv")) # write out a csv
  }
  
  heat <- merged_dat
    row.names(heat) <- heat$symbol #pheatmap must have gene names in rownames
    heat$log2FoldChange <- as.numeric(heat$log2FoldChange) #select only l2FC values
    heat <- heat %>%
      dplyr::select(log2FoldChange) %>% #select only symbols
      na.omit()
    
    p <- pheatmap(heat, 
                  kmeans_k = NA, 
                  breaks = seq(-3, 3, length.out = 101), 
                  scale = "none", 
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  clustering_distance_rows = "maximum",  # Corrected parameter name
                  clustering_method = "complete",
                  color = red2blue,
                  show_colnames = FALSE,
                  fontsize = 14)
    
    ggsave(filename = paste0(out_path,"scRNAFiltered_DEGs_002.svg"), plot = p[[4]], width = 10, height = (nrow(heat) * 0.5), units = "cm", limitsize = FALSE)
  
  
  names(heatmaps) <- names(heat)
  
  pheno_overlap <- Reduce(function(x, y) merge(x, y, by = "symbol", all = TRUE), pheno_list)
  pheno_overlap <- pheno_overlap %>% # as the reduce function generates replicate names,
    setNames(make.unique(names(.)))  # we have to manually change them to be unique
  
  heat_overlap <- pheno_overlap %>%
    dplyr::select(3,18,33,48)  %>%
    rename_with(~ c("ALL","TNFa", "Thrombin", "TGFb2")) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix
  
  rownames(heat_overlap) <-pheno_overlap$symbol
  
  p_over <- pheatmap(heat_overlap, 
                     kmeans_k = NA, 
                     breaks = seq(-3, 3, length.out = 101), 
                     scale = "none", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     clustering_distance_rows = "euclidean",  # Corrected parameter name
                     clustering_method = "median",
                     color = red2blue,
                     show_colnames = TRUE,
                     fontsize = 14)
  
  ggsave(filename = paste0(out_path, cytokine,"_scRNAFiltered_DEGs_001.svg"), plot = p_over[[4]], width = 15, height = (3 + (nrow(heat[[cyto]]) * 0.5)), units = "cm")
  #write.csv(pheno_overlap, file = paste0(out_path,cytokine,"_",phenotype,"_filtered_genes_1L2FC.csv"))
}
```