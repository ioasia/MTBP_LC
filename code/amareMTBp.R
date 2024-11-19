#-------------------------based on xCell markers -------------------------------
# load the required libraries
library(openxlsx)
library(HGNChelper)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(GSVA)
library(ggpubr)
library(cowplot)
library(tidyestimate)
library(Cairo)

# all form of prepossessing, harmonization in the gene symbol across datasets, and normalization for timsTOF are done already.

# step 1: ---------------read the data------------------------------------------
# read the meta data of the early stage cohort, 2021
meta_data <- readxl::read_excel("./43018_2021_259_MOESM3_ESM.xlsx", sheet = "Table1-EarlyStageCohortMetadat")
colnames(meta_data)[colnames(meta_data) == 'Proteome Subtype'] <- 'Subtypes'

# read the xcell 2 markers data 
cells2markers <- read.table('./xCell2markers.txt', header = TRUE, sep = '\t')

# read the timsTOF proteomics data of the early stage cohort, 2024
proteomics_timsTOF <- read.csv('./LCproteomics_timsTOF.csv', row.names = 1)

# read the transcriptomics data of the early stage cohort, 2021
transcriptomics <- read.csv('./LCtranscriptomics.csv', row.names = 1)
anyNA(transcriptomics)

# step 2: ---------------analysis-----------------------------------------------
# for each cell type, count the number of overlapping genes between the data sets
signatures <- unique(cells2markers$signature)
overlap_counts <- lapply(signatures, function(i) {
  # Filter genes for the current signature
  genes_for_signature <- cells2markers %>%
    filter(signature == i) %>%
    pull(GeneSymbol)
  
  # Total genes in cells2markers for the current signature
  total_genes_cells2markers <- length(genes_for_signature)
  
  # Find overlaps with each dataset individually
  overlap_transcriptomic <- sum(genes_for_signature %in% rownames(transcriptomics))
  overlap_proteomics <- sum(genes_for_signature %in% rownames(proteomics_timsTOF))
  
  # Find genes that are in all three datasets
  genes_in_all_three <- genes_for_signature[
    genes_for_signature %in% rownames(proteomics_timsTOF) &
      genes_for_signature %in% rownames(transcriptomics)]
  
  # Count the genes overlapping in all three
  overlap_all_three <- length(genes_in_all_three)
  
  # Combine counts in a data frame
  data.frame(
    signature = i,  # Use `i` to represent the current signature
    Xcell = total_genes_cells2markers,
    TimsTOF = overlap_proteomics,
    mRNA = overlap_transcriptomic,
    Overlap_Three = overlap_all_three
  )
})


# Combine all counts into one data frame
genes_included_dt <- do.call(rbind, overlap_counts)
genes_included_dt$signature <- factor(genes_included_dt$signature, levels = genes_included_dt$signature[order(genes_included_dt$Xcell)])


# step 3: ---------------plot-------------------------------------------------
# plot the number of genes per cell type for each data type
genes_included_plot <- reshape2::melt(genes_included_dt, id.var = 'signature')
genes_included_plot$variable <- factor(genes_included_plot$variable, 
                                       levels = c('Xcell', 'mRNA', 'TimsTOF', "Overlap_Three"), 
                                       labels = c('Xcell', 'transcriptomics', 'proteomics(timsTOF)', 'Overlap(In three)'))



# Use ggplot2 to plot the data
ggplot(genes_included_plot, mapping = aes(x = signature, y = value, group = variable)) + 
  geom_hline(yintercept = 0) + 
  geom_line() + 
  geom_point() + 
  labs(y = '# markers') + 
  facet_grid(.~variable) + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
        panel.grid.major = element_blank())

# step 4: mRNA - protein correlation-------------------------------------------
signatures <- unique(cells2markers$signature)

# Data-set pair for correlation calculation
dt1 <- transcriptomics  # mRNA data
dt2 <- proteomics_timsTOF  # timsTOF data

# Get common samples between dt1 and dt2
common_samples <- intersect(colnames(dt1), colnames(dt2))

# Calculate correlations for each cell type
dt_marker_cell_types_cor <- do.call(rbind, lapply(signatures, function(i) {
  
  # Get marker genes for the current cell type
  markers <- cells2markers$GeneSymbol[cells2markers$signature == i]
  
  # Compute correlations for each marker gene
  marker_cor <- do.call(rbind, lapply(markers, function(j) {
    
    # Extract values for the gene across common samples
    val1 <- as.numeric(dt1[j, common_samples])
    val2 <- as.numeric(dt2[j, common_samples])
    
    # Count non-NA values and overlapping non-NA values
    n1 <- sum(!is.na(val1))
    n2 <- sum(!is.na(val2))
    n_overlap <- sum(!is.na(val1) & !is.na(val2))
    
    # Perform correlation test with error handling
    res_cor <- tryCatch(
      {
        res <- cor.test(val1, val2)
        c(res$estimate, res$p.value)
      }, 
      error = function(cond) {
        return(c(NA, NA))
      }
    )
    
    # Store results for each marker gene in a data frame
    res_dt <- data.frame(
      'signatures' = i,
      'marker' = j,
      'cor_mRNA_timsTOF' = res_cor[1],
      'pval_mRNA_timsTOF' = res_cor[2],
      'n1_mRNA' = n1,
      'n2_timsTOF' = n2,
      'n_overlap' = n_overlap
    )
    
    return(res_dt)
  }))
  
  return(marker_cor)
}))


# change the row names
rownames(dt_marker_cell_types_cor) <- 1:nrow(dt_marker_cell_types_cor)


# step 5: significant correlations----------------------------------------------
ggplot(dt_marker_cell_types_cor, aes(x = cor_mRNA_timsTOF)) +
  geom_histogram(binwidth = 0.05) + 
  theme_bw()

# Filter the positive correlation > 0.5 on the samples where at least 30% of the samples have values
cormRNA_timsTOF <- dt_marker_cell_types_cor %>% 
  dplyr::filter(cor_mRNA_timsTOF > 0.5, n2_timsTOF > length(common_samples) * 0.3)

# Pairwise correlation based on timeTOF quantification for each signature 
correlation_results <- lapply(signatures, function(sig) {
  
  # Get genes associated with the current signature
  genes_for_signature <- cormRNA_timsTOF %>%
    filter(signatures == sig) %>%
    pull(marker)
  
  # Check if these genes are in the proteomics dataset
  proteomics_genes <- intersect(genes_for_signature, rownames(proteomics_timsTOF))
  
  if (length(proteomics_genes) > 1) {  # Ensure that at least two genes are available for correlation
    # Extract corresponding quantification from proteomics_timsTOF
    proteomics_data <- proteomics_timsTOF[proteomics_genes, intersect(colnames(proteomics_timsTOF), colnames(transcriptomics))]
    
    # Compute pairwise correlation matrix (Pearson correlation)
    cor_matrix <- cor(t(proteomics_data), method = "pearson", use = "pairwise.complete.obs")
    
    # Return the result
    list(Signature = sig, CorrelationMatrix = cor_matrix)
  } else {
    return(NULL)  
  }
})

# Filter out no correlations 
correlation_results <- Filter(Negate(is.null), correlation_results)

# step 6: plot the correlation heatmap and expression heatmaps-----------------
SampleData <- subset(meta_data, SampleID %in% intersect(colnames(proteomics_timsTOF), colnames(transcriptomics)))
SampleData <- as.data.frame(SampleData)
rownames(SampleData) <- SampleData$SampleID

# color for the correlation heatmap
color_function <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Initialize lists to store the plots
list_plot <- list()
list_plot2 <- list()

# run the loop for each signature
if (length(correlation_results) > 0) {
  for (i in 1:length(correlation_results)) {
    corSignature <- correlation_results[[i]]
    print(corSignature[["Signature"]])
    
    # Generate correlation heatmap 
    cor_heatmap <- ComplexHeatmap::Heatmap(
      corSignature$CorrelationMatrix,
      name = "Cor",
      col = color_function,  
      show_row_names = TRUE,
      show_column_names = TRUE,
      show_column_dend = FALSE,
      show_row_dend = FALSE,
      cluster_rows = TRUE,  
      cluster_columns = TRUE,  
      clustering_distance_rows = "pearson",
      clustering_distance_columns = "pearson",
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      border = TRUE,
      row_names_gp = gpar(fontsize = 8, fontface = "bold"),
      column_names_gp = gpar(fontsize = 8, fontface = "bold"),
      column_title = paste(corSignature[["Signature"]]),
      column_title_gp = gpar(fontsize = 15, fontface = "bold")
    )
    
    # Draw the heatmap and get the row order for the expression heatmap
    cor_heatmap <- draw(cor_heatmap)
    ordered_genes_indices <- row_order(cor_heatmap)
    ordered_genes <- rownames(corSignature$CorrelationMatrix)[ordered_genes_indices]
    
    
    # order to match the correlation heatmap
    proteomics_data <- proteomics_timsTOF[ordered_genes, intersect(colnames(proteomics_timsTOF), colnames(transcriptomics))]
    proteomics_data <- proteomics_data[, match(rownames(SampleData), colnames(proteomics_data))]
    scaleData <- as.matrix(t(scale(t(proteomics_data))))
    
    # Generate the expression heatmap 
    exp_heatmap <- Heatmap(
      scaleData,
      name = "exp",
      row_order = ordered_genes,  
      column_split = SampleData$Subtypes, 
      cluster_rows = FALSE,  
      cluster_columns = FALSE,  
      show_column_names = FALSE,
      column_title_gp = gpar(fontsize = 15, fontface = "bold"),
      border = TRUE,
      row_names_gp = gpar(fontsize = 8, fontface = "bold"),
      column_title = paste(corSignature[["Signature"]]),
      top_annotation = HeatmapAnnotation(
        Subtypes = SampleData$Subtypes,
        which = "column",
        col = list(
          Subtypes = c(
            "1" = "deepskyblue4",
            "2" = "lightseagreen",
            "3" = "burlywood3",
            "4" = "orange",
            "5" = "pink",
            "6" = "red"
          )
        ),
        annotation_legend_param = list(Subtypes = list(title = "Subtypes", nrow = 3))
      )
    )
    
    # Add to the lists
    list_plot[[i]] <- cor_heatmap
    list_plot2[[i]] <- exp_heatmap
  }
} else {
  print("No correlation matrix could be computed.")
}


# the heatmaps side by side
for (i in 1:length(list_plot)) {
  CairoPDF(width = 10, height = 6, file =
             paste0("./xCell_LC_", correlation_results[[i]]$Signature, "_correlation_exp_heatmap.pdf"))
  grid.arrange(
    grid.grabExpr(print(list_plot2[[i]])),
    grid.grabExpr(print(list_plot[[i]])), ncol = 2)
  dev.off()
  
}

# step 7: GSVA analysis-----------------------------------------
# ssGSEA based on the full transcriptomics data
transcriptomics2ssgseaFull <-transcriptomics
transcriptomics2ssgseaFull <- t(scale(t(transcriptomics2ssgseaFull)))

# siginture list
cell_list <- lapply(unique(cells2markers$signature), function(i) {
  cells2markers$GeneSymbol[cells2markers$signature == i]
})
names(cell_list) <- unique(cells2markers$signature)

# build GSVA parameter object
ssgseaObj <- GSVA::ssgseaParam(
  exprData = as.matrix(transcriptomics2ssgseaFull),
  geneSets = cell_list,
  minSize = 2,
  maxSize = Inf,
  alpha = 0.25,
  # normalizing the scores(ssGSEA) by the absolute difference between the minimum and the maximum
  normalize = TRUE)

mRNA_ssGSEAscores = GSVA::gsva(param = ssgseaObj, verbose = TRUE)
mRNA_ssGSEAscores <- as.matrix(mRNA_ssGSEAscores)

# calculate the median ssGSEA score for each sample
Sample_median_mRNA_ssGSEAscores <- apply(mRNA_ssGSEAscores, 2, median)

# save the important data for further analysis and visualization
rm(list = ls())
gc()

#------I did the same analysis for single cell atlas----------------------------

