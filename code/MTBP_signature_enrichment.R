### MTBP signature analysis

## Libraries
library(openxlsx)
library(HGNChelper)
library(ggplot2)


# Function
# Update gene symbols
current_genes <- getCurrentHumanMap()

update.genes <- function(genes, map = current_genes) {
  
  checked_genes_dt <- checkGeneSymbols(genes, map = current_genes)
  checked_genes_dt$new <- ifelse(is.na(checked_genes_dt$Suggested.Symbol), checked_genes_dt$x,  checked_genes_dt$Suggested.Symbol)
  
  # Collapse multi-names
  checked_genes_dt$new <- sapply(strsplit(x = checked_genes_dt$new, split = '\\/'), function(mult) {
    return(trimws(mult[[1]]))
    
  })
  
  # Check doublets
  freq_genes <- sort(table(checked_genes_dt$new), decreasing = TRUE)
  dup_genes <- names(freq_genes)[freq_genes > 1]
  
  if(length(dup_genes) > 0) {
    
    todiscard_all <- sapply(dup_genes, function(dup) {
      
      dup_dt <- checked_genes_dt[checked_genes_dt$new %in% dup, ]
      todiscard <- dup_dt$x[dup_dt$Approved == FALSE]
    })
    
    checked_genes_dt <- checked_genes_dt[!checked_genes_dt$x %in% unlist(todiscard_all), ]
  }
  
  res <- checked_genes_dt$new[match(genes, checked_genes_dt$x)]
  
  
 
  
  return(res)
}


# # Plot function
# p_cor_rna_protein <- lapply(genes_toplot, function(i) {
#   
#   print(i)
#   dt <- data.frame('RNA' = as.numeric(transcriptomics_retro_quant[i, ]),
#                    'Protein' = as.numeric(proteomics_retro_quant[i, ]),
#                    'Subtype' = proteomics_retro_metadata$Subtype)
#   res_cor <- cor.test(dt$RNA, dt$Protein, method = 'pearson')
#   
#   p <- ggplot(dt, aes(x = RNA, y = Protein, col = Subtype, group = 1)) + 
#     geom_point() +
#     geom_smooth(method = 'lm', se = FALSE) + 
#     scale_color_manual(values = cluster_cols) +  
#     annotate(geom = 'text', x = Inf, y = Inf, 
#              label = paste0('cor. = ', round(res_cor$estimate,2), ', p = ', format(res_cor$p.value, scientific = TRUE, digits = 2), ', n = ', sum(complete.cases(dt))), vjust = 1.1, hjust = 1) + 
#     labs(x = 'Transcriptomics (log2-TPM)', y = 'Proteomics (log2 ratio)', title = i) + 
#     theme_classic() + 
#     theme(plot.title = element_text(hjust = 0.5)) + 
#     guides(shape = guide_legend(override.aes = list(size = 3)),
#            col = guide_legend(override.aes = list(linetype = 0, fill = NA)),
#            fill = guide_legend(override.aes = list(linetype = 0)))
#   
#   

# Paths 
path_file <- '/Users/ioannis/Downloads/MTBP_project/data/'
path_save <- '/Users/ioannis/Downloads/MTBP_project/code/figures/'

# cell types
cells <- read.xlsx(paste0(path_file, 'High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer_markers_mmc3.xlsx'), startRow = 2)


# getSheetNames('~/Downloads/LungCancer_project/data/43018_2021_259_MOESM3_ESM.xlsx')
# meta data
meta_data <- read.xlsx('~/Downloads/LungCancer_project/data/43018_2021_259_MOESM3_ESM.xlsx', sheet = "Table1-EarlyStageCohortMetadat")

# Proteomics
proteomics_tmt <- read.xlsx('~/Downloads/LungCancer_project/data/43018_2021_259_MOESM3_ESM.xlsx', sheet = "Table2-DDAProteomics")
proteomics_tmt$GeneSymbol <- update.genes(genes = proteomics_tmt$GeneSymbol, map = current_genes)

# Transcriptomics
transcriptomics <- read.xlsx('~/Downloads/LungCancer_project/data/43018_2021_259_MOESM3_ESM.xlsx', 
                             sheet = "Table2-Transcriptomics")
transcriptomics$GeneSymbol <- update.genes(genes = transcriptomics$GeneSymbol, map = current_genes)


# DIA data
proteomics_dia <-  read.xlsx('~/Downloads/LungCancer_project/data/43018_2021_259_MOESM3_ESM.xlsx', sheet = "Table2-DIAProteomics")
proteomics_dia_MS2 <- proteomics_dia[, c(1, grep('PG', colnames(proteomics_dia)))]
proteomics_dia_MS2$GeneSymbol <- update.genes(genes = proteomics_dia_MS2$GeneSymbol, map = current_genes)


# Remove NA
sum(is.na(proteomics_tmt$GeneSymbol))
proteomics_tmt <- proteomics_tmt[!is.na(proteomics_tmt$GeneSymbol), ]
rownames(proteomics_tmt) <- proteomics_tmt$GeneSymbol
proteomics_tmt$GeneSymbol <- NULL

sum(is.na(transcriptomics$GeneSymbol))
transcriptomics <- transcriptomics[!is.na(transcriptomics$GeneSymbol), ]
rownames(transcriptomics) <- transcriptomics$GeneSymbol
transcriptomics$GeneSymbol <- NULL

sum(is.na(proteomics_dia_MS2$GeneSymbol))
proteomics_dia_MS2 <- proteomics_dia_MS2[!is.na(proteomics_dia_MS2$GeneSymbol), ]
rownames(proteomics_dia_MS2) <- proteomics_dia_MS2$GeneSymbol
proteomics_dia_MS2$GeneSymbol <- NULL
colnames(proteomics_dia_MS2) <- gsub('_PG.Quantity', '', colnames(proteomics_dia_MS2))

# Normalize
proteomics_dia_MS2 <- log10(proteomics_dia_MS2) 
proteomics_dia_MS2 <- sweep(proteomics_dia_MS2, 1, apply(proteomics_dia_MS2, 1, mean, na.rm = TRUE), '-')
proteomics_dia_MS2 <- sweep(proteomics_dia_MS2, 2, apply(proteomics_dia_MS2, 2, mean, na.rm = TRUE), '-')



cells$GeneSymbol <-  update.genes(genes = cells$gene_symbol, map = current_genes)
sum(is.na(cells$GeneSymbol))
cells <- cells[!is.na(cells$GeneSymbol), ]


## Overlap
genes_per_cell <-  as.data.frame(table(cells$signature))
genes_included <- lapply(list(cells$GeneSymbol, 
                              rownames(proteomics_tmt),
                              rownames(proteomics_dia_MS2),
                              rownames(transcriptomics)), function(i) {
                          
                                cells_included <- cells[cells$GeneSymbol %in% i, ]
                                genes_per_cell_included <- as.data.frame(table(cells_included$signature))   
                         })

merge.all <- function(i,j) merge(i,j, by = 'Var1', all = TRUE)

genes_included_dt <- Reduce(f = merge.all, genes_included)
genes_included_dt[is.na(genes_included_dt)] <- 0
colnames(genes_included_dt) <- c('Cell_type', 'Init', 'TMT', 'DIA_MS2', 'RNA')
write.table(genes_included_dt, paste0(path_file, 'LC_marker_overlap.txt'), sep = '\t', row.names = FALSE)


# Plot
genes_included_dt$Cell_type <- factor(genes_included_dt$Cell_type, levels = genes_included_dt$Cell_type[order(genes_included_dt$Init)])
genes_included_plot <- reshape2::melt(genes_included_dt, id.var = 'Cell_type')

p_overl <- ggplot(genes_included_plot, mapping = aes(x = Cell_type, y = value, group = variable)) + 
  geom_hline(yintercept = 0) + 
  geom_line() + 
  geom_point() + 
  labs(y = '# markers') + 
  facet_grid(.~variable) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank()) 

pdf(paste0(path_save, 'LC_marker_overlap.pdf'), width = 15, height = 5)
plot(p_overl)
dev.off()

# mRNA - protein correlation
cells_tocheck <- unique(cells$signature)


# Known marker genes for lung cancer
dt_list <- list('mRNA_tmt' = list('mRNA' = transcriptomics, 
                       'tmt' = proteomics_tmt),
                'tmt_dia' = list('tmt' = proteomics_tmt, 
                                  'dia' = proteomics_dia_MS2))

dt_marker_cell_types_cor <- lapply(names(dt_list), function(k) {
  
  dts <- dt_list[[k]]
  dt1 <- dts[[1]]
  dt2 <- dts[[2]]
  
  common_samples <- intersect(colnames(dt1), colnames(dt2))
  
  marker_cell_types_cor <- do.call(rbind, lapply(cells_tocheck, function(i) {

    # print(i)
    markers <- cells$GeneSymbol[cells$signature %in% i]
  
    marker_cor <- do.call(rbind, lapply(markers, function(j) {
    
      # print(j)
      # mRNA - protein
      val1 <- as.numeric(dt1[j, common_samples])
      val2 <-  as.numeric(dt2[j, common_samples])
  
      n1 <-  sum(!is.na(val1))
      n2 <-  sum(!is.na(val2))
      n_overlap <-  sum(!is.na(val1) & !is.na(val2))
    
      # plot(val1,  val2)
   
      res_cor <-  tryCatch(
        {
          res <- cor.test(val1, val2)
          c(res$estimate, res$p.value)
  
        }, error = function(cond) {
        
          return(c(NA, NA))
        }
      )
      
      res_dt <- data.frame(
        'cell_type' = i,
        'marker' = j,
        'cor' = res_cor[1], 
        'pval' = res_cor[2], 
        'n1' = n1,
        'n2' = n2,
        'n_overlap' = n_overlap)
      
      return(res_dt)
         }))
    }))
  colnames(marker_cell_types_cor)[-c(1,2)] <- paste0(k, '_', colnames(marker_cell_types_cor)[-c(1,2)])
  
  return(marker_cell_types_cor)
  })
  

merge.dt.cor <- function(i,j) merge(i,j, by = c( "cell_type", "marker"))

dt_marker_cell_types_cor_all <- Reduce(f = merge.dt.cor, dt_marker_cell_types_cor)



