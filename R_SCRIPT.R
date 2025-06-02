#Step 1: Data Download and Loading

BiocManager::install("STRINGdb")
# Load required libraries
library(GEOquery)
library(DESeq2)
library(limma)
library(edgeR)
library(tidyverse)
library(biomaRt)
library(dplyr)
library(stringr)

library(AnnotationDbi)
library(org.Hs.eg.db)

gse <- getGEO("GSE50760", GSEMatrix = TRUE, getGPL = FALSE)
pheno_data <- pData(gse[[1]])
expr_data <- read_tsv("D:/R programs/RNA/GSE50760_norm_counts_FPKM_GRCh38.p13_NCBI.tsv")

cat("Expression matrix dimensions:", dim(expr_data), "\n")
cat("Phenotype data dimensions:", dim(pheno_data), "\n")

# Find samples in expr_data not in pheno_data
mismatch <- setdiff(colnames(expr_data), rownames(pheno_data))
print(mismatch)

# Remove the unmatched column from expr_data
expr_data_clean <- expr_data[, colnames(expr_data) %in% rownames(pheno_data)]

# Reorder columns of expr_data_clean to match the rows of pheno_data
expr_data_clean <- expr_data_clean[, rownames(pheno_data)]

cat("Cleaned expression matrix dimensions:", dim(expr_data_clean), "\n")
cat("Phenotype data dimensions:", dim(pheno_data), "\n")

# Preview data
head(expr_data[, 1:5])
head(pheno_data[, c("title", "characteristics_ch1")])

pheno_data_df <- as.data.frame(pheno_data)

#Step 2: Data Preprocessing and Quality Control
# Extract sample information
sample_info <- pheno_data_df %>%
  dplyr::select(title, characteristics_ch1) %>%
  dplyr::mutate(
    sample_type = ifelse(grepl("normal|Normal", characteristics_ch1), 
                         "Normal", "Tumor"),
    patient_id = str_extract(title, "\\d+")
  )

# Create design matrix
design_matrix <- data.frame(
  sample_id = rownames(sample_info),
  condition = factor(sample_info$sample_type, levels = c("Normal", "Tumor")),
  patient = sample_info$patient_id
)


# Quality control plots
library(ggplot2)
library(pheatmap)

# 1. Distribution of expression values
expr_long <- expr_data %>%
  as.data.frame() %>%
  rownames_to_column("probe") %>%
  pivot_longer(-probe, names_to = "sample", values_to = "expression")

# Expression distribution plot
p1 <- ggplot(expr_long, aes(x = expression, color = sample)) +
  geom_density(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Expression Distribution Across Samples") +
  theme(legend.position = "none")

# 2. Sample correlation heatmap
cor_matrix <- cor(expr_data_clean, use = "complete.obs")
pheatmap(cor_matrix, 
         annotation_col = data.frame(
           Condition = design_matrix$condition,
           row.names = colnames(expr_data_clean)
         ),
         main = "Sample Correlation Heatmap")

# 3. Principal Component Analysis
pca_data <- prcomp(t(expr_data_clean), scale = TRUE)

# Remove genes with zero variance
non_zero_var_genes <- apply(expr_data_clean, 1, function(x) var(x) != 0)
expr_data_pca_ready <- expr_data_clean[non_zero_var_genes, ]

# Perform PCA
pca_data <- prcomp(t(expr_data_pca_ready), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Condition = design_matrix$condition
)

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Samples",
       x = paste0("PC1 (", round(summary(pca_data)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_data)$importance[2,2]*100, 1), "%)"))

print(p1)
print(p2)


'''Step 3: Differential Expression Analysis
3.1 DESeq2 Analysis
# Convert expression data to integer counts (simulate RNA-seq data)
# In real RNA-seq, you would start with raw count data
count_data <- round(2^expr_data_pca_ready)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = design_matrix,
  design = ~ condition
)
'''


#3.2 limma Analysis
# Design matrix for limma
design_limma <- model.matrix(~ condition, data = design_matrix)

# Fit linear model
fit <- lmFit(expr_data_pca_ready, design_limma)
fit <- eBayes(fit)

# Extract results
limma_results <- topTable(fit, coef = "conditionTumor", 
                          number = Inf, adjust.method = "BH")

# Add fold change and significance flags
limma_df <- limma_results %>%
  rownames_to_column("probe_id") %>%
  mutate(
    log2FoldChange = logFC,
    padj = adj.P.Val,
    gene_symbol = mapIds(
      org.Hs.eg.db,
      keys = probe_id,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  ) %>%
  arrange(padj)

# Significant genes (limma)
limma_sig <- limma_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

cat("limma: Significant DEGs:", nrow(limma_sig), "\n")


# Volcano plots
create_volcano_plot <- function(data, title, fc_col, pval_col) {
  data$significance <- "Not Significant"
  data$significance[data[[pval_col]] < 0.05 & abs(data[[fc_col]]) > 1] <- "Significant"
  
  fc_col_sym <- sym(fc_col)
  pval_col_sym <- sym(pval_col)
  
  ggplot(data, aes(
    x = !!fc_col_sym,
    y = -log10(!!pval_col_sym),
    color = significance
  )) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)
}
# Create volcano plots

p_volcano_limma <- create_volcano_plot(limma_df, "limma Volcano Plot", 
                                       "log2FoldChange", "padj")
print(p_volcano_limma)


library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(enrichplot)

# Filter out rows with missing gene symbols
limma_sig_filtered <- limma_sig[!is.na(limma_sig$gene_symbol), ]

# Convert SYMBOLs to ENTREZ IDs for pathway analysis

limma_entrez <- bitr(
  limma_sig_filtered$gene_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
limma_annotated <- merge(limma_sig_filtered, limma_entrez, 
                         by.x = "gene_symbol", by.y = "SYMBOL")

# Gene Ontology enrichment analysis
go_bp <- enrichGO(gene = limma_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

go_mf <- enrichGO(gene = limma_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

go_cc <- enrichGO(gene = limma_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
# KEGG pathway analysis
kegg_pathways <- enrichKEGG(gene = limma_entrez$ENTREZID,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)

# Reactome pathway analysis
reactome_pathways <- enrichPathway(gene = limma_entrez$ENTREZID,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   readable = TRUE)

# Visualization of enrichment results
# GO Biological Process
p_go_bp <- dotplot(go_bp, showCategory = 20) + 
  ggtitle("GO Biological Process Enrichment")
p_go_bp

# KEGG pathways
p_kegg <- dotplot(kegg_pathways, showCategory = 15) + 
  ggtitle("KEGG Pathway Enrichment")
p_kegg

# Reactome pathways
p_reactome <- dotplot(reactome_pathways, showCategory = 15) + 
  ggtitle("Reactome Pathway Enrichment")
p_reactome

# Gene-concept network
if(nrow(as.data.frame(kegg_pathways)) > 0) {
  p_cnet <- cnetplot(kegg_pathways, categorySize="pvalue", foldChange = NULL)
  print(p_cnet)
}

#Step 6: Protein-Protein Interaction (PPI) Network Analysis
library(STRINGdb)

# Initialize STRING database
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400)

# Map genes to STRING IDs
limma_mapped <- string_db$map(limma_sig_filtered, 
                              "gene_symbol", 
                              removeUnmappedRows = TRUE)

# Get PPI network
ppi_network <- string_db$get_interactions(limma_mapped$STRING_id)

# Network statistics
cat("PPI Network Statistics:\n")
cat("Nodes:", length(unique(c(ppi_network$from, ppi_network$to))), "\n")
cat("Edges:", nrow(ppi_network), "\n")


# Find hub genes (genes with many interactions)
node_degrees <- table(c(ppi_network$from, ppi_network$to))
hub_genes <- names(sort(node_degrees, decreasing = TRUE))[1:10]
str(hub_genes)

# Map back to gene symbols
# Match STRING IDs (hub_genes) to gene symbols
hub_gene_symbols <- limma_mapped[limma_mapped$STRING_id %in% hub_genes, c("gene_symbol", "STRING_id")]
print("Top 10 Hub Genes:")
print(hub_gene_symbols)


# Plot PPI network (simplified)
if(nrow(ppi_network) > 0) {
  string_db$plot_network(limma_mapped$STRING_id[1:50])  # Plot subset for clarity
}

#Step 7: Results Summary and Export
# Create comprehensive results summary


results_summary <- data.frame(
  Method = ("limma"),
  Total_DEGs =  nrow(limma_sig_filtered),
  Upregulated = sum(limma_sig_filtered$log2FoldChange > 0),
  Downregulated = sum(limma_sig_filtered$log2FoldChange < 0))

print("Results Summary:")
print(results_summary)

# Key pathways related to cancer
cancer_pathways <- c("Cell cycle", "DNA repair", "Apoptosis", "Cell adhesion", 
                     "Immune response", "Angiogenesis", "Metabolic pathways")

if(nrow(as.data.frame(go_bp)) > 0) {
  cancer_related_go <- as.data.frame(go_bp) %>%
    filter(grepl(paste(cancer_pathways, collapse = "|"), Description, ignore.case = TRUE))
  
  cat("\nCancer-related GO terms found:", nrow(cancer_related_go), "\n")
  if(nrow(cancer_related_go) > 0) {
    print(cancer_related_go[1:5, c("Description", "pvalue", "qvalue")])
  }
}



# Export results
write.csv(limma_sig_filtered, "significantgenes.csv", row.names = FALSE)
write.csv(as.data.frame(go_bp), "GO_biological_process_enrichment.csv")
write.csv(as.data.frame(kegg_pathways), "KEGG_pathway_enrichment.csv")
write.csv(results_summary, "analysis_summary.csv")

cat("\nAnalysis completed! Key files exported:\n")
cat("- significantgenes.csv: List of differentially expressed genes\n")
cat("- GO_biological_process_enrichment.csv: GO enrichment results\n")
cat("- KEGG_pathway_enrichment.csv: KEGG pathway enrichment results\n")
cat("- analysis_summary.csv: Summary of all methods\n")



