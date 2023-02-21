library(DESeq2)
library(tidyverse)
library(here)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
source("src/rules/scripts/functions.R")

results_dir <- here("results/R_analysis")

dds_obj <- readRDS(file.path(results_dir, "objs/dds_obj.rda"))
vst_obj <- readRDS(file.path(results_dir, "objs/vsd_obj.rda"))

dds_results <- read.csv(file.path(results_dir,
                                  "DE_files", 
                                  "PTPN2_KO_NoTreatment_vs_PTPN2_KO_Cytokine_Treatment.csv"))

dds_res_wt <- read.csv(file.path(results_dir,
                                 "DE_files",
                                 "WT_NoTreatment_vs_WT_Cytokine_Treatment.csv"))

dds_results$log2FoldChange <- -1 * dds_results$log2FoldChange

dds_res_wt$log2FoldChange <- -1 * dds_res_wt$log2FoldChange

dds_results_short <- dds_results %>%
  dplyr::top_n(n = -30, wt = padj) %>%
  dplyr::mutate(full_name = paste(gene_name, ens_id, sep = "_"))

# Heatmap ----------------------------------------------------------------------
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                         "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

heatmap <- make_heatmap(dds = dds_obj,
                        vsd = vst_obj,
                        de_genes = dds_results_short$full_name,
                        treatment = "PTPN2_KO_Cytokine_Treatment",
                        control = "PTPN2_KO_NoTreatment",
                        group = "group", 
                        print_genenames = TRUE, gene_identifier = "Gene_ID",
                        cluster_cols = FALSE, save_heatmap = FALSE,
                        output_dir = "results", plot_groups = "all",
                        color_test = NULL, save_name = NULL)

# Volcano ----------------------------------------------------------------------
dds_results_short <- dds_results %>%
  dplyr::top_n(n = 10, wt = abs(log2FoldChange)) %>%
  dplyr::mutate(full_name = paste(gene_name, ens_id, sep = "_"))

volcano_plot_1 <- EnhancedVolcano(dds_results,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  selectLab = dds_results_short$gene_name,
                                  lab = dds_results$gene_name,
                                  labSize = 5,
                                  colAlpha = 1,
                                  drawConnectors = TRUE,
                                  title = "PTPN2KO",
                                  subtitle = "")


volcano_plot_2 <- EnhancedVolcano(dds_res_wt,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  selectLab = dds_results_short$gene_name,
                                  lab = dds_res_wt$gene_name,
                                  labSize = 5,
                                  colAlpha = 1,
                                  drawConnectors = TRUE,
                                  title = "WT",
                                  subtitle = "")
pdf(file.path(results_dir, "images", "ptpn2_cyto_genes.pdf"))
print(heatmap)
print(volcano_plot_1)
print(volcano_plot_2)
dev.off()

# specific treatment -----------------------------------------------------------

specific_trt <- read.csv(file.path(results_dir, "DE_files",
                                   "PTPN2_specific_trt_all.csv"))


specific_trt_short <- specific_trt %>%
  dplyr::top_n(n = -30, wt = padj) %>%
  dplyr::mutate(full_name = paste(gene_name, ens_id, sep = "_"))

# Heatmap ----------------------------------------------------------------------
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                         "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
                         
heatmap <- make_heatmap(dds = dds_obj,
                        vsd = vst_obj,
                        de_genes = specific_trt_short$full_name,
                        treatment = "PTPN2_KO_Cytokine_Treatment",
                        control = "PTPN2_KO_NoTreatment",
                        group = "group", 
                        print_genenames = TRUE, gene_identifier = "Gene_ID",
                        cluster_cols = FALSE, save_heatmap = FALSE,
                        output_dir = "results", plot_groups = "all",
                        color_test = NULL, save_name = NULL)

# Volcano ----------------------------------------------------------------------
dds_results_short <- dds_results %>%
  dplyr::top_n(n = 10, wt = abs(log2FoldChange)) %>%
  dplyr::mutate(full_name = paste(gene_name, ens_id, sep = "_"))

volcano_plot_1 <- EnhancedVolcano(dds_results,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  selectLab = specific_trt_short$gene_name,
                                  lab = dds_results$gene_name,
                                  labSize = 5,
                                  colAlpha = 1,
                                  drawConnectors = TRUE,
                                  title = "PTPN2KO",
                                  subtitle = "")


volcano_plot_2 <- EnhancedVolcano(dds_res_wt,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  selectLab = specific_trt_short$gene_name,
                                  lab = dds_res_wt$gene_name,
                                  labSize = 5,
                                  colAlpha = 1,
                                  drawConnectors = TRUE,
                                  title = "WT",
                                  subtitle = "")
pdf(file.path(results_dir, "images", "ptpn2_specific_cyto_genes.pdf"))
print(heatmap)
print(volcano_plot_1)
print(volcano_plot_2)
dev.off()