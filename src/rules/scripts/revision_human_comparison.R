# Check overlaps of DE genes
# Venn diagram of each pairwise comparison?
# 1. WT vs PTPN no cytokine
# 2. WT vs PTPN cytokine
# 3. WT cytokine vs no cytokine
# 4. PTPN cytokine vs no cytokine

# The "all" files haven't been subset to only signficant genes. Read these
# in so you have more control.

# Make a function
# Takes the path to the mouse file
# Takes the path to the excel file
# Takes the name of the excel tab
# Takes a mapping of human to mouse genes
# Computes overlaps

# Compute a table of human ens ids to mouse gene names only once. Save this
# table.
# biomart kept giving errors: Error: biomaRt has encountered an
# unexpected server error.
# So I'm just using the table I saved

library(here)
library(tidyverse)
library(openxlsx)
library(VennDiagram)
library(DESeq2)
library(RColorBrewer)
library(gprofiler2)

# It would probably be beneficial to eventually use ens ids
# Name mapping files
gene_conversion_dir <- "/Users/wellskr/Documents/Analysis/references"
all_genes <- "human_mouse_gene_mapping_20221003.csv"
one2one <- "human_mouse_gene_mapping_20221003_one2one.csv"
human_conversion_all <- "human_gene_mapping_ens_20230609.csv"
human_one2one_conversion <- "human_gene_one2one_mapping_ens_20230609.csv"

# Data dirs
mouse_res <- readRDS(here("results/R_analysis/objs/dds_obj.rda"))
## excel file of human data 
## Downloaded https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
excel_file <- here("files/GSE172148_PTPN2_DoubleBatchEffect_Limma.xlsx")
all_names <- openxlsx::getSheetNames(excel_file)
## Russ PTPN2 data
## Downloaded from the supplement https://www.mdpi.com/2073-4409/11/23/3845
russ_file <- here("files/ptpn2_sBC/de_genes.xlsx")

# Set up dirs
save_dir <- here("results/R_analysis/images/human_overlap")

# Results dir
results_dir <- here("results/R_analysis/DE_files")

ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

# Read in gene mapping
all_gene_list <- read.csv(file.path(gene_conversion_dir, all_genes))
ens_mapping <- read.csv(file.path(gene_conversion_dir, human_conversion_all))

colnames(all_gene_list) <- c("X", "hgnc_symbol", "mgi_symbol")

all_gene_list <- merge(all_gene_list, ens_mapping, by = "hgnc_symbol") 
all_gene_list <- all_gene_list %>%
  dplyr::select(hgnc_symbol, ensembl_gene_id, mgi_symbol) %>%
  dplyr::distinct()

name_mapping <- list("PTPN2_WT_no_trt" = 
                       c("siControlvssiPTPN2-1417", "PTPN2_WT_noTrt_all"),
                     "WT_trt_IFNg" = 
                       c("siControlvssiControlIFNg-3853", "WT_trt_all"),
                     "PTPN2_trt_IFNg" = 
                       c("siPTPN2vssiPTPN2IFNg-6435", "PTPN2_trt_all"),
                     "WT_trt_IFNa" = 
                       c("siControlvssiControlIFNa-104", "WT_trt_all"),
                     "PTPN2_trt_IFNa" = 
                       c("siPTPN2vssiPTPN2IFNa-437", "PTPN2_trt_all"),
                     "PTPN2_WT_trt_IFNg" = 
                       c("siControlIFNgvssiPTPN2IFNg-720", "PTPN2_WT_trt_all"),
                     "PTPN2_WT_trt_IFNa" = 
                       c("siControlIFNavssiPTPN2IFNa-596", "PTPN2_WT_trt_all"))

russ_name_mapping <- list("russ_PTPN2_WT_no_trt" = 
                            c("PTPN_WT", "PTPN2_WT_noTrt_all"))

color_list <- RColorBrewer::brewer.pal(n = 3, name = "Set1")
color_list <- color_list[1:2]
names(color_list) <- c("human", "mouse")

get_intersecting_genes <- function(test_name, mapping_list, all_gene_list,
                                   excel_file, results_dir, type = "all"){
  human_de_name <- mapping_list[[test_name]][1]
  mouse_de_name <- mapping_list[[test_name]][2]
  print(human_de_name)
  human_file <- openxlsx::readWorkbook(excel_file, sheet = human_de_name)
  mouse_file <- read.csv(file.path(results_dir, paste0(mouse_de_name, ".csv")))
  
  # Upregulated genes:
  human_file_up <- human_file %>%
    dplyr::filter(LogFoldChange > 0)
  
  mouse_file_up <- mouse_file %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::filter(padj < 0.05)
  
  if(type == "russ"){
    human_gene_list_up <- all_gene_list %>%
      dplyr::filter(hgnc_symbol %in% human_file_up$`Gene-ID`)    
  } else {
    human_gene_list_up <- all_gene_list %>%
      dplyr::filter(ensembl_gene_id %in% human_file_up$`Gene-ID`)    
  }

  mouse_gene_list_up <- all_gene_list %>%
    dplyr::filter(mgi_symbol %in% mouse_file_up$gene_name)
  
  # Save list of overlapping genes and unique genes
  
  # Make venn diagrams
  venn.diagram(x = list("human" = unique(human_gene_list_up$mgi_symbol),
                        "mouse" = unique(mouse_gene_list_up$mgi_symbol)),
               fill = color_list,
               filename = file.path(save_dir, paste0(test_name, "_up", ".tiff")))
  
  
  # Down regulated genes:
  human_file_down <- human_file %>%
    dplyr::filter(LogFoldChange < 0)
  
  mouse_file_down <- mouse_file %>%
    dplyr::filter(log2FoldChange < 0) %>%
    dplyr::filter(padj < 0.05)
  
  if(type == "russ"){
    human_gene_list_down <- all_gene_list %>%
      dplyr::filter(hgnc_symbol %in% human_file_down$`Gene-ID`)    
  } else {
    human_gene_list_down <- all_gene_list %>%
      dplyr::filter(ensembl_gene_id %in% human_file_down$`Gene-ID`)    
  }

  mouse_gene_list_down <- all_gene_list %>%
    dplyr::filter(mgi_symbol %in% mouse_file_down$gene_name)
  
  # Save list of overlapping genes and unique genes
  
  # Make venn diagrams
  venn.diagram(x = list("human" = unique(human_gene_list_down$mgi_symbol),
                        "mouse" = unique(mouse_gene_list_down$mgi_symbol)),
               fill = color_list,
               filename = file.path(save_dir, paste0(test_name, "_down", ".tiff")))
  
  # All genes
  if(type == "russ"){
    human_gene_list <- all_gene_list %>%
      dplyr::filter(hgnc_symbol %in% human_file$`Gene-ID`)    
  } else {
    human_gene_list <- all_gene_list %>%
      dplyr::filter(ensembl_gene_id %in% human_file$`Gene-ID`)    
  }
  
  mouse_file_all <- mouse_file %>%
    dplyr::filter(padj < 0.05)
  
  mouse_gene_list <- all_gene_list %>%
    dplyr::filter(mgi_symbol %in% mouse_file_all$gene_name)
  
  # Save list of overlapping genes and unique genes
  
  # Make venn diagrams
  venn.diagram(x = list("human" = unique(human_gene_list$mgi_symbol),
                        "mouse" = unique(mouse_gene_list$mgi_symbol)),
               fill = color_list,
               filename = file.path(save_dir, paste0(test_name, "_all", ".tiff")))
  
  
  # Hypergeometric tests
  up_mouse <- length(unique(mouse_gene_list_up$mgi_symbol))
  up_human <- length(unique(human_gene_list_up$mgi_symbol))
  down_mouse <- length(unique(mouse_gene_list_down$mgi_symbol))
  down_human <- length(unique(human_gene_list_down$mgi_symbol))
  all_mouse <- length(unique(mouse_gene_list$mgi_symbol))
  all_human <- length(unique(human_gene_list$mgi_symbol))
  all_genes <- length(unique(all_gene_list$mgi_symbol))
  
  # x is the number of overlaps
  # k is the largest possible number of overlaps
  up_intersects <- length(intersect(unique(mouse_gene_list_up$mgi_symbol),
                                    unique(human_gene_list_up$mgi_symbol)))
  
  down_intersects <- length(intersect(unique(mouse_gene_list_down$mgi_symbol),
                                      unique(human_gene_list_down$mgi_symbol)))
  
  all_intersects <- length(intersect(unique(mouse_gene_list$mgi_symbol),
                                     unique(human_gene_list$mgi_symbol)))
  
  gene_ontology_genes <- intersect(unique(mouse_gene_list$mgi_symbol),
                                   unique(human_gene_list$mgi_symbol))
  
  save_dir_gost <- file.path(save_dir, test_name)
  ifelse(!dir.exists(save_dir_gost), dir.create(save_dir_gost), FALSE)
  gene_ontology_res <- gost(query = list("intersecting" = gene_ontology_genes),
                            organism = "mmusculus")
  
  save_gost(gost_output = gene_ontology_res, 
            save_dir_text = save_dir_gost,
            save_text = TRUE, save_excel = TRUE)
  # x is number of intersects to the total number of possible overlaps
  # m is the total genes in the human list (these would be the successes)
  # n is the total genes not in the human list - gene list - human list
  # k is the total tests or total mouse de genes
  
  
  hyper_up <- sum(dhyper(x = up_intersects:up_mouse, m = up_human, 
                         n = all_genes - up_human, k = up_mouse))
  
  hyper_down <- sum(dhyper(x = down_intersects:down_mouse, m = down_human, 
                           n = all_genes - down_human, k = down_mouse))
  
  hyper_all <- sum(dhyper(x = all_intersects:all_mouse, m = all_human, 
                          n = all_genes - all_human, k = all_mouse))
  
  # Calculated expected number of genes
  expected_num_up <- (up_human*up_mouse)/all_genes
  expected_num_down <- (down_human*down_mouse)/all_genes
  expected_num_all <- (all_human*all_mouse)/all_genes
  
  # Calcluate the representation factor
  representation_up <- up_intersects/expected_num_up
  representation_down <- down_intersects/expected_num_down
  representation_all <- all_intersects/expected_num_all
  
  return_df <- data.frame(test = test_name, gene_list = c("up", "down", "all"),
                          mouse_genes = c(up_mouse, down_mouse, all_mouse),
                          human_genes = c(up_human, down_human, all_human),
                          intersecting_genes = c(up_intersects, down_intersects,
                                                 all_intersects),
                          p_value = c(hyper_up, hyper_down, hyper_all),
                          expected_vals = c(expected_num_up, expected_num_down,
                                            expected_num_all),
                          overrepresentation = c(representation_up,
                                                 representation_down,
                                                 representation_all))
  
  return(return_df)
}

save_gost <- function(gost_output, save_dir_text,
                      save_text = TRUE, save_excel = TRUE){

  gost_text <- gost_output$result
  gost_text$parents <- NULL
  # Save all results to csv
  if(save_text){
    # Create output directory
    save_dir_text %>%
      dir.create(showWarnings = F)
    utils::write.csv(gost_text,
                     file.path(save_dir_text, "all_GSE_results.csv"))
  }
  if(save_excel){
    # Create output directory
    save_dir_text %>%
      dir.create(showWarnings = F)
    invisible(lapply(unique(gost_text$query), function(gost_query){
      print(gost_query)
      
      # Create excel wb
      gene_wb <- openxlsx::createWorkbook()
      
      # Write to excel wb
      full_list <- lapply(unique(gost_text$source), function(gost_source){
        new_df <- gost_text %>%
          dplyr::filter(source == gost_source & query == gost_query)
        gost_source_write <- sub(":", "_", gost_source)
        openxlsx::addWorksheet(gene_wb, gost_source_write)
        openxlsx::writeData(gene_wb, gost_source_write, new_df)
      })
      
      ## Save workbook to working directory
      openxlsx::saveWorkbook(gene_wb,
                             file = file.path(save_dir_text, paste0(gost_query,
                                                                    "GSE_results.xlsx")),
                             overwrite = TRUE)
    }))
  }
}


all_res <- lapply(names(name_mapping), function(test_name){
  save_df <- get_intersecting_genes(test_name = test_name,
                                    mapping_list = name_mapping,
                                    all_gene_list = all_gene_list,
                                    excel_file = excel_file,
                                    results_dir = results_dir)
  return(save_df)
})


all_res <- do.call(rbind, all_res)

new_res <- get_intersecting_genes(test_name = "russ_PTPN2_WT_no_trt",
                                  mapping_list = russ_name_mapping,
                                  all_gene_list = all_gene_list,
                                  excel_file = russ_file,
                                  results_dir = results_dir,
                                  type = "russ")

all_res <- rbind(all_res, new_res)

all_res$padj <- p.adjust(all_res$p_value, method = "fdr")

sig_res <- all_res %>%
  dplyr::filter(padj < 0.05)

save_data <- openxlsx::createWorkbook()
openxlsx::addWorksheet(save_data, sheetName = "all_overlaps")
openxlsx::writeData(wb = save_data, sheet = "all_overlaps", x = all_res)

openxlsx::addWorksheet(save_data, sheetName = "sig_overlaps")
openxlsx::writeData(wb = save_data, sheet = "sig_overlaps", x = sig_res)

openxlsx::saveWorkbook(wb = save_data, 
                       file = here("results/R_analysis", 
                                   "DE_files/human_hypergometric.xlsx"),
                       overwrite = TRUE)

# Quick check that the directionality is the same ------------------------------
# Looks like a positive log fold change is higher in PTPN2 than WT,
# Negative is higher in WT than PTPN2

## Check theirs ----------------------------------------------------------------
human_de_name <- name_mapping[[1]][1]
human_file <- openxlsx::readWorkbook(excel_file, sheet = human_de_name)

human_file %>% 
  dplyr::select(dplyr::contains("Ctl"), LogFoldChange) %>%
  dplyr::arrange(LogFoldChange) %>% 
  head

only_cpm <- human_file %>%
  dplyr::select(dplyr::contains("CPM"), LogFoldChange, `Gene-ID`)

col_mapping <- c("1-EndoC-E1-SiQ-Ctl.1.CPM" = "control_control_1",
                  "7-EndoC-E1-Val-SiQ-Ctl.7.CPM" = "control_control_2",
                  "13-EndoC-E2-Val-SiQ-Ctl.13.CPM" = "control_control_3",
                  "2-EndoC-E1-SiQ-IFNa.2.CPM" = "control_infa_1",
                  "8-EndoC-E1-Val-SiQ-IFNa.8.CPM" = "control_infa_2",
                  "14-EndoC-E2-Val-SiQ-IFNa.14.CPM" = "control_infa_3",
                  "3-EndoC-E1-SiQ-IFNY.3.CPM" = "control_infg_1",
                  "9-EndoC-E1-Val-SiQ-IFNY.9.CPM" = "control_infg_2",
                  "15-EndoC-E2-Val-SiQ-IFNY.15.CPM" = "control_infg_3",
                  "4-EndoC-E1-10-Ctl.4.CPM" = "ptpn2_control_1",
                  "10-EndoC-E1-Val-SiPTPN2-Ctl.10.CPM" = "ptpn2_control_2",
                  "16-EndoC-E2-Val-SiPTPN2-Ctl-.16.CPM" = "ptpn2_control_3",
                  "5-EndoC-E1-10-IFNa.5.CPM" = "ptpn2_infa_1",
                  "11-EndoC-E1-Val-SiPTPN2-IFNa.11.CPM" = "ptpn2_infa_2",
                  "17-EndoC-E2-Val-SiPTPN2-IFNa.17.CPM" = "ptpn2_infa_3",
                  "6-EndoC-E1-10-IFNY.6.CPM" = "ptpn2_infg_1",
                  "12-EndoC-E1-Val-SiPTPN2-IFNY.12.CPM" = "ptpn2_infg_2",
                  "18-EndoC-E2-Val-SiPTPN2-IFNY.18.CPM" = "ptpn2_infg_3",
                  "LogFoldChange" = "LogFoldChange",
                 "Gene-ID" = "Gene-ID")

new_names <- col_mapping[colnames(only_cpm)]
colnames(only_cpm) <- new_names

only_cpm %>%
  dplyr::arrange(LogFoldChange) %>%
  head

only_cpm %>%
  dplyr::arrange(desc(LogFoldChange)) %>%
  head

plot_df <- only_cpm %>%
  tibble::column_to_rownames("Gene-ID") %>%
  dplyr::mutate(LogFoldChange = NULL) %>%
  t %>%
  data.frame %>%
  dplyr::mutate(genotype = gsub("_.*", "", rownames(.)),
                treatment = gsub("^.*_([[:alnum:]]+)_.*$", "\\1", rownames(.)),
                genotype_treatment = paste(genotype, treatment, sep = "_"))

ggplot2::ggplot(plot_df, ggplot2::aes(x = genotype_treatment,
                                      y = ENSG00000091137,
                                      color = genotype_treatment)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_brewer(palette = "Set1") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Check ours ------------------------------------------------------------------
mouse_de_name <- name_mapping[[1]][2]
mouse_file <- read.csv(file.path(results_dir, paste0(mouse_de_name, ".csv")))

mouse_file %>%
  dplyr::arrange(log2FoldChange) %>%
  head

mouse_file %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  head

plotCounts(mouse_res, gene = "Serpina3n_ENSMUSG00000021091", intgroup = "group")
