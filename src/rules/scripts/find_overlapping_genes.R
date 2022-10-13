library(here)
library(gprofiler2)
library(tidyverse)
library(openxlsx)

de_path <- here("results/R_analysis/DE_files")

source("src/rules/scripts/functions.R")

count_file <- read.table(here("results/PTPN_KO_countTable_cutadapt_trim.txt"),
                         sep = "\t", header = TRUE)

gene_mapping <- count_file %>%
  select(gene) %>%
  mutate(ens_id = gsub(".*_ENSMUSG", "ENSMUSG", gene),
         gene_id = gsub("_ENSMUSG.*", "", gene))

gene_mapping$gene <- NULL

gost_res <- readRDS(file = here("results/R_analysis/GSE_files/PTPN2_WT_Trt.rds"))

results <- gost_res$result

results_short <- results %>%
  dplyr::select(!c(intersection, evidence_codes))

results_go_cc <- results_short %>%
  dplyr::filter(source == "GO:CC" & p_value < 0.05)

organelle_terms <- results_go_cc$term_name[grepl("organelle",
                                                 results_go_cc$term_name)]



gene_overlaps <- openxlsx::createWorkbook()

organelle_genes <- lapply(organelle_terms, function(x){
  single_term <- results %>%
    dplyr::filter(source == "GO:CC" & term_name == x)

  term_id <- gsub(":", "_", single_term$term_id)
  #print("")
  #print(x)
  ens_ids <- str_split(single_term$intersection, pattern = ",")[[1]]
  #print(length(ens_ids))
  #print(ens_ids)
  intersecting_genes <- data.frame(ens_id = ens_ids, term_id = term_id,
                                   term_name = x) %>%
    dplyr::left_join(gene_mapping, by = "ens_id")
  
  #print(nrow(intersecting_genes))
  #print(identical(ens_ids, intersecting_genes$ens_id))
  
  openxlsx::addWorksheet(wb = gene_overlaps, sheetName = term_id)
  openxlsx::writeData(wb = gene_overlaps, sheet = term_id,
                      x = intersecting_genes)
  
  return(intersecting_genes)
})

openxlsx::saveWorkbook(wb = gene_overlaps,
                       file = file.path(here("results/R_analysis/GSE_files",
                       "PTPN2_WT_trt_overlapping_genes.xlsx")),
                       overwrite = TRUE)
