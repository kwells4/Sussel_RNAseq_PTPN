library(here)
library(DESeq2)
library(gprofiler2)
library(tidyverse)

source("src/rules/scripts/functions.R")

dds_obj <- readRDS(here("results/R_analysis/objs/dds_obj.rda"))
de_res <- read.table(here("results/R_analysis/DE_files/PTPN2_WT_Trt.csv"),
                     header = TRUE, sep = ",")

count_file <- read.table(here("results/PTPN_KO_countTable_cutadapt_trim.txt"),
                         sep = "\t", header = TRUE)

gene_mapping <- count_file %>%
  select(gene) %>%
  mutate(ens_id = gsub(".*_ENSMUSG", "ENSMUSG", gene),
         gene_id = gsub("_ENSMUSG.*", "", gene))

gene_mapping$gene <- NULL


gost_res_trt <- readRDS(
  here("results/R_analysis/GSE_files_separated/PTPN2_KO_Cytokine_Treatment_from_PTPN2_KO_Cytokine_Treatment_vs_WT_Cytokine_Treatment.rds"))
gost_res_ctl <- readRDS(
  here("results/R_analysis/GSE_files_separated/WT_Cytokine_Treatment_from_PTPN2_KO_Cytokine_Treatment_vs_WT_Cytokine_Treatment.rds"))

all_gene_ens <- sub(".*_ENSMUSG", "ENSMUSG", rownames(dds_obj))

treatment <- "PTPN2_KO_Cytokine_Treatment"
control <- "WT_Cytokine_Treatment"

result <- run_gprofiler(gene_table = de_res,
                        pos_name = treatment,
                        neg_name = control,
                        custom_bg = all_gene_ens,
                        evcodes = TRUE)

plots <- result$plots
gost_table_trt <- result[[treatment]]$result
gost_table_ctl <- result[[control]]$result

gost_table_trt_orig <- gost_res_trt$result
gost_table_ctl_orig <- gost_res_ctl$result

gost_table_trt_short <- gost_table_trt %>%
  select(!c(evidence_codes, intersection))

gost_table_ctl_short <- gost_table_ctl %>%
  select(!c(evidence_codes, intersection))

# Make sure this is the same result from yesterday
identical(gost_table_trt_orig, gost_table_trt_short)
identical(gost_table_ctl_orig, gost_table_ctl_short)

# Find overlapping genes
results_go_wp <- gost_table_trt_short %>%
  dplyr::filter(source == "WP" & p_value < 0.05)

interesting_terms <- c("Oxidative stress and redox pathway")



gene_overlaps <- openxlsx::createWorkbook()

organelle_genes <- lapply(interesting_terms, function(x){
  single_term <- gost_table_trt %>%
    dplyr::filter(source == "WP" & term_name == x)
  
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
                       file = file.path(here("results/R_analysis/GSE_files_separated/",
                                             "PTPN2_WT_trt_WP_overlapping_genes.xlsx")),
                       overwrite = TRUE)

