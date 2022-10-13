library(here)
library(gprofiler2)
library(tidyverse)
library(openxlsx)
library(viridis)

de_path <- here("results/R_analysis/DE_files")

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source("src/rules/scripts/functions.R")

save_dir <- here("results/R_analysis/images/GSE")

count_file <- read.table(here("results/PTPN_KO_countTable_cutadapt_trim.txt"),
                         sep = "\t", header = TRUE)

gene_mapping <- count_file %>%
  dplyr::select(gene) %>%
  dplyr::mutate(ens_id = gsub(".*_ENSMUSG", "ENSMUSG", gene),
                gene_id = gsub("_ENSMUSG.*", "", gene))

gene_mapping$gene <- NULL

gost_res <- readRDS(file = here("results/R_analysis/GSE_files/PTPN2_WT_Trt.rds"))

results <- gost_res$result


results_short <- results %>%
  dplyr::select(!c(intersection, evidence_codes))

results_short <- results_short %>%
  dplyr::mutate(expected_num = (term_size * query_size) /
                  effective_domain_size) %>%
  dplyr::mutate(overrepresentation = intersection_size / expected_num)


p_val_res <- lapply(1:nrow(results_short), function(x){
  one_row <- results_short[x,]
  p_val <- sum(dhyper(x = one_row$intersection_size:one_row$query_size,
                      m = one_row$term_size,
                      n = one_row$effective_domain_size - one_row$term_size,
                      k = one_row$query_size))
  one_row$kristen_p <- p_val
  
  return(one_row)
  
  
})

p_val_res <- do.call(rbind, p_val_res)

p_val_res$kristen_adj_p <- p.adjust(p_val_res$kristen_p,
                                    method = "fdr")

p_val_res <- p_val_res %>%
  dplyr::mutate(log_p = -log10(p_value))

# The p-values here are close, but I'm currently only doing the correction on 
# significant hits. I may want to eventually double check this against all
pdf(file.path(save_dir, "gene_ontology_example.pdf"))

invisible(lapply(unique(p_val_res$source), function(x){
  res_small <- p_val_res %>%
    dplyr::filter(source == x & p_value < 0.05) %>%
    dplyr::top_n(10, wt = log_p) %>%
    dplyr::arrange(log_p)
  
  res_small$term_name <- factor(res_small$term_name,
                                levels = res_small$term_name)
  
  
  plot1 <- ggplot2::ggplot(res_small, ggplot2::aes(x = log_p,
                                              y = term_name,
                                              size = intersection_size,
                                              color = overrepresentation)) +
    ggplot2::geom_point() +
    scale_color_continuous(low="blue", high="red",
                           guide=guide_colorbar(reverse=FALSE)) +
    ggplot2::ggtitle(paste0(x, "_x_log"))
  
  
  
  res_small <- res_small %>%
    dplyr::arrange(overrepresentation)
  
  res_small$term_name <- factor(res_small$term_name,
                                    levels = res_small$term_name)
  
  
  plot2 <- ggplot2::ggplot(res_small, ggplot2::aes(color = log_p,
                                              y = term_name,
                                              size = intersection_size,
                                              x = overrepresentation)) +
    ggplot2::geom_point() +
    scale_color_continuous(low="blue", high="red",
                           guide=guide_colorbar(reverse=FALSE)) +
    ggplot2::ggtitle(paste0(x, "_x_overrepresentation"))
  
  print(plot1)
  print(plot2)
}))

dev.off()
