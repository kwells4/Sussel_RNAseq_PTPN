library(DESeq2)
library(openxlsx)

base_dir =
  "/Users/wellskr/Documents/Analysis/Lori_sussel/Young_Kim/Sussel_RNASeq_PTPN/"

cytokine_treated_gene_list <- "files/WT_vs_KO_cytokine_RNA_Seq_categorized.xlsx"
untreated_gene_list <- "files/WT_vs_KO_Non_treated_RNA_Seq_categorized.xlsx"

# Function to make a heatmap
make_heatmap <- function(dds, vsd, de_genes, treatment, control, group, 
                         print_genenames = FALSE, gene_identifier = "Gene_ID",
                         cluster_cols = FALSE, save_heatmap = TRUE,
                         output_dir = "/results/", plot_groups = "all",
                         color_test = NULL){
  # First make a "col data" data frame. This is just telling
  # pheatmap what you want to use to color the columns.
  df <- as.data.frame(colData(dds)[,c(group)])
  sample_info <- colData(dds)
  colnames(df) <- group
  rownames(df) <- rownames(colData(dds))
  
  # grab genes
  if(isS4(de_genes)){
    genes <- rownames(de_genes)
  } else if (is.character(de_genes)){
    genes <- de_genes
  } else {
    stop("de_genes must be a S4 output from DESeq2 (res$DE_genes) or a list of genes")
  }
  if(is.null(color_test)){
    color_test <- brewer.pal(length(levels(dds[[group]])), "Set1")
    names(color_test) <- levels(dds[[group]])
  }
  coloring <- list(color_test)
  names(coloring) <- group
  heatmap_df <- assay(vsd)[genes,]
  if(!("all" %in% plot_groups)){
    sample_info <- colData(dds)
    sample_plot <- sample_info[sample_info[[group]] %in% plot_groups, ]
    heatmap_df <- heatmap_df[ , rownames(sample_plot)]
  } else {
    sample_plot <- sample_info
  }
  # Heatmaps are always mean cenetered. Meaning the value shown is actually the 
  # expression value of the sample minus the mean. This command centers the
  # dataframe so it can be plotted. When it isn't centered it is very hard to
  # see any trends because all genes have very different expression levels.
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  palOut <- colorRampPalette(blueYellow)(256)
  if(!cluster_cols){
    if(!identical(colnames(heatmap_scale), rownames(sample_plot))){
      heatmap_scale <- heatmap_scale[ , rownames(sample_plot)]
    }
    # Order based on the comparison
    col_order <- c(
      grep(control, sample_plot[[group]]),
      grep(treatment, sample_plot[[group]]),
      grep((paste0(treatment, "|", control)), sample_plot[[group]],
           invert = TRUE)
    )
    heatmap_scale <- heatmap_scale[ , col_order]
  }
  if(print_genenames){
    gene_names <- rownames(heatmap_scale)
    if(gene_identifier == "Gene_ID"){
      gene_names <- sub("_ENSMUSG[0-9]*", "", gene_names)
    } else if (gene_identifier == "ENS_ID"){
      gene_names <- sub("*_ENSMUSG", "ENSMUSG", gene_names)
    }
    rownames(heatmap_scale) <- gene_names
  }
  # Here we make the heatmap
  heatmap <- pheatmap(heatmap_scale, cluster_rows = TRUE,
                      cluster_cols = cluster_cols,
                      show_rownames = print_genenames,
                      show_colnames = TRUE, annotation_col = df,
                      annotation_colors = coloring, color = blueYellow,
                      border_color = NA, clustering_method = "complete",
                      silent = TRUE)
  
  if(save_heatmap){
    comparison <- paste0(control, "_vs_", treatment, ".pdf")
    pdf(file.path(output_dir, comparison), width = 10, height = 10)
    print(heatmap)
    dev.off()
  }
  
  return(heatmap)
}

dds <- readRDS(paste0(base_dir, "results/R_analysis/objs/dds_obj.rda"))
vsd <- readRDS(paste0(base_dir, "results/R_analysis/objs/vsd_obj.rda"))

control <- "WT_NoTreatment"
treatment <- "PTPN2_KO_NoTreatment"

gene_lists <- openxlsx::getSheetNames(paste0(base_dir,
                                             untreated_gene_list))
gene_file <- paste0(base_dir, untreated_gene_list)

gene_list <- "Immune response related genes"
gene_table <- read.xlsx(gene_file, sheet = gene_list)

gene_ids <- gene_table$ens_id
all_gene_ids <- lapply(gene_table$ens_id, function(x){
  if(!is.na(x)){
    gene_id <- rownames(assay(vsd))[grepl(x, rownames(assay(vsd)))]
  }
})

all_gene_ids <- unlist(all_gene_ids)

heatmap_1 <- make_heatmap(dds = dds,
                          vsd = vsd,
                          de_genes = all_gene_ids,
                          treatment = treatment,
                          control = control,
                          group = params$sample_column,
                          print_genenames = TRUE,
                          cluster_cols = FALSE,
                          save_heatmap = FALSE,
                          color_test = sample_colors)
print(heatmap_1)


heatmap_1 <- make_heatmap(dds = dds,
                          vsd = vsd,
                          de_genes = all_gene_ids,
                          treatment = treatment,
                          control = control,
                          group = "group",
                          print_genenames = TRUE,
                          cluster_cols = TRUE,
                          save_heatmap = FALSE,
                          color_test = sample_colors)
