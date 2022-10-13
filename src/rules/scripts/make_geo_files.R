library(DESeq2)
library(here)

dds_obj <- readRDS(here("results/R_analysis/objs/dds_obj.rda"))

all_counts <- read.table(here("results/PTPN_KO_countTable_cutadapt_trim.txt"))
filtered_counts <- DESeq2::counts(dds_obj)

all_metadata <- colData(dds_obj)

write.table(filtered_counts, file = here("manuscript/filtered_counts.tsv"),
            sep = "\t")

write.table(all_metadata, file = here("manuscript/metadata.tsv"),
            sep = "\t")
