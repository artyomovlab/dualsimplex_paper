library(TCGAbiolinks)
library(dplyr)

# TCGA-HNSC, NCICCR-DLBCL, TCGA-KIRC, TCGA-PAAD
project <- "TCGA-HNSC"
subproject_only <- str_split(project, "-")[[1]][2]
counts_colname <- "tpm_unstranded"

genex_query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
)

res <- getResults(genex_query)
GDCdownload(genex_query)

base_path <- file.path(
  getwd(),
  "GDCdata",
  genex_query$project[[1]],
  "harmonized",
  "Transcriptome_Profiling",
  "Gene_Expression_Quantification"
)


all_counts <- list()

for (i in 1:nrow(res)) {
  fp <- file.path(base_path, res[i, "id"], res[i, "file_name"])
  
  # Read table
  df <- as.data.frame(data.table::fread(
    fp,
    header = T,
    select = c("gene_name", counts_colname, "gene_type")
  )[-c(1:4),])
  # TODO: keep max expressed of duplicated gene_names
  
  # Record into final list
  counts <- df[, counts_colname]
  names(counts) <- df$gene_name
  all_counts[[res[i, "cases"]]] <- counts
}

all_counts_matrix <- do.call(cbind, all_counts)

# cat(rownames(all_counts_matrix), file = "gene_annotation_lists/CODING.txt", sep=' ')

saveRDS(all_counts_matrix, file.path(subproject_only, "gene_counts.Rda"))

