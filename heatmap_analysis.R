library(pheatmap)
library(RColorBrewer)

#Input and output files
data_dir <- "C:/Users/katza/Uni/R/project"
out_dir  <- "C:/Users/katza/Uni/R/project/image"

#Define file paths and group labels for easy classification in plot
file_info <- data.frame(
  filename = c("0.csv",
               "3_rep1.csv", "3_rep2.csv",
               "7_rep1.csv", "7_rep2.csv",
               "14_rep1_low.csv", "14_rep2_low.csv",
               "14_rep1_med.csv", "14_rep2_med.csv",
               "14_rep1_high.csv", "14_rep2_high.csv"),
  
  label = c("D0",
            "D3 R1", "D3 R2",
            "D7 R1", "D7 R2",
            "D14 Cyc R1", "D14 Cyc R2",
            "D14 Mod R1", "D14 Mod R2",
            "D14 Non-cyc R1", "D14 Non-cyc R2"),
  
  timepoint = c("Day 0",
                "Day 3", "Day 3",
                "Day 7", "Day 7",
                "Day 14", "Day 14",
                "Day 14", "Day 14",
                "Day 14", "Day 14"),
  
  avg_group = c("Day 0",
                "Day 3", "Day 3",
                "Day 7", "Day 7",
                "Day 14\nCycling", "Day 14\nCycling",
                "Day 14\nModerate", "Day 14\nModerate",
                "Day 14\nNon-cycling", "Day 14\nNon-cycling"),
  
  stringsAsFactors = FALSE
)


#Define the 15 target genes, grouped by category
gene_groups <- list(
  "Cell Cycle"      = c("MKI67", "TOP2A", "CDK1", "CCNB1", "PCNA"),
  "Stress/Survival" = c("HSPA1A", "HSP90AA1", "VIM", "HMGA1"),
  "Housekeeping"    = c("ACTB", "GAPDH"),
  "Apoptosis"       = c("BAX", "BCL2"),
  "Drug Resistance" = c("EGFR", "AXL")
)

target_genes <- unlist(gene_groups, use.names = FALSE)# all of the genes in one vector

#Make a small dataframe to link each gene to his category, to make the classification in the plot easier
gene_categories <- data.frame(
  Category = rep(names(gene_groups), lengths(gene_groups)),
  row.names = target_genes
)

row_gaps <- cumsum(lengths(gene_groups))
row_gaps <- row_gaps[-length(row_gaps)]

mean_expr_list <- list()

#Read each one of the files and print wich file is in process
for (i in seq_len(nrow(file_info))) {
  fpath <- file.path(data_dir, file_info$filename[i])
  cat("  Processing:", file_info$label[i], "...\n")
  
  #Read only the header to process only the relevant genes
  header <- read.csv(fpath, nrows = 0, check.names = FALSE)
  all_genes <- colnames(header)

  gene_cols <- which(all_genes %in% target_genes)
  col_classes <- rep("NULL", length(all_genes))
  col_classes[gene_cols] <- "numeric"
  col_classes[1] <- "character"

  df <- read.csv(fpath, colClasses = col_classes, check.names = FALSE)
  gene_data <- df[, -1, drop = FALSE]

  #Calculate the mean appearance of each genes for the heatmap 
  mean_expr_list[[file_info$label[i]]] <- colMeans(gene_data, na.rm = TRUE)
}

#Build a matrix with the value
expr_matrix <- do.call(cbind, mean_expr_list)
expr_matrix <- expr_matrix[target_genes, ]

#log the value to avoid the difference between value to be huge
log_matrix <- log1p(expr_matrix)



#Build averaged matrix for the days that have multiple sample
avg_groups <- unique(file_info$avg_group)
avg_matrix <- matrix(NA, nrow = length(target_genes), ncol = length(avg_groups),
                     dimnames = list(target_genes, avg_groups))

for (g in avg_groups) {
  idx <- which(file_info$avg_group == g)
  cols <- file_info$label[idx]
  if (length(cols) == 1) {
    avg_matrix[, g] <- log_matrix[, cols]
  } else {
    avg_matrix[, g] <- rowMeans(log_matrix[, cols])
  }
}

#Creation of the color palete for the heatmap
heatmap_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)#color gradient

category_colors <- list(
  Category  = c("Cell Cycle"      = "#E41A1C",
                "Stress/Survival" = "#FF7F00",
                "Housekeeping"    = "#4DAF4A",
                "Apoptosis"       = "#377EB8",
                "Drug Resistance" = "#984EA3"),
  
  Timepoint = c("Day 0"  = "#FFFFB2",
                "Day 3"  = "#FED976",
                "Day 7"  = "#FE9666",
                "Day 14" = "#E31A1C")
)

#Column annotation: timepoint only (no status bar)
avg_annot <- data.frame(
  Timepoint = c("Day 0", "Day 3", "Day 7", "Day 14", "Day 14", "Day 14"),
  row.names = avg_groups
)

col_gaps_avg <- c(1, 2, 3)  #gaps after Day 0, Day 3, Day 7

#Create heatmap directly in file
out_avg <- file.path(out_dir, "heatmap_averaged.png")
png(out_avg, width = 9, height = 8, units = "in", res = 300)

#Heatmap settings
pheatmap(
  avg_matrix,
  main              = "PC9 Gene Expression Under Osimertinib",
  color             = heatmap_colors,
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  gaps_row          = row_gaps,
  gaps_col          = col_gaps_avg,
  annotation_row    = gene_categories,
  annotation_col    = avg_annot,
  annotation_colors = category_colors,
  fontsize          = 11,
  fontsize_row      = 12,
  fontsize_col      = 11,
  angle_col         = 45,
  border_color      = NA,
  cellwidth         = 45,
  cellheight        = 22
)

#Finish the program
dev.off()
cat("Heatmap saved to:", out_avg, "\n")
