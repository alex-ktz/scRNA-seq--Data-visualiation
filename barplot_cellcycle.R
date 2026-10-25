library(ggplot2)
library(RColorBrewer)

#file path
data_dir <- "C:/Users/katza/Uni/R/project"
out_dir  <- "C:/Users/katza/Uni/R/project/image"

#Define file paths and group labels 
file_info <- data.frame(
  filename = c("0.csv",
               "3_rep1.csv", "3_rep2.csv",
               "7_rep1.csv", "7_rep2.csv",
               "14_rep1_low.csv", "14_rep2_low.csv",
               "14_rep1_med.csv", "14_rep2_med.csv",
               "14_rep1_high.csv", "14_rep2_high.csv"),
  
  group = c("Day0",
            "Day3_rep1", "Day3_rep2",
            "Day7_rep1", "Day7_rep2",
            "Day14_cycling_rep1", "Day14_cycling_rep2",
            "Day14_mod_rep1", "Day14_mod_rep2",
            "Day14_noncycling_rep1", "Day14_noncycling_rep2"),
  
  avg_group = c("Day 0",
                "Day 3", "Day 3",
                "Day 7", "Day 7",
                "Day 14\nCycling", "Day 14\nCycling",
                "Day 14\nModerate", "Day 14\nModerate",
                "Day 14\nNon-cycling", "Day 14\nNon-cycling"),
  
  stringsAsFactors = FALSE
)


#Define the 5 cell cycle target genes
target_genes <- c("MKI67", "TOP2A", "CDK1", "CCNB1", "PCNA")


mean_expr_list <- list()

#Read each file one by one
for (i in seq_len(nrow(file_info))) {
  fpath <- file.path(data_dir, file_info$filename[i])
  cat("  Processing:", file_info$group[i], "...\n")

  #take only the header
  header <- read.csv(fpath, nrows = 0, check.names = FALSE)
  all_genes <- colnames(header)
  
  #take only the relevant genes
  gene_cols <- which(all_genes %in% target_genes)
  col_classes <- rep("NULL", length(all_genes))
  col_classes[gene_cols] <- "numeric"
  col_classes[1] <- "character"

  #Make a df table only with the relevant gene
  df <- read.csv(fpath, colClasses = col_classes, check.names = FALSE)
  gene_data <- df[, -1, drop = FALSE]

  #Mean UMI count per gene, then log to have value on the same scale
  gene_means <- log1p(colMeans(gene_data, na.rm = TRUE))

  #insert the mean in a list
  mean_expr_list[[file_info$group[i]]] <- gene_means
}


#Making a small matrix from the data
expr_matrix <- do.call(cbind, mean_expr_list)
expr_matrix <- expr_matrix[target_genes, ]

#Average replicates into 6 groups
avg_groups <- unique(file_info$avg_group)
avg_matrix <- matrix(NA, nrow = length(target_genes), ncol = length(avg_groups),
                     dimnames = list(target_genes, avg_groups))

#calculate the mean for each one of the genes
for (g in avg_groups) {
  idx <- which(file_info$avg_group == g)
  cols <- file_info$group[idx]
  if (length(cols) == 1) {
    avg_matrix[, g] <- expr_matrix[, cols]
  } else {
    avg_matrix[, g] <- rowMeans(expr_matrix[, cols])
  }
}

#Convert to long format
plot_data <- data.frame(
  gene       = rep(rownames(avg_matrix), times = ncol(avg_matrix)),
  group      = rep(colnames(avg_matrix), each = nrow(avg_matrix)),
  expression = as.vector(avg_matrix),
  stringsAsFactors = FALSE
)

# Set group as an ordered factor
plot_data$group <- factor(plot_data$group, levels = avg_groups)

#Choosing the color palette
gene_colors <- brewer.pal(5, "Set1")

#Making the ggplot
p <- ggplot(plot_data, aes(x = group, y = expression, fill = gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = gene_colors) +
  labs(
    title = "Cell Cycle Gene Expression Across Different Days",
    x     = "Treatment Group",
    y     = "Mean Expression (log)",
    fill  = "Gene"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    plot.title          = element_text(face = "bold", size = 14),
    legend.position     = "right",
    panel.grid.major.x  = element_blank()
  )

#saving as a png
out_file <- file.path(out_dir, "explanatory.png")

ggsave(out_file, plot = p, width = 12, height = 6, dpi = 300)

cat("\nBar plot saved to:", out_file, "\n")
