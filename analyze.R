# Load libraries
library(limma)
library(ggplot2)
library(readr)
library(dplyr)

# --- USER-DEFINED VARIABLES ---
# IMPORTANT: Edit these to match your experimental design!
PROTEIN_FILE <- "peptideshaker_output/protein_reports.tsv"
# Define which columns belong to which group (by column number)
CONTROL_COLS <- c(7, 8)      # e.g., columns 7 and 8 are control replicates
TREATMENT_COLS <- c(9, 10)   # e.g., columns 9 and 10 are treatment replicates

# --- 1. Read and Prepare Data ---
data <- read_tsv(PROTEIN_FILE)

# Extract quantitative columns and log2 transform
quant_cols <- colnames(data)[c(CONTROL_COLS, TREATMENT_COLS)]
log_data <- data[, quant_cols] %>% log2()

# Replace -Inf from log2(0) with NA for imputation
log_data[is.infinite(log_data)] <- NA

# Simple imputation: replace missing values with a small value
min_val <- min(log_data, na.rm = TRUE)
log_data[is.na(log_data)] <- min_val - 1

# --- 2. Statistical Analysis with limma ---
# Create a design matrix based on your groups
groups <- factor(c(rep("Control", length(CONTROL_COLS)), rep("Treatment", length(TREATMENT_COLS))))
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit the linear model and calculate statistics
fit <- lmFit(log_data, design)
contrast.matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract all results
results <- topTable(fit2, adjust = "fdr", number = Inf)
# Add protein names back to the results
results$Protein <- data$Protein 
write_csv(results, "significant_proteins.csv")

# --- 3. Create Volcano Plot ---
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.6, aes(color = (adj.P.Val < 0.05 & abs(logFC) > 1))) +
  scale_color_manual(values = c("grey70", "red"), name = "Significant") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot – Treatment vs Control",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  )

# Save the plot
ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)
cat("Analysis complete. Results saved to significant_proteins.csv and volcano_plot.png\n")
