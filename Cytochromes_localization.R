#!/usr/bin/env Rscript

# ============================================================================
# Numerical Heatmap Generator for Bicarbonate Proteomics Data
# Uses Excel file with ProteinID, Description, Localization, Control, Test columns
# ============================================================================

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,      # Data manipulation
  readxl,         # Read Excel files
  ComplexHeatmap, # Primary heatmap package
  circlize,       # Color functions
  RColorBrewer,   # Color palettes
  grid            # For drawing separators
)

# ===== CONFIGURATION =========================================================

INPUT_FILE <- '/Users/sysoevm/Library/CloudStorage/Dropbox/Pascal Lab/Desulfuromonas/Bicarbonate Chapter/Stats DIA-NN/heatmaps/Cytochromes_regulation.xlsx'
OUTPUT_DIR <- "/Users/sysoevm/Library/CloudStorage/Dropbox/Pascal Lab/Desulfuromonas/Bicarbonate Chapter/Stats DIA-NN/heatmaps"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Color scheme - matching your original aesthetic
REGULATION_COLORS <- c(
  "Control" = "#377eb8",  # Blue - with acetate+bicarbonate
  "Test" = "#e41a1c"       # Red - without acetate
)

# Localization order (top to bottom in heatmap)
LOCALIZATION_ORDER <- c(
  "Cytoplasmic",
  "Cytoplasmic membrane", 
  "Periplasmic",
  "Outer membrane",
  "Extracellular"
)

# Localization abbreviations for row titles
LOCALIZATION_ABBREV <- c(
  "Cytoplasmic" = "C",
  "Cytoplasmic membrane" = "CM", 
  "Periplasmic" = "P",
  "Outer membrane" = "OM",
  "Extracellular" = "EC"
)

# ===== HELPER FUNCTIONS ======================================================

#' Load and prepare data from Excel file
load_regulation_data <- function(file_path) {
  message("Loading regulation data from Excel...")
  
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path)
  }
  
  # Read data
  data <- read_excel(file_path)
  
  # Check for required columns - mean log2 intensities for Control and Test
  required_cols <- c("ProteinID", "Description", "Localization", "mean_control", "mean_test")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  message(sprintf("  Loaded %d proteins", nrow(data)))
  
  # Standardize localization names (case-insensitive matching)
  data <- data %>%
    mutate(
      Localization = case_when(
        grepl("cytoplasmic membrane", Localization, ignore.case = TRUE) ~ "Cytoplasmic membrane",
        grepl("outer membrane", Localization, ignore.case = TRUE) ~ "Outer membrane",
        grepl("cytoplasm", Localization, ignore.case = TRUE) ~ "Cytoplasmic",
        grepl("periplasm", Localization, ignore.case = TRUE) ~ "Periplasmic",
        grepl("extracellular", Localization, ignore.case = TRUE) ~ "Extracellular",
        TRUE ~ Localization
      ),
      # Ensure localization is a factor with proper order
      Localization = factor(Localization, levels = LOCALIZATION_ORDER)
    ) %>%
    # Remove any proteins with unrecognized localization
    filter(!is.na(Localization))
  
  # Center each protein to its own mean for visualization
  # This shows relative differences between Control and Test for each protein
  data <- data %>%
    mutate(
      protein_mean = (mean_control + mean_test) / 2,
      Control_centered = mean_control - protein_mean,
      Test_centered = mean_test - protein_mean,
      logFC = mean_test - mean_control  # For reference
    )
  
  # Check localization distribution
  message("  Localization distribution:")
  loc_table <- table(data$Localization)
  for (loc in names(loc_table)) {
    message(sprintf("    %s: %d proteins", loc, loc_table[loc]))
  }
  
  # Check data range
  message("  Mean expression range (log2 scale):")
  message(sprintf("    Control: %.2f to %.2f", min(data$mean_control, na.rm = TRUE), max(data$mean_control, na.rm = TRUE)))
  message(sprintf("    Test: %.2f to %.2f", min(data$mean_test, na.rm = TRUE), max(data$mean_test, na.rm = TRUE)))
  message(sprintf("    logFC range: %.2f to %.2f", min(data$logFC, na.rm = TRUE), max(data$logFC, na.rm = TRUE)))
  
  # Count positive vs negative logFC
  n_control <- sum(data$logFC < 0, na.rm = TRUE)
  n_test <- sum(data$logFC > 0, na.rm = TRUE)
  message(sprintf("    %d proteins more abundant in Control (negative logFC)", n_control))
  message(sprintf("    %d proteins more abundant in Test (positive logFC)", n_test))
  
  return(data)
}

#' Create row labels with descriptions
create_row_labels <- function(protein_id, description, max_desc_length = 50) {
  # Find the longest protein ID to determine padding width
  max_id_length <- max(nchar(protein_id))
  
  # Pad protein IDs with spaces on the right for alignment
  padded_ids <- sprintf(paste0("%-", max_id_length, "s"), protein_id)
  
  # Truncate descriptions if needed
  desc_clean <- ifelse(
    nchar(description) > max_desc_length,
    paste0(substr(description, 1, max_desc_length - 3), "..."),
    description
  )
  
  # Combine with consistent dash position
  paste(padded_ids, "-", desc_clean)
}

#' Create localization-organized heatmap with two columns (centered)
create_localization_heatmap <- function(data, 
                                        output_file = "heatmap_localization.pdf") {
  
  message("Creating localization-organized heatmap...")
  
  # Sort data by localization only
  data_sorted <- data %>%
    arrange(Localization) %>%
    mutate(
      RowLabel = create_row_labels(ProteinID, Description)
    )
  
  # Create matrix using CENTERED values (relative to each protein's mean)
  expr_matrix <- as.matrix(data_sorted[, c("Control_centered", "Test_centered")])
  rownames(expr_matrix) <- data_sorted$RowLabel
  colnames(expr_matrix) <- c("Acetate", "Bicarbonate")
  
  message(sprintf("  Matrix dimensions: %d proteins x 2 columns", nrow(expr_matrix)))
  message(sprintf("  Centered data range: %.2f to %.2f", min(expr_matrix, na.rm = TRUE), max(expr_matrix, na.rm = TRUE)))
  
  # Create localization annotation for rows
  localization_vector <- setNames(
    as.character(data_sorted$Localization),
    data_sorted$RowLabel
  )
  
  # Color scheme: Blue-White-Red diverging, centered at 0
  # Negative (blue) = below protein's mean = relatively lower
  # Zero (white) = at protein's mean
  # Positive (red) = above protein's mean = relatively higher
  max_abs_val <- max(abs(range(expr_matrix, na.rm = TRUE)))
  
  col_fun <- colorRamp2(
    breaks = c(-max_abs_val, 0, max_abs_val),
    colors = c("#377eb8", "white", "#e41a1c")  # Blue - White - Red
  )
  
  # Calculate split points for localization groups
  loc_splits <- data_sorted$Localization
  
  # Create abbreviated row titles
  loc_abbrev_vector <- LOCALIZATION_ABBREV[as.character(loc_splits)]
  
  # Create heatmap
  ht <- Heatmap(
    expr_matrix,
    name = "Relative to\nprotein mean",
    col = col_fun,
    
    # Row settings
    cluster_rows = FALSE,  # Keep localization order
    show_row_dend = FALSE,
    row_split = loc_splits,
    row_title = loc_abbrev_vector[!duplicated(loc_splits)],  # Use abbreviations
    row_title_rot = 0,
    row_title_side = "left",  # Localization labels on the left
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_gap = unit(3, "mm"),  # Space between localization groups
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "right",  # Protein labels on the right
    
    # Column settings
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_side = "top",  # Move labels to top
    column_names_rot = 0,
    column_names_centered = TRUE,
    show_column_names = TRUE,
    column_gap = unit(5, "mm"),  # Add space between columns
    
    # Overall appearance
    border = FALSE,
    width = unit(4, "cm"),  # Increased from 3cm to give more space for column labels
    
    # Legend
    heatmap_legend_param = list(
      title = "log2 deviation\nfrom protein mean",
      title_position = "topcenter",
      direction = "horizontal",
      legend_width = unit(5, "cm"),
      labels_gp = gpar(fontsize = 9),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    )
  )
  
  # Calculate appropriate height based on number of proteins
  plot_height <- min(30, max(8, 0.15 * nrow(expr_matrix)))
  
  # Plot width for two columns
  plot_width <- 7.5
  
  # Save PDF
  pdf(file.path(OUTPUT_DIR, output_file), 
      width = plot_width, 
      height = plot_height)
  
  # Draw heatmap
  draw_result <- ComplexHeatmap::draw(
    ht, 
    heatmap_legend_side = "bottom",
    padding = unit(c(2, 2, 2, 2), "mm")
  )
  
  dev.off()
  message("  Saved PDF to: ", file.path(OUTPUT_DIR, output_file))
  
  # Save PNG
  png_file <- sub("\\.pdf$", ".png", output_file)
  png(file.path(OUTPUT_DIR, png_file), 
      width = plot_width, 
      height = plot_height,
      units = "in",
      res = 300)
  
  draw_result <- ComplexHeatmap::draw(
    ht, 
    heatmap_legend_side = "bottom",
    padding = unit(c(2, 2, 2, 2), "mm")
  )
  
  dev.off()
  message("  Saved PNG to: ", file.path(OUTPUT_DIR, png_file))
  
  return(ht)
}

#' Create summary statistics
create_summary_stats <- function(data) {
  message("Creating summary statistics...")
  
  # Overall summary
  overall <- data %>%
    summarise(
      N_proteins = n(),
      Control_mean = mean(mean_control, na.rm = TRUE),
      Control_sd = sd(mean_control, na.rm = TRUE),
      Test_mean = mean(mean_test, na.rm = TRUE),
      Test_sd = sd(mean_test, na.rm = TRUE),
      Mean_logFC = mean(logFC, na.rm = TRUE),
      SD_logFC = sd(logFC, na.rm = TRUE),
      N_higher_in_control = sum(logFC < 0, na.rm = TRUE),
      N_higher_in_test = sum(logFC > 0, na.rm = TRUE)
    )
  
  # By localization
  by_localization <- data %>%
    group_by(Localization) %>%
    summarise(
      N = n(),
      Control_mean = mean(mean_control, na.rm = TRUE),
      Test_mean = mean(mean_test, na.rm = TRUE),
      Mean_logFC = mean(logFC, na.rm = TRUE),
      N_higher_in_control = sum(logFC < 0, na.rm = TRUE),
      N_higher_in_test = sum(logFC > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(Localization, LOCALIZATION_ORDER))
  
  message("Overall statistics:")
  print(overall)
  
  message("\nBy localization:")
  print(by_localization)
  
  return(list(
    overall = overall,
    by_localization = by_localization
  ))
}

# ===== MAIN EXECUTION ========================================================

message("\n=== Starting Numerical Heatmap Generation ===\n")

# Load data
protein_data <- load_regulation_data(INPUT_FILE)

# Create summary statistics
summary_stats <- create_summary_stats(protein_data)

# Generate timestamp for output files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

# Create heatmap
heatmap_result <- create_localization_heatmap(
  protein_data,
  output_file = paste0("heatmap_localization_twocol_", timestamp, ".pdf")
)

# Save summary statistics
summary_output <- file.path(OUTPUT_DIR, paste0("summary_stats_", timestamp, ".txt"))
sink(summary_output)
cat("=== Protein Regulation Summary ===\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Total proteins:", nrow(protein_data), "\n\n")
cat("Overall Statistics:\n")
print(summary_stats$overall)
cat("\n\nBy Localization:\n")
print(summary_stats$by_localization)
sink()
message("Summary statistics saved to: ", summary_output)

message("\n=== Heatmap Generation Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("Total proteins plotted: ", nrow(protein_data))