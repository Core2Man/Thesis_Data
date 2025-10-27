#!/usr/bin/env Rscript

# ============================================================================
# Standalone Heatmap Generator for Bicarbonate Proteomics Data
# Uses cleaned MSnSet objects from main bicarbonate analysis
# ============================================================================

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,      # Data manipulation
  readxl,         # Read Excel files
  ComplexHeatmap, # Primary heatmap package
  pheatmap,       # Alternative heatmap package
  circlize,       # Color functions
  RColorBrewer,   # Color palettes
  writexl,        # Save results
  MSnbase         # For MSnSet objects
)

# ===== CONFIGURATION =========================================================

ANNOTATION_FILE <- '/Users/sysoevm/Library/CloudStorage/Dropbox/Pascal Lab/Desulfuromonas/Bicarbonate Chapter/Stats DIA-NN/FinalTP/tables/ALL_RESULTS_COMBINED.xlsx'
OUTPUT_DIR <- "/Users/sysoevm/Library/CloudStorage/Dropbox/Pascal Lab/Desulfuromonas/Bicarbonate Chapter/Stats DIA-NN/heatmaps"
PROTEIN_ID_COLUMN <- "uniprot_id"      # Primary identifier
DISPLAY_ID_COLUMN <- "refseq_id"       # For row labels
DESCRIPTION_COLUMN <- "uniprot_description"  # For protein descriptions

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
SIGNIFICANCE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.1  # Matching your main script threshold
TOP_N_PROTEINS <- 50
MIN_PROTEINS_FOR_HEATMAP <- 5

# Color scheme for groups - Bicarbonate experiment (Control vs Test)
GROUP_COLORS <- c(
  Control = "#377eb8",  # Blue - with acetate+bicarbonate
  Test = "#e41a1c"      # Red - without acetate
)

# ===== HELPER FUNCTIONS ======================================================

#' Load expression data from MSnSet object
load_expression_from_msnset <- function(msn_object = NULL, msnset_file = NULL) {
  message("Loading expression data from MSnSet...")
  
  # Load from file if provided
  if (!is.null(msnset_file)) {
    if (!file.exists(msnset_file)) {
      stop("MSnSet file not found: ", msnset_file)
    }
    msn <- readRDS(msnset_file)
    message(sprintf("  Loaded from file: %s", basename(msnset_file)))
  } else if (!is.null(msn_object)) {
    msn <- msn_object
    message("  Using MSnSet from environment")
  } else {
    stop("Must provide either msn_object or msnset_file")
  }
  
  # Extract expression matrix
  expr_matrix <- exprs(msn)
  
  message(sprintf("  Matrix: %d proteins x %d samples", 
                  nrow(expr_matrix), ncol(expr_matrix)))
  
  # Final validation - should not be needed if data is truly cleaned!
  invalid_rows <- apply(expr_matrix, 1, function(x) any(!is.finite(x)))
  if (any(invalid_rows)) {
    warning(sprintf("Found %d proteins with invalid values - removing", sum(invalid_rows)))
    expr_matrix <- expr_matrix[!invalid_rows, , drop = FALSE]
  }
  
  return(expr_matrix)
}

#' Load annotations with full support for uniprot_id, refseq_id, and descriptions
load_annotations <- function(file_path, protein_ids) {
  message("Loading protein annotations...")
  
  if (!file.exists(file_path)) {
    warning("Annotation file not found - will use protein IDs only")
    return(list(
      annotations = data.frame(uniprot_id = protein_ids),
      descriptions = setNames(rep("", length(protein_ids)), protein_ids),
      id_map = setNames(protein_ids, protein_ids)
    ))
  }
  
  annot_data <- read_excel(file_path)
  
  # Check for required columns
  required_cols <- c(PROTEIN_ID_COLUMN)
  missing_cols <- setdiff(required_cols, names(annot_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for optional columns
  has_display_id <- DISPLAY_ID_COLUMN %in% names(annot_data)
  has_description <- DESCRIPTION_COLUMN %in% names(annot_data)
  
  if (!has_display_id) {
    warning(sprintf("%s column not found - will use %s", 
                    DISPLAY_ID_COLUMN, PROTEIN_ID_COLUMN))
  }
  if (!has_description) {
    warning(sprintf("%s column not found - no descriptions will be shown", 
                    DESCRIPTION_COLUMN))
  }
  
  # Extract annotations for the proteins we're plotting
  select_cols <- c(PROTEIN_ID_COLUMN)
  if (has_display_id) select_cols <- c(select_cols, DISPLAY_ID_COLUMN)
  if (has_description) select_cols <- c(select_cols, DESCRIPTION_COLUMN)
  
  annotations <- annot_data %>%
    dplyr::select(all_of(select_cols)) %>%
    distinct(.data[[PROTEIN_ID_COLUMN]], .keep_all = TRUE) %>%
    filter(.data[[PROTEIN_ID_COLUMN]] %in% protein_ids)
  
  # Create ID mapping for row labels
  if (has_display_id) {
    id_map <- setNames(
      annotations[[DISPLAY_ID_COLUMN]],
      annotations[[PROTEIN_ID_COLUMN]]
    )
    message(sprintf("  Using %s for row labels", DISPLAY_ID_COLUMN))
  } else {
    id_map <- setNames(
      annotations[[PROTEIN_ID_COLUMN]],
      annotations[[PROTEIN_ID_COLUMN]]
    )
    message(sprintf("  Using %s for row labels", PROTEIN_ID_COLUMN))
  }
  
  # Create descriptions mapping
  if (has_description) {
    descriptions <- setNames(
      annotations[[DESCRIPTION_COLUMN]],
      annotations[[PROTEIN_ID_COLUMN]]
    )
  } else {
    descriptions <- setNames(
      rep("", nrow(annotations)),
      annotations[[PROTEIN_ID_COLUMN]]
    )
  }
  
  message(sprintf("  Loaded annotations for %d proteins", nrow(annotations)))
  
  return(list(
    annotations = annotations,
    descriptions = descriptions,
    id_map = id_map
  ))
}

#' Get significant proteins from results
get_significant_proteins <- function(results_file, comparison = NULL,
                                     p_cutoff = 0.05, fc_cutoff = 1) {
  message("Identifying significant proteins...")
  
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }
  
  results <- read_excel(results_file)
  
  # Debug: show column names
  message("  Available columns: ", paste(names(results), collapse = ", "))
  
  # Filter by comparison if specified
  if (!is.null(comparison) && "Comparison" %in% names(results)) {
    results <- filter(results, Comparison == comparison)
    message(sprintf("  Filtered to comparison: %s", comparison))
  }
  
  # Get significant proteins (matching your main script logic)
  sig_proteins <- results %>%
    filter(adj.P.Val <= p_cutoff, abs(logFC) >= fc_cutoff) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    pull(!!sym(PROTEIN_ID_COLUMN)) %>%
    unique()
  
  message(sprintf("  Found %d significant proteins", length(sig_proteins)))
  
  return(sig_proteins)
}

#' Create heatmap using ComplexHeatmap
create_complex_heatmap <- function(expr_matrix, annot_info, 
                                   output_file = "heatmap_complex.pdf",
                                   title = "Protein Expression Heatmap",
                                   cluster_rows = TRUE,
                                   cluster_cols = FALSE,
                                   use_zscore = FALSE) {
  
  message("Creating ComplexHeatmap...")
  
  # Identify Control samples (R1F, R2F, R3_2F)
  control_samples <- c("R1F", "R2F", "R3_2F")
  control_cols <- which(colnames(expr_matrix) %in% control_samples)
  
  if (length(control_cols) == 0) {
    stop("No control samples found (R1F, R2F, R3_2F)")
  }
  
  message(sprintf("  Using %d control samples: %s", 
                  length(control_cols), 
                  paste(colnames(expr_matrix)[control_cols], collapse = ", ")))
  
  # Calculate baseline (mean of control samples)
  baseline <- rowMeans(expr_matrix[, control_cols, drop = FALSE])
  baseline[baseline == 0] <- min(baseline[baseline > 0], na.rm = TRUE) * 0.01
  
  # Calculate log2FC relative to control
  expr_scaled <- expr_matrix - baseline
  
  if (use_zscore) {
    message("Applying z-score transformation (row-wise)...")
    expr_scaled <- t(scale(t(expr_scaled)))  # z-score across samples for each protein
    scale_label <- "Z-score"
    color_breaks <- seq(-2, 2, length.out = 9)
  } else {
    scale_label <- "log2FC vs Control"
    color_breaks <- seq(-3, 3, length.out = 9)
  }
  
  message(sprintf("Data range: %.3f to %.3f", 
                  min(expr_scaled, na.rm=TRUE), 
                  max(expr_scaled, na.rm=TRUE)))
  
  # Remove any invalid rows
  invalid_rows <- apply(expr_scaled, 1, function(x) any(!is.finite(x)))
  if (any(invalid_rows)) {
    message(sprintf("  Removing %d proteins with invalid values", sum(invalid_rows)))
    expr_scaled <- expr_scaled[!invalid_rows, , drop = FALSE]
  }
  
  # Determine groups for each sample
  groups <- character(ncol(expr_matrix))
  for (i in 1:ncol(expr_matrix)) {
    sample_name <- colnames(expr_matrix)[i]
    if (sample_name %in% control_samples) {
      groups[i] <- "Control"
    } else {
      groups[i] <- "Test"
    }
  }
  
  message(sprintf("  Sample grouping: %s", 
                  paste(unique(groups), table(groups)[unique(groups)], 
                        sep = "=", collapse = ", ")))
  
  # Column annotation
  #col_ha <- HeatmapAnnotation(
  #  Group = factor(groups, levels = c("Control", "Test")),
  #  col = list(Group = GROUP_COLORS),
  #  annotation_name_side = "left",
  #  annotation_legend_param = list(
  #    Group = list(
  #      title = "Condition",
  #      labels = c("Control\n(acetate+bicarbonate)", "Test\n(no acetate)")
  #    )
  #  )
  #)
  
  # Row labels with refseq_id and descriptions (matching Mg script logic)
  if (!is.null(annot_info)) {
    row_labels <- sapply(rownames(expr_scaled), function(uniprot_id) {
      # Get refseq ID (or display ID)
      refseq <- annot_info$id_map[uniprot_id]
      if (is.na(refseq)) refseq <- uniprot_id
      
      # Get description
      desc <- annot_info$descriptions[uniprot_id]
      if (!is.na(desc) && nchar(desc) > 0) {
        desc_short <- substr(desc, 1, 40)
        if (nchar(desc) > 40) desc_short <- paste0(desc_short, "...")
        return(paste(refseq, "-", desc_short))
      } else {
        return(refseq)
      }
    })
    rownames(expr_scaled) <- row_labels
  }
  
  # Color scheme - Blue-White-Red for log2FC
  col_fun <- colorRamp2(
    breaks = seq(-3, 3, length.out = 9),
    colors = rev(RColorBrewer::brewer.pal(9, "RdBu"))
  )
  
  # Create heatmap
  ht <- Heatmap(
    expr_scaled,
    name = scale_label,
    col = col_fun,
    #top_annotation = col_ha,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_dend = cluster_rows,
    show_column_dend = cluster_cols,
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 11),
    column_names_rot = 45,
    column_split = factor(groups, levels = c("Control", "Test")),
    column_title = c("Acetate+\nbicarbonate", 
                     "Only bicarbonate \n(without acetate)"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_gap = unit(2, "mm"),
    border = TRUE,
    width = unit(10, "cm"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter",
      legend_width = unit(6, "cm")
    )
  )
  
  # Save
  pdf(file.path(OUTPUT_DIR, output_file), 
      width = 10, 
      height = min(20, 5 + nrow(expr_scaled) * 0.15))
  ComplexHeatmap::draw(ht, 
                       heatmap_legend_side = "bottom",
                       annotation_legend_side = "bottom",
                       merge_legend = TRUE) 
  dev.off()
  
  message("  Saved to: ", file.path(OUTPUT_DIR, output_file))
  png_file <- sub("\\.pdf$", ".png", output_file)
  png(file.path(OUTPUT_DIR, png_file), 
      width = 10, 
      height = min(20, 5 + nrow(expr_scaled) * 0.15),
      units = "in",
      res = 300)
  ComplexHeatmap::draw(ht, 
                       heatmap_legend_side = "bottom",
                       annotation_legend_side = "bottom",
                       merge_legend = TRUE) 
  dev.off()
  
  message("  Saved PNG to: ", file.path(OUTPUT_DIR, png_file))
  return(ht)
}

#' Main function to generate heatmaps
generate_heatmaps <- function(msn_object = NULL,
                              msnset_file = NULL,
                              specific_proteins = NULL, 
                              comparison = NULL,
                              top_n = TOP_N_PROTEINS) {
  
  # Load expression data from MSnSet
  expr_matrix <- load_expression_from_msnset(
    msn_object = msn_object,
    msnset_file = msnset_file
  )
  
  # Get significant proteins if not specified
  if (is.null(specific_proteins)) {
    sig_proteins <- get_significant_proteins(
      ANNOTATION_FILE, 
      comparison = comparison,
      p_cutoff = SIGNIFICANCE_CUTOFF,
      fc_cutoff = LOGFC_CUTOFF
    )
    proteins_to_plot <- head(sig_proteins, top_n)
  } else {
    proteins_to_plot <- specific_proteins
  }
  
  # Check minimum proteins
  if (length(proteins_to_plot) < MIN_PROTEINS_FOR_HEATMAP) {
    stop(sprintf("Not enough proteins for heatmap (%d < %d)", 
                 length(proteins_to_plot), MIN_PROTEINS_FOR_HEATMAP))
  }
  
  # Filter to proteins present in expression data
  proteins_to_plot <- proteins_to_plot[proteins_to_plot %in% rownames(expr_matrix)]
  expr_subset <- expr_matrix[proteins_to_plot, , drop = FALSE]
  
  message(sprintf("Creating heatmap for %d proteins", nrow(expr_subset)))
  
  # Load annotations
  annot_info <- load_annotations(ANNOTATION_FILE, proteins_to_plot)
  
  # Generate output filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  if (!is.null(comparison)) {
    comp_clean <- gsub("[^A-Za-z0-9]", "_", comparison)
    output_base <- paste0("heatmap_bicarbonate_", comp_clean, "_", timestamp)
  } else {
    output_base <- paste0("heatmap_bicarbonate_all_sig_", timestamp)
  }
  
  # Create heatmap (log2FC version)
  create_complex_heatmap(
    expr_subset, annot_info,
    output_file = paste0(output_base, "_log2FC.pdf"),
    title = ifelse(!is.null(comparison), comparison, "Significant Proteins"),
    use_zscore = FALSE
  )
  
  # Save protein list with full annotations
  protein_list <- annot_info$annotations %>%
    filter(.data[[PROTEIN_ID_COLUMN]] %in% proteins_to_plot)
  
  write_xlsx(protein_list, 
             file.path(OUTPUT_DIR, paste0(output_base, "_protein_list.xlsx")))
  
  message("Heatmap generation complete!")
  
  return(list(
    proteins = proteins_to_plot,
    annotations = annot_info$annotations,
    expression = expr_subset
  ))
}

# ===== RUN HEATMAP GENERATION ================================================

message("\n=== Starting Bicarbonate Heatmap Generation ===\n")

# Check if msn_set exists in environment
if (!exists("msn_set")) {
  stop("Error: msn_set not found in environment. Please run the main bicarbonate analysis script first!")
}

# Option 1: Use MSnSet from environment (after running main analysis)
results <- generate_heatmaps(
  msn_object = msn_set,  # Use the cleaned MSnSet from your bicarbonate analysis
  top_n = 50
)

message("\n=== Heatmap Generation Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("Number of proteins plotted: ", length(results$proteins))