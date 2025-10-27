#!/usr/bin/env Rscript

## --------------------------------------------------------------
## Bicarbonate Proteomics Analysis - Final Clean Version
## Based on Mg analysis best practices
## --------------------------------------------------------------

## Packages ------------------------------------------------------
libs <- c("MSnbase", "impute", "limma", "readr", "dplyr", "tibble",
          "ComplexHeatmap", "ggplot2", "clusterProfiler", "tidyr",
          "GO.db", "forcats", "yaml", "ggrepel", "showtext", 
          "AnnotationDbi", "writexl", "uwot", "patchwork", "readxl",
          "circlize", "enrichplot", "ggfortify")

# Enhanced package installer for CRAN + Bioconductor
load_or_install_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "... Installing", pkg, "\n")
      if (pkg %in% rownames(available.packages())) {
        install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
      } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cran.rstudio.com/")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}
load_or_install_packages(libs)

font_add_google("Lato", "lato")
showtext_auto()

## Utility Functions ---------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x  # Coalesce operator

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "...", ..., "\n")
warn <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "!!! WARNING:", ..., "\n")

monitor_memory <- function(label = "") {
  gc_result <- gc(reset = TRUE)
  mem_used <- sum(gc_result[, 2])  # Memory in Mb
  msg(sprintf("Memory check %s: %.1f Mb used", label, mem_used))
  return(mem_used)
}

## Configuration -------------------------------------------------
cfg_file <- "/Users/sysoevm/Library/CloudStorage/Dropbox/Pascal Lab/Desulfuromonas/Bicarbonate Chapter/Stats DIA-NN/bicarb_config.yml"
if (!file.exists(cfg_file)) stop("Config file not found: ", cfg_file)

cfg <- yaml::read_yaml(cfg_file)
out_dir      <- cfg$out_dir
raw_file     <- cfg$raw_file
go_file      <- cfg$go_file
sig_adjP     <- cfg$sig_adjP %||% 0.05
sig_logFC    <- cfg$sig_logFC %||% 1
go_min_score <- cfg$go_min_score %||% 0.5

# Add timestamp for unique file naming
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

# Create output directories
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(out_dir, "plots")
tables_dir <- file.path(out_dir, "tables")
dir.create(plots_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)

# Validate inputs
for (f in c(raw_file, go_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

# Test Excel file readability
tryCatch({
  test <- readxl::read_excel(raw_file, n_max = 2)
}, error = function(e) stop("Failed to read Excel file: ", e$message))

msg("Configuration loaded and validated")
msg("Configuration Summary:")
msg(sprintf("  Output directory: %s", out_dir))
msg(sprintf("  Significance thresholds: adj.P <= %.3f, |logFC| >= %.2f", sig_adjP, sig_logFC))
msg(sprintf("  GO score threshold: >= %.2f", go_min_score))
msg(sprintf("  Timestamp: %s", timestamp))

# Function to extract protein accessions robustly
id_regex <- "^(?:sp|tr)\\|([^|]+)\\|.*$"
acc_from <- function(x) {
  out <- sub(id_regex, "\\1", x)
  need_fallback <- out == x
  out[need_fallback] <- sub("^([^|]+).*$", "\\1", x[need_fallback])
  out
}

## Data Loading --------------------------------------------------
msg("Reading normalized data from:", raw_file)
normalized_data_raw <- readxl::read_excel(raw_file) %>%
  dplyr::select(ProteinID, everything())

# Convert to matrix format for MSnSet
tab <- normalized_data_raw %>%
  column_to_rownames(var = "ProteinID") %>%
  as.matrix()

# Apply robust ID extraction
rownames(tab) <- acc_from(rownames(tab))
monitor_memory("after data load")

msg(sprintf("Loaded data: %d proteins x %d samples", nrow(tab), ncol(tab)))

# Check data scale
data_range <- range(tab, na.rm = TRUE)
msg(sprintf("Data range: [%.2f, %.2f]", data_range[1], data_range[2]))
if (data_range[1] < 0 || data_range[2] <= 50) {
  msg("Data likely log-scale (negatives or small dynamic range)")
} else if (data_range[2] >= 1e4) {
  warn("Data likely raw (very large values) - log2 transform recommended before limma")
}

# Missing data summary
msg("Missing data summary:")
msg(sprintf("  Total NA values: %d (%.2f%%)", 
            sum(is.na(tab)), 
            100 * sum(is.na(tab)) / length(tab)))
msg(sprintf("  Proteins with any NA: %d (%.2f%%)", 
            sum(rowSums(is.na(tab)) > 0),
            100 * sum(rowSums(is.na(tab)) > 0) / nrow(tab)))

## Build MSnSet Object -------------------------------------------
msg("Building MSnSet object")
msn_set <- MSnSet(exprs = tab)

# Define sample groups
control_samples <- c("R1F", "R2F", "R3_2F")
test_samples <- c("R4F", "R5_2F", "R6_2F", "R7F")

# Annotate groups
annotate_groups <- function(msn) {
  rn <- sampleNames(msn)
  
  # Check if all samples are accounted for
  all_samples <- c(control_samples, test_samples)
  if (!all(rn %in% all_samples)) {
    warn("Some samples not assigned to groups:", 
         paste(rn[!rn %in% all_samples], collapse = ", "))
  }
  
  Group <- factor(
    ifelse(rn %in% control_samples, "Control", 
           ifelse(rn %in% test_samples, "Test", "Unknown")),
    levels = c("Control", "Test")
  )
  
  pData(msn)$Group <- Group
  pData(msn)$SampleName <- rn
  
  return(msn)
}

msn_set <- annotate_groups(msn_set)

msg("Sample distribution:")
print(table(pData(msn_set)$Group))

# Check for minimum replicates
if (any(table(pData(msn_set)$Group) < 3)) {
  warn("Some groups have <3 replicates - statistical power may be limited")
}

## Advanced Imputation Function ----------------------------------
impute_msn_data <- function(msn, tag) {
  msg(sprintf("Starting mixed imputation for %s", tag))
  
  expr_data <- exprs(msn)
  expr_data_original <- expr_data  # Keep original for validation
  groups <- pData(msn)$Group
  present_levels <- levels(groups)[table(groups) > 0]
  
  # Calculate missingness per protein per group
  missing_pattern <- matrix(0, nrow = nrow(expr_data), 
                            ncol = length(present_levels))
  colnames(missing_pattern) <- present_levels
  
  for (i in seq_along(present_levels)) {
    g <- present_levels[i]
    idx <- which(groups == g)
    missing_pattern[, i] <- rowSums(is.na(expr_data[, idx, drop = FALSE]))
  }
  
  # Get actual number of samples per group
  n_per_group <- table(groups)
  
  # Classify proteins based on missingness pattern
  completely_missing <- apply(missing_pattern, 1, function(x) {
    any(x == n_per_group[colnames(missing_pattern)])
  })
  
  partially_missing <- rowSums(missing_pattern > 0) > 0 & !completely_missing
  
  msg(sprintf("Proteins with complete group absence (MNAR): %d", sum(completely_missing)))
  msg(sprintf("Proteins with partial missingness (MAR): %d", sum(partially_missing)))
  msg(sprintf("Complete proteins: %d", sum(!completely_missing & !partially_missing)))
  
  # MAR imputation using kNN
  if (sum(partially_missing) > 0) {
    partial_rows <- which(partially_missing)
    expr_partial <- expr_data[partial_rows, , drop = FALSE]
    
    # Use impute.knn with appropriate k
    k_value <- min(5, floor(ncol(expr_partial) / 2))
    if (k_value > 0) {
      expr_partial_imputed <- impute.knn(expr_partial, k = k_value)$data
      expr_data[partial_rows, ] <- expr_partial_imputed
      msg(sprintf("kNN imputation completed for %d MAR proteins (k=%d)", 
                  sum(partially_missing), k_value))
    }
  }
  
  # MNAR imputation - Perseus-like method
  if (sum(completely_missing) > 0) {
    msg("Performing MNAR imputation using Perseus-like method")
    
    # Calculate parameters from the entire dataset
    all_values <- as.vector(expr_data)
    global_mean <- mean(all_values, na.rm = TRUE)
    global_sd <- sd(all_values, na.rm = TRUE)
    
    # Down-shifted distribution for imputation
    impute_mean <- global_mean - 1.8 * global_sd
    impute_sd <- global_sd * 0.3
    
    # Impute MNAR values
    complete_idx <- which(completely_missing)
    for (p in complete_idx) {
      for (g_idx in split(seq_len(ncol(expr_data)), groups)) {
        vals <- expr_data[p, g_idx]
        if (all(is.na(vals))) {
          set.seed(42 + p)  # For reproducibility
          expr_data[p, g_idx] <- rnorm(length(vals), 
                                       mean = impute_mean, 
                                       sd = impute_sd)
        }
      }
    }
    msg(sprintf("MNAR imputation completed for %d proteins", sum(completely_missing)))
  }
  
  # Validation checks
  post_var <- apply(expr_data, 1, var, na.rm = TRUE)
  pre_var <- apply(expr_data_original, 1, var, na.rm = TRUE)
  
  had_missing <- rowSums(is.na(expr_data_original)) > 0
  if (sum(had_missing) > 0) {
    var_ratio <- post_var[had_missing] / pre_var[had_missing]
    var_preserved <- median(var_ratio, na.rm = TRUE)
    msg(sprintf("Median variance preserved after imputation: %.1f%%", 
                var_preserved * 100))
    
    if (var_preserved < 0.8) {
      warn("Imputation reduced variance by >20% - consider alternative methods")
    }
    
    # Correlation check
    mask_not_na <- !is.na(expr_data_original)
    if (sum(mask_not_na) > 0) {
      cor_check <- cor(as.vector(expr_data_original[mask_not_na]), 
                       as.vector(expr_data[mask_not_na]))
      msg(sprintf("Correlation between original and post-imputation: %.3f", 
                  cor_check))
    }
  }
  
  # Update expression data
  exprs(msn) <- expr_data
  
  # Remove proteins with remaining NAs
  remaining_na <- rowSums(is.na(expr_data)) > 0
  if (any(remaining_na)) {
    msg(sprintf("Removing %d proteins with remaining NAs", sum(remaining_na)))
    msn <- msn[!remaining_na, ]
  }
  
  msg(sprintf("Imputation complete for %s: %d proteins retained", 
              tag, nrow(msn)))
  return(msn)
}

# Apply imputation
set.seed(42)
msn_set <- impute_msn_data(msn_set, "Bicarbonate_Analysis")
gc()  # Force garbage collection
monitor_memory("after imputation")

# Save imputed data
imputed_data_export <- as.data.frame(exprs(msn_set)) %>%
  rownames_to_column("ProteinID")
write_xlsx(imputed_data_export, 
           path = file.path(tables_dir, 
                            paste0("Final_Imputed_Data_", timestamp, ".xlsx")))

## Enhanced Limma Analysis Function ------------------------------
run_limma_analysis <- function(msn, tag, timestamp = format(Sys.time(), "%Y%m%d_%H%M")) {
  msg("Running enhanced limma analysis on", tag, "data")
  
  # Get sample lists for each group
  sample_lists <- lapply(levels(pData(msn)$Group), function(g) {
    sampleNames(msn)[pData(msn)$Group == g]
  })
  names(sample_lists) <- levels(pData(msn)$Group)
  
  # Set up design matrix
  pData(msn)$Group <- relevel(pData(msn)$Group, ref = "Control")
  design <- model.matrix(~ Group, data = pData(msn))
  
  # Fit linear model
  fit <- lmFit(exprs(msn), design)
  fit2 <- eBayes(fit)
  
  all_results <- list()
  
  # Process each comparison (for now just Test vs Control)
  for (coef_name in colnames(design)[-1]) {
    comp_label <- gsub("^Group", "", coef_name)
    comp_clean <- paste(comp_label, "vs Control")
    
    msg(sprintf("Processing comparison: %s", comp_clean))
    
    # Extract results
    tt <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
    
    # Ensure protein IDs match
    common_ids <- intersect(rownames(tt), rownames(exprs(msn)))
    if (length(common_ids) == 0) {
      warn("No overlapping protein IDs for", comp_label)
      next
    }
    
    # Calculate mean intensities
    mean_group <- if (comp_label %in% names(sample_lists) && 
                      length(sample_lists[[comp_label]]) > 0) {
      rowMeans(exprs(msn)[common_ids, sample_lists[[comp_label]], drop = FALSE], 
               na.rm = TRUE)
    } else {
      rep(NA, length(common_ids))
    }
    mean_control <- rowMeans(exprs(msn)[common_ids, sample_lists[["Control"]], drop = FALSE], 
                             na.rm = TRUE)
    
    # Enhance results with metadata
    tt_enhanced <- tt[common_ids, , drop = FALSE] %>%
      rownames_to_column("ProteinID") %>%
      mutate(
        Direction = case_when(
          logFC > 0 ~ paste("Up in", comp_label),
          logFC < 0 ~ "Up in Control",
          TRUE ~ "No change"
        ),
        Regulation = case_when(
          adj.P.Val <= sig_adjP & logFC >= sig_logFC ~ paste("Significantly up in", comp_label),
          adj.P.Val <= sig_adjP & logFC <= -sig_logFC ~ "Significantly up in Control",
          adj.P.Val <= sig_adjP & abs(logFC) < sig_logFC ~ "Significant but low fold change",
          TRUE ~ "Not significant"
        ),
        FoldChange_linear = 2^logFC,
        abs_logFC = abs(logFC),
        !!paste0("Mean_", comp_label) := mean_group,
        Mean_Control = mean_control,
        Comparison = comp_clean,
        Analysis = tag
      ) %>%
      dplyr::select(
        ProteinID, logFC, FoldChange_linear, Direction, AveExpr,
        starts_with("Mean_"), t, P.Value, adj.P.Val, B,
        Regulation, abs_logFC, Comparison, Analysis
      ) %>%
      arrange(adj.P.Val, desc(abs_logFC))
    
    # Save detailed results
    output_file <- file.path(tables_dir, 
                             paste0("Results_", comp_label, "_vs_Control_", 
                                    tag, "_", timestamp, ".xlsx"))
    write_xlsx(tt_enhanced, output_file)
    msg(sprintf("Saved results to: %s", basename(output_file)))
    
    # Create summary statistics
    summary_stats <- tt_enhanced %>%
      group_by(Regulation) %>%
      summarise(
        Count = n(),
        Mean_logFC = mean(logFC),
        Median_logFC = median(logFC),
        Max_abs_logFC = max(abs(logFC)),
        .groups = 'drop'
      )
    
    write_xlsx(summary_stats, 
               file.path(tables_dir, 
                         paste0("Summary_", comp_label, "_vs_Control_", 
                                tag, "_", timestamp, ".xlsx")))
    
    # Print summary
    msg(strrep("=", 50))
    msg(sprintf("SUMMARY: %s (%s)", comp_clean, tag))
    msg(sprintf("Total proteins analyzed: %d", nrow(tt_enhanced)))
    msg(sprintf("Significantly up in %s: %d", comp_label, 
                sum(tt_enhanced$Regulation == paste("Significantly up in", comp_label))))
    msg(sprintf("Significantly up in Control: %d", 
                sum(tt_enhanced$Regulation == "Significantly up in Control")))
    msg(sprintf("Significance thresholds: adj.P.Val <= %.3f, |logFC| >= %.2f", 
                sig_adjP, sig_logFC))
    msg(strrep("=", 50))
    
    all_results[[comp_label]] <- tt_enhanced
    
    # QC Plots
    pdf(file.path(plots_dir, 
                  paste0("QC_", comp_label, "_", tag, "_", timestamp, ".pdf")), 
        width = 10, height = 5)
    par(mfrow = c(1,2))
    plotSA(fit2, main = paste("Mean-variance trend -", tag))
    hist(fit2$p.value[, coef_name], 50, xlim = c(0,1), 
         main = paste("p-value distribution (", comp_label, " vs Control)"), 
         xlab = "p-value")
    dev.off()
    
    # Enhanced Volcano Plot
    volc_colors <- c(
      "Not significant" = "grey70",
      "Significant but low fold change" = "#feb24c",
      "Significantly up in Control" = "#377eb8"
    )
    volc_colors[paste("Significantly up in", comp_label)] <- "#e41a1c"
    
    # Calculate label positions
    x_right <- quantile(tt_enhanced$logFC, 0.98, na.rm = TRUE)
    x_left  <- quantile(tt_enhanced$logFC, 0.02, na.rm = TRUE)
    y_top   <- quantile(-log10(tt_enhanced$adj.P.Val), 0.98, na.rm = TRUE)
    
    # Identify top proteins for labeling (optional)
    top_proteins <- tt_enhanced %>%
      filter(adj.P.Val < 0.01) %>%
      arrange(desc(abs_logFC)) %>%
      head(10)
    
    p_volcano <- ggplot(tt_enhanced, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(colour = Regulation), size = 2, alpha = 0.7) +
      scale_colour_manual(values = volc_colors) +
      geom_vline(xintercept = c(-sig_logFC, sig_logFC), 
                 linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(sig_adjP), 
                 linetype = "dashed", alpha = 0.5) +
      annotate("text", x = x_right, y = y_top, 
               label = paste("More abundant\nin", comp_label), 
               hjust = 1, size = 3.5, color = "#e41a1c") +
      annotate("text", x = x_left, y = y_top, 
               label = "More abundant\nin Control", 
               hjust = 0, size = 3.5, color = "#377eb8") +
      labs(
        title = paste("Volcano Plot:", comp_clean),
        subtitle = paste("Analysis:", tag, "| Positive logFC = Higher in", comp_label),
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"Adjusted p-value")
      ) +
      theme_classic(base_size = 12) + 
      theme(legend.position = "bottom")
    
    # Add labels for top proteins if needed
    if (nrow(top_proteins) > 0) {
      p_volcano <- p_volcano +
        geom_text_repel(
          data = top_proteins,
          aes(label = ProteinID),
          size = 3, 
          box.padding = 0.5,
          max.overlaps = 20
        )
    }
    
    ggsave(file.path(plots_dir, 
                     paste0("Volcano_", comp_label, "_vs_Control_", 
                            tag, "_", timestamp, ".pdf")),
           plot = p_volcano, width = 9, height = 7)
    
    # Distribution bar plot
    regulation_summary <- tt_enhanced %>%
      filter(Regulation %in% c(paste("Significantly up in", comp_label), 
                               "Significantly up in Control")) %>%
      group_by(Direction) %>%
      summarise(Count = n(), .groups = 'drop')
    
    if (nrow(regulation_summary) > 0) {
      p_bar <- ggplot(regulation_summary, aes(x = Direction, y = Count, fill = Direction)) +
        geom_bar(stat = "identity", width = 0.6) +
        scale_fill_manual(values = c(
          setNames("#e41a1c", paste("Up in", comp_label)),
          "Up in Control" = "#377eb8"
        )) +
        geom_text(aes(label = Count), vjust = -0.5, size = 5) +
        theme_classic(base_size = 12) +
        labs(
          title = "Distribution of Significantly Differentially Expressed Proteins",
          subtitle = sprintf("Total significant: %d proteins", sum(regulation_summary$Count)),
          x = "",
          y = "Number of Proteins"
        ) +
        theme(legend.position = "none")
      
      ggsave(file.path(plots_dir, 
                       paste0("Significant_proteins_distribution_", 
                              comp_label, "_", timestamp, ".pdf")), 
             plot = p_bar, width = 6, height = 5)
    }
    
    # Enhanced Heatmap
    sig_proteins <- tt_enhanced %>% 
      filter(grepl("Significantly up", Regulation)) %>% 
      pull(ProteinID)
    
    min_proteins <- 2
    max_proteins <- 100
    
    if (length(sig_proteins) >= min_proteins) {
      top_sig <- head(sig_proteins, max_proteins)
      expr_sig <- exprs(msn)[top_sig, , drop = FALSE]
      
      # Row annotation
      row_anno_data <- tt_enhanced %>%
        filter(ProteinID %in% top_sig) %>%
        dplyr::select(ProteinID, Direction, logFC)
      
      row_ha <- rowAnnotation(
        Direction = row_anno_data$Direction[match(rownames(expr_sig), 
                                                  row_anno_data$ProteinID)],
        logFC = row_anno_data$logFC[match(rownames(expr_sig), 
                                          row_anno_data$ProteinID)],
        col = list(
          Direction = c(setNames("#e41a1c", paste("Up in", comp_label)), 
                        setNames("#377eb8", "Up in Control")),
          logFC = circlize::colorRamp2(
            c(-max(abs(row_anno_data$logFC)), 0, max(abs(row_anno_data$logFC))), 
            c("#377eb8", "white", "#e41a1c"))
        )
      )
      
      # Column annotation
      col_ha <- HeatmapAnnotation(
        Group = pData(msn)$Group,
        col = list(Group = c(Control = "#377eb8", Test = "#e41a1c"))
      )
      
      pdf(file.path(plots_dir, 
                    paste0("Heatmap_", comp_label, "_vs_Control_", 
                           tag, "_", timestamp, ".pdf")),
          width = 10, height = min(15, 4 + length(top_sig) * 0.2))
      print(ComplexHeatmap::Heatmap(
        expr_sig,
        name = "log2 intensity",
        top_annotation = col_ha,
        right_annotation = row_ha,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        column_split = pData(msn)$Group,
        row_title = sprintf("Top %d Significant Proteins", length(top_sig)),
        column_title = paste(comp_clean, "-", tag)
      ))
      dev.off()
    } else {
      warn(sprintf("Too few significant proteins (%d) for heatmap in %s", 
                   length(sig_proteins), comp_clean))
    }
  }
  
  return(list(fit = fit2, results = all_results))
}

## Run Analysis --------------------------------------------------
msg("\n", strrep("=", 60))
msg("ANALYZING BICARBONATE DATA")
msg(strrep("=", 60))
results <- run_limma_analysis(msn_set, "Bicarbonate", timestamp)
gc()
monitor_memory("after limma analysis")

## Enhanced GO Enrichment Analysis --------------------------------
run_go_enrichment <- function(results_list, msn_current, tag_str, timestamp) {
  msg(sprintf("\nPerforming GO enrichment analysis (%s)", tag_str))
  
  # Read and parse GO annotations
  go_raw <- readr::read_tsv(go_file, col_names = FALSE, na = c("", "NA"), 
                            show_col_types = FALSE) %>%
    rename(ProteinID = X1) %>%
    tidyr::pivot_longer(-ProteinID, names_to = "col", values_to = "GO_score", 
                        values_drop_na = TRUE) %>%
    dplyr::select(-col) %>%
    tidyr::separate(GO_score, into = c("GO_ID", "score"), sep = "\\|", 
                    convert = TRUE, fill = "right") %>%
    filter(!is.na(GO_ID), score >= go_min_score) %>%
    mutate(Acc = acc_from(ProteinID))
  
  if (nrow(go_raw) == 0) { 
    warn("No valid GO annotations parsed - check format and score threshold")
    return(invisible(NULL))
  }
  
  # Filter for Biological Process terms
  bp_ids <- AnnotationDbi::select(GO.db, keys = unique(go_raw$GO_ID), 
                                  columns = "ONTOLOGY", keytype = "GOID") %>%
    filter(ONTOLOGY == "BP") %>% 
    pull(GOID)
  go_bp <- go_raw %>% 
    filter(GO_ID %in% bp_ids)
  
  # Prepare TERM2GENE and TERM2NAME
  TERM2GENE <- distinct(go_bp, GO_ID, Acc)
  colnames(TERM2GENE) <- c("term", "gene")
  TERM2GENE$term <- as.character(TERM2GENE$term)
  TERM2GENE$gene <- as.character(TERM2GENE$gene)
  
  TERM2NAME <- AnnotationDbi::select(GO.db, keys = unique(TERM2GENE$term), 
                                     columns = "TERM", keytype = "GOID") %>% 
    distinct()
  colnames(TERM2NAME) <- c("term", "name")
  TERM2NAME$term <- as.character(TERM2NAME$term)
  TERM2NAME$name <- as.character(TERM2NAME$name)
  
  # Define universe (all tested proteins with GO annotations)
  universe_tp <- intersect(acc_from(rownames(exprs(msn_current))), 
                           unique(TERM2GENE$gene))
  universe_tp <- unique(as.character(universe_tp))
  msg(sprintf("GO annotation summary (%s): %d proteins in universe", 
              tag_str, length(universe_tp)))
  
  # Safe enricher function with error handling
  safe_enricher <- function(genes, universe_arg, direction_lab) {
    genes <- intersect(unique(as.character(genes)), universe_arg)
    if (length(genes) < 5) {
      msg(sprintf("Too few genes (%d) for %s enrichment", 
                  length(genes), direction_lab))
      return(NULL)
    }
    tryCatch(
      clusterProfiler::enricher(
        gene = genes, 
        universe = universe_arg,
        TERM2GENE = TERM2GENE, 
        TERM2NAME = TERM2NAME,
        pAdjustMethod = "BH", 
        pvalueCutoff = 1, 
        qvalueCutoff = 1,
        minGSSize = 5, 
        maxGSSize = 5000
      ),
      error = function(e) {
        warn(sprintf("enricher() failed for %s (%s): %s", 
                     direction_lab, tag_str, e$message))
        NULL
      }
    )
  }
  
  # Store ORA file paths for comparison
  ora_paths <- list()
  
  # Process each comparison
  for (comp_name in names(results_list)) {
    msg(sprintf("\nRunning GO for %s vs Control (%s)", comp_name, tag_str))
    df <- results_list[[comp_name]]
    
    # Get tested IDs and create contrast-specific universe
    tested_ids <- acc_from(df$ProteinID)
    universe_contrast <- intersect(tested_ids, unique(TERM2GENE$gene))
    
    if (length(universe_contrast) < 10) {
      warn(sprintf("Universe too small for %s (%s): %d proteins", 
                   comp_name, tag_str, length(universe_contrast)))
      next
    }
    
    msg(sprintf("Contrast-specific universe: %d proteins", length(universe_contrast)))
    
    # Get UP-regulated proteins
    hits_up <- df %>% 
      filter(Regulation == paste("Significantly up in", comp_name)) %>%
      transmute(Acc = acc_from(ProteinID)) %>% 
      pull(Acc) %>% 
      intersect(universe_contrast)
    
    # Get DOWN-regulated proteins (up in Control)
    hits_dn <- df %>% 
      filter(Regulation == "Significantly up in Control") %>%
      transmute(Acc = acc_from(ProteinID)) %>% 
      pull(Acc) %>% 
      intersect(universe_contrast)
    
    msg(sprintf("Proteins for enrichment: %d UP, %d DOWN", 
                length(hits_up), length(hits_dn)))
    
    # ORA for UP-regulated proteins
    up_path <- NA_character_
    ora_up <- safe_enricher(hits_up, universe_contrast, 
                            sprintf("UP_in_%s", comp_name))
    if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
      up_path <- file.path(tables_dir, 
                           sprintf("ORA_BP_UP_in_%s_vs_Control_%s_%s.csv", 
                                   comp_name, tag_str, timestamp))
      write.csv(ora_up@result, up_path, row.names = FALSE)
      
      # Create bar plot for UP
      ora_plot_data <- ora_up@result %>% 
        filter(p.adjust < 0.05) %>% 
        slice_min(p.adjust, n = 20) %>%
        mutate(Description = fct_reorder(Description, -log10(p.adjust)))
      
      if (nrow(ora_plot_data) > 0) {
        p_up <- ggplot(ora_plot_data, aes(x = -log10(p.adjust), y = Description)) +
          geom_col(fill = "#e41a1c", alpha = 0.8) +
          labs(
            title = sprintf("ORA: UP in %s (vs Control, %s)", comp_name, tag_str),
            x = expression(-log[10]~"Adjusted p-value"),
            y = NULL
          ) +
          theme_classic(base_size = 11)
        
        ggsave(file.path(plots_dir, 
                         sprintf("ORA_barplot_UP_in_%s_%s_%s.pdf", 
                                 comp_name, tag_str, timestamp)), 
               plot = p_up, width = 10, height = 8)
      }
    }
    
    # ORA for DOWN-regulated proteins
    down_path <- NA_character_
    ora_dn <- safe_enricher(hits_dn, universe_contrast, 
                            sprintf("DOWN_in_%s", comp_name))
    if (!is.null(ora_dn) && nrow(ora_dn@result) > 0) {
      down_path <- file.path(tables_dir, 
                             sprintf("ORA_BP_DOWN_in_%s_vs_Control_%s_%s.csv", 
                                     comp_name, tag_str, timestamp))
      write.csv(ora_dn@result, down_path, row.names = FALSE)
      
      # Create bar plot for DOWN
      ora_plot_data <- ora_dn@result %>% 
        filter(p.adjust < 0.05) %>% 
        slice_min(p.adjust, n = 20) %>%
        mutate(Description = fct_reorder(Description, -log10(p.adjust)))
      
      if (nrow(ora_plot_data) > 0) {
        p_dn <- ggplot(ora_plot_data, aes(x = -log10(p.adjust), y = Description)) +
          geom_col(fill = "#377eb8", alpha = 0.8) +
          labs(
            title = sprintf("ORA: DOWN in %s (vs Control, %s)", comp_name, tag_str),
            x = expression(-log[10]~"Adjusted p-value"),
            y = NULL
          ) +
          theme_classic(base_size = 11)
        
        ggsave(file.path(plots_dir, 
                         sprintf("ORA_barplot_DOWN_in_%s_%s_%s.pdf", 
                                 comp_name, tag_str, timestamp)), 
               plot = p_dn, width = 10, height = 8)
      }
    }
    
    ora_paths[[comp_name]] <- list(up = up_path, down = down_path)
    
    # GSEA (ranked by t-statistic)
    ranked_df <- df %>%
      mutate(Acc = acc_from(ProteinID)) %>%
      filter(Acc %in% universe_contrast, is.finite(t)) %>%
      group_by(Acc) %>% 
      summarise(t = mean(t), .groups = "drop") %>%
      arrange(desc(t))
    
    if (nrow(ranked_df) >= 10) {
      geneList <- ranked_df$t
      names(geneList) <- ranked_df$Acc
      
      gsea <- tryCatch(
        clusterProfiler::GSEA(geneList, 
                              TERM2GENE = TERM2GENE, 
                              TERM2NAME = TERM2NAME,
                              pvalueCutoff = 1, 
                              verbose = FALSE, 
                              minGSSize = 5, 
                              maxGSSize = 5000),
        error = function(e) {
          warn(sprintf("GSEA failed for %s (%s): %s", 
                       comp_name, tag_str, e$message))
          NULL
        }
      )
      
      if (!is.null(gsea) && nrow(gsea@result) > 0) {
        # Save GSEA results
        write.csv(gsea@result, 
                  file.path(tables_dir, 
                            sprintf("GSEA_BP_%s_vs_Control_%s_%s.csv", 
                                    comp_name, tag_str, timestamp)), 
                  row.names = FALSE)
        
        # Create GSEA dotplot
        p_gsea_dot <- enrichplot::dotplot(gsea, showCategory = 20) + 
          ggtitle(sprintf("GSEA: %s vs Control (%s)", comp_name, tag_str))
        
        ggsave(file.path(plots_dir, 
                         sprintf("GSEA_dotplot_%s_vs_Control_%s_%s.pdf", 
                                 comp_name, tag_str, timestamp)),
               plot = p_gsea_dot, width = 10, height = 8)
        
        # Create NES bar plot
        gsea_df <- gsea@result %>%
          arrange(desc(NES))
        
        # Select top positive and negative terms
        top_n <- 10
        top_terms <- bind_rows(
          gsea_df %>% filter(NES > 0) %>% head(top_n),
          gsea_df %>% filter(NES < 0) %>% tail(top_n)
        )
        
        if (nrow(top_terms) > 0) {
          p_gsea_bar <- ggplot(top_terms, aes(x = NES, y = fct_reorder(Description, NES))) +
            geom_col(aes(fill = NES > 0)) +
            scale_fill_manual(
              values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
              name = "Enrichment Direction",
              labels = c("FALSE" = "Enriched in Control", 
                         "TRUE" = paste("Enriched in", comp_name))
            ) +
            labs(
              title = "GSEA: Top Enriched Pathways by NES",
              subtitle = sprintf("%s vs Control (%s)", comp_name, tag_str),
              x = "Normalized Enrichment Score (NES)",
              y = "GO Term"
            ) +
            theme_classic(base_size = 12) +
            theme(legend.position = "bottom")
          
          ggsave(file.path(plots_dir, 
                           sprintf("GSEA_barplot_NES_%s_vs_Control_%s_%s.pdf", 
                                   comp_name, tag_str, timestamp)),
                 plot = p_gsea_bar, width = 10, height = 8)
        }
      }
    } else {
      warn(sprintf("Too few proteins for GSEA in %s", comp_name))
    }
  }
  
  # Create comparison table if both UP and DOWN results exist
  for (comp_name in names(ora_paths)) {
    paths <- ora_paths[[comp_name]]
    if (is.na(paths$up) || is.na(paths$down)) next
    
    ora_up_df <- tryCatch(read.csv(paths$up), error = function(e) NULL)
    ora_dn_df <- tryCatch(read.csv(paths$down), error = function(e) NULL)
    
    if (!is.null(ora_up_df) && !is.null(ora_dn_df)) {
      # Create side-by-side comparison
      comparison_df <- data.frame(
        UP_in_Treatment = c(utils::head(ora_up_df$Description, 10), 
                            rep("", max(0, 10 - nrow(ora_up_df)))),
        UP_pvalue = c(utils::head(ora_up_df$p.adjust, 10), 
                      rep(NA, max(0, 10 - nrow(ora_up_df)))),
        DOWN_in_Treatment = c(utils::head(ora_dn_df$Description, 10), 
                              rep("", max(0, 10 - nrow(ora_dn_df)))),
        DOWN_pvalue = c(utils::head(ora_dn_df$p.adjust, 10), 
                        rep(NA, max(0, 10 - nrow(ora_dn_df))))
      )
      
      write.csv(comparison_df, 
                file.path(tables_dir, 
                          sprintf("ORA_Comparison_%s_%s_%s.csv", 
                                  comp_name, tag_str, timestamp)), 
                row.names = FALSE)
    }
  }
  
  invisible(NULL)
}

# Run GO enrichment
run_go_enrichment(results$results, msn_set, "Bicarbonate", timestamp)
## Enhanced GO Analysis with Filtering and Umbrella Terms -------
run_enhanced_go_analysis <- function(results_list, msn_current, tag_str, timestamp) {
  msg("\n", strrep("=", 60))
  msg(sprintf("ENHANCED GO ANALYSIS WITH CLUSTERING + UMBRELLAS (%s)", tag_str))
  msg(strrep("=", 60))
  
  # Define negative filter list for eukaryotic terms (irrelevant for prokaryotes)
  eukaryotic_filters <- c(
    "organelle", "cilium", "cilia", "flagell", "nucleus", "nuclear",
    "chromatin", "chromosome", "histone", "nucleosome", "mitochondri",
    "endoplasmic reticulum", "golgi", "peroxisome", "lysosome",
    "endosome", "cytoskeleton", "microtubule", "actin", "tubulin",
    "spindle", "centrosome", "centriole", "kinetochore", "telomere",
    "spliceosome", "proteasome", "ubiquitin", "sumo", "autophagy",
    "endocytosis", "exocytosis", "vesicle", "clathrin", "adaptin",
    "kinesin", "dynein", "myosin", "lamina", "heterochromatin",
    "euchromatin", "nucleolus", "nucleolar", "ribosome biogenesis",
    "rRNA processing", "mRNA splicing", "intron", "exon"
  )
  
  # Function to check if a GO term contains eukaryotic keywords
  is_eukaryotic_term <- function(term_descriptions) {
    pattern <- paste(eukaryotic_filters, collapse = "|")
    grepl(pattern, term_descriptions, ignore.case = TRUE)
  }
  
  # Define umbrella terms and their children - adapted for bicarbonate metabolism
  umbrella_mappings <- list(
    # Transport processes
    "transmembrane transport" = c(
      "transport",
      "transmembrane transport",
      "monoatomic ion transport",
      "monoatomic cation transport", 
      "monoatomic anion transport",
      "monoatomic ion transmembrane transport",
      "monoatomic cation transmembrane transport",
      "inorganic ion transmembrane transport",
      "inorganic cation transmembrane transport",
      "inorganic anion transport",
      "metal ion transport",
      "ion transmembrane transport",
      "cation transmembrane transport",
      "anion transmembrane transport",
      "organic anion transport",
      "transition metal ion transport",
      "establishment of localization",
      "cellular localization",
      "carboxylic acid transmembrane transport",
      "organic acid transport",
      "protein transport",
      "protein transmembrane transport",
      "macromolecule localization",
      "localization",
      "establishment of localization"
    ),
    
    # Carbon metabolism - relevant for bicarbonate
    "carbon metabolism" = c(
      "carbon utilization",
      "one-carbon metabolic process",
      "carbon fixation",
      "organic acid metabolic process",
      "carboxylic acid metabolic process",
      "organic acid biosynthetic process",
      "carboxylic acid biosynthetic process"
    ),
    
    # Energy metabolism
    "energy generation" = c(
      "generation of precursor metabolites and energy",
      "nucleoside triphosphate metabolic process",
      "ribonucleotide metabolic process",
      "purine ribonucleotide metabolic process",
      "ribose phosphate metabolic process",
      "energy derivation by oxidation of organic compounds",
      "anaerobic respiration",
      "cellular respiration",
      "ATP metabolic process",
      "ATP biosynthetic process"
    ),
    
    # Electron transport
    "electron transport and respiration" = c(
      "respiratory electron transport chain",
      "electron transport chain",
      "cellular respiration",
      "oxidation-reduction process"
    ),
    
    # Homeostasis
    "homeostasis processes" = c(
      "homeostatic process",
      "chemical homeostasis",
      "cellular homeostasis",
      "intracellular chemical homeostasis",
      "pH regulation",
      "cellular pH homeostasis"
    ),
    
    # Lipid metabolism
    "lipid metabolism" = c(
      "lipid biosynthetic process",
      "lipid transport",
      "lipid localization",
      "phospholipid metabolic process",
      "phospholipid biosynthetic process",
      "lipid metabolic process",
      "fatty acid metabolic process",
      "fatty acid biosynthetic process"
    ),
    
    # Transcription regulation
    "transcription regulation" = c(
      "regulation of cellular process",
      "regulation of biological process",
      "positive regulation of cellular process",
      "positive regulation of metabolic process",
      "positive regulation of macromolecule biosynthetic process",
      "positive regulation of biosynthetic process",
      "negative regulation of metabolic process",
      "regulation of RNA metabolic process",
      "regulation of DNA-templated transcription",
      "regulation of RNA biosynthetic process",
      "positive regulation of macromolecule metabolic process",
      "regulation of gene expression",
      "regulation of transcription, DNA-templated"
    ),
    
    # Protein synthesis
    "protein synthesis" = c(
      "translation",
      "peptide biosynthetic process",
      "amide biosynthetic process",
      "peptide metabolic process",
      "cellular amide metabolic process",
      "L-amino acid metabolic process"
    ),
    
    # Cell motility
    "cell motility metabolism" = c(
      "cell projection organization",
      "bacterial-type flagellum-dependent cell motility",
      "archaeal or bacterial-type flagellum-dependent cell motility",
      "cell motility",
      "locomotion",
      "membraneless organelle assembly",
      "cell projection assembly",
      "cellular component organization",
      "cellular component biogenesis",
      "cellular component assembly"
    ),
    
    # Response to stimulus
    "response to stimulus" = c(
      "taxis",
      "chemotaxis",
      "response to chemical",
      "response to external stimulus",
      "cellular response to stimulus"
    ),
    
    # Small molecule metabolism
    "small molecule metabolic process" = c(
      "small molecule biosynthetic process",
      "small molecule catabolic process",
      "alcohol metabolic process",
      "amino acid metabolic process",
      "amino acid biosynthetic process",
      "nucleobase-containing small molecule metabolic process"
    ),
    
    # Cofactor metabolism
    "cofactor metabolic process" = c(
      "cofactor biosynthetic process",
      "coenzyme metabolic process",
      "coenzyme biosynthetic process",
      "vitamin metabolic process",
      "vitamin biosynthetic process",
      "water-soluble vitamin metabolic process"
    )
  )
  
  # Function to apply umbrella term grouping - handles both ORA and GSEA
  apply_umbrella_terms <- function(go_df, method = "GSEA") {
    if (nrow(go_df) == 0) return(go_df)
    
    # Add column to track original terms
    go_df$original_term <- go_df$Description
    go_df$was_grouped <- FALSE
    
    # Apply each umbrella mapping
    for (umbrella in names(umbrella_mappings)) {
      children <- umbrella_mappings[[umbrella]]
      
      # Find matches for child terms
      for (child in children) {
        # Use exact matching
        matches <- grepl(paste0("^", child, "$"), go_df$Description, ignore.case = TRUE)
        if (any(matches)) {
          go_df$Description[matches] <- umbrella
          go_df$was_grouped[matches] <- TRUE
        }
      }
    }
    
    # Merge duplicate entries based on method
    if (method == "GSEA") {
      go_df_merged <- go_df %>%
        group_by(Description) %>%
        summarise(
          ID = dplyr::first(ID),
          setSize = max(setSize),  # Take the largest set size
          enrichmentScore = dplyr::first(enrichmentScore),
          NES = dplyr::first(NES[which.max(abs(NES))]),  # Take the most extreme NES
          pvalue = min(pvalue),  # Take the best p-value
          p.adjust = min(p.adjust),  # Take the best adjusted p-value
          qvalue = min(qvalue),
          rank = dplyr::first(rank),
          leading_edge = dplyr::first(leading_edge),
          core_enrichment = dplyr::first(core_enrichment),
          original_term = paste(unique(original_term), collapse = " | "),
          was_grouped = any(was_grouped),
          n_merged = n(),
          .groups = 'drop'
        )
    } else if (method == "ORA") {
      go_df_merged <- go_df %>%
        group_by(Description) %>%
        summarise(
          ID = dplyr::first(ID),
          GeneRatio = dplyr::first(GeneRatio),
          BgRatio = dplyr::first(BgRatio),
          pvalue = min(pvalue),
          p.adjust = min(p.adjust),
          qvalue = min(qvalue),
          geneID = paste(unique(unlist(strsplit(geneID, "/"))), collapse = "/"),
          Count = sum(Count),
          original_term = paste(unique(original_term), collapse = " | "),
          was_grouped = any(was_grouped),
          n_merged = n(),
          .groups = 'drop'
        )
    }
    
    return(go_df_merged)
  }
  
  # Function to process GO results with filtering and grouping
  process_go_results <- function(go_result, method = "GSEA", nes_threshold = 1) {
    if (is.null(go_result) || nrow(go_result@result) == 0) {
      return(NULL)
    }
    
    go_df <- go_result@result
    
    # Step 1: Filter by NES if GSEA, by p.adjust if ORA
    if (method == "GSEA") {
      go_filtered <- go_df %>%
        filter(abs(NES) >= nes_threshold)
      msg(sprintf("  Filtered by |NES|≥%.1f: %d -> %d terms", 
                  nes_threshold, nrow(go_df), nrow(go_filtered)))
    } else {
      # For ORA, filter by p.adjust
      go_filtered <- go_df %>%
        filter(p.adjust < 0.05)
      msg(sprintf("  Filtered by p.adjust<0.05: %d -> %d terms", 
                  nrow(go_df), nrow(go_filtered)))
    }
    
    # Step 2: Remove eukaryotic terms
    go_filtered <- go_filtered %>%
      filter(!is_eukaryotic_term(Description))
    msg(sprintf("  After removing eukaryotic terms: %d terms", nrow(go_filtered)))
    
    # Step 3: Apply umbrella grouping
    if (nrow(go_filtered) > 0) {
      go_grouped <- apply_umbrella_terms(go_filtered, method = method)
      msg(sprintf("  After umbrella grouping: %d terms", nrow(go_grouped)))
      return(go_grouped)
    } else {
      return(go_filtered)
    }
  }
  
  # Load GO annotations
  go_raw <- readr::read_tsv(go_file, col_names = FALSE, na = c("", "NA"), 
                            show_col_types = FALSE) %>%
    rename(ProteinID = X1) %>%
    tidyr::pivot_longer(-ProteinID, names_to = "col", values_to = "GO_score", 
                        values_drop_na = TRUE) %>%
    dplyr::select(-col) %>%
    tidyr::separate(GO_score, into = c("GO_ID", "score"), sep = "\\|", 
                    convert = TRUE, fill = "right") %>%
    filter(!is.na(GO_ID), score >= go_min_score) %>%
    mutate(Acc = acc_from(ProteinID))
  
  # Get biological process terms
  bp_ids <- AnnotationDbi::select(GO.db, keys = unique(go_raw$GO_ID), 
                                  columns = "ONTOLOGY", keytype = "GOID") %>%
    filter(ONTOLOGY == "BP") %>% 
    pull(GOID)
  go_bp <- go_raw %>% 
    filter(GO_ID %in% bp_ids)
  
  # Prepare TERM2GENE and TERM2NAME
  TERM2GENE <- distinct(go_bp, GO_ID, Acc)
  colnames(TERM2GENE) <- c("term", "gene")
  TERM2GENE$term <- as.character(TERM2GENE$term)
  TERM2GENE$gene <- as.character(TERM2GENE$gene)
  
  TERM2NAME <- AnnotationDbi::select(GO.db, keys = unique(TERM2GENE$term), 
                                     columns = "TERM", keytype = "GOID") %>% 
    distinct()
  colnames(TERM2NAME) <- c("term", "name")
  TERM2NAME$term <- as.character(TERM2NAME$term)
  TERM2NAME$name <- as.character(TERM2NAME$name)
  
  # Process each comparison
  for (comp_name in names(results_list)) {
    msg(sprintf("\nProcessing GO for %s vs Control (%s)", comp_name, tag_str))
    df <- results_list[[comp_name]]
    
    # Prepare universe
    tested_ids <- acc_from(df$ProteinID)
    universe_contrast <- intersect(tested_ids, unique(TERM2GENE$gene))
    
    if (length(universe_contrast) < 10) {
      warn(sprintf("Universe too small for %s: %d proteins", comp_name, length(universe_contrast)))
      next
    }
    
    msg(sprintf("Contrast-specific universe: %d proteins", length(universe_contrast)))
    
    ## ========== ORA ANALYSIS ==========
    msg("\n--- Running ORA ---")
    
    # Get UP-regulated proteins
    hits_up <- df %>% 
      filter(Regulation == paste("Significantly up in", comp_name)) %>%
      transmute(Acc = acc_from(ProteinID)) %>% 
      pull(Acc) %>% 
      intersect(universe_contrast)
    
    # Get DOWN-regulated proteins
    hits_dn <- df %>% 
      filter(Regulation == "Significantly up in Control") %>%
      transmute(Acc = acc_from(ProteinID)) %>% 
      pull(Acc) %>% 
      intersect(universe_contrast)
    
    msg(sprintf("Proteins for ORA: %d UP, %d DOWN", length(hits_up), length(hits_dn)))
    
    # Run enhanced ORA for UP
    if (length(hits_up) >= 5) {
      ora_up_result <- tryCatch(
        clusterProfiler::enricher(
          gene = hits_up, universe = universe_contrast,
          TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
          minGSSize = 5, maxGSSize = 5000
        ),
        error = function(e) { warn(sprintf("ORA UP failed: %s", e$message)); NULL }
      )
      
      if (!is.null(ora_up_result) && nrow(ora_up_result@result) > 0) {
        ora_up_enhanced <- process_go_results(ora_up_result, method = "ORA")
        
        if (!is.null(ora_up_enhanced) && nrow(ora_up_enhanced) > 0) {
          write.csv(ora_up_enhanced, 
                    file.path(tables_dir, sprintf("ORA_ENHANCED_UP_%s_%s_%s.csv", 
                                                  comp_name, tag_str, timestamp)), 
                    row.names = FALSE)
          msg(sprintf("  Saved enhanced ORA UP results: %d terms", nrow(ora_up_enhanced)))
        }
      }
    }
    
    # Run enhanced ORA for DOWN
    if (length(hits_dn) >= 5) {
      ora_dn_result <- tryCatch(
        clusterProfiler::enricher(
          gene = hits_dn, universe = universe_contrast,
          TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
          minGSSize = 5, maxGSSize = 5000
        ),
        error = function(e) { warn(sprintf("ORA DOWN failed: %s", e$message)); NULL }
      )
      
      if (!is.null(ora_dn_result) && nrow(ora_dn_result@result) > 0) {
        ora_dn_enhanced <- process_go_results(ora_dn_result, method = "ORA")
        
        if (!is.null(ora_dn_enhanced) && nrow(ora_dn_enhanced) > 0) {
          write.csv(ora_dn_enhanced, 
                    file.path(tables_dir, sprintf("ORA_ENHANCED_DOWN_%s_%s_%s.csv", 
                                                  comp_name, tag_str, timestamp)), 
                    row.names = FALSE)
          msg(sprintf("  Saved enhanced ORA DOWN results: %d terms", nrow(ora_dn_enhanced)))
        }
      }
    }
    
    
    ## ========== GSEA ANALYSIS ==========
    msg("\n--- Running GSEA ---")
    
    # Prepare ranked list
    ranked_df <- df %>%
      mutate(Acc = acc_from(ProteinID)) %>%
      filter(Acc %in% universe_contrast, is.finite(t)) %>%
      group_by(Acc) %>% 
      summarise(t = mean(t), .groups = "drop") %>%
      arrange(desc(t))
    
    if (nrow(ranked_df) >= 10) {
      geneList <- ranked_df$t
      names(geneList) <- ranked_df$Acc
      
      gsea_result <- tryCatch(
        clusterProfiler::GSEA(geneList, 
                              TERM2GENE = TERM2GENE, 
                              TERM2NAME = TERM2NAME,
                              pvalueCutoff = 1, 
                              verbose = FALSE, 
                              minGSSize = 5, 
                              maxGSSize = 5000),
        error = function(e) { 
          warn(sprintf("GSEA failed: %s", e$message))
          NULL 
        }
      )
      
      if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
        gsea_enhanced <- process_go_results(gsea_result, method = "GSEA", nes_threshold = 1)
        
        if (!is.null(gsea_enhanced) && nrow(gsea_enhanced) > 0) {
          # Save enhanced GSEA results
          write.csv(gsea_enhanced, 
                    file.path(tables_dir, sprintf("GSEA_ENHANCED_%s_%s_%s.csv", 
                                                  comp_name, tag_str, timestamp)), 
                    row.names = FALSE)
          msg(sprintf("  Saved enhanced GSEA results: %d terms", nrow(gsea_enhanced)))
          
          # Save grouping report if applicable
          if (any(gsea_enhanced$was_grouped)) {
            grouping_report <- gsea_enhanced %>%
              filter(was_grouped) %>%
              dplyr::select(Description, original_term, NES, p.adjust, n_merged)
            
            write.csv(grouping_report,
                      file.path(tables_dir, sprintf("GO_grouping_report_%s_%s_%s.csv",
                                                    comp_name, tag_str, timestamp)),
                      row.names = FALSE)
          }
          
          # Create visualization
          top_n <- 10
          top_terms <- bind_rows(
            gsea_enhanced %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_n),
            gsea_enhanced %>% filter(NES < 0) %>% arrange(NES) %>% head(top_n)
          )
          
          if (nrow(top_terms) > 0) {
            top_terms <- top_terms %>%
              mutate(
                Description_short = ifelse(nchar(Description) > 60,
                                           paste0(substr(Description, 1, 57), "..."),
                                           Description),
                Description_short = stringr::str_to_title(Description_short)
              )
            
            p_enhanced <- ggplot(top_terms, 
                                 aes(x = NES, y = fct_reorder(Description_short, NES))) +
              geom_col(aes(fill = NES > 0), width = 0.7) +
              scale_fill_manual(
                values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
                name = element_blank(),
                labels = c("FALSE" = "Enriched in Control (with acetate+bicarbonate)", 
                           "TRUE" = "Enriched in Test (without acetate)")
              ) +
              scale_x_continuous(
                limits = c(min(top_terms$NES) - 0.5, max(top_terms$NES) + 0.5)
              ) +
              labs(
                title = "Enhanced GO Analysis: Filtered & Grouped Terms",
                subtitle = sprintf("%s vs Control (%s) | GSEA with |NES|≥1", 
                                   comp_name, tag_str),
                x = "Normalized Enrichment Score (NES)",
                y = NULL
              ) +
              theme_classic(base_size = 14) +
              theme(
                legend.position = "bottom",
                legend.box.margin = margin(t = -10),
                axis.text.y = element_text(size = 12, colour = "black"),
                panel.grid.major.x = element_line(color = "grey90", size = 0.2),
                plot.title = element_text(face = "bold", size = 16),
                plot.subtitle = element_text(size = 12, color = "grey30")
              ) +
              geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5)
            
            ggsave(
              file.path(plots_dir, sprintf("GO_ENHANCED_%s_%s_%s.pdf", 
                                           comp_name, tag_str, timestamp)),
              plot = p_enhanced,
              width = 12,
              height = 10
            )
          }
        }
      }
    }
  }
  
  invisible(NULL)
}

# Add this function call after the standard GO enrichment in your script:
run_enhanced_go_analysis(results$results, msn_set, "Bicarbonate", timestamp)

## Dimensionality Reduction (PCA and UMAP) ----------------------
run_dim_reduction <- function(msn, tag, timestamp) {
  msg("Running PCA and UMAP on", tag)
  
  data_for_dim_red <- t(exprs(msn))
  sample_info <- pData(msn)
  
  # Define consistent color scheme
  group_colors <- c(Control = "#377eb8", Test = "#e41a1c")
  
  # Create sample labels
  sample_labels <- sub("^(R[0-9]+).*", "\\1", rownames(data_for_dim_red))
  
  # PCA
  pca_result <- prcomp(data_for_dim_red, scale. = TRUE, center = TRUE)
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  pca_data$Group <- sample_info$Group
  pca_data$Label <- sample_labels
  
  # UMAP
  set.seed(42)
  n_neighbors <- min(10, nrow(data_for_dim_red) - 1)
  umap_result <- umap(data_for_dim_red, 
                      n_neighbors = n_neighbors, 
                      min_dist = 0.1, 
                      n_components = 2,
                      metric = "euclidean",
                      n_epochs = 200)
  
  umap_data <- as.data.frame(umap_result)
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  umap_data$Group <- sample_info$Group
  umap_data$Label <- sample_labels
  
  # PCA Plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Group), size = 6, alpha = 0.8) +
    geom_text_repel(
      aes(label = Label),
      size = 4,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.3,
      segment.alpha = 0.5,
      show.legend = FALSE,
      seed = 42
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      legend.position = "bottom",
      aspect.ratio = 1
    ) +
    labs(
      title = paste("PCA —", tag),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    scale_color_manual(values = group_colors) +
    coord_fixed(ratio = 1)
  
  # UMAP Plot
  umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Group), size = 6, alpha = 0.8) +
    geom_text_repel(
      aes(label = Label),
      size = 4,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.3,
      segment.alpha = 0.5,
      show.legend = FALSE,
      seed = 42
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      legend.position = "bottom",
      aspect.ratio = 1
    ) +
    labs(
      title = paste("UMAP —", tag),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    scale_color_manual(values = group_colors) +
    coord_fixed(ratio = 1)
  
  # Save individual plots
  ggsave(file.path(plots_dir, paste0("PCA_", tag, "_", timestamp, ".pdf")), 
         plot = pca_plot, width = 8, height = 7)
  ggsave(file.path(plots_dir, paste0("UMAP_", tag, "_", timestamp, ".pdf")), 
         plot = umap_plot, width = 8, height = 7)
  
  # Combined plot
  combined_plot <- pca_plot + umap_plot + 
    plot_layout(guides = "collect", ncol = 2) & 
    theme(legend.position = 'bottom')
  
  ggsave(file.path(plots_dir, paste0("PCA_UMAP_", tag, "_", timestamp, ".png")), 
         plot = combined_plot, width = 14, height = 7, dpi = 300, bg = "white")
  
  msg("Dimensionality reduction plots created")
}

# Run dimensionality reduction
run_dim_reduction(msn_set, "Bicarbonate", timestamp)

## Combine All Results -------------------------------------------
all_results_combined <- dplyr::bind_rows(
  lapply(names(results$results), function(x) results$results[[x]])
)

write_xlsx(all_results_combined, 
           file.path(tables_dir, paste0("ALL_RESULTS_COMBINED_", timestamp, ".xlsx")))

msg(sprintf("Combined results saved: %d proteins", nrow(all_results_combined)))

## Clean up large objects ----------------------------------------
cleanup_large_objects <- function(threshold_mb = 100) {
  obj_sizes <- sapply(ls(envir = .GlobalEnv), function(x) {
    object.size(get(x, envir = .GlobalEnv)) / 1024^2
  })
  large_objs <- names(obj_sizes[obj_sizes > threshold_mb])
  if (length(large_objs) > 0) {
    msg(sprintf("Large objects (>%d Mb): %s", threshold_mb, 
                paste(large_objs, round(obj_sizes[large_objs]), 
                      sep = "=", collapse = ", ")))
  }
}

cleanup_large_objects()

# Force garbage collection
gc()
final_memory <- monitor_memory("final")

## Session Info --------------------------------------------------
sink(file.path(out_dir, paste0("session_info_", timestamp, ".txt")))
cat("Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Configuration:\n")
cat("  Output directory:", out_dir, "\n")
cat("  Raw data file:", raw_file, "\n")
cat("  GO annotation file:", go_file, "\n")
cat("  Significance thresholds: adj.P <=", sig_adjP, ", |logFC| >=", sig_logFC, "\n")
cat("  GO score threshold:", go_min_score, "\n")
cat("  Final memory usage:", round(final_memory, 1), "Mb\n\n")
print(sessionInfo())
sink()

msg("")
msg(strrep("=", 60))
msg("ANALYSIS COMPLETE!")
msg(sprintf("All results saved to: %s", out_dir))
msg(strrep("=", 60))

