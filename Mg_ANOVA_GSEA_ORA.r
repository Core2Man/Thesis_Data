#!/usr/bin/env Rscript

## --------------------------------------------------------------
## Mg Proteomics Analysis – Final Version (No Combine)
## Experimental Design:
## - R1-3: Low Mg (both F and MP)
## - R4-6: No Mg at Final (F), Normal Mg at Midpoint (MP)
## - R7-9: Control with Normal Mg (both F and MP)
## --------------------------------------------------------------

## Packages ------------------------------------------------------
libs <- c("MSnbase", "impute", "limma", "readr", "dplyr", "tibble",
          "ComplexHeatmap", "ggplot2", "clusterProfiler", "tidyr",
          "GO.db", "forcats", "yaml", "ggrepel", "showtext", "foreach", "doParallel",
          "AnnotationDbi", "writexl", "uwot", "patchwork", "readxl",
          "circlize", "enrichplot", "fgsea", "BiocParallel", "parallel", "future", "furrr")

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
`%||%` <- function(x, y) if (is.null(x)) y else x  # Coalesce

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "...", ..., "\n")
warn <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "!!! WARNING:", ..., "\n")

## Configuration -------------------------------------------------
cfg_file <- "/Users/sysoevm/Library/CloudStorage/Dropbox/Pascal Lab/Desulfuromonas/Mg Chapter/Stat analysis/config.yml"
if (!file.exists(cfg_file)) stop("Config file not found: ", cfg_file)

cfg <- yaml::read_yaml(cfg_file)
out_dir      <- cfg$out_dir
raw_file     <- cfg$raw_file
go_file      <- cfg$go_file
sig_adjP     <- cfg$sig_adjP %||% 0.05
#sig_logFC    <- cfg$sig_logFC %||% 1
go_min_score <- cfg$go_min_score %||% 0.5
#id_regex     <- cfg$id_extraction_regex %||% "^([^|]+).*$"
#id_regex <- "^tr\\|([^|]+)\\|.*$"
id_regex <- "^(?:sp|tr)\\|([^|]+)\\|.*$"
acc_from <- function(x) {
  out <- sub(id_regex, "\\1", x)
  need_fallback <- out == x
  out[need_fallback] <- sub("^([^|]+).*$", "\\1", x[need_fallback])
  out
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(out_dir, "plots")
tables_dir <- file.path(out_dir, "tables")
dir.create(plots_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)

# Validate inputs
for (f in c(raw_file, go_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

tryCatch({
  test <- readxl::read_excel(raw_file, n_max = 2)
}, error = function(e) stop("Failed to read Excel file: ", e$message))

msg("Configuration loaded and validated")
monitor_memory <- function(label = "") {
  gc_result <- gc(reset = TRUE)
  mem_used <- sum(gc_result[, 2])  # Memory in Mb
  msg(sprintf("Memory check %s: %.1f Mb used", label, mem_used))
  return(mem_used)
}
configure_parallel <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- min(detectCores() - 1, 4)  # Leave one core free, max 4
  }
  
  # For BiocParallel (used by many Bioconductor packages)
  if (Sys.info()["sysname"] == "Windows") {
    register(SnowParam(workers = n_cores))
  } else {
    register(MulticoreParam(workers = n_cores))
  }
  
  # For future/furrr (modern tidyverse parallel)
  plan(multisession, workers = n_cores)
  
  msg(sprintf("Parallel processing configured with %d cores", n_cores))
  return(n_cores)
}

n_cores <- configure_parallel()
# Function to clean up large objects
cleanup_large_objects <- function(threshold_mb = 100) {
  obj_sizes <- sapply(ls(envir = .GlobalEnv), function(x) {
    object.size(get(x, envir = .GlobalEnv)) / 1024^2
  })
  large_objs <- names(obj_sizes[obj_sizes > threshold_mb])
  if (length(large_objs) > 0) {
    msg(sprintf("Large objects (>%d Mb): %s", threshold_mb, 
                paste(large_objs, round(obj_sizes[large_objs]), sep = "=", collapse = ", ")))
  }
}
## Data Loading and Preparation ----------------------------------
msg("Reading intensity matrix from:", raw_file)
tab <- readxl::read_excel(raw_file) %>%
  column_to_rownames("Accession") %>%
  as.matrix()
rownames(tab) <- acc_from(rownames(tab))
monitor_memory("after data load")
msg(sprintf("Loaded data: %d proteins x %d samples", nrow(tab), ncol(tab)))

data_range <- range(tab, na.rm = TRUE)
msg(sprintf("Data range: [%.2f, %.2f]", data_range[1], data_range[2]))
if (data_range[1] < 0 || data_range[2] <= 50) {
  msg("Data likely log-scale (negatives or small dynamic range)")
} else if (data_range[2] >= 1e4) {
  warn("Data likely raw (very large values) — log2 transform recommended before limma")
}

F_cols  <- grep("F$", colnames(tab), value = TRUE)
MP_cols <- grep("MP$", colnames(tab), value = TRUE)

if (length(F_cols) == 0 || length(MP_cols) == 0) {
  stop("Could not find columns ending in 'F' or 'MP' – check naming convention.")
}

msg(sprintf("Found %d Final samples and %d Midpoint samples", length(F_cols), length(MP_cols)))

exprs_f  <- tab[, F_cols, drop = FALSE]
exprs_mp <- tab[, MP_cols, drop = FALSE]

## Build MSnSet Objects ------------------------------------------
msg("Building MSnSet objects")
msn_f  <- MSnSet(exprs = exprs_f)
msn_mp <- MSnSet(exprs = exprs_mp)
mg_cols <- c(Control = "#0073C2FF", NormalMg = "#EFC000FF", LowMg = "#CD534CFF", NoMg = "#EFC000FF")

annotate_Mg_groups <- function(msn, timepoint) {
  rn <- sampleNames(msn)
  reactor_nums <- as.integer(gsub("^R([0-9]+)[FM]P?$", "\\1", rn))
  
  if (any(is.na(reactor_nums))) {
    stop("Could not extract reactor numbers from sample names: ", paste(rn[is.na(reactor_nums)], collapse = ", "))
  }
  
  Mg_group <- character(length(reactor_nums))
  
  for (i in seq_along(reactor_nums)) {
    r <- reactor_nums[i]
    if (r %in% 1:3) {
      Mg_group[i] <- "LowMg"
    } else if (r %in% 4:6) {
      Mg_group[i] <- if (timepoint == "Final") "NoMg" else "NormalMg"
    } else if (r %in% 7:9) {
      Mg_group[i] <- "Control"
    } else {
      stop("Unexpected reactor number: ", r)
    }
  }
  
  pData(msn)$Mg <- factor(Mg_group, levels = c("Control", "NormalMg", "LowMg", "NoMg"))
  pData(msn)$Reactor <- factor(reactor_nums)
  pData(msn)$Timepoint <- timepoint
  
  return(msn)
}

msn_f  <- annotate_Mg_groups(msn_f, "Final")
msn_mp <- annotate_Mg_groups(msn_mp, "Midpoint")

msg("Sample distribution (Final):")
print(table(pData(msn_f)$Mg))
msg("Sample distribution (Midpoint):")
print(table(pData(msn_mp)$Mg))

## DATA IMPUTATION AND FILTERING ----------------------------------
msg("Starting data imputation and filtering")

# Function to impute data for a given MSnSet

impute_msn_data <- function(msn, tag) {
  msg(sprintf("Starting mixed imputation for %s", tag))
  
  expr_data <- exprs(msn)
  expr_data_original <- expr_data  # Keep original for variance check
  groups <- pData(msn)$Mg
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
  
  n_per_group <- 3
  
  # Classify proteins
  completely_missing <- rowSums(missing_pattern == n_per_group) > 0
  partially_missing <- rowSums(missing_pattern > 0 & missing_pattern < n_per_group) > 0 &
    !completely_missing
  
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
  
  # MNAR imputation - using a global, down-shifted distribution
  if (sum(completely_missing) > 0) {
    msg("Performing MNAR imputation using a global Perseus-like method.")
    
    # 1. Calculate parameters from the ENTIRE dataset ONCE
    all_values <- as.vector(expr_data)
    global_mean <- mean(all_values, na.rm = TRUE)
    global_sd <- sd(all_values, na.rm = TRUE)
    
    # 2. Define the down-shifted distribution for imputation
    impute_mean <- global_mean - 1.8 * global_sd
    impute_sd <- global_sd * 0.3
    
    # 3. Impute all MNAR values by drawing from this single distribution
    complete_idx <- which(completely_missing)
    for (p in complete_idx) {
      for (g_idx in split(seq_len(ncol(expr_data)), groups)) {
        vals <- expr_data[p, g_idx]
        if (all(is.na(vals))) {
          set.seed(42 + p) # Seed for reproducibility
          expr_data[p, g_idx] <- rnorm(length(vals), mean = impute_mean, sd = impute_sd)
        }
      }
    }
    msg(sprintf("MNAR imputation completed for %d proteins", sum(completely_missing)))
  } # <-- This brace correctly closes the MNAR if-statement.
  
  # Variance check
  post_var <- apply(expr_data, 1, var, na.rm = TRUE)
  pre_var <- apply(expr_data_original, 1, var, na.rm = TRUE)
  
  # Only check proteins that had missing values
  had_missing <- rowSums(is.na(expr_data_original)) > 0
  if (sum(had_missing) > 0) {
    var_ratio <- post_var[had_missing] / pre_var[had_missing]
    var_preserved <- median(var_ratio, na.rm = TRUE)
    msg(sprintf("Median variance preserved after imputation: %.1f%%", var_preserved * 100))
    
    if (var_preserved < 0.8) {
      warn("Imputation reduced variance by >20% - consider alternative methods")
    }
  }
  
  # Update expression data
  exprs(msn) <- expr_data
  
  # Remove proteins still with any NA
  remaining_na <- rowSums(is.na(expr_data)) > 0
  if (any(remaining_na)) {
    msg(sprintf("Removing %d proteins with remaining NAs", sum(remaining_na)))
    msn <- msn[!remaining_na, ]
  }
  
  msg(sprintf("Imputation complete for %s: %d proteins retained", tag, nrow(msn)))
  return(msn)
  
}

# Apply imputation to both datasets
set.seed(42)
msn_f  <- impute_msn_data(msn_f,  "Final")
msn_mp <- impute_msn_data(msn_mp, "Midpoint")
gc()  # Force garbage collection
monitor_memory("after imputation")


# Check if phenoData exists before trying to drop levels
if (!is.null(pData(msn_f)) && "Mg" %in% colnames(pData(msn_f))) {
  pDat <- pData(msn_f)
  pDat$Mg <- droplevels(pDat$Mg)
  pData(msn_f) <- pDat
  
  # Verify
  msg(sprintf("Final Mg levels: %s", paste(levels(pData(msn_f)$Mg), collapse=", ")))
} else {
  stop("phenoData was lost during imputation - check the impute_msn_data function")
}

# Do the same for midpoint if needed
if (!is.null(pData(msn_mp)) && "Mg" %in% colnames(pData(msn_mp))) {
  pDat <- pData(msn_mp)
  pDat$Mg <- droplevels(pDat$Mg)
  pData(msn_mp) <- pDat
  msg(sprintf("Midpoint Mg levels: %s", paste(levels(pData(msn_mp)$Mg), collapse=", ")))
}

# Verify
levels(pData(msn_f)$Mg) 

## Limma Analysis Function ------------------------------
run_limma_analysis <- function(msn, tag, timestamp = format(Sys.time(), "%Y%m%d_%H%M")) {
  msg("Running enhanced limma analysis on", tag, "data")
  
  sample_lists <- lapply(levels(pData(msn)$Mg), function(g) {
    sampleNames(msn)[pData(msn)$Mg == g]
  })
  names(sample_lists) <- levels(pData(msn)$Mg)
  
  pData(msn)$Mg <- relevel(pData(msn)$Mg, ref = "Control")
  design <- model.matrix(~ Mg, data = pData(msn))
  
  fit <- lmFit(exprs(msn), design)
  fit2 <- eBayes(fit)
  
  all_results <- list()
  
  for (coef_name in colnames(design)[-1]) {
    comp_label <- gsub("^Mg", "", coef_name)
    comp_clean <- paste(comp_label, "vs Control")
    
    msg(sprintf("Processing comparison: %s", comp_clean))
    
    tt <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
    common_ids <- intersect(rownames(tt), rownames(exprs(msn)))
    if (length(common_ids) == 0) {
      warn("No overlapping protein IDs for", comp_label)
      next
    }
    
    mean_group <- if (comp_label %in% names(sample_lists) && length(sample_lists[[comp_label]]) > 0) {
      rowMeans(exprs(msn)[common_ids, sample_lists[[comp_label]], drop = FALSE], na.rm = TRUE)
    } else {
      rep(NA, length(common_ids))
    }
    mean_control <- rowMeans(exprs(msn)[common_ids, sample_lists[["Control"]], drop = FALSE], na.rm = TRUE)
    
    tt_enhanced <- tt[common_ids, , drop = FALSE] %>%
      rownames_to_column("ProteinID") %>%
      mutate(
        Direction = case_when(
          logFC > 0 ~ paste("Up in", comp_label),
          logFC < 0 ~ "Up in Control",
          TRUE ~ "No change"
        ),
        Regulation = case_when(
          adj.P.Val <= sig_adjP & logFC > 0 ~ paste("Significantly up in", comp_label),
          adj.P.Val <= sig_adjP & logFC < 0 ~ "Significantly up in Control",
          TRUE ~ "Not significant"
        ),
        FoldChange_linear = 2^logFC,
        abs_logFC = abs(logFC),
        !!paste0("Mean_", comp_label) := mean_group,
        Mean_Control = mean_control,
        Comparison = comp_clean,
        Timepoint = tag
      ) %>%
      dplyr::select(
        ProteinID, logFC, FoldChange_linear, Direction, AveExpr,
        starts_with("Mean_"), t, P.Value, adj.P.Val, B,
        Regulation, abs_logFC, Comparison, Timepoint
      ) %>%
      arrange(adj.P.Val, desc(abs_logFC))
    
    output_file <- file.path(tables_dir, paste0("Results_", comp_label, "_vs_Control_", tag, "_", timestamp, ".xlsx"))
    write_xlsx(tt_enhanced, output_file)
    msg(sprintf("Saved results to: %s", basename(output_file)))
    
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
               file.path(tables_dir, paste0("Summary_", comp_label, "_vs_Control_", tag, "_", timestamp, ".xlsx")))
    
    msg(strrep("=", 50))
    msg(sprintf("SUMMARY: %s (%s)", comp_clean, tag))
    msg(sprintf("Total proteins analyzed: %d", nrow(tt_enhanced)))
    msg(sprintf("Significantly up in %s: %d", comp_label, 
                sum(tt_enhanced$Regulation == paste("Significantly up in", comp_label))))
    msg(sprintf("Significantly up in Control: %d", 
                sum(tt_enhanced$Regulation == "Significantly up in Control")))
    msg(sprintf("Significance threshold: adj.P.Val <= %.3f", sig_adjP))
    msg(strrep("=", 50))
    all_results[[comp_label]] <- tt_enhanced
    
    
    # QC Plots
    pdf(file.path(plots_dir, paste0("QC_", comp_label, "_", tag, "_", timestamp, ".pdf")), width = 10, height = 5)
    par(mfrow = c(1,2))
    plotSA(fit2, main = paste("Mean-variance trend –", tag))
    hist(fit2$p.value[, coef_name], 50, xlim = c(0,1), main = paste("p-value distribution (", comp_label, " vs Control)"), xlab = "p-value")
    dev.off()
    
    # Volcano Plot
    volc_colors <- c(
      "Not significant" = "grey70",
      "Significantly up in Control" = "#377eb8"
    )
    volc_colors[paste("Significantly up in", comp_label)] <- "#e41a1c"
    x_right <- quantile(tt_enhanced$logFC, 0.98, na.rm = TRUE)
    x_left  <- quantile(tt_enhanced$logFC, 0.02, na.rm = TRUE)
    y_top   <- quantile(-log10(tt_enhanced$adj.P.Val), 0.98, na.rm = TRUE)
    
    p_volcano <- ggplot(tt_enhanced, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(colour = Regulation), size = 2, alpha = 0.7) +
      scale_colour_manual(values = volc_colors) +
      geom_hline(yintercept = -log10(sig_adjP), linetype = "dashed", alpha = 0.5) +
      annotate("text", x = x_right, y = y_top, 
               label = paste("More abundant \nin", comp_label), 
               hjust = 1, size = 3.5, color = "#e41a1c") +
      annotate("text", x = x_left, y = y_top, 
               label = "More abundant\nin Control", 
               hjust = 0, size = 3.5, color = "#377eb8") +
      labs(
        title = paste("Volcano Plot:", comp_clean),
        subtitle = paste("Timepoint:", tag, "| Positive logFC = Higher in", comp_label),
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"Adjusted p-value")
      ) +
      theme_classic(base_size = 12) + theme(legend.position = "bottom")
    
    ggsave(file.path(plots_dir, paste0("Volcano_", comp_label, "_vs_Control_", tag, "_", timestamp, ".pdf")),
           plot = p_volcano, width = 9, height = 7)
    
    # Heatmap for significant proteins
    sig_proteins <- tt_enhanced %>% 
      filter(grepl("Significantly up", Regulation)) %>% 
      pull(ProteinID)
    
    min_proteins <- cfg$min_sig_proteins_heatmap %||% 2
    max_proteins <- cfg$max_heatmap_proteins %||% 50
    
    if (length(sig_proteins) >= min_proteins) {
      n_to_show <- 100
      top_sig <- head(sig_proteins, n_to_show)
      expr_sig <- exprs(msn)[top_sig, , drop = FALSE]
      
      # Calculate control-relative z-scores
      control_cols <- sample_lists[["Control"]]
      control_mean <- rowMeans(expr_sig[, control_cols, drop = FALSE], na.rm = TRUE)
      expr_sig_lfc <- expr_sig - control_mean
      color_fun <- circlize::colorRamp2(
        breaks = c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3),
        colors = c("#053061", "#2166ac", "#4393c3", "#92c5de", 
                   "#f7f7f7", 
                   "#fddbc7", "#f4a582", "#d6604d", "#67001f")
      )
      
      # Define annotations
      row_anno_data <- tt_enhanced %>%
        filter(ProteinID %in% top_sig) %>%
        dplyr::select(ProteinID, Direction, logFC)
      
      col_ha <- HeatmapAnnotation(
        Group = pData(msn)$Mg,
        col = list(Group = mg_cols[levels(pData(msn)$Mg)])
      )
      
      # Calculate appropriate figure height
      fig_height <- max(8, min(20, 4 + n_to_show * 0.2))
      
      pdf(file.path(plots_dir, paste0("Heatmap_", comp_label, "_vs_Control_", tag, "_", timestamp, ".pdf")),
          width = 12, height = fig_height)
      
      print(ComplexHeatmap::Heatmap(
        expr_sig_lfc,
        name = "Log2FC\nvs Control",  # Clear label
        col = color_fun,
        top_annotation = col_ha,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        cluster_columns = FALSE,
        column_split = pData(msn)$Mg,
        row_title = sprintf("Top %d Significant Proteins", n_to_show),
        column_title = paste(comp_clean, "-", tag)
      ))
      dev.off()
    } else {
      warn(sprintf("Too few significant proteins (%d) for heatmap in %s", length(sig_proteins), comp_clean))
    }
  }
  if (tag == "Final" && "NoMg" %in% levels(pData(msn)$Mg) && "LowMg" %in% levels(pData(msn)$Mg)) {
    msg("Adding special comparison: NoMg vs LowMg (R4-6 vs R1-3)")
    
    # Create a subset with just NoMg and LowMg
    samples_subset <- sampleNames(msn)[pData(msn)$Mg %in% c("NoMg", "LowMg")]
    msn_subset <- msn[, samples_subset]
    
    # Re-factor to have just these two levels
    pData(msn_subset)$Mg <- factor(pData(msn_subset)$Mg, levels = c("LowMg", "NoMg"))
    
    # Design matrix for this comparison
    design_special <- model.matrix(~ Mg, data = pData(msn_subset))
    
    # Run limma
    fit_special <- lmFit(exprs(msn_subset), design_special)
    fit2_special <- eBayes(fit_special)
    
    # Extract results
    tt_special <- topTable(fit2_special, coef = "MgNoMg", number = Inf, adjust.method = "BH")
    
    # Process and save similar to other comparisons
    tt_enhanced_special <- tt_special %>%
      rownames_to_column("ProteinID") %>%
      mutate(
        Direction = case_when(
          logFC > 0 ~ "Up in NoMg (R4-6)",
          logFC < 0 ~ "Up in LowMg (R1-3)",
          TRUE ~ "No change"
        ),
        Regulation = case_when(
          adj.P.Val <= sig_adjP & logFC > 0 ~ paste("Significantly up in", comp_label),
          adj.P.Val <= sig_adjP & logFC < 0 ~ "Significantly up in Control",
          TRUE ~ "Not significant"
        ),
        FoldChange_linear = 2^logFC,
        abs_logFC = abs(logFC),
        Comparison = "NoMg vs LowMg (0 g/L vs 0 g/L different times)",
        Timepoint = tag
      ) %>%
      arrange(adj.P.Val, desc(abs_logFC))
    
    # Save this comparison
    write_xlsx(tt_enhanced_special, 
               file.path(tables_dir, paste0("Results_NoMg_vs_LowMg_Final_", timestamp, ".xlsx")))
    
    msg(sprintf("NoMg vs LowMg comparison:"))
    msg(sprintf("  Significantly up in NoMg (R4-6): %d", 
                sum(tt_enhanced_special$Regulation == "Significantly up in NoMg")))
    msg(sprintf("  Significantly up in LowMg (R1-3): %d", 
                sum(tt_enhanced_special$Regulation == "Significantly up in LowMg")))
    
    all_results[["NoMg_vs_LowMg"]] <- tt_enhanced_special
  }
  return(list(fit = fit2, results = all_results))
}

## Run Analyses --------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

msg("\n", strrep("=", 60))
msg("ANALYZING FINAL TIMEPOINT")
msg(strrep("=", 60))
res_F <- run_limma_analysis(msn_f, "Final", timestamp)
gc()
monitor_memory("after Final analysis")

msg("\n", strrep("=", 60))
msg("ANALYZING MIDPOINT TIMEPOINT")
msg(strrep("=", 60))
res_MP <- run_limma_analysis(msn_mp, "Midpoint", timestamp)

## Combine All Results -------------------------------------------
all_results_combined <- dplyr::bind_rows(
  lapply(names(res_F$results), function(x) res_F$results[[x]]),
  lapply(names(res_MP$results), function(x) res_MP$results[[x]])
)
write_xlsx(all_results_combined, file.path(tables_dir, paste0("ALL_RESULTS_COMBINED_", timestamp, ".xlsx")))

## PCA and UMAP (Separate by Timepoint) -------------------------
run_dim_reduction <- function(msn, tag, timestamp) {
  msg("Running PCA and UMAP on", tag)
  
  data_for_dim_red <- t(exprs(msn))
  sample_info <- pData(msn)
  
  # PCA
  pca_result <- prcomp(data_for_dim_red, scale. = TRUE, center = TRUE)
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  pca_data$Mg <- sample_info$Mg
  pca_data$Label <- paste0("R", sample_info$Reactor)
  
  # UMAP
  set.seed(42)
  umap_result <- umap(data_for_dim_red, n_neighbors = min(10, nrow(data_for_dim_red) - 1), min_dist = 0.1, n_components = 2)
  umap_data <- as.data.frame(umap_result)
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  umap_data$Mg <- sample_info$Mg
  umap_data$Label <- paste0("R", sample_info$Reactor)
  
  # Color scheme
  group_colors <- c(Control = "#0073C2FF", NormalMg = "#EFC000FF", LowMg = "#CD534CFF", NoMg = "#EFC000FF")
  
  # PCA Plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Mg), size = 6, alpha = 0.8) +
    geom_text_repel(aes(label = Label), size = 4, box.padding = 0.5, point.padding = 0.3, segment.size = 0, show.legend = FALSE, seed = 42) +
    theme_classic(base_size = 14) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1), legend.position = "bottom", aspect.ratio = 1) +
    labs(title = paste("PCA —", tag), x = sprintf("PC1 (%.1f%%)", var_explained[1]), y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
    scale_color_manual(values = group_colors) + coord_fixed(ratio = 1)
  
  # UMAP Plot
  umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Mg), size = 6, alpha = 0.8) +
    geom_text_repel(aes(label = Label), size = 4, box.padding = 0.5, point.padding = 0.3, segment.size = 0, show.legend = FALSE, seed = 42) +
    theme_classic(base_size = 14) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1), legend.position = "bottom", aspect.ratio = 1) +
    labs(title = paste("UMAP —", tag), x = "UMAP1", y = "UMAP2") +
    scale_color_manual(values = group_colors) + coord_fixed(ratio = 1)
  
  # Save
  ggsave(file.path(plots_dir, paste0("PCA_", tag, "_", timestamp, ".pdf")), plot = pca_plot, width = 8, height = 7)
  ggsave(file.path(plots_dir, paste0("UMAP_", tag, "_", timestamp, ".pdf")), plot = umap_plot, width = 8, height = 7)
  ggsave(file.path(plots_dir, paste0("PCA_UMAP_", tag, "_", timestamp, ".png")), plot = (pca_plot + umap_plot), width = 14, height = 7, dpi = 300, bg = "white")
}

run_dim_reduction(msn_f, "Final", timestamp)
run_dim_reduction(msn_mp, "Midpoint", timestamp)
rm(list = c("tab", "exprs_f", "exprs_mp"))  # Remove raw data no longer needed
gc()

## GO Enrichment Analysis ----------------------------------------
msg("\nPerforming GO enrichment analysis")
run_go_enrichment <- function(results_list, msn_current, tag_str, timestamp) {
  msg(sprintf("\nPerforming GO enrichment analysis (%s)", tag_str))
  
  go_raw <- readr::read_tsv(go_file, col_names = FALSE, na = c("", "NA"), show_col_types = FALSE) %>%
    rename(ProteinID = X1) %>%
    tidyr::pivot_longer(-ProteinID, names_to = "col", values_to = "GO_score", values_drop_na = TRUE) %>%
    dplyr::select(-col) %>%
    tidyr::separate(GO_score, into = c("GO_ID", "score"), sep = "\\|", convert = TRUE, fill = "right") %>%
    filter(!is.na(GO_ID), score >= go_min_score) %>%
    mutate(Acc = acc_from(ProteinID))
  
  if (nrow(go_raw) == 0) { warn("No valid GO annotations parsed – check format and score threshold"); return(invisible(NULL)) }
  
  bp_ids <- AnnotationDbi::select(GO.db, keys = unique(go_raw$GO_ID), columns = "ONTOLOGY", keytype = "GOID") %>%
    filter(ONTOLOGY == "BP") %>% pull(GOID)
  go_bp <- go_raw %>% filter(GO_ID %in% bp_ids)
  
  TERM2GENE <- distinct(go_bp, GO_ID, Acc); colnames(TERM2GENE) <- c("term","gene")
  TERM2GENE$term <- as.character(TERM2GENE$term); TERM2GENE$gene <- as.character(TERM2GENE$gene)
  TERM2NAME <- AnnotationDbi::select(GO.db, keys = unique(TERM2GENE$term), columns = "TERM", keytype = "GOID") %>% distinct()
  colnames(TERM2NAME) <- c("term","name"); TERM2NAME$term <- as.character(TERM2NAME$term); TERM2NAME$name <- as.character(TERM2NAME$name)
  
  # timepoint-wide universe (used only to sanity-check), not for ORA
  universe_tp <- intersect(acc_from(rownames(exprs(msn_current))), unique(TERM2GENE$gene))
  universe_tp <- unique(as.character(universe_tp))
  msg(sprintf("GO annotation summary (%s): %d proteins in timepoint universe", tag_str, length(universe_tp)))
  
  safe_enricher <- function(genes, universe_arg, direction_lab) {
    genes <- intersect(unique(as.character(genes)), universe_arg)
    if (length(genes) < 5) return(NULL)
    tryCatch(
      clusterProfiler::enricher(
        gene=genes, universe=universe_arg,
        TERM2GENE=TERM2GENE, TERM2NAME=TERM2NAME,
        pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
        minGSSize=5, maxGSSize=5000
        # Remove: by = "fgsea", eps = 1e-20, nPermSimple = 10000
      ),
      error=function(e){
        warn(sprintf("enricher() failed for %s (%s): %s", direction_lab, tag_str, e$message))
        NULL
      }
    )
  }
  
  # keep track of files for the optional comparison table
  ora_paths <- list()
  
  for (comp_name in names(results_list)) {
    msg(sprintf("\nRunning GO for %s vs Control (%s)", comp_name, tag_str))
    df <- results_list[[comp_name]]
    
    tested_ids <- acc_from(df$ProteinID)
    universe_contrast <- intersect(tested_ids, unique(TERM2GENE$gene))
    if (length(universe_contrast) < 10) { warn(sprintf("Tiny universe for %s (%s)", comp_name, tag_str)); next }
    
    hits_up <- df %>% filter(Regulation == paste("Significantly up in", comp_name)) %>%
      transmute(Acc = acc_from(ProteinID)) %>% pull(Acc) %>% intersect(universe_contrast)
    hits_dn <- df %>% filter(Regulation == "Significantly up in Control") %>%
      transmute(Acc = acc_from(ProteinID)) %>% pull(Acc) %>% intersect(universe_contrast)
    
    # ORA UP
    up_path <- NA_character_
    ora_up <- safe_enricher(hits_up, universe_contrast, sprintf("UP_in_%s", comp_name))
    if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
      up_path <- file.path(tables_dir, sprintf("ORA_BP_UP_in_%s_vs_Control_%s_%s.csv", comp_name, tag_str, timestamp))
      write.csv(ora_up@result, up_path, row.names = FALSE)
      p_up <- barplot(ora_up, showCategory = 20) + 
        ggtitle(sprintf("ORA: UP in %s (vs Control, %s)", comp_name, tag_str))
      ggsave(file.path(plots_dir, sprintf("ORA_barplot_UP_in_%s_%s_%s.pdf", comp_name, tag_str, timestamp)), plot = p_up, width = 10, height = 8)
    }
    
    # ORA DOWN
    down_path <- NA_character_
    ora_dn <- safe_enricher(hits_dn, universe_contrast, sprintf("DOWN_in_%s", comp_name))
    if (!is.null(ora_dn) && nrow(ora_dn@result) > 0) {
      down_path <- file.path(tables_dir, sprintf("ORA_BP_DOWN_in_%s_vs_Control_%s_%s.csv", comp_name, tag_str, timestamp))
      write.csv(ora_dn@result, down_path, row.names = FALSE)
      p_dn <- barplot(ora_dn, showCategory = 20) + 
        ggtitle(sprintf("ORA: DOWN in %s (vs Control, %s)", comp_name, tag_str))
      ggsave(file.path(plots_dir, sprintf("ORA_barplot_DOWN_in_%s_%s_%s.pdf", comp_name, tag_str, timestamp)), plot = p_dn, width = 10, height = 8)
    }
    
    ora_paths[[comp_name]] <- list(up = up_path, down = down_path)
    
    # GSEA (ranked by t)
    ranked_df <- df %>%
      mutate(Acc = acc_from(ProteinID)) %>%
      filter(Acc %in% universe_contrast, is.finite(t)) %>%
      group_by(Acc) %>% summarise(t = mean(t), .groups = "drop") %>%
      arrange(desc(t))
    if (nrow(ranked_df) >= 10) {
      geneList <- ranked_df$t; names(geneList) <- ranked_df$Acc
      gsea <- tryCatch(
        clusterProfiler::GSEA(geneList, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                              pvalueCutoff = 1, verbose = FALSE, minGSSize = 5, maxGSSize = 5000),
        error = function(e) { warn(sprintf("GSEA failed for %s (%s): %s", comp_name, tag_str, e$message)); NULL }
      )
      if (!is.null(gsea) && nrow(gsea@result) > 0) {
        write.csv(gsea@result, file.path(tables_dir, sprintf("GSEA_BP_%s_vs_Control_%s_%s.csv", comp_name, tag_str, timestamp)), row.names = FALSE)
        ggsave(file.path(plots_dir, sprintf("GSEA_dotplot_%s_vs_Control_%s_%s.pdf", comp_name, tag_str, timestamp)),
               plot = enrichplot::dotplot(gsea, showCategory = 20) + ggtitle(sprintf("GSEA: %s vs Control (%s)", comp_name, tag_str)),
               width = 10, height = 8)
        gsea_df <- gsea@result %>%
          arrange(desc(NES))
        
        # Select top N positive and negative terms for a clean plot
        top_n <- 10
        top_terms <- bind_rows(
          gsea_df %>% filter(NES > 0) %>% head(top_n),
          gsea_df %>% filter(NES < 0) %>% tail(top_n)
        )
        
        if (nrow(top_terms) > 0) {
          treatment_label <- if (comp_name == "LowMg") {
            "Enriched in 0 g/L"
          } else {
            sprintf("Enriched in %s", comp_name)
          }
          p_gsea_bar <- ggplot(top_terms, aes(x = NES, y = fct_reorder(Description, NES))) +
            geom_col(aes(fill = NES > 0)) +
            scale_fill_manual(
              values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
              name = "Enrichment Direction",
              # 2. Use the new labels
              labels = c("Enriched in 0.4 g/L", treatment_label)
            ) +
            labs(
              title = "GSEA: Top Enriched Pathways vs. NES",
              subtitle = sprintf("%s vs Control (%s)", comp_name, tag_str),
              x = "Normalized Enrichment Score (NES)",
              y = "GO Term"
            ) +
            theme_classic(base_size = 12) +
            theme(legend.position = "bottom")
          
          ggsave(
            file.path(plots_dir, sprintf("GSEA_barplot_NES_%s_vs_Control_%s_%s.pdf", comp_name, tag_str, timestamp)),
            plot = p_gsea_bar,
            width = 10,
            height = 8
          )
        }
      }
    }
  }
  
  # Optional: UP/DOWN side-by-side summary if both exist (now inside the function)
  for (comp_name in names(ora_paths)) {
    paths <- ora_paths[[comp_name]]
    if (is.na(paths$up) || is.na(paths$down)) next
    ora_up_df <- tryCatch(read.csv(paths$up),  error = function(e) NULL)
    ora_dn_df <- tryCatch(read.csv(paths$down), error = function(e) NULL)
    if (is.null(ora_up_df) || is.null(ora_dn_df)) next
    comparison_df <- data.frame(
      UP_in_Treatment   = utils::head(ora_up_df$Description, 10),
      UP_pvalue         = utils::head(ora_up_df$p.adjust, 10),
      DOWN_in_Treatment = utils::head(ora_dn_df$Description, 10),
      DOWN_pvalue       = utils::head(ora_dn_df$p.adjust, 10)
    )
    write.csv(comparison_df, file.path(tables_dir, sprintf("ORA_Comparison_%s_%s_%s.csv", comp_name, tag_str, timestamp)), row.names = FALSE)
  }
  
  invisible(NULL)
}

# Run GO/GSEA for Final and Midpoint using per-contrast universes
#run_go_enrichment(res_F$results, msn_f,  "Final",    timestamp)
#run_go_enrichment(res_MP$results, msn_mp, "Midpoint", timestamp)

run_enhanced_go_analysis <- function(results_list, msn_current, tag_str, timestamp) {
  msg("\n", strrep("=", 60))
  msg(sprintf("ENHANCED GO ANALYSIS WITH CLUSTERING + UMBRELLAS (%s)", tag_str))
  msg(strrep("=", 60))
  
  # Define negative filter list for eukaryotic terms
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
  
  # Define umbrella terms
  umbrella_mappings <- list(
    "transmembrane transport" = c(
      "transport", "transmembrane transport", "monoatomic ion transport",
      "monoatomic cation transport", "monoatomic anion transport",
      "monoatomic ion transmembrane transport", "monoatomic cation transmembrane transport",
      "inorganic ion transmembrane transport", "inorganic cation transmembrane transport",
      "inorganic anion transport", "metal ion transport", "ion transmembrane transport",
      "cation transmembrane transport", "anion transmembrane transport",
      "organic anion transport", "transition metal ion transport",
      "establishment of localization", "cellular localization",
      "carboxylic acid transmembrane transport", "organic acid transport",
      "protein transport", "protein transmembrane transport",
      "macromolecule localization", "localization", "magnesium ion transport",
      "divalent metal ion transport", "divalent inorganic cation transport"
    ),
    "electron transport and respiration" = c(
      "respiratory electron transport chain", "electron transport chain",
      "cellular respiration", "oxidation-reduction process",
      "electron transfer activity", "anaerobic respiration",
      "respiratory chain complex", "electron carrier activity",
      "heme-copper terminal oxidase activity"
    ),
    "energy generation" = c(
      "generation of precursor metabolites and energy",
      "nucleoside triphosphate metabolic process", "ribonucleotide metabolic process",
      "purine ribonucleotide metabolic process", "ribose phosphate metabolic process",
      "energy derivation by oxidation of organic compounds",
      "ATP metabolic process", "ATP biosynthetic process",
      "proton-transporting ATP synthase activity", "ATPase activity",
      "phosphorylation", "oxidative phosphorylation"
    ),
    "sulfur metabolism" = c(
      "sulfur compound metabolic process", "sulfate assimilation",
      "sulfur compound biosynthetic process", "cysteine metabolic process",
      "cysteine biosynthetic process", "sulfate metabolic process",
      "sulfate transport", "sulfur amino acid metabolic process",
      "sulfur amino acid biosynthetic process"
    ),
    "stress response" = c(
      "response to stress", "cellular response to stress",
      "response to oxidative stress", "response to osmotic stress",
      "response to nutrient levels", "response to starvation",
      "cellular response to nutrient levels", "response to metal ion",
      "cellular response to metal ion", "response to magnesium ion"
    ),
    "cell motility and attachment" = c(
      "cell projection organization", "bacterial-type flagellum-dependent cell motility",
      "pilus organization", "pilus assembly", "type IV pilus assembly",
      "cell adhesion", "biofilm formation", "surface attachment",
      "twitching motility", "cell motility", "locomotion"
    ),
    "protein synthesis and folding" = c(
      "translation", "peptide biosynthetic process", "amide biosynthetic process",
      "protein folding", "chaperone-mediated protein folding",
      "response to unfolded protein", "protein stabilization"
    )
  )
  
  # Function to apply umbrella term grouping
  apply_umbrella_terms <- function(go_df, method = "GSEA") {
    if (nrow(go_df) == 0) return(go_df)
    
    go_df$original_term <- go_df$Description
    go_df$was_grouped <- FALSE
    
    for (umbrella in names(umbrella_mappings)) {
      children <- umbrella_mappings[[umbrella]]
      for (child in children) {
        matches <- grepl(paste0("^", child, "$"), go_df$Description, ignore.case = TRUE)
        if (any(matches)) {
          go_df$Description[matches] <- umbrella
          go_df$was_grouped[matches] <- TRUE
        }
      }
    }
    
    if (method == "GSEA") {
      go_df_merged <- go_df %>%
        group_by(Description) %>%
        summarise(
          ID = dplyr::first(ID), # Explicitly use dplyr's first()
          setSize = max(setSize),
          enrichmentScore = dplyr::first(enrichmentScore),
          NES = dplyr::first(NES[which.max(abs(NES))]),
          pvalue = min(pvalue),
          p.adjust = min(p.adjust),
          qvalue = min(qvalue),
          rank = dplyr::first(rank),
          leading_edge = dplyr::first(leading_edge),
          core_enrichment = dplyr::first(core_enrichment),
          original_term = paste(unique(original_term), collapse = " | "),
          was_grouped = any(was_grouped),
          n_merged = n(),
          .groups = 'drop'
        )
    } else { # ORA method
      go_df_merged <- go_df %>%
        group_by(Description) %>%
        dplyr::summarise(
          ID = dplyr::first(ID), # Explicitly use dplyr's first()
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
  
  # Process GO results with filtering and grouping
  process_go_results <- function(go_result, method = "GSEA", nes_threshold = 1.5) {
    if (is.null(go_result) || nrow(go_result@result) == 0) {
      return(NULL)
    }
    
    go_df <- go_result@result
    
    # Filter by significance
    if (method == "GSEA") {
      go_filtered <- go_df %>%
        filter(abs(NES) >= nes_threshold, p.adjust < 0.05)
      msg(sprintf("  Filtered by |NES|≥%.1f & p.adj<0.05: %d -> %d terms", 
                  nes_threshold, nrow(go_df), nrow(go_filtered)))
    } else {
      go_filtered <- go_df %>%
        filter(p.adjust < 0.05)
      msg(sprintf("  Filtered by p.adjust<0.05: %d -> %d terms", 
                  nrow(go_df), nrow(go_filtered)))
    }
    
    # Remove eukaryotic terms
    go_filtered <- go_filtered %>%
      filter(!is_eukaryotic_term(Description))
    msg(sprintf("  After removing eukaryotic terms: %d terms", nrow(go_filtered)))
    
    # Apply umbrella grouping
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
  
  if (nrow(go_raw) == 0) { 
    warn("No valid GO annotations for enhanced analysis")
    return(invisible(NULL))
  }
  
  # Get biological process terms
  bp_ids <- AnnotationDbi::select(GO.db, keys = unique(go_raw$GO_ID), 
                                  columns = "ONTOLOGY", keytype = "GOID") %>%
    filter(ONTOLOGY == "BP") %>% 
    pull(GOID)
  go_bp <- go_raw %>% filter(GO_ID %in% bp_ids)
  
  # Prepare TERM2GENE and TERM2NAME
  TERM2GENE <- distinct(go_bp, GO_ID, Acc)
  colnames(TERM2GENE) <- c("term", "gene")
  TERM2GENE$term <- as.character(TERM2GENE$term)
  TERM2GENE$gene <- as.character(TERM2GENE$gene)
  
  TERM2NAME <- AnnotationDbi::select(GO.db, keys = unique(TERM2GENE$term), 
                                     columns = "TERM", keytype = "GOID") %>% 
    distinct()
  colnames(TERM2NAME) <- c("term", "name")
  
  # Process each comparison
  for (comp_name in names(results_list)) {
    msg(sprintf("\nProcessing enhanced GO for %s (%s)", comp_name, tag_str))
    df <- results_list[[comp_name]]
    
    # Prepare universe
    tested_ids <- acc_from(df$ProteinID)
    universe_contrast <- intersect(tested_ids, unique(TERM2GENE$gene))
    
    if (length(universe_contrast) < 10) {
      warn(sprintf("Universe too small for %s: %d proteins", comp_name, length(universe_contrast)))
      next
    }
    
    msg(sprintf("Contrast-specific universe: %d proteins", length(universe_contrast)))
    
    # ========== GSEA SECTION ==========
    msg("--- Running enhanced GSEA ---")
    
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
        clusterProfiler::GSEA(geneList,  # ADD THIS
                              TERM2GENE = TERM2GENE,  # ADD THIS
                              TERM2NAME = TERM2NAME,  # ADD THIS
                              pvalueCutoff = 1,  # ADD THIS
                              verbose = FALSE,  # ADD THIS
                              minGSSize = 5,  # ADD THIS
                              maxGSSize = 5000),  # ADD THIS
        error = function(e) { 
          warn(sprintf("GSEA failed: %s", e$message))  # FIX THIS
          NULL 
        }
      )
      
      if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
        gsea_enhanced <- process_go_results(gsea_result, method = "GSEA", nes_threshold = 1.5)
        
        if (!is.null(gsea_enhanced) && nrow(gsea_enhanced) > 0) {
          # Save GSEA results
          output_name <- if(comp_name == "NoMg_vs_LowMg") {
            sprintf("GSEA_ENHANCED_NoMg_vs_LowMg_%s_%s.csv", tag_str, timestamp)
          } else {
            sprintf("GSEA_ENHANCED_%s_vs_Control_%s_%s.csv", comp_name, tag_str, timestamp)
          }
          write.csv(gsea_enhanced, file.path(tables_dir, output_name), row.names = FALSE)
          msg(sprintf("  Saved enhanced GSEA results: %d terms", nrow(gsea_enhanced)))
          
          # ========== ADD THIS PLOTTING CODE ==========
          # Create GSEA barplot
          gsea_plot_data <- gsea_enhanced %>%
            arrange(desc(NES))
          
          # Select top N positive and negative terms
          top_n <- 10
          top_gsea_terms <- bind_rows(
            gsea_plot_data %>% filter(NES > 0) %>% head(top_n),
            gsea_plot_data %>% filter(NES < 0) %>% tail(top_n)
          )
          
          if (nrow(top_gsea_terms) > 0) {
            treatment_label <- if (comp_name == "LowMg") {
              "Enriched in 0 g/L"
            } else if (comp_name == "NoMg_vs_LowMg") {
              "Enriched in NoMg"
            } else {
              sprintf("Enriched in %s", comp_name)
            }
            
            control_label <- if (comp_name == "NoMg_vs_LowMg") {
              "Enriched in LowMg"
            } else {
              "Enriched in 0.4 g/L"
            }
            
            p_gsea_enhanced <- ggplot(top_gsea_terms, aes(x = NES, y = fct_reorder(Description, NES))) +
              geom_col(aes(fill = NES > 0)) +
              scale_fill_manual(
                values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
                name = "Enrichment Direction",
                labels = c(control_label, treatment_label)
              ) +
              labs(
                title = "GSEA Enhanced: Top Enriched Pathways",
                subtitle = sprintf("%s (%s) - Filtered & Grouped", comp_name, tag_str),
                x = "Normalized Enrichment Score (NES)",
                y = "GO Term (Umbrella Categories)"
              ) +
              theme_classic(base_size = 12) +
              theme(
                legend.position = "bottom",
                axis.text.y = element_text(size = 10)
              )
            
            plot_name <- if(comp_name == "NoMg_vs_LowMg") {
              sprintf("GSEA_ENHANCED_barplot_NoMg_vs_LowMg_%s_%s.pdf", tag_str, timestamp)
            } else {
              sprintf("GSEA_ENHANCED_barplot_%s_vs_Control_%s_%s.pdf", comp_name, tag_str, timestamp)
            }
            
            ggsave(
              file.path(plots_dir, plot_name),
              plot = p_gsea_enhanced,
              width = 10,
              height = 8
            )
            
            msg(sprintf("  Saved enhanced GSEA barplot: %s", plot_name))
          }
          # ========== END OF PLOTTING CODE ==========
        }
      }
    }
    
    msg("--- Running enhanced ORA ---")
    
    # Get significantly up/down proteins FOR THIS COMPARISON
    hits_up <- df %>% 
      filter(grepl("Significantly up in", Regulation) & 
               grepl(comp_name, Regulation)) %>%
      transmute(Acc = acc_from(ProteinID)) %>% 
      pull(Acc) %>% 
      intersect(universe_contrast)
    
    hits_dn <- df %>% 
      filter(Regulation == "Significantly up in Control") %>%
      transmute(Acc = acc_from(ProteinID)) %>% 
      pull(Acc) %>% 
      intersect(universe_contrast)
    
    # ORA for UP-regulated proteins
    if (length(hits_up) >= 5) {
      ora_up <- tryCatch(
        clusterProfiler::enricher(
          gene = hits_up, 
          universe = universe_contrast,
          TERM2GENE = TERM2GENE, 
          TERM2NAME = TERM2NAME,
          pAdjustMethod = "BH", 
          pvalueCutoff = 1, 
          qvalueCutoff = 1,
          minGSSize = 5, 
          maxGSSize = 5000
        ),
        error = function(e) { 
          warn(sprintf("ORA UP failed: %s", e$message))  # FIX THIS
          NULL 
        }
      )
      
      if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
        ora_up_enhanced <- process_go_results(ora_up, method = "ORA")
        
        if (!is.null(ora_up_enhanced) && nrow(ora_up_enhanced) > 0) {
          # Save and plot ORA UP results
          output_name <- if(comp_name == "NoMg_vs_LowMg") {
            sprintf("ORA_ENHANCED_UP_NoMg_vs_LowMg_%s_%s.csv", tag_str, timestamp)
          } else {
            sprintf("ORA_ENHANCED_UP_%s_vs_Control_%s_%s.csv", comp_name, tag_str, timestamp)
          }
          write.csv(ora_up_enhanced, file.path(tables_dir, output_name), row.names = FALSE)
          msg(sprintf("  Saved enhanced ORA UP results: %d terms", nrow(ora_up_enhanced)))
          
          # Create barplot [your plotting code]
          top_terms <- ora_up_enhanced %>% arrange(p.adjust) %>% head(20)
          if (nrow(top_terms) > 0) {
            p_ora_up <- ggplot(top_terms, aes(x = Count, y = fct_reorder(Description, Count))) +
              geom_col(fill = "#e41a1c") +
              labs(
                title = sprintf("ORA: Enriched in %s", comp_name),
                subtitle = sprintf("%s vs Control (%s)", comp_name, tag_str),
                x = "Gene Count",
                y = "GO Term"
              ) +
              theme_classic(base_size = 12) +
              theme(axis.text.y = element_text(size = 10))
            
            ggsave(
              file.path(plots_dir, sprintf("ORA_ENHANCED_barplot_UP_%s_%s_%s.pdf", 
                                           comp_name, tag_str, timestamp)),
              plot = p_ora_up,
              width = 10,
              height = 8
            )
          }
        }
      }
    } else {
      msg(sprintf("  Too few UP proteins for ORA: %d", length(hits_up)))
    }
    
    # ORA for DOWN-regulated proteins
    if (length(hits_dn) >= 5) {
      ora_dn <- tryCatch(
        clusterProfiler::enricher(
          gene = hits_dn,  # ADD THIS
          universe = universe_contrast,  # ADD THIS
          TERM2GENE = TERM2GENE,  # ADD THIS
          TERM2NAME = TERM2NAME,  # ADD THIS
          pAdjustMethod = "BH",  # ADD THIS
          pvalueCutoff = 1,  # ADD THIS
          qvalueCutoff = 1,  # ADD THIS
          minGSSize = 5,  # ADD THIS
          maxGSSize = 5000  # ADD THIS
        ),
        error = function(e) { 
          warn(sprintf("ORA DOWN failed: %s", e$message))  # FIX THIS
          NULL 
        }
      )
      
      if (!is.null(ora_dn) && nrow(ora_dn@result) > 0) {
        ora_dn_enhanced <- process_go_results(ora_dn, method = "ORA")
        
        if (!is.null(ora_dn_enhanced) && nrow(ora_dn_enhanced) > 0) {
          # Save results
          output_name <- if(comp_name == "NoMg_vs_LowMg") {
            sprintf("ORA_ENHANCED_DOWN_NoMg_vs_LowMg_%s_%s.csv", tag_str, timestamp)
          } else {
            sprintf("ORA_ENHANCED_DOWN_%s_vs_Control_%s_%s.csv", comp_name, tag_str, timestamp)
          }
          
          write.csv(ora_dn_enhanced, file.path(tables_dir, output_name), row.names = FALSE)
          msg(sprintf("  Saved enhanced ORA DOWN results: %d terms", nrow(ora_dn_enhanced)))
          
          # Create barplot
          top_terms <- ora_dn_enhanced %>%
            arrange(p.adjust) %>%
            head(20)
          
          if (nrow(top_terms) > 0) {
            p_ora_dn <- ggplot(top_terms, aes(x = Count, y = fct_reorder(Description, Count))) +
              geom_col(fill = "#377eb8") +
              labs(
                title = "ORA: Enriched in Control",
                subtitle = sprintf("%s vs Control (%s)", comp_name, tag_str),
                x = "Gene Count",
                y = "GO Term"
              ) +
              theme_classic(base_size = 12) +
              theme(axis.text.y = element_text(size = 10))
            
            ggsave(
              file.path(plots_dir, sprintf("ORA_ENHANCED_barplot_DOWN_%s_%s_%s.pdf", 
                                           comp_name, tag_str, timestamp)),
              plot = p_ora_dn,
              width = 10,
              height = 8
            )
          }
        }
      } 
    } else {
      msg(sprintf("  Too few DOWN proteins for ORA: %d", length(hits_dn)))
    }
    
  }
  
  invisible(NULL)
}

# Call this after your regular GO enrichment
# For Final timepoint
run_enhanced_go_analysis(res_F$results, msn_f, "Final", timestamp)
run_enhanced_go_analysis(res_MP$results, msn_mp, "Midpoint", timestamp)

# For timepoint comparisons
if (length(timepoint_results) > 0) {
  # You'll need to adapt the timepoint results format to match
  # the expected structure for GO analysis
  run_enhanced_go_analysis(timepoint_results, msn_f, "Timepoint_Comparisons", timestamp)
}

## Cross-Timepoint Analysis ---------------------------------------
run_timepoint_comparison <- function(msn_f, msn_mp, timestamp) {
  msg("\n", strrep("=", 60))
  msg("ANALYZING ACROSS TIMEPOINTS")
  msg(strrep("=", 60))
  
  # Get common proteins between timepoints
  common_proteins <- intersect(rownames(exprs(msn_f)), rownames(exprs(msn_mp)))
  msg(sprintf("Common proteins between timepoints: %d", length(common_proteins)))
  
  all_timepoint_results <- list()
  
  # 1. Standard comparisons for groups that are the same across timepoints
  # (Control and LowMg - R1-3 and R7-9)
  stable_groups <- c("Control", "LowMg")
  
  for (group in stable_groups) {
    # Pick samples for this group at each timepoint
    samples_f  <- sampleNames(msn_f)[pData(msn_f)$Mg == group]
    samples_mp <- sampleNames(msn_mp)[pData(msn_mp)$Mg == group]
    if (length(samples_f) == 0 || length(samples_mp) == 0) {
      warn(sprintf("Skipping %s: not present in both timepoints", group))
      next
    }
    
    # Combine expression data
    expr_combined <- cbind(
      exprs(msn_f)[common_proteins, samples_f, drop = FALSE],
      exprs(msn_mp)[common_proteins, samples_mp, drop = FALSE]
    )
    
    # Build metadata for pairing by Reactor
    meta_f  <- pData(msn_f)[samples_f, c("Reactor"), drop = FALSE]
    meta_mp <- pData(msn_mp)[samples_mp, c("Reactor"), drop = FALSE]
    reactor_block <- factor(c(as.character(meta_f$Reactor), as.character(meta_mp$Reactor)))
    
    time <- factor(c(rep("Final", length(samples_f)), rep("Midpoint", length(samples_mp))),
                   levels = c("Midpoint", "Final"))
    design <- model.matrix(~ time)
    
    # Paired limma with duplicateCorrelation
    corfit <- duplicateCorrelation(expr_combined, design, block = reactor_block)
    fit    <- lmFit(expr_combined, design, block = reactor_block, correlation = corfit$consensus)
    fit2   <- eBayes(fit)
    
    # Extract results
    tt <- topTable(fit2, coef = "timeFinal", number = Inf, adjust.method = "BH")
    
    # Enhance results
    tt_enhanced <- tt %>%
      rownames_to_column("ProteinID") %>%
      mutate(
        Direction = case_when(
          logFC > 0 ~ "Higher at Final",
          logFC < 0 ~ "Higher at Midpoint",
          TRUE ~ "No change"
        ),
        Regulation = case_when(
          adj.P.Val <= sig_adjP & logFC > 0 ~ "Significantly increased after Mg removal",
          adj.P.Val <= sig_adjP & logFC < 0 ~ "Significantly decreased after Mg removal",
          TRUE ~ "Not significant"
        ),
        FoldChange_linear = 2^logFC,
        abs_logFC = abs(logFC),
        Mean_Final = rowMeans(expr_combined[, samples_f, drop = FALSE], na.rm = TRUE),
        Mean_Midpoint = rowMeans(expr_combined[, samples_mp, drop = FALSE], na.rm = TRUE),
        Group = group,
        Comparison = paste(group, "Final vs Midpoint")
      ) %>%
      arrange(adj.P.Val, desc(abs_logFC))
    
    # Save results
    output_file <- file.path(tables_dir, 
                             paste0("Timepoint_Comparison_", group, "_Final_vs_Midpoint_", 
                                    timestamp, ".xlsx"))
    write_xlsx(tt_enhanced, output_file)
    
    msg(sprintf("Results for %s:", group))
    msg(sprintf("  Total proteins: %d", nrow(tt_enhanced)))
    msg(sprintf("  Significantly higher at Final: %d", 
                sum(tt_enhanced$Regulation == "Significantly higher at Final")))
    msg(sprintf("  Significantly higher at Midpoint: %d", 
                sum(tt_enhanced$Regulation == "Significantly higher at Midpoint")))
    
    all_timepoint_results[[group]] <- tt_enhanced
    
    # Create volcano plot
    p_volcano <- ggplot(tt_enhanced, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(colour = Regulation), size = 2, alpha = 0.7) +
      scale_colour_manual(values = c(
        "Not significant" = "grey70",
        "Significantly higher at Final" = "#e41a1c",
        "Significantly higher at Midpoint" = "#377eb8"
      )) +
      #geom_vline(xintercept = c(-sig_logFC, sig_logFC), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(sig_adjP), linetype = "dashed", alpha = 0.5) +
      labs(
        title = paste("Timepoint Comparison:", group),
        subtitle = "Final vs Midpoint | Positive logFC = Higher at Final",
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"Adjusted p-value")
      ) +
      theme_classic(base_size = 12) + 
      theme(legend.position = "bottom")
    
    ggsave(file.path(plots_dir, 
                     paste0("Volcano_Timepoint_", group, "_Final_vs_Midpoint_", 
                            timestamp, ".pdf")),
           plot = p_volcano, width = 9, height = 7)
  }
  
  # 2. SPECIAL COMPARISON: R4-6 (NoMg at Final vs NormalMg at Midpoint)
  # This shows the effect of Mg removal after midpoint
  msg("\n", strrep("-", 50))
  msg("SPECIAL COMPARISON: R4-6 Reactors")
  msg("Comparing NoMg (Final) vs NormalMg (Midpoint)")
  msg("This shows the effect of Mg removal after midpoint")
  msg(strrep("-", 50))
  
  # Get R4-6 samples
  samples_f_r456 <- sampleNames(msn_f)[pData(msn_f)$Mg == "NoMg"]
  samples_mp_r456 <- sampleNames(msn_mp)[pData(msn_mp)$Mg == "NormalMg"]
  reactors_f <- gsub("^R([456])F$", "\\1", samples_f_r456)
  reactors_mp <- gsub("^R([456])MP$", "\\1", samples_mp_r456)
  common_reactors <- intersect(reactors_f, reactors_mp)
  if (length(common_reactors) >= 2) {  # Need at least 2 pairs
    msg(sprintf("Found %d paired reactors: R%s", 
                length(common_reactors), 
                paste(common_reactors, collapse = ", R")))
    
    # Order samples to ensure pairing
    samples_f_paired <- paste0("R", sort(common_reactors), "F")
    samples_mp_paired <- paste0("R", sort(common_reactors), "MP")
    
    # Extract paired expression data
    expr_f_paired <- exprs(msn_f)[common_proteins, samples_f_paired, drop = FALSE]
    expr_mp_paired <- exprs(msn_mp)[common_proteins, samples_mp_paired, drop = FALSE]
    
    # Method 1: Paired limma using duplicateCorrelation
    expr_combined <- cbind(expr_f_paired, expr_mp_paired)
    
    # Create design with blocking
    condition <- factor(c(rep("NoMg_Final", length(samples_f_paired)),
                          rep("NormalMg_Midpoint", length(samples_mp_paired))),
                        levels = c("NormalMg_Midpoint", "NoMg_Final"))
    
    reactor_block <- factor(c(common_reactors, common_reactors))
    
    design <- model.matrix(~ condition)
    
    # Calculate correlation within reactors
    corfit <- duplicateCorrelation(expr_combined, design, block = reactor_block)
    msg(sprintf("Consensus correlation within reactors: %.3f", corfit$consensus))
    
    # Fit model with blocking
    fit <- lmFit(expr_combined, design, block = reactor_block, 
                 correlation = corfit$consensus)
    fit2 <- eBayes(fit)
    
    # Extract results
    tt_r456 <- topTable(fit2, coef = "conditionNoMg_Final", 
                        number = Inf, adjust.method = "BH")
    
    # Method 2 (Alternative): Direct paired differences
    # Calculate within-reactor differences
    expr_diff <- expr_f_paired - expr_mp_paired
    colnames(expr_diff) <- paste0("R", common_reactors, "_diff")
    
    # Test if mean difference != 0
    design_diff <- matrix(1, ncol = 1, nrow = ncol(expr_diff))
    fit_diff <- lmFit(expr_diff, design_diff)
    fit2_diff <- eBayes(fit_diff)
    tt_diff <- topTable(fit2_diff, coef = 1, number = Inf, adjust.method = "BH")
    
    # Use the blocked analysis (Method 1) for main results
    tt_r456_enhanced <- tt_r456 %>%
      rownames_to_column("ProteinID") %>%
      mutate(
        Direction = case_when(
          logFC > 0 ~ "Higher after Mg removal (Final)",
          logFC < 0 ~ "Higher with Mg present (Midpoint)",
          TRUE ~ "No change"
        ),
        Regulation = case_when(
          adj.P.Val <= sig_adjP & logFC > 0 ~ "Significantly increased after Mg removal",
          adj.P.Val <= sig_adjP & logFC < 0 ~ "Significantly decreased after Mg removal",
          TRUE ~ "Not significant"
        ),
        FoldChange_linear = 2^logFC,
        abs_logFC = abs(logFC),
        Mean_NoMg_Final = rowMeans(expr_f_paired, na.rm = TRUE),
        Mean_NormalMg_Midpoint = rowMeans(expr_mp_paired, na.rm = TRUE),
        Mean_Difference = rowMeans(expr_diff, na.rm = TRUE),
        Comparison = "R4-6: NoMg(Final) vs NormalMg(Midpoint) [PAIRED]",
        Description = "Paired analysis - Effect of Mg removal after midpoint"
      ) %>%
      arrange(adj.P.Val, desc(abs_logFC))
    
    # Add reactor-specific changes for transparency
    for (r in common_reactors) {
      col_name <- paste0("Reactor_R", r, "_diff")
      r_diff <- expr_f_paired[, paste0("R", r, "F")] - 
        expr_mp_paired[, paste0("R", r, "MP")]
      tt_r456_enhanced[[col_name]] <- r_diff[tt_r456_enhanced$ProteinID]
    }
    
    # Save results
    output_file <- file.path(tables_dir, 
                             paste0("R456_MgRemoval_PAIRED_Final_vs_Midpoint_", 
                                    timestamp, ".xlsx"))
    write_xlsx(tt_r456_enhanced, output_file)
    
    msg(sprintf("PAIRED Results for R4-6 (Mg removal effect):"))
    msg(sprintf("  Total proteins: %d", nrow(tt_r456_enhanced)))
    msg(sprintf("  Significantly increased after Mg removal: %d", 
                sum(tt_r456_enhanced$Regulation == "Significantly increased after Mg removal")))
    msg(sprintf("  Significantly decreased after Mg removal: %d", 
                sum(tt_r456_enhanced$Regulation == "Significantly decreased after Mg removal")))
    
    # Create comparison with unpaired analysis (optional)
    msg("\nComparison note: Paired analysis typically yields more significant results")
    msg("due to controlling for reactor-specific variation.")
    
    all_timepoint_results[["R456_MgRemoval_Paired"]] <- tt_r456_enhanced
    
    x_right <- quantile(tt_r456_enhanced$logFC, 0.98, na.rm = TRUE)
    x_left  <- quantile(tt_r456_enhanced$logFC, 0.02, na.rm = TRUE)
    y_top   <- quantile(-log10(tt_r456_enhanced$adj.P.Val), 0.98, na.rm = TRUE)
    
    # Create volcano plot with special labeling
    p_volcano_r456 <- ggplot(tt_r456_enhanced, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(colour = Regulation), size = 2, alpha = 0.7) +
      scale_colour_manual(values = c(
        "Not significant" = "grey70",
        "Significant but low fold change" = "#feb24c",
        "Significantly increased after Mg removal" = "#e41a1c",
        "Significantly decreased after Mg removal" = "#377eb8"
      )) +
      #geom_vline(xintercept = c(-sig_logFC, sig_logFC), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(sig_adjP), linetype = "dashed", alpha = 0.5) +
      annotate("text", x = x_right, y = y_top, label = "Increased after\nMg removal",
               hjust = 1, size = 3.5, color = "#e41a1c") +
      annotate("text", x = x_left,  y = y_top, label = "Decreased after\nMg removal",
               hjust = 0, size = 3.5, color = "#377eb8") +
      labs(
        title = "R4-6 Reactors: Effect of Mg Removal",
        subtitle = "NoMg (Final) vs NormalMg (Midpoint) | Same reactors, different conditions",
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"Adjusted p-value")
      ) +
      theme_classic(base_size = 12) + 
      theme(legend.position = "bottom",
            legend.box.margin = margin(t = -10, r = -2, b = 0, l = 0))
    
    ggsave(file.path(plots_dir, 
                     paste0("Volcano_R456_MgRemoval_Effect_", timestamp, ".pdf")),
           plot = p_volcano_r456, width = 10, height = 8)
  } else {
    warn("Could not find R4-6 samples for special comparison")
  }
  
  return(all_timepoint_results)
}

timepoint_results <- run_timepoint_comparison(msn_f, msn_mp, timestamp)
if (length(timepoint_results) > 0) {
  # Run GO for the R4-6 comparison specifically
  if ("R456_MgRemoval_Paired" %in% names(timepoint_results)) {
    # Create a pseudo results_list for GO enrichment
    r456_results <- list(R456_MgRemoval = timepoint_results[["R456_MgRemoval_Paired"]])
    
    # Need to modify the regulation column to match expected format
    r456_results$R456_MgRemoval <- r456_results$R456_MgRemoval %>%
      mutate(
        Regulation = case_when(
          Regulation == "Significantly increased after Mg removal" ~ "Significantly up in NoMg_Final",
          Regulation == "Significantly decreased after Mg removal" ~ "Significantly up in NormalMg_Midpoint",
          TRUE ~ "Not significant"
        )
      )
    
    # Run GO enrichment - you'll need to determine which msn to use for universe
    run_go_enrichment(r456_results, msn_f, "R456_Timepoint", timestamp)
  }
}

# Combine all timepoint comparison results
if (length(timepoint_results) > 0) {
  all_timepoint_combined <- dplyr::bind_rows(
    lapply(names(timepoint_results), function(x) timepoint_results[[x]])
  )
  write_xlsx(all_timepoint_combined, 
             file.path(tables_dir, paste0("ALL_TIMEPOINT_COMPARISONS_", timestamp, ".xlsx")))
  
  msg(sprintf("\nSaved %d timepoint comparison results", length(timepoint_results)))
}

## Session Info --------------------------------------------------
sink(file.path(out_dir, paste0("session_info_", timestamp, ".txt")))
cat("Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

msg("")
msg(strrep("=", 60))
msg("ANALYSIS COMPLETE!")
msg(sprintf("All results saved to: %s", out_dir))
msg(strrep("=", 60))