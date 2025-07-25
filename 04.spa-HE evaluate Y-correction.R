### 1. Load and preprocess ST data and SpaHE prediction results
### 2. Aggregate data into spatial grids
### 3. Calculate grid-level Pearson correlation coefficients (PCC)
### 4. Generate visualization results to display spatial correlations

evaluate_spahe_accuracy <- function(st_rds_path, 
                                    spahe_pred_path, 
                                    output_prefix = "spahe_evaluation",
                                    grid_size = 100,  
                                    min_counts = 1,  
                                    cell_types = NULL,
                                    min_cells_per_grid = 1,  
                                    min_variance_threshold = 0.001,
                                    spahe_col_mapping = list(
                                      cell_id = "cell_id",
                                      x = "x",
                                      y = "y",
                                      cell_type = "predicted_celltype"
                                    ),
                                    scale_spahe_coords = TRUE,  
                                    st_x_range = NULL,  
                                    st_y_range = NULL  
) {
  # Load required packages
  required_packages <- c("ggplot2", "dplyr", "tidyr", "purrr", "patchwork", "Seurat", "viridis", "scales")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing_packages, update = FALSE)
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
  
  # Cell type color scheme
  celltype_colors <- c(
    "Tumor cells" = "#FF6F61", "Normal epithelial cells" = "#FFD800", 
    "CD8 T cells" = "#7FDBFF", "CD4 T cells" = "#1E90FF", 
    "Macrophages" = "#9370DB", "Neutrophils" = "#00CED1", 
    "NK cells" = "#DA70D6", "Tregs" = "#CD853F", 
    "Fibroblasts" = "#4CAF50", "Endothelial cells" = "#9ACD32", 
    "B cells" = "#F8C3CD", "Dendritic cells" = "#20B2AA"
  )
  
  # Universal theme (fixed margin parameters)
  theme_uniform <- theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      aspect.ratio = 1,  # Force 1:1 aspect ratio
      axis.line = element_line(colour = "black", linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
      axis.title = element_text(face = "bold", size = 10),  
      axis.text = element_text(size = 9),  
      plot.margin = unit(c(5, 5, 5, 5), "pt") 
    )
  
  # Custom theme for first row plots (enforce consistent aspect ratio)
  theme_first_row <- theme_uniform +
    theme(
      aspect.ratio = 1,
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    )
  
  # --------------------------
  # 1. Data loading and basic validation
  # --------------------------
  message("=== Loading data ===")
  if (!file.exists(st_rds_path)) stop("ST file does not exist: ", st_rds_path)
  st_data <- readRDS(st_rds_path)
  if (!"auc_celltype" %in% colnames(st_data@meta.data)) {
    stop("'auc_celltype' column not found in ST data meta.data")
  }
  
  if (!file.exists(spahe_pred_path)) stop("SpaHE file does not exist: ", spahe_pred_path)
  load(spahe_pred_path)  # Assuming result_df is loaded
  
  # Validate SpaHE column mapping
  message("Checking SpaHE data column mapping...")
  available_cols <- colnames(result_df)
  for (col_key in names(spahe_col_mapping)) {
    expected_col <- spahe_col_mapping[[col_key]]
    if (!expected_col %in% available_cols) {
      matched_cols <- grep(expected_col, available_cols, ignore.case = TRUE, value = TRUE)
      if (length(matched_cols) > 0) {
        message(paste0("Auto-mapping: '", col_key, "' -> '", matched_cols[1], "'"))
        spahe_col_mapping[[col_key]] <- matched_cols[1]
      } else {
        stop(paste0("Column matching '", expected_col, "' not found in SpaHE data"))
      }
    }
  }
  
  message("ST data: ", ncol(st_data), " cells, number of cell types: ", length(unique(st_data$auc_celltype)))
  message("SpaHE predictions: ", nrow(result_df), " cells")
  
  # --------------------------
  # 2. Data preprocessing
  # --------------------------
  message("\n=== Data preprocessing ===")
  # ST data coordinate extraction and grid calculation
  st_coords <- GetTissueCoordinates(st_data)
  
  # Check coordinate data integrity
  valid_cells <- colnames(st_data) %in% rownames(st_coords)
  valid_cell_count <- sum(valid_cells)
  missing_cell_count <- sum(!valid_cells)
  
  if (missing_cell_count > 0) {
    warning("Warning: ", missing_cell_count, " cells in ST data lack coordinate information and will be filtered out")
    st_data <- st_data[, valid_cells]
    message("Retained cells with coordinates: ", valid_cell_count)
  }
  
  # Re-extract coordinates to ensure match with cell count
  st_coords <- GetTissueCoordinates(st_data)
  
  st_meta <- data.frame(
    cell_id = colnames(st_data),
    x = st_coords[, 1],
    y = st_coords[, 2],
    cell_type = st_data$auc_celltype,
    row.names = NULL,
    check.names = FALSE
  ) %>% filter(!is.na(cell_type))
  st_meta <- st_meta %>%
    mutate(
      grid_x = round(x / grid_size) * grid_size,
      grid_y = round(y / grid_size) * grid_size
    )
  
  # SpaHE data coordinate extraction and grid calculation
  spahe_meta <- result_df[, c(spahe_col_mapping$cell_id, spahe_col_mapping$x, spahe_col_mapping$y, spahe_col_mapping$cell_type)]
  colnames(spahe_meta) <- c("cell_id", "x", "y", "cell_type")
  spahe_meta <- spahe_meta %>% filter(!is.na(cell_type))
  # Coordinate scaling
  st_x_range <- if (!is.null(st_x_range)) st_x_range else range(st_meta$x, na.rm = TRUE)
  st_y_range <- if (!is.null(st_y_range)) st_y_range else range(st_meta$y, na.rm = TRUE)
  if (scale_spahe_coords) {
    message("Scaling SpaHE coordinates to ST range: x[", st_x_range[1], "~", st_x_range[2], "], y[", st_y_range[1], "~", st_y_range[2],"]")
    spahe_meta$x <- scales::rescale(spahe_meta$x, to = st_x_range)
    spahe_meta$y <- scales::rescale(spahe_meta$y, to = st_y_range)
  }
  spahe_meta <- spahe_meta %>%
    mutate(
      grid_x = round(x / grid_size) * grid_size,
      grid_y = round(y / grid_size) * grid_size
    )
  
  # Coordinate range validation
  message("ST grid coordinate range: ")
  message("  x: ", min(st_meta$grid_x), "~", max(st_meta$grid_x), 
          " (number of grids: ", length(unique(st_meta$grid_x)), ")")
  message("  y: ", min(st_meta$grid_y), "~", max(st_meta$grid_y), 
          " (number of grids: ", length(unique(st_meta$grid_y)), ")")
  message("SpaHE grid coordinate range: ")
  message("  x: ", min(spahe_meta$grid_x), "~", max(spahe_meta$grid_x), 
          " (number of grids: ", length(unique(spahe_meta$grid_x)), ")")
  message("  y: ", min(spahe_meta$grid_y), "~", max(spahe_meta$grid_y), 
          " (number of grids: ", length(unique(spahe_meta$grid_y)), ")")
  
  # Overlapping grid calculation
  st_grid_ids <- paste(st_grid_wide$grid_x, st_grid_wide$grid_y, sep = "_")
  spahe_grid_ids <- paste(spahe_grid_wide$grid_x, spahe_grid_wide$grid_y, sep = "_")
  overlap_grids <- intersect(st_grid_ids, spahe_grid_ids)
  message("Number of overlapping grids: ", length(overlap_grids))
  
  if (length(overlap_grids) < 500) {
    warning("⚠️ Insufficient overlapping grids (less than 500)! Consider reducing grid_size (current: ", grid_size, ")")
  }
  
  # --------------------------
  # 3. Spatial grid aggregation
  # --------------------------
  message("\n=== Spatial grid aggregation ===")
  # ST grid aggregation
  st_grid <- st_meta %>%
    group_by(grid_x, grid_y, cell_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(grid_x, grid_y) %>%
    mutate(
      total = sum(count),
      proportion = count / total
    ) %>%
    filter(total >= min_cells_per_grid)
  
  # Check for required columns
  required_cols <- c("grid_x", "grid_y", "cell_type", "proportion")
  missing_cols <- setdiff(required_cols, colnames(st_grid))
  if (length(missing_cols) > 0) {
    stop(paste("ST grid data missing required columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  st_grid_wide <- st_grid %>%
    dplyr::select(all_of(required_cols)) %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = proportion, values_fill = 0)
  message("Number of ST grids (after aggregation): ", nrow(st_grid_wide))
  
  # SpaHE grid aggregation
  spahe_grid <- spahe_meta %>%
    group_by(grid_x, grid_y, cell_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(grid_x, grid_y) %>%
    mutate(
      total = sum(count),
      proportion = count / total
    ) %>%
    filter(total >= min_cells_per_grid)
  
  missing_cols <- setdiff(required_cols, colnames(spahe_grid))
  if (length(missing_cols) > 0) {
    stop(paste("SpaHE grid data missing required columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  spahe_grid_wide <- spahe_grid %>%
    dplyr::select(all_of(required_cols)) %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = proportion, values_fill = 0)
  message("Number of SpaHE grids (after aggregation): ", nrow(spahe_grid_wide))
  
  # --------------------------
  # 4. Merging grid data
  # --------------------------
  message("\n=== Merging grid data ===")
  if (nrow(st_grid_wide) == 0 || nrow(spahe_grid_wide) == 0) {
    stop("ST or SpaHE grid data is empty! Please check parameters or data")
  }
  
  merged_grid <- merge(
    st_grid_wide, 
    spahe_grid_wide, 
    by = c("grid_x", "grid_y"), 
    suffixes = c("_st", "_spahe"),
    all = FALSE
  )
  
  all_celltypes <- union(
    setdiff(colnames(st_grid_wide), c("grid_x", "grid_y")),
    setdiff(colnames(spahe_grid_wide), c("grid_x", "grid_y"))
  )
  if (is.null(cell_types)) cell_types <- all_celltypes
  for (ct in cell_types) {
    if (!paste0(ct, "_st") %in% colnames(merged_grid)) merged_grid[[paste0(ct, "_st")]] <- 0
    if (!paste0(ct, "_spahe") %in% colnames(merged_grid)) merged_grid[[paste0(ct, "_spahe")]] <- 0
  }
  message("Number of merged grids: ", nrow(merged_grid), " (target: around 1000)")
  
  # --------------------------
  # 5. PCC calculation
  # --------------------------
  message("\n=== Calculating Pearson correlation coefficients ===")
  st_counts <- st_meta %>%
    group_by(grid_x, grid_y) %>%
    summarise(st_total = n(), .groups = "drop")
  spahe_counts <- spahe_meta %>%
    group_by(grid_x, grid_y) %>%
    summarise(spahe_total = n(), .groups = "drop")
  merged_grid <- merged_grid %>%
    left_join(st_counts, by = c("grid_x", "grid_y")) %>%
    left_join(spahe_counts, by = c("grid_x", "grid_y")) %>%
    mutate(
      st_total = ifelse(is.na(st_total), 0, st_total),
      spahe_total = ifelse(is.na(spahe_total), 0, spahe_total),
      total_counts = st_total + spahe_total,
      low_quality = total_counts < min_counts
    )
  
  # Calculate overall PCC
  calculate_pcc <- function(grid_row, cell_types) {
    if (grid_row["low_quality"]) return(NA)
    st_vals <- as.numeric(grid_row[paste0(cell_types, "_st")])
    spahe_vals <- as.numeric(grid_row[paste0(cell_types, "_spahe")])
    if (all(st_vals == 0) && all(spahe_vals == 0)) return(NA)
    tryCatch(cor(st_vals, spahe_vals, method = "pearson"), error = function(e) NA)
  }
  merged_grid$pcc <- apply(merged_grid, 1, function(row) calculate_pcc(row, cell_types))
  
  valid_pcc <- sum(!is.na(merged_grid$pcc))
  message("Number of valid PCC grids: ", valid_pcc, "/", nrow(merged_grid))
  
  # Calculate individual PCC for each cell type
  message("\n=== Calculating individual PCC for each cell type ===")
  celltype_stats <- list()
  celltype_data <- list()
  
  for (ct in cell_types) {
    st_col <- paste0(ct, "_st")
    spahe_col <- paste0(ct, "_spahe")
    
    st_values <- merged_grid[[st_col]][!merged_grid$low_quality]
    spahe_values <- merged_grid[[spahe_col]][!merged_grid$low_quality]
    
    st_nonzero <- sum(st_values > 0, na.rm = TRUE)
    spahe_nonzero <- sum(spahe_values > 0, na.rm = TRUE)
    
    if (st_nonzero == 0) {
      message(sprintf("Skipping cell type '%s': ST non-zero grid count = 0", ct))
      celltype_stats[[ct]] <- list(
        st_nonzero = st_nonzero,
        spahe_nonzero = spahe_nonzero,
        pcc = NA,
        valid = FALSE
      )
      next
    }
    
    valid_grids <- (st_values > 0) | (spahe_values > 0)
    valid_st <- st_values[valid_grids]
    valid_spahe <- spahe_values[valid_grids]
    n_valid <- sum(valid_grids)
    pcc_value <- NA
    
    if (n_valid >= 2) {
      if (sd(valid_st) > min_variance_threshold && sd(valid_spahe) > min_variance_threshold) {
        pcc_value <- tryCatch(
          cor(valid_st, valid_spahe, method = "pearson"),
          error = function(e) {
            message(sprintf("Error calculating PCC for cell type '%s': %s", ct, e$message))
            return(NA)
          }
        )
      }
    }
    
    merged_grid[[paste0("pcc_", ct)]] <- ifelse(valid_grids, pcc_value, NA)
    celltype_stats[[ct]] <- list(
      st_nonzero = st_nonzero,
      spahe_nonzero = spahe_nonzero,
      pcc = pcc_value,
      valid = TRUE
    )
    celltype_data[[ct]] <- data.frame(
      st = st_values,
      spahe = spahe_values,
      valid = valid_grids,
      pcc = pcc_value
    )
    
    message(sprintf("Cell type '%s': ST non-zero grids = %d, SpaHE non-zero grids = %d, PCC = %.4f", 
                    ct, st_nonzero, spahe_nonzero, pcc_value))
  }
  
  # Output cell type statistical summary
  celltype_stats_df <- do.call(rbind, lapply(names(celltype_stats), function(ct) {
    data.frame(
      cell_type = ct,
      st_nonzero_grids = celltype_stats[[ct]]$st_nonzero,
      spahe_nonzero_grids = celltype_stats[[ct]]$spahe_nonzero,
      pcc = celltype_stats[[ct]]$pcc,
      valid = celltype_stats[[ct]]$valid
    )
  })) %>% filter(valid)
  
  message("\nCell type statistical summary:")
  print(celltype_stats_df)
  
  # --------------------------
  # 6. Visualization
  # --------------------------
  message("\n=== Generating visualization results ===")
  
  # Helper function: check for valid PCC data
  has_valid_pcc <- function(pcc_values) {
    sum(!is.na(pcc_values)) > 0
  }
  
  # PCC distribution map
  if (has_valid_pcc(merged_grid$pcc)) {
    p_pcc_map <- ggplot(merged_grid, aes(x = grid_x, y = grid_y, color = pcc)) +
      geom_point(size = 1.5) +  # Increase point size
      scale_color_gradient2(
        low = "#2c7fb8", mid = "#ffffcc", high = "#e41a1c",
        na.value = "#aaaaaa", limits = c(-1, 1)
      ) +
      scale_y_reverse() +
      labs(title = "Overall PCC Distribution", x = "X", y = "Y") +
      theme_first_row
  } else {
    p_pcc_map <- ggplot() + 
      geom_text(aes(0, 0), label = "No Valid PCC Data", color = "red", size = 5) +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = "Overall PCC Distribution") +
      theme_first_row
  }
  
  # Overall PCC density plot
  valid_pcc_data <- merged_grid$pcc[!is.na(merged_grid$pcc)]
  if (length(valid_pcc_data) > 0) {
    p_density_overall <- ggplot(data.frame(pcc = valid_pcc_data), aes(x = pcc)) +
      geom_density(fill = "#4D8066", alpha = 0.7, linewidth = 0.8) +  # Thicken lines
      xlim(-1, 1) + 
      labs(title = "Overall PCC Density", x = "PCC", y = "Density") +
      theme_first_row
  } else {
    p_density_overall <- ggplot() + 
      geom_text(aes(0, 0), label = "No Valid Overall PCC", color = "red", size = 5) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = "Overall PCC Density") +
      theme_first_row
  }
  
  # Cell type PCC plots
  if (nrow(celltype_stats_df) > 0) {
    # Create PCC bar plot
    p_celltype_bar <- ggplot(celltype_stats_df, aes(x = reorder(cell_type, -pcc), y = pcc, fill = cell_type)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = sprintf("%.2f", pcc)), vjust = -0.5, size = 4) +  # Increase text size
      scale_fill_manual(values = celltype_colors) +
      labs(title = "Cell Type - Pearson Correlation Coefficient", x = "Cell Type", y = "PCC") +
      theme_first_row +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),  # Increase axis labels
        legend.position = "none"
      )
    
    # Create list of scatter plots (enlarge subplots)
    scatter_plots <- lapply(names(celltype_data), function(ct) {
      if (!celltype_stats[[ct]]$valid) return(NULL)
      
      ct_data <- celltype_data[[ct]]
      pcc_val <- celltype_stats[[ct]]$pcc
      ct_color <- ifelse(ct %in% names(celltype_colors), celltype_colors[ct], "#000000")
      
      valid_points <- ct_data[!is.na(ct_data$st) & !is.na(ct_data$spahe), ]
      if (nrow(valid_points) == 0) {
        p <- ggplot() + 
          geom_text(aes(0, 0), label = "No Valid Data", color = "red", size = 5) +
          xlim(0, 1) + ylim(0, 1) +
          labs(title = ct) +
          theme_uniform
        return(p)
      }
      
      p <- ggplot(ct_data, aes(x = st, y = spahe)) +
        geom_point(aes(color = valid), alpha = 0.7, size = 2) +  # Increase point size
        scale_color_manual(values = c("gray50", ct_color)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.6) +  # Thicken lines
        annotate("text", x = max(ct_data$st, na.rm = TRUE) * 0.1, 
                 y = max(ct_data$spahe, na.rm = TRUE) * 0.9, 
                 label = sprintf("PCC = %.2f", pcc_val), size = 4) +  # Increase text size
        labs(title = ct, x = "ST Proportion", y = "SpaHE Proportion") +
        theme_uniform +
        theme(legend.position = "none")
      
      return(p)
    })
    
    # Filter out NULL plots
    scatter_plots <- scatter_plots[!sapply(scatter_plots, is.null)]
    
    # Create list of density plots (enlarge subplots)
    density_plots <- lapply(names(celltype_data), function(ct) {
      if (!celltype_stats[[ct]]$valid) return(NULL)
      
      pcc_val <- celltype_stats[[ct]]$pcc
      if (is.na(pcc_val)) {
        p <- ggplot() + 
          geom_text(aes(0, 0), label = "NA", color = "red", size = 5) +
          xlim(-1, 1) + ylim(0, 1) +
          labs(title = ct) +
          theme_uniform
        return(p)
      }
      
      ct_color <- ifelse(ct %in% names(celltype_colors), celltype_colors[ct], "#000000")
      
      p <- ggplot() +
        geom_density(aes(x = rep(pcc_val, 1000)), fill = ct_color, alpha = 0.7, linewidth = 0.8) +  # Thicken lines
        geom_vline(xintercept = pcc_val, linetype = "dashed", color = "black", size = 0.6) +  # Thicken lines
        annotate("text", x = pcc_val, y = 0, label = sprintf("%.2f", pcc_val), 
                 vjust = -0.5, hjust = 0.5, size = 4) +  # Increase text size
        xlim(-1, 1) +
        labs(title = ct, x = "PCC", y = "Density") +
        theme_uniform +
        theme(plot.title = element_text(hjust = 0.5, size = 12))  # Increase title size
      
      return(p)
    })
    
    # Filter out NULL plots
    density_plots <- density_plots[!sapply(density_plots, is.null)]
    
    # Combine scatter plots and density plots (reduce columns, increase individual subplot size)
    if (length(scatter_plots) > 0 && length(density_plots) > 0) {
      if (length(scatter_plots) == length(density_plots)) {
        combined_plots <- list()
        for (i in seq_along(scatter_plots)) {
          combined_plots[[i]] <- scatter_plots[[i]] + density_plots[[i]] + plot_layout(widths = c(1, 1))
        }
        # Key adjustment: reduce to 2 columns (from 3) to increase individual subplot size
        p_scatter_density <- wrap_plots(combined_plots, ncol = 2) + 
          plot_annotation(title = "ST vs SpaHE Proportions & PCC Density by Cell Type", 
                          theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
      } else {
        p_scatter <- wrap_plots(scatter_plots, ncol = 2) +  # Reduce columns
          plot_annotation(title = "ST vs SpaHE Proportions by Cell Type", 
                          theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
        p_density <- wrap_plots(density_plots, ncol = 2) +  # Reduce columns
          plot_annotation(title = "PCC Density by Cell Type", 
                          theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
        p_scatter_density <- p_scatter / p_density
      }
    } else {
      p_scatter_density <- ggplot() + 
        geom_text(aes(0, 0), label = "No Valid Cell Type Data", color = "red", size = 6) +
        labs(title = "Cell Type Analysis") +
        theme_uniform
    }
  } else {
    p_celltype_bar <- ggplot() + 
      geom_text(aes(0, 0), label = "No Valid Cell Type Data", color = "red", size = 5) +
      labs(title = "Cell Type - Pearson Correlation Coefficient") +
      theme_first_row
    
    p_scatter_density <- ggplot() + 
      geom_text(aes(0, 0), label = "No Valid Cell Type Data", color = "red", size = 6) +
      labs(title = "Cell Type Analysis") +
      theme_uniform
  }
  
  # Integrate into a single PDF page (increase overall size)
  p_top <- p_pcc_map + p_density_overall + p_celltype_bar + plot_layout(widths = c(1, 1, 1))
  p_bottom <- p_scatter_density
  p_combined <- p_top / p_bottom + plot_layout(heights = c(1, 2.5))  # Increase bottom height ratio
  
  # Save as larger PDF
  ggsave(
    paste0(output_prefix, "_combined_plots.pdf"), 
    p_combined, 
    width = 22,  
    height = 28,  
    dpi = 300,
    device = "pdf"
  )
  
  message("Integrated chart saved: ", paste0(output_prefix, "_combined_plots.pdf"))
  
  # Return results
  list(
    grid_data = merged_grid, 
    combined_plot = p_combined,
    celltype_pcc_means = celltype_stats_df,
    celltype_data = celltype_data
  )
}