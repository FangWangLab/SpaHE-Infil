# Construction of training set based on spatial transcriptome annotation and HE image features
# Input: rds format (containing spatial transcriptome data and HE images)
# Output: Merged data statistics, prediction model (rds) with HE image features (nuclear/cytoplasmic morphology, cell density), and CSV with model and feature information

# Load necessary packages
library(Seurat)
library(EBImage)
library(OpenImageR)
library(randomForest)
library(imager)
library(caret)
library(tidyverse)
library(FNN)
library(dplyr)


# Define %||% operator
`%||%` <- function(x, y) if (is.null(x)) y else x


# Define TME cell type list
TME_celltypes <- c("CD8 T cells", "CD4 T cells", "Macrophages", 
                   "Neutrophils", "NK cells", "Tregs", "Fibroblasts", 
                   "B cells", "Dendritic cells")


# Core function: Merge RDS files and output cell type statistics
assess_celltype <- function(rds_files, output_stats = "celltype_combined_stats.csv") {
  # 1. Merge cell type annotations from all RDS files
  celltype_data <- purrr::map_dfr(rds_files, function(file) {
    dat <- readRDS(file)
    tibble(cell_type = dat$auc_celltype, source_file = basename(file))
  })
  
  # 2. Calculate overall counts and proportions
  combined_stats <- celltype_data %>%
    count(cell_type, name = "total_count") %>%
    mutate(
      proportion = total_count / sum(total_count),
      proportion_percent = proportion * 100,
      is_TME = cell_type %in% TME_celltypes
    ) %>%
    arrange(desc(total_count))
  
  # 3. Calculate cell type distribution per file
  file_stats <- celltype_data %>%
    count(source_file, cell_type, name = "file_count") %>%
    pivot_wider(names_from = source_file, values_from = file_count, values_fill = 0)
  
  # 4. Merge statistics and save results
  final_stats <- combined_stats %>%
    left_join(file_stats, by = "cell_type") %>%
    select(cell_type, total_count, proportion_percent, is_TME, everything())
  
  write.csv(final_stats, output_stats, row.names = FALSE)
  cat("===== Merged cell type statistics completed =====\n")
  cat("Statistics saved to:", output_stats, "\n", sep = "")
  print(final_stats %>% select(cell_type, total_count, proportion_percent, is_TME))
  
  return(final_stats)
}


# Function: Calculate texture features (keep original logic)
compute_features <- function(img) {
  if (length(dim(img)) > 2) {
    img <- img[,,1]
  }
  
  features <- tryCatch({
    EBImage::computeFeatures.haralick(
      img,
      ref = img,
      haralick.nbins = 32,
      haralick.scales = c(1, 2)
    )
  }, error = function(e) {
    message("Warning: computeFeatures.haralick failed - ", e$message)
    matrix(NA, nrow = 1, 
           ncol = 13,
           dimnames = list(NULL, c("h.1.contrast", "h.1.correlation", "h.1.energy", 
                                   "h.1.homogeneity", "h.2.contrast", "h.2.correlation",
                                   "h.2.energy", "h.2.homogeneity", "h.entropy",
                                   "h.variance", "h.sumAverage", "h.sumVariance",
                                   "h.sumEntropy")))
  })
  
  if (all(is.na(features))) {
    return(list(
      contrast = 0,
      correlation = 0,
      energy = 0,
      homogeneity = 0
    ))
  }
  
  col_patterns <- c(
    contrast = "contrast",
    correlation = "correlation",
    energy = "energy",
    homogeneity = "homogeneity"
  )
  
  feature_values <- sapply(col_patterns, function(pattern) {
    matched_col <- grep(pattern, colnames(features), value = TRUE)
    if (length(matched_col) > 0) {
      val <- mean(features[, matched_col], na.rm = TRUE)
      ifelse(is.na(val), 0, val)
    } else {
      0
    }
  })
  
  as.list(feature_values)
}


# Function: Calculate morphology features (reused for nuclei/cytoplasm)
compute_morphology_features <- function(bin_img) {
  # 1. Check input image validity
  if (!inherits(bin_img, "Image")) {
    bin_img <- EBImage::Image(bin_img, colormode = "Grayscale")
  }
  img_dims <- dim(bin_img)
  if (length(img_dims) < 2 || img_dims[1] < 3 || img_dims[2] < 3) {
    message("Warning: Image dimensions too small (", paste(img_dims, collapse="x"), "), cannot compute morphology features")
    return(list(area=0, perimeter=0, circularity=0, aspect_ratio=0, solidity=0))
  }
  
  # 2. Calculate connected components and morphology features
  morph_features <- tryCatch({
    connected <- EBImage::bwlabel(bin_img)  # Label connected components
    if (max(connected) == 0) {  # No valid regions
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    
    # 3. Verify connected component sizes
    region_sizes <- table(connected[connected > 0])
    if (length(region_sizes) == 0) {
      message("Warning: No valid connected components")
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    
    # 4. Select largest connected component
    main_region <- as.integer(names(which.max(region_sizes)))
    if (is.na(main_region) || main_region < 1 || main_region > max(connected)) {
      message("Warning: Invalid largest component index (", main_region, ")")
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    main_mask <- connected == main_region  # Extract largest region mask
    
    # 5. Check mask dimensions
    mask_dims <- dim(main_mask)
    if (mask_dims[1] != img_dims[1] || mask_dims[2] != img_dims[2]) {
      message("Warning: Mask dimensions don't match original image, cannot compute morphology features")
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    
    # 6. Calculate morphology features
    shape_stats <- EBImage::computeFeatures.shape(EBImage::Image(main_mask, colormode = "Grayscale"))
    
    # 7. Extract features (compatible with different version outputs)
    data.frame(
      area = ifelse("area" %in% rownames(shape_stats), shape_stats["area", ], 0),
      perimeter = ifelse("perimeter" %in% rownames(shape_stats), shape_stats["perimeter", ], 0),
      circularity = ifelse("circularity" %in% rownames(shape_stats), shape_stats["circularity", ], 0),
      aspect_ratio = ifelse("aspect.ratio" %in% rownames(shape_stats), shape_stats["aspect.ratio", ], 0),
      solidity = ifelse("solidity" %in% rownames(shape_stats), shape_stats["solidity", ], 0)
    )
  }, error = function(e) {
    message("Warning: Morphology feature calculation failed - ", e$message)
    data.frame(area=0, perimeter=0, circularity=0, aspect_ratio=0, solidity=0)
  })
  
  # 8. Final NA handling
  morph_features[is.na(morph_features)] <- 0
  as.list(morph_features)
}


# New function: Segment nuclei and cytoplasm based on HE staining characteristics
segment_he_nucleus_cytoplasm <- function(he_region) {
  # HE staining characteristics: Nuclei (hematoxylin) blue/dark, cytoplasm (eosin) red/light
  if (!inherits(he_region, "Image")) {
    he_region <- EBImage::Image(he_region, colormode = "Color")
  }
  
  # Extract RGB channels (HE images typically in RGB format)
  red <- he_region[,,1]
  blue <- he_region[,,3]  # Blue channel better for nuclei segmentation
  
  # Nuclei segmentation (enhance nuclear signal using blue channel)
  nucleus_bin <- tryCatch({
    # Normalize and threshold (nuclei brighter in blue channel)
    blue_norm <- EBImage::normalize(blue)
    EBImage::thresh(blue_norm, w = 5, h = 5, offset = 0.3)  # High threshold keeps strong signals (nuclei)
  }, error = function(e) {
    message("Warning: Nuclei segmentation failed - ", e$message)
    NULL
  })
  
  # Cytoplasm segmentation (cell region minus nuclei)
  cytoplasm_bin <- tryCatch({
    # Use difference between red and blue channels to extract cell regions
    he_gray <- EBImage::channel(he_region, "gray")
    gray_norm <- EBImage::normalize(he_gray)
    cell_bin <- EBImage::thresh(gray_norm, w = 7, h = 7, offset = 0.1)  # Low threshold keeps cell regions
    
    if (!is.null(nucleus_bin)) {
      # Cytoplasm = cell region - nuclear region
      cell_bin <- cell_bin - nucleus_bin
      cell_bin[cell_bin < 0] <- 0  # Ensure non-negative
    }
    cell_bin
  }, error = function(e) {
    message("Warning: Cytoplasm segmentation failed - ", e$message)
    NULL
  })
  
  list(nucleus = nucleus_bin, cytoplasm = cytoplasm_bin)
}


# New function: Calculate cell distribution density features
compute_density_features <- function(coords, index, density_radius = 20) {
  # coords: All cell coordinates; index: Current cell index
  if (nrow(coords) < 2) {
    return(list(
      cell_density = 0,          # Cells per unit area
      nearest_neighbor_dist = 0, # Nearest neighbor distance
      neighbor_dist_sd = 0       # Neighborhood distance SD (uniformity)
    ))
  }
  
  # Extract current cell coordinates
  current_x <- coords$imagecol[index]
  current_y <- coords$imagerow[index]
  
  # Calculate Euclidean distance to all other cells
  distances <- sqrt(
    (coords$imagecol[-index] - current_x)^2 + 
      (coords$imagerow[-index] - current_y)^2
  )
  
  # 1. Cell density (cells within radius / circle area)
  in_radius <- sum(distances <= density_radius)
  area <- pi * density_radius^2
  cell_density <- ifelse(area > 0, in_radius / area, 0)
  
  # 2. Nearest neighbor distance
  nearest_neighbor_dist <- if (length(distances) > 0) min(distances) else 0
  
  # 3. Neighborhood distance SD (reflects distribution uniformity)
  neighbor_dist_sd <- if (length(distances) >= 2) sd(distances) else 0
  
  list(
    cell_density = cell_density,
    nearest_neighbor_dist = nearest_neighbor_dist,
    neighbor_dist_sd = neighbor_dist_sd
  )
}


# Optimized: Extract training data (integrate HE image nuclear/cytoplasmic features and density features)
extract_training_data <- function(spatial_rds, density_radius = 20) {
  spatial_data <- readRDS(spatial_rds)
  
  # Verify HE image data exists (typically stored in slice1)
  if (is.null(spatial_data@images$slice1)) {
    stop("HE image data not found, please check input file")
  }
  he_img <- spatial_data@images$slice1@image  # HE image
  he_array <- as.array(he_img)
  he_eb <- EBImage::Image(he_array, colormode = "Color")  # Convert to EBImage format
  
  # Get cell coordinates and annotated types
  coords <- GetTissueCoordinates(spatial_data)
  celltypes <- spatial_data$auc_celltype
  
  # Loop to extract features for each cell
  features <- purrr::map_dfr(1:nrow(coords), function(i) {
    x <- coords$imagecol[i]
    y <- coords$imagerow[i]
    radius <- 5  # Cell region extraction radius
    
    # 1. Calculate region boundaries (avoid exceeding image dimensions)
    img_rows <- dim(he_eb)[1]
    img_cols <- dim(he_eb)[2]
    x_min <- max(1, x - radius)
    x_max <- min(img_rows, x + radius)
    y_min <- max(1, y - radius)
    y_max <- min(img_cols, y + radius)
    
    # 2. Check region size validity
    region_width <- x_max - x_min + 1
    region_height <- y_max - y_min + 1
    if (region_width < 5 || region_height < 5) {
      message("Warning: Cell region too small (", region_width, "x", region_height, "), using default features")
      return(dplyr::tibble(
        cell_id = rownames(coords)[i], cell_type = celltypes[i],
        # Original color features
        red_mean=0, red_sd=0, green_mean=0, green_sd=0, blue_mean=0, blue_sd=0,
        # Intensity and texture features
        intensity_mean=0, intensity_sd=0, contrast=0, correlation=0, energy=0, homogeneity=0,
        # Nuclear morphology features
        nucleus_area=0, nucleus_perimeter=0, nucleus_circularity=0, nucleus_aspect_ratio=0, nucleus_solidity=0,
        # Cytoplasmic morphology features
        cytoplasm_area=0, cytoplasm_perimeter=0, cytoplasm_circularity=0, cytoplasm_aspect_ratio=0, cytoplasm_solidity=0,
        # Cell distribution density features
        cell_density=0, nearest_neighbor_dist=0, neighbor_dist_sd=0
      ))
    }
    
    # 3. Extract cell region from HE image
    region <- tryCatch({
      he_eb[x_min:x_max, y_min:y_max, , drop=FALSE]  # Maintain 3D structure
    }, error = function(e) {
      message("Warning: HE image region extraction failed - ", e$message)
      NULL
    })
    if (is.null(region)) {
      return(dplyr::tibble(
        cell_id = rownames(coords)[i], cell_type = celltypes[i],
        red_mean=0, red_sd=0, green_mean=0, green_sd=0, blue_mean=0, blue_sd=0,
        intensity_mean=0, intensity_sd=0, contrast=0, correlation=0, energy=0, homogeneity=0,
        nucleus_area=0, nucleus_perimeter=0, nucleus_circularity=0, nucleus_aspect_ratio=0, nucleus_solidity=0,
        cytoplasm_area=0, cytoplasm_perimeter=0, cytoplasm_circularity=0, cytoplasm_aspect_ratio=0, cytoplasm_solidity=0,
        cell_density=0, nearest_neighbor_dist=0, neighbor_dist_sd=0
      ))
    }
    
    # 4. Extract color features from HE image (RGB channels)
    red_channel <- region[,,1]
    green_channel <- region[,,2]
    blue_channel <- region[,,3]
    red_mean <- mean(red_channel, na.rm=TRUE) %||% 0
    red_sd <- sd(red_channel, na.rm=TRUE) %||% 0
    green_mean <- mean(green_channel, na.rm=TRUE) %||% 0
    green_sd <- sd(green_channel, na.rm=TRUE) %||% 0
    blue_mean <- mean(blue_channel, na.rm=TRUE) %||% 0
    blue_sd <- sd(blue_channel, na.rm=TRUE) %||% 0
    
    # 5. Extract intensity and texture features
    gray_img <- EBImage::channel(region, "gray")
    texture_features <- compute_features(gray_img)
    intensity_mean <- mean(gray_img, na.rm=TRUE) %||% 0
    intensity_sd <- sd(gray_img, na.rm=TRUE) %||% 0
    
    # 6. Segment and calculate nuclear and cytoplasmic morphology features
    segmented <- segment_he_nucleus_cytoplasm(region)  # Segmentation based on HE characteristics
    
    nucleus_features <- if (!is.null(segmented$nucleus)) {
      compute_morphology_features(segmented$nucleus)
    } else {
      list(area=0, perimeter=0, circularity=0, aspect_ratio=0, solidity=0)
    }
    
    cytoplasm_features <- if (!is.null(segmented$cytoplasm)) {
      compute_morphology_features(segmented$cytoplasm)
    } else {
      list(area=0, perimeter=0, circularity=0, aspect_ratio=0, solidity=0)
    }
    
    # 7. Calculate cell distribution density features
    density_features <- compute_density_features(coords, i, density_radius = density_radius)
    
    # 8. Combine all features and return
    dplyr::tibble(
      cell_id = rownames(coords)[i], cell_type = celltypes[i],
      # Color features
      red_mean=red_mean, red_sd=red_sd, green_mean=green_mean, green_sd=green_sd,
      blue_mean=blue_mean, blue_sd=blue_sd,
      # Intensity and texture features
      intensity_mean=intensity_mean, intensity_sd=intensity_sd,
      contrast=texture_features$contrast, correlation=texture_features$correlation,
      energy=texture_features$energy, homogeneity=texture_features$homogeneity,
      # Nuclear morphology features
      nucleus_area=nucleus_features$area, 
      nucleus_perimeter=nucleus_features$perimeter, 
      nucleus_circularity=nucleus_features$circularity, 
      nucleus_aspect_ratio=nucleus_features$aspect_ratio, 
      nucleus_solidity=nucleus_features$solidity,
      # Cytoplasmic morphology features
      cytoplasm_area=cytoplasm_features$area, 
      cytoplasm_perimeter=cytoplasm_features$perimeter, 
      cytoplasm_circularity=cytoplasm_features$circularity, 
      cytoplasm_aspect_ratio=cytoplasm_features$aspect_ratio, 
      cytoplasm_solidity=cytoplasm_features$solidity,
      # Cell distribution density features
      cell_density=density_features$cell_density,
      nearest_neighbor_dist=density_features$nearest_neighbor_dist,
      neighbor_dist_sd=density_features$neighbor_dist_sd
    )
  }, .id = "point_id")
  
  # Filter rows with missing values
  features[complete.cases(features), ]
}


# Main training function (supports TME_adjust=FALSE option, no calibration)
train_celltype_model <- function(rds_files, output_model = "he_based_celltype_model.rds", 
                                 min_samples_per_type = 1, min_total_samples = 5, allow_unknown = FALSE,
                                 density_radius = 20,
                                 TME_adjust = FALSE) {  # New FALSE option, no calibration by default
  
  # Validate TME_adjust parameter (supports FALSE or 0-1 numeric value)
  if (is.logical(TME_adjust)) {
    if (TME_adjust) {
      stop("When TME_adjust=TRUE, please specify a value between 0-1 (0=TME proportion 5%, 1=TME proportion 50%)")
    } else {
      # TME_adjust=FALSE: No calibration, use original proportions
      cat("TME_adjust=FALSE, will use original training set TME cell proportion without calibration\n")
      target_tme_ratio <- NULL  # Mark as no calibration
    }
  } else if (is.numeric(TME_adjust)) {
    if (TME_adjust < 0 || TME_adjust > 1) {
      stop("TME_adjust must be between 0-1 (0=TME proportion 5%, 1=TME proportion 50%)")
    }
    target_tme_ratio <- 0.05 + 0.45 * TME_adjust  # Calculate target TME proportion
    cat("TME adjustment parameter set:", TME_adjust, ", target training set TME cell proportion:", 
        round(target_tme_ratio*100, 1), "%\n", sep="")
  } else {
    stop("Invalid TME_adjust parameter, please enter FALSE or a value between 0-1")
  }
  
  # Extract features from all samples (including new HE image features)
  all_features <- purrr::map_dfr(rds_files, function(file) {
    extract_training_data(file, density_radius = density_radius)
  })
  
  # Data filtering and validation
  cat("\nOriginal data statistics:\n")
  print(table(all_features$cell_type, useNA = "always"))
  
  all_features <- all_features[!is.na(all_features$cell_type), ]
  type_counts <- table(all_features$cell_type)
  valid_types <- names(type_counts[type_counts >= min_samples_per_type])
  
  if (length(valid_types) == 0) {
    warning("All cell types have insufficient samples, keeping the most abundant type")
    valid_types <- names(which.max(type_counts))
  }
  
  # Keep Unknown type (if needed)
  if (allow_unknown && "Unknown" %in% all_features$cell_type) {
    valid_types <- unique(c(valid_types, "Unknown"))
    cat("Kept 'Unknown' type for training\n")
  }
  
  filtered_features <- all_features %>% dplyr::filter(cell_type %in% valid_types)
  
  # Mark TME cell types (already includes Dendritic cells)
  filtered_features <- filtered_features %>%
    dplyr::mutate(is_TME = cell_type %in% TME_celltypes)
  
  cat("\nFiltered data statistics (with TME markers):\n")
  tme_summary <- filtered_features %>%
    dplyr::group_by(cell_type, is_TME) %>%
    dplyr::count() %>%
    dplyr::arrange(dplyr::desc(is_TME), dplyr::desc(n))
  print(tme_summary)
  
  # Calculate original TME proportion
  original_tme_count <- sum(filtered_features$is_TME)
  original_tme_ratio <- original_tme_count / nrow(filtered_features)
  cat("\nOriginal training set TME cell proportion:", round(original_tme_ratio*100, 1), "% (", 
      original_tme_count, "/", nrow(filtered_features), ")\n", sep="")
  
  # Only adjust proportions if TME_adjust is numeric
  if (is.numeric(TME_adjust)) {
    # Separate TME and non-TME cells
    tme_data <- filtered_features %>% dplyr::filter(is_TME)
    non_tme_data <- filtered_features %>% dplyr::filter(!is_TME)
    
    # Handle edge cases (error if no TME cells)
    if (nrow(tme_data) == 0) {
      stop("No TME cell types detected in training data, cannot apply TME adjustment")
    }
    if (nrow(non_tme_data) == 0) {
      warning("Training data contains only TME cells, no adjustment needed")
    } else {
      # Calculate target sample size (maintain total samples, adjust TME/non-TME ratio)
      total_samples <- nrow(filtered_features)
      target_tme_count <- round(target_tme_ratio * total_samples)
      target_non_tme_count <- total_samples - target_tme_count
      
      # Ensure target counts meet minimum samples (avoid extreme adjustments)
      target_tme_count <- max(target_tme_count, min_samples_per_type)
      target_non_tme_count <- max(target_non_tme_count, min_samples_per_type)
      
      # Oversample TME cells (if insufficient) or undersample (if excessive)
      if (nrow(tme_data) < target_tme_count) {
        # Oversample: Sample with replacement to increase TME samples
        tme_adjusted <- tme_data[sample(1:nrow(tme_data), size = target_tme_count, replace = TRUE), ]
        cat("TME cell oversampling: from", nrow(tme_data), "to", target_tme_count, "\n")
      } else {
        # Undersample: Sample without replacement to reduce TME samples
        tme_adjusted <- tme_data[sample(1:nrow(tme_data), size = target_tme_count, replace = FALSE), ]
        cat("TME cell undersampling: from", nrow(tme_data), "to", target_tme_count, "\n")
      }
      
      # Oversample non-TME cells (if insufficient) or undersample (if excessive)
      if (nrow(non_tme_data) < target_non_tme_count) {
        non_tme_adjusted <- non_tme_data[sample(1:nrow(non_tme_data), size = target_non_tme_count, replace = TRUE), ]
        cat("Non-TME cell oversampling: from", nrow(non_tme_data), "to", target_non_tme_count, "\n")
      } else {
        non_tme_adjusted <- non_tme_data[sample(1:nrow(non_tme_data), size = target_non_tme_count, replace = FALSE), ]
        cat("Non-TME cell undersampling: from", nrow(non_tme_data), "to", target_non_tme_count, "\n")
      }
      
      # Merge adjusted samples and shuffle
      adjusted_features <- dplyr::bind_rows(tme_adjusted, non_tme_adjusted) %>%
        dplyr::sample_frac(1)  # Random shuffle
      
      # Verify adjusted proportion
      adjusted_tme_ratio <- sum(adjusted_features$is_TME) / nrow(adjusted_features)
      cat("Adjusted training set TME cell proportion:", round(adjusted_tme_ratio*100, 1), "% (", 
          sum(adjusted_features$is_TME), "/", nrow(adjusted_features), ")\n", sep="")
      
      # Use adjusted feature data
      filtered_features <- adjusted_features
    }
  } else {
    # When TME_adjust=FALSE, use original data directly
    cat("TME_adjust=FALSE, no proportion calibration, using original data for training\n")
  }
  
  # Final filtering (ensure sample size meets requirements)
  if (nrow(filtered_features) < min_total_samples) {
    warning(paste("Extremely small sample size (", nrow(filtered_features), " samples), model may not be meaningful"))
  }
  
  # Prepare training data
  train_data <- filtered_features %>%
    dplyr::select(-cell_id, -point_id, -is_TME) %>%  # Exclude non-feature columns
    dplyr::mutate(cell_type = factor(cell_type)) %>%
    na.omit()
  
  cat("\nProcessed training data dimensions:", dim(train_data), "\n")
  cat("Processed cell type sample counts:\n")
  print(table(train_data$cell_type))
  
  if(nrow(train_data) < 5) {
    warning("Processed data too small, model may be unstable")
  }
  
  # Set cross-validation strategy
  if (nrow(train_data) <= 5) {
    train_control <- trainControl(method = "none", savePredictions = "all", classProbs = TRUE)
    cat("\nSample size too small, no cross-validation\n")
  } else if (nrow(train_data) < 10) {
    train_control <- trainControl(method = "LOOCV", savePredictions = "final")
    cat("\nUsing leave-one-out cross-validation (LOOCV)\n")
  } else {
    train_control <- trainControl(method = "cv", number = min(3, floor(nrow(train_data)/2)), savePredictions = "final")
    cat("\nUsing", train_control$number, "-fold cross-validation\n")
  }
  
  # Train random forest model (automatically includes all HE-related features)
  set.seed(123)
  model <- train(
    cell_type ~ ., data = train_data, method = "rf",
    trControl = train_control, importance = TRUE,
    ntree = min(200, nrow(train_data)*5),  # Adaptive tree count based on sample size
    nodesize = max(1, floor(nrow(train_data)/10))  # Minimum samples per leaf node
  )
  
  # Save model (includes TME adjustment parameter information)
  model$TME_adjust <- TME_adjust
  model$target_tme_ratio <- if(is.numeric(TME_adjust)) target_tme_ratio else NULL
  model$original_tme_ratio <- original_tme_ratio  # Record original proportion
  saveRDS(model, output_model)
  cat("\nModel saved to:", output_model, "\n", sep="")
  
  # Generate and save cell type model feature information CSV
  # 1. Collect basic information
  basic_info <- data.frame(
    category = "basic_info",
    parameter = c(
      "input_files",
      "total_samples",
      "total_features",
      "TME_adjust",
      "original_TME_ratio_percent",
      "adjusted_TME_ratio_percent",
      "cv_method",
      "cv_folds",
      "model_method",
      "random_seed"
    ),
    value = c(
      paste(rds_files, collapse = "; "),
      as.character(nrow(train_data)),
      as.character(ncol(train_data) - 1),  # Subtract cell_type column
      as.character(TME_adjust),
      as.character(round(original_tme_ratio * 100, 2)),
      ifelse(exists("adjusted_tme_ratio"), as.character(round(adjusted_tme_ratio * 100, 2)), "Not adjusted"),
      train_control$method,
      ifelse(train_control$method == "cv", as.character(train_control$number), 
             ifelse(train_control$method == "LOOCV", "LOOCV", "none")),
      "rf",
      "123"
    ),
    stringsAsFactors = FALSE
  )
  
  # 2. Collect cell type distribution information
  cell_type_counts <- table(train_data$cell_type)
  cell_type_dist <- data.frame(
    category = "cell_type_distribution",
    parameter = names(cell_type_counts),
    value = as.character(cell_type_counts),
    proportion_percent = as.character(round(cell_type_counts / sum(cell_type_counts) * 100, 2)),
    stringsAsFactors = FALSE
  )
  
  # 3. Merge and save
  combined_features <- dplyr::bind_rows(basic_info, cell_type_dist)
  write.csv(combined_features, file = "celltype_model_features.csv", row.names = FALSE, quote = FALSE)
  cat("Cell type model feature information saved to: celltype_model_features.csv\n")
  
  # Output model information
  print(model)
  if ("pred" %in% names(model)) {
    print(confusionMatrix(model$pred$pred, model$pred$obs))
  }
  
  # Output feature importance (including nuclear/cytoplasmic morphology and density features)
  if (!is.null(model$finalModel)) {
    cat("\nFeature importance (top 15, including HE-related features):\n")
    importance <- varImp(model)$importance
    print(head(importance[order(-importance[,1]), , drop = FALSE], 15))
  }
  
  return(model)
}