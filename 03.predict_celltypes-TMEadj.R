# Cell Type Prediction Script Based on HE Images
# Input: 1) celltype_model.rds 2) HE image 3) TME_adjust(FALSE/0-1) 4) step_size(1-1000) 
#       5) point_size(0.1-5) 6) point_alpha(0.1-1)
# Output: Calibrated prediction results and visualization files

library(Seurat)
library(EBImage)
library(OpenImageR)
library(randomForest)
library(imager)
library(caret)
library(tidyverse)
library(FNN)
library(patchwork)  # Ensure the patchwork package is loaded for combining plots

# Define %||% operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# TME cell type list (Dendritic cells added)
TME_celltypes <- c("CD8 T cells", "CD4 T cells", "Macrophages", 
                   "Neutrophils", "NK cells", "Tregs", "Fibroblasts",
                   "B cells", "Dendritic cells")  # Added Dendritic cells

# Cell type color scheme (supplemented with Dendritic cells color)
celltype_colors <- c(
  "Tumor cells" = "#FF6F61",                 # Pink frost (main visual color, most prominent)
  "Normal epithelial cells" = "#FFD800",     # Weiss gold (highly saturated gold)
  "CD8 T cells" = "#7FDBFF",                 # Mint blue (cool tone)
  "CD4 T cells" = "#1E90FF",                 # Sea green (deep cool green, contrasting with tumor pink)
  "Macrophages" = "#9370DB",                 # Medium violet (medium saturation)
  "Neutrophils" = "#00CED1",                 # Deep sapphire blue (highly saturated cool color)
  "NK cells" = "#DA70D6",                    # Orchid purple (reddish purple, distinct from macrophages)
  "Tregs" = "#CD853F",                       # Peru brown (warm brown)
  "Fibroblasts" = "#4CAF50",                 # Tomato red (highly saturated warm color)
  "Endothelial cells" = "#9ACD32",           # Yellow-green (transitional color)
  "B cells" = "#F8C3CD",                     # Bright coral (highly saturated pink)
  "Dendritic cells" = "#20B2AA"              # Light sea green (cool tone ending)
)


### 1. Texture feature calculation function
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


### 2. Morphological feature calculation function
compute_morphology_features <- function(bin_img) {
  if (!inherits(bin_img, "Image")) {
    bin_img <- EBImage::Image(bin_img, colormode = "Grayscale")
  }
  img_dims <- dim(bin_img)
  if (length(img_dims) < 2 || img_dims[1] < 3 || img_dims[2] < 3) {
    message("Warning: Image size is too small (", paste(img_dims, collapse="x"), "), cannot compute morphological features")
    return(list(area=0, perimeter=0, circularity=0, aspect_ratio=0, solidity=0))
  }
  
  morph_features <- tryCatch({
    connected <- EBImage::bwlabel(bin_img)
    if (max(connected) == 0) {
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    
    region_sizes <- table(connected[connected > 0])
    if (length(region_sizes) == 0) {
      message("Warning: No valid connected regions")
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    
    main_region <- as.integer(names(which.max(region_sizes)))
    if (is.na(main_region) || main_region < 1 || main_region > max(connected)) {
      message("Warning: Invalid index for the largest connected region (", main_region, ")")
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    main_mask <- connected == main_region
    
    mask_dims <- dim(main_mask)
    if (mask_dims[1] != img_dims[1] || mask_dims[2] != img_dims[2]) {
      message("Warning: Mask dimensions do not match the original image, cannot compute morphological features")
      return(data.frame(
        area = 0, perimeter = 0, circularity = 0, aspect_ratio = 0, solidity = 0
      ))
    }
    
    shape_stats <- EBImage::computeFeatures.shape(EBImage::Image(main_mask, colormode = "Grayscale"))
    
    data.frame(
      area = ifelse("area" %in% rownames(shape_stats), shape_stats["area", ], 0),
      perimeter = ifelse("perimeter" %in% rownames(shape_stats), shape_stats["perimeter", ], 0),
      circularity = ifelse("circularity" %in% rownames(shape_stats), shape_stats["circularity", ], 0),
      aspect_ratio = ifelse("aspect.ratio" %in% rownames(shape_stats), shape_stats["aspect.ratio", ], 0),
      solidity = ifelse("solidity" %in% rownames(shape_stats), shape_stats["solidity", ], 0)
    )
  }, error = function(e) {
    message("Warning: Skipping morphological feature calculation - ", e$message)
    data.frame(area=0, perimeter=0, circularity=0, aspect_ratio=0, solidity=0)
  })
  
  morph_features[is.na(morph_features)] <- 0
  as.list(morph_features)
}


### 3. Segment nucleus and cytoplasm based on HE staining
segment_he_nucleus_cytoplasm <- function(he_region) {
  if (!inherits(he_region, "Image")) {
    he_region <- EBImage::Image(he_region, colormode = "Color")
  }
  
  red <- he_region[,,1]
  blue <- he_region[,,3]
  
  # Nucleus segmentation (blue channel enhancement)
  nucleus_bin <- tryCatch({
    blue_norm <- EBImage::normalize(blue)
    EBImage::thresh(blue_norm, w = 5, h = 5, offset = 0.3)
  }, error = function(e) {
    message("Warning: Skipping nucleus segmentation - ", e$message)
    NULL
  })
  
  # Cytoplasm segmentation (cell region minus nucleus)
  cytoplasm_bin <- tryCatch({
    he_gray <- EBImage::channel(he_region, "gray")
    gray_norm <- EBImage::normalize(he_gray)
    cell_bin <- EBImage::thresh(gray_norm, w = 7, h = 7, offset = 0.1)
    
    if (!is.null(nucleus_bin)) {
      cell_bin <- cell_bin - nucleus_bin
      cell_bin[cell_bin < 0] <- 0
    }
    cell_bin
  }, error = function(e) {
    message("Warning: Skipping cytoplasm segmentation - ", e$message)
    NULL
  })
  
  list(nucleus = nucleus_bin, cytoplasm = cytoplasm_bin)
}


### 4. Calculate cell distribution density features
compute_density_features <- function(coords, index, density_radius = 20) {
  if (nrow(coords) < 2) {
    return(list(
      cell_density = 0,
      nearest_neighbor_dist = 0,
      neighbor_dist_sd = 0
    ))
  }
  
  current_x <- coords$imagecol[index]
  current_y <- coords$imagerow[index]
  
  distances <- sqrt(
    (coords$imagecol[-index] - current_x)^2 + 
      (coords$imagerow[-index] - current_y)^2
  )
  
  in_radius <- sum(distances <= density_radius)
  area <- pi * density_radius^2
  cell_density <- ifelse(area > 0, in_radius / area, 0)
  
  nearest_neighbor_dist <- if (length(distances) > 0) min(distances) else 0
  neighbor_dist_sd <- if (length(distances) >= 2) sd(distances) else 0
  
  list(
    cell_density = cell_density,
    nearest_neighbor_dist = nearest_neighbor_dist,
    neighbor_dist_sd = neighbor_dist_sd
  )
}


### 5. Extract features for prediction
extract_prediction_features <- function(he_image_path, coords = NULL, 
                                        bg_threshold = 230/255,
                                        radius = 5,
                                        step_size = 10,
                                        density_radius = 20) {
  if (!file.exists(he_image_path)) {
    stop("HE image does not exist: ", he_image_path)
  }
  he_img <- EBImage::readImage(he_image_path)
  if (length(dim(he_img)) != 3) {
    stop("HE image must be in RGB three-channel format")
  }
  cat("HE image dimensions: ", dim(he_img)[1:2], "(width x height)\n")
  cat("Prediction point step size setting: ", step_size, "pixels\n")
  
  he_img_255 <- he_img * 255
  
  if (is.null(coords)) {
    cat("Automatically generating grid coordinates based on step size...\n")
    x_seq <- seq(radius + 1, dim(he_img)[1] - radius, by = step_size)
    y_seq <- seq(radius + 1, dim(he_img)[2] - radius, by = step_size)
    
    if (tail(x_seq, 1) < dim(he_img)[1] - radius) {
      x_seq <- c(x_seq, dim(he_img)[1] - radius)
    }
    if (tail(y_seq, 1) < dim(he_img)[2] - radius) {
      y_seq <- c(y_seq, dim(he_img)[2] - radius)
    }
    
    coords <- expand.grid(imagecol = x_seq, imagerow = y_seq)
    rownames(coords) <- paste0("pred_", 1:nrow(coords))
    cat("Number of automatically generated prediction points: ", nrow(coords), "\n")
  } else {
    if (!all(c("imagecol", "imagerow") %in% colnames(coords))) {
      stop("Coordinate data frame must contain 'imagecol' and 'imagerow' columns")
    }
    cat("Using custom coordinates, total of", nrow(coords), "prediction points\n")
  }
  
  density_features_list <- if (nrow(coords) >= 1) {
    purrr::map(1:nrow(coords), function(i) {
      compute_density_features(coords, i, density_radius = density_radius)
    })
  } else {
    list()
  }
  
  pred_features <- purrr::map_dfr(1:nrow(coords), function(i) {
    x <- coords$imagecol[i]
    y <- coords$imagerow[i]
    cell_id <- rownames(coords)[i]
    
    img_rows <- dim(he_img)[1]
    img_cols <- dim(he_img)[2]
    x_min <- max(1, x - radius)
    x_max <- min(img_rows, x + radius)
    y_min <- max(1, y - radius)
    y_max <- min(img_cols, y + radius)
    
    region_width <- x_max - x_min + 1
    region_height <- y_max - y_min + 1
    if (region_width < 5 || region_height < 5) {
      message("Warning: Region size is too small (", region_width, "x", region_height, "), using default features (", cell_id, ")")
      return(dplyr::tibble(
        cell_id = cell_id,
        x = x, y = y,
        red_mean = 0, red_sd = 0,
        green_mean = 0, green_sd = 0,
        blue_mean = 0, blue_sd = 0,
        intensity_mean = 0, intensity_sd = 0,
        contrast = 0, correlation = 0,
        energy = 0, homogeneity = 0,
        nucleus_area = 0, nucleus_perimeter = 0, nucleus_circularity = 0,
        nucleus_aspect_ratio = 0, nucleus_solidity = 0,
        cytoplasm_area = 0, cytoplasm_perimeter = 0, cytoplasm_circularity = 0,
        cytoplasm_aspect_ratio = 0, cytoplasm_solidity = 0,
        cell_density = 0, nearest_neighbor_dist = 0, neighbor_dist_sd = 0,
        tme_feature_score = 0  # TME feature score
      ))
    }
    
    region <- tryCatch({
      he_img_255[x_min:x_max, y_min:y_max, , drop = FALSE]
    }, error = function(e) {
      message("Warning: Skipping region extraction (", cell_id, ") - ", e$message)
      NULL
    })
    
    is_background <- FALSE
    if (!is.null(region) && length(dim(region)) == 3) {
      red_mean <- mean(region[,,1], na.rm = TRUE)
      green_mean <- mean(region[,,2], na.rm = TRUE)
      blue_mean <- mean(region[,,3], na.rm = TRUE)
      
      if (red_mean >= 230 && green_mean >= 230 && blue_mean >= 230) {
        message("Info: Blank background region, skipping (", cell_id, ")")
        is_background <- TRUE
      }
    }
    if (is_background) {
      return(NULL)
    }
    
    if (is.null(region)) {
      return(dplyr::tibble(
        cell_id = cell_id,
        x = x, y = y,
        red_mean = 0, red_sd = 0,
        green_mean = 0, green_sd = 0,
        blue_mean = 0, blue_sd = 0,
        intensity_mean = 0, intensity_sd = 0,
        contrast = 0, correlation = 0,
        energy = 0, homogeneity = 0,
        nucleus_area = 0, nucleus_perimeter = 0, nucleus_circularity = 0,
        nucleus_aspect_ratio = 0, nucleus_solidity = 0,
        cytoplasm_area = 0, cytoplasm_perimeter = 0, cytoplasm_circularity = 0,
        cytoplasm_aspect_ratio = 0, cytoplasm_solidity = 0,
        cell_density = 0, nearest_neighbor_dist = 0, neighbor_dist_sd = 0,
        tme_feature_score = 0  # TME feature score
      ))
    }
    
    region_01 <- region / 255
    
    red_channel <- region_01[,,1]
    green_channel <- region_01[,,2]
    blue_channel <- region_01[,,3]
    red_mean <- mean(red_channel, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    red_sd <- sd(red_channel, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    green_mean <- mean(green_channel, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    green_sd <- sd(green_channel, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    blue_mean <- mean(blue_channel, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    blue_sd <- sd(blue_channel, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    
    gray_img <- tryCatch({
      EBImage::channel(EBImage::Image(region_01, colormode = "Color"), "gray")
    }, error = function(e) {
      message("Warning: Skipping grayscale conversion (", cell_id, ")")
      matrix(0, nrow = 5, ncol = 5)
    })
    texture_features <- compute_features(gray_img)
    intensity_mean <- mean(gray_img, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    intensity_sd <- sd(gray_img, na.rm = TRUE) %>% ifelse(is.na(.), 0, .)
    
    segmented <- segment_he_nucleus_cytoplasm(region_01)
    
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
    
    density_features <- if (i <= length(density_features_list)) {
      density_features_list[[i]]
    } else {
      list(cell_density=0, nearest_neighbor_dist=0, neighbor_dist_sd=0)
    }
    
    # Calculate TME feature score
    tme_feature_score <- (nucleus_features$circularity * 0.3) +  # Nuclear circularity
      (density_features$cell_density * 0.4) +   # Cell density
      (blue_mean * 0.3)                         # Blue channel intensity
    
    dplyr::tibble(
      cell_id = cell_id,
      x = x, y = y,
      red_mean = red_mean, red_sd = red_sd,
      green_mean = green_mean, green_sd = green_sd,
      blue_mean = blue_mean, blue_sd = blue_sd,
      intensity_mean = intensity_mean, intensity_sd = intensity_sd,
      contrast = texture_features$contrast %||% 0,
      correlation = texture_features$correlation %||% 0,
      energy = texture_features$energy %||% 0,
      homogeneity = texture_features$homogeneity %||% 0,
      nucleus_area = nucleus_features$area %||% 0,
      nucleus_perimeter = nucleus_features$perimeter %||% 0,
      nucleus_circularity = nucleus_features$circularity %||% 0,
      nucleus_aspect_ratio = nucleus_features$aspect_ratio %||% 0,
      nucleus_solidity = nucleus_features$solidity %||% 0,
      cytoplasm_area = cytoplasm_features$area %||% 0,
      cytoplasm_perimeter = cytoplasm_features$perimeter %||% 0,
      cytoplasm_circularity = cytoplasm_features$circularity %||% 0,
      cytoplasm_aspect_ratio = cytoplasm_features$aspect_ratio %||% 0,
      cytoplasm_solidity = cytoplasm_features$solidity %||% 0,
      cell_density = density_features$cell_density %||% 0,
      nearest_neighbor_dist = density_features$nearest_neighbor_dist %||% 0,
      neighbor_dist_sd = density_features$neighbor_dist_sd %||% 0,
      tme_feature_score = tme_feature_score  # TME feature score
    )
  })
  
  required_features <- c("red_mean", "red_sd", "green_mean", "green_sd", 
                         "blue_mean", "blue_sd", "intensity_mean", "intensity_sd",
                         "contrast", "correlation", "energy", "homogeneity",
                         "nucleus_area", "nucleus_perimeter", "nucleus_circularity",
                         "nucleus_aspect_ratio", "nucleus_solidity",
                         "cytoplasm_area", "cytoplasm_perimeter", "cytoplasm_circularity",
                         "cytoplasm_aspect_ratio", "cytoplasm_solidity",
                         "cell_density", "nearest_neighbor_dist", "neighbor_dist_sd")
  
  missing_features <- setdiff(required_features, colnames(pred_features))
  if (length(missing_features) > 0) {
    stop("Feature extraction failed, missing required features: ", paste(missing_features, collapse = ", "))
  } else {
    cat("Feature extraction successful,", nrow(pred_features), "valid region samples remaining after filtering invalid regions\n")
  }
  
  return(pred_features)
}


### 6. Main prediction function (supports TME_adjust=FALSE, using original training set proportions)
generate_predictions <- function(he_image_path, 
                                 model_path = "celltype_model_with_morphology_density.rds", 
                                 output_prefix = "he_prediction",
                                 coords = NULL,
                                 TME_adjust = FALSE,  # Standardized to logical FALSE, default no calibration
                                 step_size = 10,    # Step size (1-1000)
                                 density_radius = 20) {
  
  # Process TME_adjust parameter: only supports logical FALSE and numeric values between 0-1
  if (is.character(TME_adjust)) {
    # Compatible with old string input, automatically convert to logical
    if (tolower(TME_adjust) == "false") {
      TME_adjust <- FALSE
      cat("Warning: TME_adjust has been automatically converted from string 'false' to logical FALSE\n")
    } else {
      stop("Invalid TME_adjust parameter, only supports logical FALSE or numeric values between 0-1")
    }
  }
  
  # Parameter validity check
  if (is.logical(TME_adjust)) {
    if (TME_adjust) {
      stop("When TME_adjust is TRUE, please specify a numeric value between 0-1 (0=TME proportion 5%, 1=TME proportion 50%)")
    } else {
      # TME_adjust=FALSE: use original TME proportion from training set
      use_original_tme <- TRUE
    }
  } else if (is.numeric(TME_adjust)) {
    if (TME_adjust < 0 || TME_adjust > 1) {
      stop("TME_adjust must be between 0-1 (0=TME proportion 5%, 1=TME proportion 50%)")
    }
    use_original_tme <- FALSE
  } else {
    stop("Invalid TME_adjust parameter, please input FALSE (no calibration) or a numeric value between 0-1 (calibration proportion)")
  }
  
  # Load model and data
  if (!file.exists(model_path)) {
    stop("Model file does not exist: ", model_path)
  }
  model <- readRDS(model_path)
  train_data <- model$trainingData
  train_celltypes <- train_data$.outcome
  train_ratio <- prop.table(table(train_celltypes))
  
  # Define valid cell types
  valid_celltypes <- setdiff(names(train_ratio), c("Extracellular matrix", "Unknown", "Background"))
  train_tme_types <- intersect(valid_celltypes, TME_celltypes)  # Including newly added Dendritic cells
  train_non_tme_types <- setdiff(valid_celltypes, train_tme_types)
  
  if (length(train_tme_types) == 0) {
    stop("No TME cell types detected in training data, cannot perform TME adjustment")
  }
  
  # Standardize training proportions
  train_ratio <- train_ratio[valid_celltypes]
  train_ratio <- train_ratio / sum(train_ratio)
  cat("Standardized valid cell type proportions:\n")
  print(train_ratio)
  
  # Calculate target TME proportion (two scenarios)
  if (use_original_tme) {
    # Use original TME proportion from training set
    original_tme_ratio <- sum(train_ratio[train_tme_types])
    target_tme_ratio <- original_tme_ratio
    cat("TME_adjust=FALSE, will use original TME cell proportion from training set: ", 
        round(original_tme_ratio*100, 1), "%\n", sep="")
  } else {
    # Calculate target proportion based on TME_adjust
    target_tme_ratio <- 0.05 + 0.45 * TME_adjust
    cat("TME adjustment strength: ", TME_adjust, ", target total TME cell proportion: ", 
        round(target_tme_ratio*100, 1), "%\n", sep="")
  }
  
  # Extract prediction features
  pred_features <- extract_prediction_features(
    he_image_path = he_image_path,
    coords = coords,
    radius = 5,
    step_size = step_size,
    density_radius = density_radius
  )
  
  # Prepare prediction data and generate initial prediction probabilities
  pred_data <- pred_features %>% dplyr::select(-cell_id, -x, -y, -tme_feature_score)
  pred_probs <- predict(model, newdata = pred_data, type = "prob")
  pred_classes <- colnames(pred_probs)
  
  # Match training types
  pred_probs <- pred_probs[, intersect(pred_classes, valid_celltypes), drop = FALSE]
  pred_classes <- colnames(pred_probs)
  
  # Calculate target number of TME cells
  target_tme_count <- round(nrow(pred_features) * target_tme_ratio)
  
  # Strengthen TME proportion adjustment (only executed when not using original proportion)
  # 1. Initial prediction
  initial_predictions <- colnames(pred_probs)[max.col(pred_probs)]
  initial_tme_count <- sum(initial_predictions %in% TME_celltypes)  # Including Dendritic cells
  initial_tme_ratio <- initial_tme_count / nrow(pred_features) * 100
  cat("Initial predicted TME cell proportion: ", round(initial_tme_ratio, 1), "%\n", sep="")
  
  # 2. If not using original proportion and initial TME cells are insufficient, perform calibration allocation
  final_predictions <- initial_predictions  # Default to initial prediction
  if (!use_original_tme && initial_tme_count < target_tme_count) {
    cat("Initiating TME cell calibration allocation mechanism...\n")
    
    # Sort by TME feature score
    pred_features$predicted_celltype <- initial_predictions
    pred_features <- pred_features %>% arrange(desc(tme_feature_score))
    
    # Calculate number of TME cells needed to supplement
    need_tme_count <- target_tme_count - initial_tme_count
    
    # Select from non-TME samples with high TME feature scores
    non_tme_samples <- pred_features %>% 
      filter(!(predicted_celltype %in% TME_celltypes)) %>%
      head(need_tme_count)
    
    if (nrow(non_tme_samples) > 0) {
      # Allocate according to the proportion of TME types in the training set (including Dendritic cells)
      tme_types <- intersect(colnames(pred_probs), TME_celltypes)
      tme_ratios <- train_ratio[tme_types] / sum(train_ratio[tme_types])  # Normalize internal TME proportions
      
      # Allocate TME types to each sample needing conversion according to proportion
      for (i in 1:nrow(non_tme_samples)) {
        cell_id <- non_tme_samples$cell_id[i]
        new_type <- sample(tme_types, 1, prob = tme_ratios)
        
        # Update prediction probability matrix
        pred_probs[pred_features$cell_id == cell_id, ] <- 0
        pred_probs[pred_features$cell_id == cell_id, new_type] <- 1
      }
      
      final_predictions <- colnames(pred_probs)[max.col(pred_probs)]
      cat("Successfully allocated", nrow(non_tme_samples), "samples to TME cell types\n")
    } else {
      warning("Insufficient non-TME samples to convert to TME cell types")
    }
  }
  
  # 3. Final prediction results
  final_tme_count <- sum(final_predictions %in% TME_celltypes)  # Including Dendritic cells
  final_tme_ratio <- final_tme_count / nrow(pred_features) * 100
  cat("Adjusted actual TME cell proportion: ", round(final_tme_ratio, 1), "%\n", sep="")
  
  # 4. Integrate results (save TME_adjust parameter information)
  result_df <- dplyr::bind_cols(
    pred_features %>% dplyr::select(cell_id, x, y),
    predicted_celltype = final_predictions,
    as.data.frame(pred_probs)
  )
  
  # Add attributes for visualization
  attr(result_df, "TME_adjust") <- if (use_original_tme) "original" else TME_adjust
  attr(result_df, "step_size") <- step_size
  
  # Calculate and save proportion statistics
  ratio_df <- result_df %>%
    dplyr::count(predicted_celltype) %>%
    dplyr::mutate(
      proportion = (n / sum(n)) * 100,
      is_TME = predicted_celltype %in% TME_celltypes,  # Including Dendritic cells judgment
      train_proportion = as.numeric(train_ratio[as.character(predicted_celltype)]) * 100,
      cell_type = predicted_celltype
    ) %>%
    dplyr::select(cell_type, is_TME, count = n, predicted_proportion = proportion, 
                  train_proportion = train_proportion) %>%
    dplyr::arrange(dplyr::desc(predicted_proportion))
  
  # Calculate total cell count and add to results
  total_cells <- nrow(result_df)
  total_row <- data.frame(
    cell_type = "Total cells",
    is_TME = NA,  # TME label not applicable to total cells
    count = total_cells,
    predicted_proportion = 100,  # Total proportion is 100%
    train_proportion = NA  # No training set total proportion
  )
  
  # Add total cell row to the end of the data frame
  ratio_df <- dplyr::bind_rows(ratio_df, total_row)
  
  # Save results (prediction results and proportion statistics)
  save(result_df, ratio_df, train_ratio, file = paste0(output_prefix, "_predictions.RData"))
  write.csv(ratio_df, paste0(output_prefix, "_celltype_ratio.csv"), row.names = FALSE)
  cat("Prediction result data saved: ", paste0(output_prefix, "_predictions.RData"), "\n")
  cat("Prediction proportion statistics (including total cell count) saved: ", paste0(output_prefix, "_celltype_ratio.csv"), "\n")
  
  return(list(
    full_results = result_df,
    predicted_ratio = ratio_df,
    train_ratio = train_ratio,
    TME_adjust = if (use_original_tme) FALSE else TME_adjust,
    target_tme_ratio = target_tme_ratio,
    actual_tme_ratio = final_tme_ratio / 100,
    step_size = step_size,
    density_radius = density_radius
  ))
}


# Visualization function (optimized version, added cell type count bar plot)
generate_visualizations <- function(he_image_path,
                                    prediction_data_path,  # Path to the .RData file generated in the first section
                                    output_prefix = "he_prediction",
                                    point_size = 1.0,    # Scatter point size
                                    point_alpha = 0.8) {  # Scatter point transparency
  
  # Parameter validity check
  if (point_size <= 0 || point_size > 10) {
    stop("point_size must be a positive number and not greater than 10")
  }
  if (point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be between 0-1")
  }
  
  # Load prediction results
  if (!file.exists(prediction_data_path)) {
    stop("Prediction result file does not exist: ", prediction_data_path)
  }
  load(prediction_data_path)  # Load result_df, ratio_df, train_ratio, etc.
  
  # Extract key parameters from prediction results
  TME_adjust <- attr(result_df, "TME_adjust") %||% FALSE
  step_size <- attr(result_df, "step_size") %||% 10
  final_tme_ratio <- ratio_df %>%
    dplyr::filter(cell_type %in% TME_celltypes) %>%  # Including Dendritic cells
    dplyr::summarise(total = sum(predicted_proportion)) %>%
    dplyr::pull(total)
  
  # Generate visualization images
  # Read original HE image (not normalized)
  he_img_raw <- EBImage::readImage(he_image_path)
  
  # Process original image data for visualization (no normalization)
  img_df_raw <- reshape2::melt(he_img_raw) %>%
    dplyr::rename(x = Var1, y = Var2, channel = Var3, value = value) %>%
    dplyr::mutate(
      channel = factor(channel, labels = c("red", "green", "blue")),
      y = max(y) - y,  # Flip y-axis for correct display
      value = pmin(pmax(value, 0), 1)  # Only limit range, no stretching
    ) %>%
    tidyr::pivot_wider(names_from = channel, values_from = value)
  
  # Upper original HE image (using unnormalized data)
  p_original <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = img_df_raw %>% dplyr::sample_frac(0.5),
      ggplot2::aes(x = x, y = y, fill = rgb(red, green, blue))
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void() +
    ggplot2::ggtitle("Original HE Image") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16,
                                         margin = ggplot2::margin(b = 10))
    ) +
    ggplot2::coord_fixed()
  
  # Middle annotated image (including Dendritic cells color)
  he_img <- EBImage::readImage(he_image_path)
  img_df <- reshape2::melt(he_img) %>%
    dplyr::rename(x = Var1, y = Var2, channel = Var3, value = value) %>%
    dplyr::mutate(
      channel = factor(channel, labels = c("red", "green", "blue")),
      y = max(y) - y,  # Flip y-axis for correct display
      value = pmin(pmax(value, 0), 1)
    ) %>%
    tidyr::pivot_wider(names_from = channel, values_from = value)
  
  # Generate title (distinguish original proportion and adjusted proportion)
  if (TME_adjust == "original" || isFALSE(TME_adjust)) {
    plot_title <- paste0("Predicted Cell Types (Original TME Ratio, Step=", step_size, 
                         ", Actual TME=", round(final_tme_ratio, 1), "%)")
  } else {
    plot_title <- paste0("Predicted Cell Types (TME=", TME_adjust, 
                         ", Step=", step_size, ", Actual TME=", round(final_tme_ratio, 1), "%)")
  }
  
  # Create annotated HE image
  p_annotated_base <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = img_df %>% dplyr::sample_frac(0.5),
      ggplot2::aes(x = x, y = y, fill = rgb(red, green, blue))
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_point(
      data = result_df,
      ggplot2::aes(x = x, y = max(img_df$y) - y, color = predicted_celltype),
      size = point_size,
      alpha = point_alpha,
      shape = 16
    ) +
    ggplot2::scale_color_manual(values = celltype_colors) +  # Including Dendritic cells color
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",  # Remove legend as we will add bar plot
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20,
                                         margin = ggplot2::margin(b = 10, t = 10))
    ) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::coord_fixed()
  
  # Prepare bar plot data: count the number of each cell type
  bar_data <- result_df %>%
    dplyr::count(predicted_celltype) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::mutate(
      predicted_celltype = factor(predicted_celltype, levels = predicted_celltype),
      label = paste0(predicted_celltype, "\n(n=", n, ")")
    )
  
  # Create right bar plot
  p_bar <- ggplot2::ggplot(bar_data, ggplot2::aes(x = predicted_celltype, y = n, fill = predicted_celltype)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.1, size = 4) +
    ggplot2::scale_fill_manual(values = celltype_colors) +  # Use the same color scheme
    ggplot2::scale_y_reverse() +  # Reverse y-axis to match legend order
    ggplot2::coord_flip(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(l = 10, r = 10)
    ) +
    ggplot2::ylim(0, max(bar_data$n) * 1.1)  # Leave space for text
  
  # Combine annotated image and bar plot
  p_annotated <- p_annotated_base + p_bar +
    patchwork::plot_layout(widths = c(4, 1))  # Width ratio of main plot to bar plot is 4:1
  
  # Pie chart (including Dendritic cells)
  pie_data <- ratio_df %>%
    dplyr::filter(cell_type != "Total cells") %>%  # Exclude total cells row
    dplyr::select(cell_type, proportion = predicted_proportion) %>%
    dplyr::mutate(
      is_TME = ifelse(cell_type %in% TME_celltypes, "TME", "Non-TME"),  # Including Dendritic cells judgment
      label = ifelse(is_TME == "TME", 
                     paste0(cell_type, " (TME)\n", round(proportion, 1), "%"),
                     paste0(cell_type, "\n", round(proportion, 1), "%"))
    )
  
  p_pie <- ggplot2::ggplot(pie_data, 
                           ggplot2::aes(x = "", y = proportion, fill = cell_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::scale_fill_manual(values = celltype_colors) +  # Including Dendritic cells color
    ggplot2::geom_text(
      aes(label = label), 
      position = ggplot2::position_stack(vjust = 0.5), 
      size = 4,
      color = "white"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12,
                                         margin = ggplot2::margin(t = 10)),
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::ggtitle(
      paste0("Predicted Cell Type Proportions (Total: ", nrow(result_df), " cells)")
    )
  
  # Training vs prediction proportion comparison
  compare_data <- ratio_df %>%
    dplyr::filter(cell_type != "Total cells") %>%  # Exclude total cells row
    dplyr::select(cell_type, predicted = predicted_proportion, train = train_proportion) %>%
    tidyr::pivot_longer(cols = -cell_type, names_to = "type", values_to = "proportion")
  
  p_compare <- ggplot2::ggplot(compare_data,
                               ggplot2::aes(x = cell_type, y = proportion, 
                                            fill = type, group = type)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::scale_fill_manual(values = c("predicted" = "#F4A261", "train" = "#888888")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(title = "Predicted vs Training Proportions",
                  x = "Cell Type", y = "Proportion (%)", fill = "Type") +
    ggplot2::ylim(0, max(compare_data$proportion) * 1.2)
  
  # Combine images
  combined_plot <- (p_original / p_annotated / (p_pie + p_compare)) +
    patchwork::plot_layout(heights = c(3, 3, 3)) +
    patchwork::plot_annotation(
      title = if (TME_adjust == "original" || isFALSE(TME_adjust)) {
        paste0("HE Prediction (Original TME Ratio, Step=", step_size, ")")
      } else {
        paste0("HE Prediction (TME=", TME_adjust, ", Step=", step_size, ")")
      },
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 18))
    )
  
  # Save PDF image
  ggplot2::ggsave(
    paste0(output_prefix, "_report.pdf"),
    plot = combined_plot,
    width = 14, height = 18, dpi = 300, device = "pdf"
  )
  cat("Comprehensive comparison plot saved: ", paste0(output_prefix, "_report.pdf"), "\n")
  
  # Save annotated HE image (including bar plot)
  ggplot2::ggsave(
    paste0(output_prefix, "_annotated_he.tif"),
    plot = p_annotated,
    width = (dim(he_img)[1] / 300 * 1.5) * 1.25,  # Increase width to accommodate bar plot
    height = dim(he_img)[2] / 300 * 1.5,
    dpi = 300, device = "tiff", compression = "lzw"
  )
  cat("Annotated HE image saved: ", paste0(output_prefix, "_annotated_he.tif"), "\n")
  
  # Generate feature expression plot (including TME marker for Dendritic cells)
  feature_names <- c("red_mean", "green_mean", "blue_mean", 
                     "intensity_mean", "contrast", "correlation", 
                     "energy", "homogeneity",
                     "nucleus_area", "nucleus_circularity",
                     "cytoplasm_area", "cytoplasm_circularity",
                     "cell_density", "nearest_neighbor_dist",
                     "tme_feature_score")
  
  existing_features <- intersect(feature_names, colnames(result_df))
  
  if (length(existing_features) > 0) {
    # Merge feature data (requires features from original pred_features)
    pred_features_with_expr <- pred_features %>%
      dplyr::select(cell_id, dplyr::all_of(existing_features)) %>%
      dplyr::inner_join(result_df %>% dplyr::select(cell_id, predicted_celltype), by = "cell_id")
    
    expr_df <- pred_features_with_expr %>%
      dplyr::group_by(predicted_celltype) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(existing_features), mean), .groups = "drop") %>%
      dplyr::mutate(
        is_TME = predicted_celltype %in% TME_celltypes  # Including Dendritic cells judgment
      ) %>%
      tidyr::pivot_longer(cols = -c(predicted_celltype, is_TME), 
                          names_to = "feature", values_to = "value")
    
    p_expr <- ggplot2::ggplot(
      data = expr_df,
      ggplot2::aes(x = predicted_celltype, y = value, 
                   fill = predicted_celltype, color = is_TME)
    ) +
      ggplot2::geom_col(width = 0.7, size = 0.8) +
      ggplot2::facet_wrap(~feature, scales = "free_y", ncol = 3) +
      ggplot2::scale_fill_manual(values = celltype_colors) +  # Including Dendritic cells color
      ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.x = ggplot2::element_blank(),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16)
      ) +
      ggplot2::labs(title = "Mean Feature Values (TME Cells Marked in Red Border)")
    
    ggplot2::ggsave(
      paste0(output_prefix, "_feature_expression.pdf"),
      plot = p_expr,
      width = 11, height = 9, dpi = 300
    )
    cat("Feature expression plot (with TME markers) saved: ", paste0(output_prefix, "_feature_expression.pdf"), "\n")
  }
  
  return(list(
    plots = list(
      combined = combined_plot,
      annotated = p_annotated,
      pie = p_pie,
      report = p_compare,
      bar = p_bar
    )
  ))
}


# Example running code (please modify according to actual paths)
# 1. Generate prediction results (using original TME proportion from training set)
# predictions <- generate_predictions(
#   he_image_path = "path/to/your/he_image.tif",
#   model_path = "celltype_model.rds",
#   output_prefix = "he_prediction_original_tme",
#   TME_adjust = FALSE,  # Use original TME proportion, no calibration
#   step_size = 10
# )

# 2. Generate prediction results (using TME adjustment)
# predictions <- generate_predictions(
#   he_image_path = "path/to/your/he_image.tif",
#   model_path = "celltype_model.rds",
#   output_prefix = "he_prediction_adjusted_tme",
#   TME_adjust = 0.5,  # Adjust TME proportion to approximately 27.5%
#   step_size = 10
# )

# 3. Generate visualization results
# visualizations <- generate_visualizations(
#   he_image_path = "path/to/your/he_image.tif",
#   prediction_data_path = "he_prediction_original_tme_predictions.RData",  # Or the adjusted result file
#   output_prefix = "he_prediction_original_tme",
#   point_size = 1.0,
#   point_alpha = 0.8
# )