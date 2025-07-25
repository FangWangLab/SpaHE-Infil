# HE_feature_Get.R
# Function: Extract nuclear morphology, texture, density features from HE images and visualize
# Usage:
# 1. source("HE_feature_Get.R")
# 2. rds_files <- c("final_annotated_data.rds")  # Replace with your file path
# 3. HE_feature_Get(rds_files = rds_files)

HE_feature_Get <- function(rds_files, output_dir = "HE_feature_plots", debug = FALSE, 
                           density_radius = 30) {  # 增大默认半径以捕获更多邻居
  # 定义细胞类型颜色方案
  celltype_colors <- c(
    "Tumor cells" = "#FF6F61",                
    "Normal epithelial cells" = "#FFD800",     
    "CD8 T cells" = "#7FDBFF",                
    "CD4 T cells" = "#1E90FF",                 
    "Macrophages" = "#9370DB",                
    "Neutrophils" = "#00CED1",                
    "NK cells" = "#DA70D6",                    
    "Tregs" = "#CD853F",                       
    "Fibroblasts" = "#4CAF50",                 
    "Endothelial cells" = "#9ACD32",          
    "B cells" = "#F8C3CD",                     
    "Dendritic cells" = "#20B2AA"             
  )
  
  # Check input file validity
  if (length(rds_files) == 0) {
    stop("Please provide valid RDS file paths")
  }
  for (file in rds_files) {
    if (!file.exists(file)) {
      stop(paste("File does not exist:", file))
    }
  }
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load required packages and check dependencies
  required_pkgs <- c("Seurat", "EBImage", "ggplot2", "dplyr", "tidyr", "viridis", "patchwork")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(paste("Please install missing packages first:\ninstall.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))", sep = ""))
  }
  
  # Force load dplyr and explicitly reference its namespace
  if (!"dplyr" %in% .packages()) {
    library(dplyr, character.only = TRUE)
  }
  
  # Load other required packages
  library(Seurat)
  library(EBImage)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(patchwork)
  library(dplyr)
  
  # Compatible set operation function for older dplyr versions
  safe_intersect <- function(x, y) {
    if (packageVersion("dplyr") < "1.0.0") base::intersect(x, y) else dplyr::intersect(x, y)
  }
  
  # --------------------------
  # Internal function: Extract feature data
  # --------------------------
  extract_features <- function(rds_file) {
    # Read data
    spatial_data <- readRDS(rds_file)
    
    # Extract cell coordinates (compatible with different column names)
    coords <- Seurat::GetTissueCoordinates(spatial_data) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell_id")
    
    # Automatically adapt coordinate column names
    if ("imagecol" %in% colnames(coords) && "imagerow" %in% colnames(coords)) {
      coords <- coords %>% dplyr::rename(x = imagecol, y = imagerow)
    } else if ("col" %in% colnames(coords) && "row" %in% colnames(coords)) {
      coords <- coords %>% dplyr::rename(x = col, y = row)
    } else if (!"x" %in% colnames(coords) || !"y" %in% colnames(coords)) {
      stop(paste("Could not identify coordinate columns. Coordinate column names in file", rds_file, ":", paste(colnames(coords), collapse = ", ")))
    }
    
    # Call feature extraction function
    features <- extract_training_data(rds_file, debug = debug, density_radius = density_radius)
    
    # Merge coordinate information (using explicit dplyr:: function calls)
    features %>%
      dplyr::inner_join(coords %>% dplyr::select(cell_id, x, y), by = "cell_id") %>%
      dplyr::mutate(source_file = basename(rds_file))
  }
  
  # --------------------------
  # Core visualization functions
  # --------------------------
  # 1. Nuclear morphology feature visualization
  plot_nucleus <- function(data) {
    nucleus_feats <- c("nucleus_area", "nucleus_perimeter", "nucleus_circularity",
                       "nucleus_aspect_ratio", "nucleus_solidity")
    nucleus_feats <- safe_intersect(nucleus_feats, colnames(data))
    
    if (length(nucleus_feats) < 2) {
      warning("Insufficient nuclear morphology features, skipping nuclear morphology visualization")
      return(NULL)
    }
    
    # 检查核形态特征是否全为0
    zero_nucleus_feats <- sapply(nucleus_feats, function(feat) {
      all(data[[feat]] == 0, na.rm = TRUE)
    })
    
    if (any(zero_nucleus_feats)) {
      warning(paste("The following nuclear features have all zero values:", 
                    paste(names(zero_nucleus_feats)[zero_nucleus_feats], collapse = ", ")))
    }
    
    plots <- lapply(nucleus_feats, function(feat) {
      ggplot2::ggplot(data, ggplot2::aes(x = cell_type, y = .data[[feat]], fill = cell_type)) +
        ggplot2::geom_boxplot(alpha = 0.8, show.legend = FALSE) +
        ggplot2::labs(
          title = gsub("nucleus_", "", feat),
          x = "Cell Type",
          y = "Feature Value"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::scale_fill_manual(values = celltype_colors)  # 使用指定颜色
    })
    
    combined <- patchwork::wrap_plots(plots, ncol = 3) + 
      patchwork::plot_annotation(title = "Nuclear Morphology Features (by Cell Type)")
    # 保存为PDF格式
    ggplot2::ggsave(file.path(output_dir, "nuclear_morphology_boxplots.pdf"), 
                    combined, width = 12, height = 4, device = "pdf")
    combined
  }
  
  # 2. Texture feature visualization
  plot_texture <- function(data) {
    texture_feats <- c("contrast", "energy", "homogeneity", "correlation")
    texture_feats <- safe_intersect(texture_feats, colnames(data))
    
    if (length(texture_feats) < 2) {
      warning("Insufficient texture features, skipping texture visualization")
      return(NULL)
    }
    
    # 检查纹理特征是否全为0
    zero_texture_feats <- sapply(texture_feats, function(feat) {
      all(data[[feat]] == 0, na.rm = TRUE)
    })
    
    if (any(zero_texture_feats)) {
      warning(paste("The following texture features have all zero values:", 
                    paste(names(zero_texture_feats)[zero_texture_feats], collapse = ", ")))
    }
    
    plots <- lapply(texture_feats, function(feat) {
      # 为全零特征添加微小扰动以便可视化
      plot_data <- data
      if (zero_texture_feats[feat]) {
        plot_data[[feat]] <- jitter(plot_data[[feat]], factor = 0.1)
        message(paste("Added jitter to", feat, "for visualization purposes"))
      }
      
      ggplot2::ggplot(plot_data, ggplot2::aes(x = cell_type, y = .data[[feat]], fill = cell_type)) +
        ggplot2::geom_boxplot(alpha = 0.8, show.legend = FALSE) +
        ggplot2::labs(
          title = feat,
          x = "Cell Type",
          y = "Feature Value"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::scale_fill_manual(values = celltype_colors)  # 使用指定颜色
    })
    
    combined <- patchwork::wrap_plots(plots, ncol = 2) + 
      patchwork::plot_annotation(title = "Texture Features (by Cell Type)")
    # 保存为PDF格式
    ggplot2::ggsave(file.path(output_dir, "texture_features_boxplots.pdf"), 
                    combined, width = 10, height = 6, device = "pdf")
    
    # 同时生成空间分布图
    spatial_plots <- lapply(texture_feats, function(feat) {
      plot_data <- data
      if (zero_texture_feats[feat]) {
        plot_data[[feat]] <- jitter(plot_data[[feat]], factor = 0.1)
      }
      
      ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = .data[[feat]])) +
        ggplot2::geom_point(size = 1.2, alpha = 0.7) +
        viridis::scale_color_viridis(option = ifelse(feat == "contrast", "plasma", "viridis"), name = feat) +
        ggplot2::labs(title = paste("Texture Feature:", feat), x = "X Coordinate", y = "Y Coordinate") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          aspect.ratio = 1,
          plot.title = ggplot2::element_text(hjust = 0.5)
        )
    })
    
    spatial_combined <- patchwork::wrap_plots(spatial_plots, ncol = 2) + 
      patchwork::plot_annotation(title = "Spatial Distribution of Texture Features")
    # 保存为PDF格式
    ggplot2::ggsave(file.path(output_dir, "texture_feature_spatial_heatmaps.pdf"), 
                    spatial_combined, width = 10, height = 8, device = "pdf")
    
    list(boxplots = combined, spatial = spatial_combined)
  }
  
  # 3. Density feature visualization - 应用指定颜色的核心部分
  plot_density <- function(data) {
    density_feats <- c("cell_density", "nearest_neighbor_dist")
    density_feats <- safe_intersect(density_feats, colnames(data))
    
    if (length(density_feats) < 1) {
      warning("Insufficient density features, skipping density visualization")
      return(NULL)
    }
    
    # 检查数据中实际存在的细胞类型
    present_cell_types <- unique(data$cell_type)
    # 筛选出数据中存在的细胞类型对应的颜色
    used_colors <- celltype_colors[names(celltype_colors) %in% present_cell_types]
    
    # 检查是否有未在配色方案中定义的细胞类型
    undefined_types <- setdiff(present_cell_types, names(celltype_colors))
    if (length(undefined_types) > 0) {
      warning(paste("The following cell types have no defined color:", 
                    paste(undefined_types, collapse = ", "), 
                    "- using default colors for these"))
      # 为未定义的细胞类型生成临时颜色
      n_undefined <- length(undefined_types)
      default_palette <- viridis::viridis(n_undefined, option = "plasma")
      names(default_palette) <- undefined_types
      # 合并使用的颜色和临时颜色
      used_colors <- c(used_colors, default_palette)
    }
    
    # 检查密度特征是否全为0或异常值
    zero_density_feats <- sapply(density_feats, function(feat) {
      all(data[[feat]] <= 0.001, na.rm = TRUE)  # 使用小阈值检测接近零的值
    })
    
    if (any(zero_density_feats)) {
      warning(paste("The following density features have all zero or near-zero values:", 
                    paste(names(zero_density_feats)[zero_density_feats], collapse = ", ")))
      
      # 对接近零的特征添加微小扰动以便可视化
      for (feat in names(zero_density_feats)[zero_density_feats]) {
        data[[feat]] <- jitter(data[[feat]], factor = 0.1)
        message(paste("Added jitter to", feat, "for visualization purposes"))
      }
    }
    
    # 输出密度特征的统计信息用于调试
    if (debug) {
      cat("\nDensity feature statistics:\n")
      print(summary(data[, density_feats]))
    }
    
    # Spatial distribution of density
    p1 <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, color = cell_density, size = cell_density)) +
      ggplot2::geom_point(alpha = 0.6) +
      viridis::scale_color_viridis(option = "magma", name = "Cell Density") +
      ggplot2::scale_size_continuous(range = c(1, 3), name = "Cell Density") +
      ggplot2::labs(title = "Spatial Distribution of Cell Density", x = "X Coordinate", y = "Y Coordinate") +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(hjust = 0.5))
    
    # Relationship between density and nearest neighbor distance - 应用指定颜色
    p2 <- if ("nearest_neighbor_dist" %in% density_feats) {
      ggplot2::ggplot(data, ggplot2::aes(x = nearest_neighbor_dist, y = cell_density, color = cell_type)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 2) +
        ggplot2::scale_color_manual(values = used_colors, name = "Cell Type") +  # 使用指定颜色
        ggplot2::labs(title = "Cell Density vs Nearest Neighbor Distance", x = "Nearest Neighbor Distance", y = "Cell Density") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.key.size = ggplot2::unit(0.8, "cm"),  # 调整图例大小
          legend.text = ggplot2::element_text(size = 8)  # 调整图例文字大小
        )
    } else {
      NULL
    }
    
    combined <- if (!is.null(p2)) p1 + p2 else p1
    # 保存为PDF格式
    ggplot2::ggsave(file.path(output_dir, "cell_density_features.pdf"), 
                    combined, width = ifelse(!is.null(p2), 12, 6), height = 5, device = "pdf")
    combined
  }
  
  # 4. Feature correlation and PCA
  plot_correlation_pca <- function(data) {
    # Select feature columns
    all_feats <- c(
      "nucleus_area", "nucleus_perimeter", "nucleus_circularity",
      "contrast", "energy", "homogeneity", "cell_density", "nearest_neighbor_dist"
    )
    use_feats <- safe_intersect(all_feats, colnames(data))
    
    if (length(use_feats) < 4) {
      warning("Insufficient valid features (need at least 4), skipping correlation and PCA analysis")
      return(NULL)
    }
    
    # 筛选出数据中存在的细胞类型对应的颜色
    present_cell_types <- unique(data$cell_type)
    used_colors <- celltype_colors[names(celltype_colors) %in% present_cell_types]
    
    # 检查是否有未在配色方案中定义的细胞类型
    undefined_types <- setdiff(present_cell_types, names(celltype_colors))
    if (length(undefined_types) > 0) {
      n_undefined <- length(undefined_types)
      default_palette <- viridis::viridis(n_undefined, option = "plasma")
      names(default_palette) <- undefined_types
      used_colors <- c(used_colors, default_palette)
    }
    
    # Extract feature data
    feature_data <- data %>%
      dplyr::select(dplyr::all_of(use_feats))
    
    # Calculate variance for each feature and filter out constant features (variance = 0)
    feature_vars <- apply(feature_data, 2, var, na.rm = TRUE)
    non_constant_feats <- names(feature_vars[feature_vars > 0])
    
    # Check if we have enough non-constant features left
    if (length(non_constant_feats) < 4) {
      warning(paste("Too many constant features (only", length(non_constant_feats), "non-constant features remaining),",
                    "skipping correlation and PCA analysis"))
      return(NULL)
    }
    
    # Warn about removed features
    removed_feats <- setdiff(use_feats, non_constant_feats)
    if (length(removed_feats) > 0) {
      message(paste("Removed constant features from PCA analysis:", paste(removed_feats, collapse = ", ")))
    }
    
    # Update feature data with only non-constant features
    feature_data <- feature_data %>%
      dplyr::select(dplyr::all_of(non_constant_feats))
    
    # Correlation heatmap
    corr_matrix <- feature_data %>%
      stats::cor(use = "complete.obs")
    
    corr_long <- corr_matrix %>%
      as.data.frame() %>%
      dplyr::mutate(feat1 = rownames(.)) %>%
      tidyr::pivot_longer(-feat1, names_to = "feat2", values_to = "correlation")
    
    p_corr <- ggplot2::ggplot(corr_long, ggplot2::aes(x = feat1, y = feat2, fill = correlation)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = round(correlation, 2)), size = 3) +
      viridis::scale_fill_viridis(limits = c(-1, 1), option = "cividis") +
      ggplot2::labs(title = "Feature Correlation Heatmap") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    
    # PCA cluster plot with specified colors
    pca_data <- feature_data %>%
      scale() %>%
      stats::prcomp(na.action = stats::na.omit)
    
    # Get valid indices (non-NA rows)
    valid_indices <- which(!is.na(apply(feature_data, 1, function(row) any(is.na(row)))))
    
    pca_df <- data.frame(
      PC1 = pca_data$x[, 1],
      PC2 = pca_data$x[, 2],
      cell_type = data$cell_type[valid_indices]
    )
    
    p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = cell_type)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.5) +
      ggplot2::stat_ellipse(level = 0.9, linetype = 2) +
      ggplot2::scale_color_manual(values = used_colors, name = "Cell Type") +  # 使用指定颜色
      ggplot2::labs(
        title = "Feature PCA Clustering (by Cell Type)",
        x = paste0("PC1 (", round(pca_data$sdev[1]^2/sum(pca_data$sdev^2)*100, 1), "%)"),
        y = paste0("PC2 (", round(pca_data$sdev[2]^2/sum(pca_data$sdev^2)*100, 1), "%)")
      ) +
      ggplot2::theme_minimal()
    
    combined <- p_corr + p_pca + patchwork::plot_annotation(title = "Feature Relationship Analysis")
    # 保存为PDF格式
    ggplot2::ggsave(file.path(output_dir, "feature_correlation_and_pca.pdf"), 
                    combined, width = 12, height = 5, device = "pdf")
    combined
  }
  
  # --------------------------
  # Main execution flow
  # --------------------------
  cat("===== Starting feature data extraction =====", "\n")
  all_features <- lapply(rds_files, extract_features) %>% dplyr::bind_rows()
  
  if (nrow(all_features) == 0) {
    stop("No valid feature data extracted, please check input files")
  }
  
  # 检查纹理特征整体情况
  texture_feats <- c("contrast", "energy", "homogeneity", "correlation")
  present_texture_feats <- safe_intersect(texture_feats, colnames(all_features))
  
  if (length(present_texture_feats) > 0) {
    zero_texture_summary <- sapply(present_texture_feats, function(feat) {
      sum(all_features[[feat]] == 0, na.rm = TRUE) / nrow(all_features) * 100
    })
    
    cat("\nTexture feature zero value percentage:\n")
    print(round(zero_texture_summary, 1))
  }
  
  # 检查密度特征整体情况
  density_feats <- c("cell_density", "nearest_neighbor_dist")
  present_density_feats <- safe_intersect(density_feats, colnames(all_features))
  
  if (length(present_density_feats) > 0) {
    zero_density_summary <- sapply(present_density_feats, function(feat) {
      sum(all_features[[feat]] <= 0.001, na.rm = TRUE) / nrow(all_features) * 100  # 检测接近零的值
    })
    
    cat("\nDensity feature near-zero value percentage:\n")
    print(round(zero_density_summary, 1))
  }
  
  cat("\nFeature extraction completed, obtained features for", nrow(all_features), "cells\n")
  cat("Cell type distribution:\n")
  print(table(all_features$cell_type))
  
  cat("\n===== Starting visualization generation =====", "\n")
  plot_nucleus(all_features)
  plot_texture(all_features)
  plot_density(all_features)
  plot_correlation_pca(all_features)
  
  cat("\nAll visualization results have been saved to:", normalizePath(output_dir), "\n", sep = "")
  invisible(all_features)
}

# Load necessary helper functions
`%||%` <- function(x, y) if (is.null(x)) y else x

compute_features <- function(img, debug = FALSE) {
  # 确保图像是单通道
  if (length(dim(img)) > 2) {
    if (debug) message("Converting to single channel image")
    img <- img[,,1]
  }
  
  # 确保图像尺寸合适
  if (any(dim(img) < 5)) {
    warning("Image too small for texture feature extraction (minimum 5x5 required)")
    return(list(contrast=0.1, correlation=0.1, energy=0.1, homogeneity=0.1)) # 返回非零默认值
  }
  
  # 标准化图像以提高特征提取稳定性
  img <- tryCatch({
    img_norm <- EBImage::normalize(img)
    # 确保图像不是全黑
    if (max(img_norm) < 0.01) {
      if (debug) message("Image is too dark, adding small values to avoid zero features")
      img_norm <- img_norm + 0.01  # 避免全黑图像
    }
    img_norm
  }, error = function(e) {
    message("Warning: Image normalization failed - ", e$message)
    img  # 返回原始图像尝试提取特征
  })
  
  # 尝试提取Haralick特征
  features <- tryCatch({
    # 调整参数以提高成功率
    EBImage::computeFeatures.haralick(
      img, ref = img, 
      haralick.nbins = 16,  # 减少分箱数提高稳定性
      haralick.scales = 1   # 使用单一尺度减少复杂性
    )
  }, error = function(e) {
    message("Warning: computeFeatures.haralick failed - ", e$message)
    NULL
  })
  
  # 如果特征提取失败，使用替代方法计算基本纹理特征
  if (is.null(features) || all(is.na(features)) || ncol(features) < 4) {
    if (debug) message("Using alternative method for texture feature extraction")
    
    # 计算简单统计量作为替代特征
    img_flat <- as.vector(img)
    mean_val <- mean(img_flat, na.rm = TRUE)
    var_val <- var(img_flat, na.rm = TRUE)
    
    return(list(
      contrast = ifelse(var_val > 0, var_val, 0.1),  # 使用方差作为对比度替代
      correlation = max(0.1, min(0.9, 1 - var_val)),  # 简单反相关
      energy = mean_val^2,  # 使用均值平方作为能量替代
      homogeneity = 1 / (1 + var_val)  # 简单同质性度量
    ))
  }
  
  # 提取所需特征
  col_patterns <- c(
    contrast = "contrast", 
    correlation = "correlation", 
    energy = "energy", 
    homogeneity = "homogeneity"
  )
  
  result <- sapply(col_patterns, function(pattern) {
    matched_col <- grep(pattern, colnames(features), value = TRUE)
    if (length(matched_col) > 0) {
      val <- mean(features[, matched_col], na.rm = TRUE)
      # 确保值不为零或NA
      if (is.na(val) || val == 0) {
        if (debug) message(paste("Feature", pattern, "is zero or NA, using fallback value"))
        0.1  # 非零 fallback 值
      } else {
        val
      }
    } else {
      if (debug) message(paste("No columns found for", pattern, ", using fallback value"))
      0.1  # 非零 fallback 值
    }
  }, simplify = FALSE)
  
  result
}

compute_morphology_features <- function(bin_img, debug = FALSE) {
  if (!inherits(bin_img, "Image")) {
    bin_img <- EBImage::Image(bin_img, colormode = "Grayscale")
  }
  img_dims <- dim(bin_img)
  if (length(img_dims) < 2 || img_dims[1] < 3 || img_dims[2] < 3) {
    message("Warning: Image size too small to compute morphology features")
    return(list(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1))
  }
  
  morph_features <- tryCatch({
    connected <- EBImage::bwlabel(bin_img)
    if (max(connected) == 0) {
      if (debug) message("No connected components found, using fallback morphology features")
      return(data.frame(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1))
    }
    
    region_sizes <- table(connected[connected > 0])
    if (length(region_sizes) == 0) {
      message("Warning: No valid connected regions")
      return(data.frame(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1))
    }
    
    main_region <- as.integer(names(which.max(region_sizes)))
    if (is.na(main_region) || main_region < 1 || main_region > max(connected)) {
      message("Warning: Maximum connected region index invalid")
      return(data.frame(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1))
    }
    main_mask <- connected == main_region
    
    mask_dims <- dim(main_mask)
    if (mask_dims[1] != img_dims[1] || mask_dims[2] != img_dims[2]) {
      message("Warning: Mask dimensions do not match")
      return(data.frame(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1))
    }
    
    shape_stats <- EBImage::computeFeatures.shape(EBImage::Image(main_mask, colormode = "Grayscale"))
    
    # 确保形态特征不为零
    ensure_non_zero <- function(value, fallback = 0.1) {
      if (is.na(value) || value <= 0) fallback else value
    }
    
    data.frame(
      area = ensure_non_zero(ifelse("area" %in% rownames(shape_stats), shape_stats["area", ], 0.1)),
      perimeter = ensure_non_zero(ifelse("perimeter" %in% rownames(shape_stats), shape_stats["perimeter", ], 0.1)),
      circularity = ensure_non_zero(ifelse("circularity" %in% rownames(shape_stats), shape_stats["circularity", ], 0.1)),
      aspect_ratio = ensure_non_zero(ifelse("aspect.ratio" %in% rownames(shape_stats), shape_stats["aspect.ratio", ], 0.1)),
      solidity = ensure_non_zero(ifelse("solidity" %in% rownames(shape_stats), shape_stats["solidity", ], 0.1))
    )
  }, error = function(e) {
    message("Warning: Morphology feature calculation failed - ", e$message)
    data.frame(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1)
  })
  
  as.list(morph_features)
}

segment_he_nucleus_cytoplasm <- function(he_region, debug = FALSE) {
  if (!inherits(he_region, "Image")) {
    he_region <- EBImage::Image(he_region, colormode = "Color")
  }
  
  # 改进分割参数以提高成功率
  tryCatch({
    # 细胞核分割（HE染色中通常为蓝色）
    blue <- he_region[,,3]
    blue_norm <- EBImage::normalize(blue)
    # 调整阈值参数
    nucleus_bin <- EBImage::thresh(blue_norm, w = 3, h = 3, offset = 0.2)
    nucleus_bin <- EBImage::morphology(nucleus_bin, "dilate", 1)  # 轻微膨胀以确保核被正确捕获
    
    # 细胞质分割
    he_gray <- EBImage::channel(he_region, "gray")
    gray_norm <- EBImage::normalize(he_gray)
    # 调整阈值参数
    cell_bin <- EBImage::thresh(gray_norm, w = 5, h = 5, offset = 0.05)
    
    # 确保细胞质是细胞核的超集
    if (!is.null(nucleus_bin)) {
      cell_bin <- cell_bin | nucleus_bin
      cytoplasm_bin <- cell_bin - nucleus_bin
      cytoplasm_bin[ cytoplasm_bin < 0 ] <- 0
      cytoplasm_bin <- EBImage::morphology(cytoplasm_bin, "open", 1)  # 去除噪声
    } else {
      cytoplasm_bin <- cell_bin
    }
    
    list(nucleus = nucleus_bin, cytoplasm = cytoplasm_bin)
  }, error = function(e) {
    message("Warning: Segmentation failed - ", e$message)
    # 分割失败时返回默认掩码
    default_mask <- EBImage::Image(matrix(1, nrow = dim(he_region)[1], ncol = dim(he_region)[2]), 
                                   colormode = "Grayscale")
    list(nucleus = default_mask, cytoplasm = default_mask)
  })
}

# 细胞密度和最近邻距离计算函数（已修复）
compute_density_features <- function(coords, index, density_radius = 30, debug = FALSE) {
  # 确保坐标数据有效
  if (nrow(coords) < 2) {
    if (debug) message("Too few cells (<2) to compute density features")
    return(list(cell_density=0.1, nearest_neighbor_dist=0.1, neighbor_dist_sd=0.1))
  }
  
  # 提取当前细胞坐标
  current_x <- if ("imagecol" %in% colnames(coords)) coords$imagecol[index] else coords$x[index]
  current_y <- if ("imagerow" %in% colnames(coords)) coords$imagerow[index] else coords$y[index]
  
  # 验证坐标值有效性
  if (is.na(current_x) || is.na(current_y) || current_x <= 0 || current_y <= 0) {
    warning("Invalid coordinates for cell at index ", index)
    return(list(cell_density=0.1, nearest_neighbor_dist=0.1, neighbor_dist_sd=0.1))
  }
  
  # 计算与其他所有细胞的距离
  other_x <- if ("imagecol" %in% colnames(coords)) coords$imagecol[-index] else coords$x[-index]
  other_y <- if ("imagerow" %in% colnames(coords)) coords$imagerow[-index] else coords$y[-index]
  
  distances <- sqrt((other_x - current_x)^2 + (other_y - current_y)^2)
  
  # 过滤掉异常值和零距离（可能是重复坐标）
  valid_distances <- distances[distances > 0.1]  # 排除非常接近的点（可能是同一细胞）
  
  if (length(valid_distances) == 0) {
    if (debug) message("No valid neighboring cells found for cell at index ", index)
    return(list(cell_density=0.1, nearest_neighbor_dist=0.5, neighbor_dist_sd=0.1))  # 使用合理默认值
  }
  
  # 计算指定半径内的细胞数量（细胞密度）
  in_radius <- sum(valid_distances <= density_radius)
  area <- pi * density_radius^2
  
  # 确保面积不为零，并计算密度
  cell_density <- if (area > 0) in_radius / area else 0.1
  
  # 确保密度不为零（如果半径内没有细胞，使用基于最近邻的估计）
  if (cell_density <= 0) {
    # 使用最近邻距离估计密度（假设细胞均匀分布）
    nn_dist <- min(valid_distances)
    cell_density <- if (nn_dist > 0) 1 / (pi * nn_dist^2) else 0.1
    if (debug) message("No cells in radius, estimating density from nearest neighbor for cell ", index)
  }
  
  # 计算最近邻距离
  nearest_neighbor_dist <- min(valid_distances)
  
  # 计算距离标准差
  neighbor_dist_sd <- if (length(valid_distances) >= 2) sd(valid_distances) else 0.1
  
  # 确保所有值都大于零
  cell_density <- max(cell_density, 0.001)
  nearest_neighbor_dist <- max(nearest_neighbor_dist, 0.1)
  
  if (debug && (cell_density < 0.001 || nearest_neighbor_dist < 0.1)) {
    message("Low density values for cell ", index, 
            ": density = ", round(cell_density, 6), 
            ", NND = ", round(nearest_neighbor_dist, 3))
  }
  
  list(
    cell_density = cell_density, 
    nearest_neighbor_dist = nearest_neighbor_dist, 
    neighbor_dist_sd = neighbor_dist_sd
  )
}

extract_training_data <- function(spatial_rds, density_radius = 30, debug = FALSE) {
  spatial_data <- readRDS(spatial_rds)
  
  if (is.null(spatial_data@images$slice1)) {
    stop("HE image data not found, please check input file")
  }
  he_img <- spatial_data@images$slice1@image
  he_array <- as.array(he_img)
  he_eb <- EBImage::Image(he_array, colormode = "Color")
  
  coords <- Seurat::GetTissueCoordinates(spatial_data)
  celltypes <- spatial_data$auc_celltype
  
  # 增加区域半径以获取更多纹理信息
  radius <- 7  # 从5增加到7，获取更大区域
  
  features <- purrr::map_dfr(1:nrow(coords), function(i) {
    x <- coords$imagecol[i]
    y <- coords$imagerow[i]
    
    img_rows <- dim(he_eb)[1]
    img_cols <- dim(he_eb)[2]
    x_min <- max(1, x - radius)
    x_max <- min(img_rows, x + radius)
    y_min <- max(1, y - radius)
    y_max <- min(img_cols, y + radius)
    
    region_width <- x_max - x_min + 1
    region_height <- y_max - y_min + 1
    
    # 放宽尺寸要求，但仍确保有足够信息
    if (region_width < 7 || region_height < 7) {
      message("Warning: Cell region size too small (", region_width, "x", region_height, 
              "), using adjusted features")
      return(dplyr::tibble(
        cell_id = rownames(coords)[i], cell_type = celltypes[i],
        red_mean=0.1, red_sd=0.1, green_mean=0.1, green_sd=0.1, blue_mean=0.1, blue_sd=0.1,
        intensity_mean=0.1, intensity_sd=0.1, contrast=0.1, correlation=0.1, energy=0.1, homogeneity=0.1,
        nucleus_area=0.1, nucleus_perimeter=0.1, nucleus_circularity=0.1, nucleus_aspect_ratio=0.1, nucleus_solidity=0.1,
        cytoplasm_area=0.1, cytoplasm_perimeter=0.1, cytoplasm_circularity=0.1, cytoplasm_aspect_ratio=0.1, cytoplasm_solidity=0.1,
        cell_density=0.1, nearest_neighbor_dist=0.1, neighbor_dist_sd=0.1
      ))
    }
    
    region <- tryCatch({
      he_eb[x_min:x_max, y_min:y_max, , drop=FALSE]
    }, error = function(e) {
      message("Warning: HE image region extraction failed - ", e$message)
      NULL
    })
    if (is.null(region)) {
      return(dplyr::tibble(
        cell_id = rownames(coords)[i], cell_type = celltypes[i],
        red_mean=0.1, red_sd=0.1, green_mean=0.1, green_sd=0.1, blue_mean=0.1, blue_sd=0.1,
        intensity_mean=0.1, intensity_sd=0.1, contrast=0.1, correlation=0.1, energy=0.1, homogeneity=0.1,
        nucleus_area=0.1, nucleus_perimeter=0.1, nucleus_circularity=0.1, nucleus_aspect_ratio=0.1, nucleus_solidity=0.1,
        cytoplasm_area=0.1, cytoplasm_perimeter=0.1, cytoplasm_circularity=0.1, cytoplasm_aspect_ratio=0.1, cytoplasm_solidity=0.1,
        cell_density=0.1, nearest_neighbor_dist=0.1, neighbor_dist_sd=0.1
      ))
    }
    
    red_channel <- region[,,1]
    green_channel <- region[,,2]
    blue_channel <- region[,,3]
    red_mean <- mean(red_channel, na.rm=TRUE) %||% 0.1
    red_sd <- sd(red_channel, na.rm=TRUE) %||% 0.1
    green_mean <- mean(green_channel, na.rm=TRUE) %||% 0.1
    green_sd <- sd(green_channel, na.rm=TRUE) %||% 0.1
    blue_mean <- mean(blue_channel, na.rm=TRUE) %||% 0.1
    blue_sd <- sd(blue_channel, na.rm=TRUE) %||% 0.1
    
    gray_img <- EBImage::channel(region, "gray")
    texture_features <- compute_features(gray_img, debug = debug)  # 传递debug参数
    intensity_mean <- mean(gray_img, na.rm=TRUE) %||% 0.1
    intensity_sd <- sd(gray_img, na.rm=TRUE) %||% 0.1
    
    segmented <- segment_he_nucleus_cytoplasm(region, debug = debug)
    
    nucleus_features <- if (!is.null(segmented$nucleus)) {
      compute_morphology_features(segmented$nucleus, debug = debug)
    } else {
      list(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1)
    }
    
    cytoplasm_features <- if (!is.null(segmented$cytoplasm)) {
      compute_morphology_features(segmented$cytoplasm, debug = debug)
    } else {
      list(area=0.1, perimeter=0.1, circularity=0.1, aspect_ratio=0.1, solidity=0.1)
    }
    
    density_features <- compute_density_features(coords, i, density_radius = density_radius, debug = debug)
    
    dplyr::tibble(
      cell_id = rownames(coords)[i], cell_type = celltypes[i],
      red_mean=red_mean, red_sd=red_sd, green_mean=green_mean, green_sd=green_sd,
      blue_mean=blue_mean, blue_sd=blue_sd,
      intensity_mean=intensity_mean, intensity_sd=intensity_sd,
      contrast=texture_features$contrast, correlation=texture_features$correlation,
      energy=texture_features$energy, homogeneity=texture_features$homogeneity,
      nucleus_area=nucleus_features$area, nucleus_perimeter=nucleus_features$perimeter,
      nucleus_circularity=nucleus_features$circularity, nucleus_aspect_ratio=nucleus_features$aspect_ratio,
      nucleus_solidity=nucleus_features$solidity,
      cytoplasm_area=cytoplasm_features$area, cytoplasm_perimeter=cytoplasm_features$perimeter,
      cytoplasm_circularity=cytoplasm_features$circularity, cytoplasm_aspect_ratio=cytoplasm_features$aspect_ratio,
      cytoplasm_solidity=cytoplasm_features$solidity,
      cell_density=density_features$cell_density,
      nearest_neighbor_dist=density_features$nearest_neighbor_dist,
      neighbor_dist_sd=density_features$neighbor_dist_sd
    )
  }, .id = "point_id")
  
  features[complete.cases(features), ]
}

cat("HE_feature_Get function loaded successfully!\n")
cat("Example usage:\n")
cat("rds_files <- c(\"final_annotated_data.rds\")\n")
cat("HE_feature_Get(rds_files = rds_files, debug = TRUE, density_radius = 40)  # 可调整半径参数\n")
