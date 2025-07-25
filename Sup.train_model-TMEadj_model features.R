# 02.train_model-TMEadj-features.R
# Function: Read trained model file and visualize cross-validation results, confusion matrix, feature importance, etc.
# Usage: After source("02.train_model-TMEadj-features.R"), run model_features(rds_files = rds_files)

model_features <- function(rds_files, model_path = "celltype_model.rds") {
  # Load required libraries with explicit namespace calls
  required_pkgs <- c("ggplot2", "caret", "dplyr", "tidyr", "ggpubr", "RColorBrewer")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it first with install.packages('", pkg, "')", sep = ""))
    }
  }
  
  # Set safe font configuration to avoid font category errors
  set_safe_font <- function() {
    # Use base R font without specifying family to avoid platform issues
    ggplot2::theme_update(
      text = ggplot2::element_text(family = ""),
      axis.text = ggplot2::element_text(family = ""),
      plot.title = ggplot2::element_text(family = "", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(family = "", hjust = 0.5),
      # 恢复坐标线和网格线默认设置
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_line(color = "grey95"),
      axis.line = ggplot2::element_line(color = "black")
    )
  }
  set_safe_font()
  
  # Define custom color palette from user-provided colors
  custom_colors <- c("#FF6F61", "#FFD800", "#7FDBFF", "#1E90FF", "#9370DB", 
                     "#00CED1", "#DA70D6", "#CD853F", "#4CAF50", "#9ACD32", 
                     "#F8C3CD", "#20B2AA")
  
  # 细胞类型配色方案
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
  
  # Read model file
  if (!file.exists(model_path)) {
    stop("Model file does not exist. Please check the path: ", model_path)
  }
  model <- readRDS(model_path)
  cat("Successfully read model file: ", model_path, "\n", sep = "")
  
  # Create result saving directory
  if (!dir.exists("model_visualizations")) {
    dir.create("model_visualizations", showWarnings = FALSE)
  }
  
  # 1. Cross-validation results visualization (3:4)
  if (!is.null(model$resample)) {
    cv_results <- model$resample
    
    # Reshape cross-validation results
    cv_long <- tidyr::pivot_longer(cv_results, cols = c(Accuracy, Kappa), names_to = "Metric", values_to = "Value")
    
    # 调整颜色为#00CED1和#FFD800
    cv_colors <- c("#00CED1", "#FFD800")[1:length(unique(cv_long$Metric))]
    
    cv_plot <- ggplot2::ggplot(cv_long, ggplot2::aes(x = Metric, y = Value, fill = Metric)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
      ggplot2::ggtitle("Distribution of Cross-Validation Metrics", subtitle = paste("Method:", model$control$method)) +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_fill_manual(values = cv_colors) +
      # 确保坐标线可见
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95"),
        axis.line = ggplot2::element_line(color = "black")
      )
    
    print(cv_plot)
    ggplot2::ggsave("model_visualizations/01.cv_performance.pdf", cv_plot, width = 4, height = 6, dpi = 300)
    cat("Cross-validation results saved to: model_visualizations/cv_performance.pdf\n")
  } else {
    cat("No cross-validation results in model, skipping this visualization\n")
  }
  
  # 2. Confusion matrix visualization (4:4)
  if (!is.null(model$pred)) {
    # Extract prediction results
    pred_results <- model$pred
    pred_results$pred <- factor(pred_results$pred, levels = levels(pred_results$obs))
    
    # Calculate confusion matrix
    cm <- caret::confusionMatrix(pred_results$pred, pred_results$obs)
    cm_df <- as.data.frame(cm$table)
    colnames(cm_df) <- c("Predicted", "Actual", "Count")
    
    # Calculate percentages
    cm_df <- dplyr::group_by(cm_df, Actual) %>%
      dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
      dplyr::ungroup()
    
    # Plot confusion matrix heatmap with custom gradient
    cm_plot <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Predicted, y = Actual, fill = Percentage)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%d (%.1f%%)", Count, Percentage)), size = 3) +
      ggplot2::scale_fill_gradient(low = custom_colors[3], high = custom_colors[4]) +  # Light to dark blue
      ggplot2::labs(title = "Confusion Matrix", 
                    subtitle = paste("Overall Accuracy:", sprintf("%.2f%%", cm$overall["Accuracy"] * 100))) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95"),
        axis.line = ggplot2::element_line(color = "black")
      )
    
    print(cm_plot)
    ggplot2::ggsave("model_visualizations/02.confusion_matrix.pdf", cm_plot, width = 9, height = 8, dpi = 300)
    cat("Confusion matrix saved to: model_visualizations/confusion_matrix.pdf\n")
    
    # Plot class-specific performance metrics (4:3)
    # Handle cases where Accuracy column might not exist
    class_stats <- as.data.frame(cm$byClass) %>%
      dplyr::mutate(Class = rownames(.))
    
    # Check available metrics and select only existing ones
    available_metrics <- colnames(class_stats)
    desired_metrics <- c("Accuracy", "Sensitivity", "Specificity", "Precision")
    selected_metrics <- intersect(desired_metrics, available_metrics)
    
    # 为类别性能图指定特定颜色: #F8C3CD、#9ACD32、#1E90FF
    class_performance_colors <- c("#F8C3CD", "#9ACD32", "#9370DB")
    
    # If Accuracy not available, use what's available
    if (!"Accuracy" %in% selected_metrics) {
      warning("Accuracy metric not available for class-specific analysis, using available metrics instead")
    }
    
    # Only proceed if there are metrics to plot
    if (length(selected_metrics) > 0) {
      # 根据指标数量调整颜色（循环使用指定颜色）
      metric_colors <- rep(class_performance_colors, length.out = length(selected_metrics))
      
      class_stats <- class_stats %>%
        dplyr::select(dplyr::all_of(c("Class", selected_metrics))) %>%
        tidyr::pivot_longer(cols = -Class, names_to = "Metric", values_to = "Value")
      
      class_plot <- ggplot2::ggplot(class_stats, ggplot2::aes(x = Class, y = Value, fill = Metric)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Value)), 
                           position = ggplot2::position_dodge(0.8), vjust = -0.3, size = 3) +
        ggplot2::ylim(0, 1.1) +
        ggplot2::labs(title = "Class-specific Performance Metrics", x = "Cell Type", y = "Value") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.grid.major = ggplot2::element_line(color = "grey90"),
          panel.grid.minor = ggplot2::element_line(color = "grey95"),
          axis.line = ggplot2::element_line(color = "black")
        ) +
        ggplot2::scale_fill_manual(values = metric_colors)
      
      print(class_plot)
      ggplot2::ggsave("model_visualizations/03.class_performance.pdf", class_plot, width = 7, height = 7, dpi = 300)  # 4:3
      cat("Class performance metrics saved to: model_visualizations/class_performance.pdf\n")
    } else {
      cat("No valid performance metrics available for class-specific visualization\n")
    }
  } else {
    cat("No prediction results in model, cannot generate confusion matrix visualization\n")
  }
  
  # 3. Feature importance visualization (4:4)
  if (!is.null(model$finalModel) && !is.null(caret::varImp(model)$importance)) {
    # Extract feature importance
    imp <- caret::varImp(model)$importance
    imp_df <- data.frame(
      Feature = rownames(imp),
      Importance = imp[, 1],
      row.names = NULL
    ) %>%
      dplyr::arrange(dplyr::desc(Importance)) %>%
      dplyr::slice_head(n = 20)  # Top 20 most important features
    
    # Plot feature importance with custom gradient
    imp_plot <- ggplot2::ggplot(imp_df, ggplot2::aes(x = reorder(Feature, Importance), y = Importance, fill = Importance)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(low = custom_colors[9], high = custom_colors[10]) +  # Green gradient
      ggplot2::labs(title = "Feature Importance Ranking (Top 20)", x = "Feature", y = "Importance Score") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95"),
        axis.line = ggplot2::element_line(color = "black")
      )
    
    print(imp_plot)
    ggplot2::ggsave("model_visualizations/04.feature_importance.pdf", imp_plot, width = 4, height = 3, dpi = 300)
    cat("Feature importance plot saved to: model_visualizations/feature_importance.pdf\n")
  } else {
    cat("No feature importance information in model, skipping this visualization\n")
  }
  
  # 4. Class distribution visualization (4:3)
  if (!is.null(model$trainingData)) {
    # Extract class distribution from training data
    class_dist <- data.frame(
      Cell_Type = model$trainingData$.outcome
    ) %>%
      dplyr::group_by(Cell_Type) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(Count)) %>%
      dplyr::ungroup()
    
    # 匹配细胞类型配色（仅使用存在的细胞类型对应的颜色）
    existing_types <- class_dist$Cell_Type
    dist_colors <- celltype_colors[match(existing_types, names(celltype_colors))]
    # 若有未匹配的细胞类型，使用默认备用色
    if (any(is.na(dist_colors))) {
      na_indices <- which(is.na(dist_colors))
      dist_colors[na_indices] <- custom_colors[1:length(na_indices)]
    }
    
    # Plot class distribution
    dist_plot <- ggplot2::ggplot(class_dist, ggplot2::aes(x = reorder(Cell_Type, -Count), y = Count, fill = Cell_Type)) +
      ggplot2::geom_bar(stat = "identity", show.legend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.3, size = 3) +
      ggplot2::labs(title = "Distribution of Sample Counts by Cell Type", x = "Cell Type", y = "Sample Count") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95"),
        axis.line = ggplot2::element_line(color = "black")
      ) +
      ggplot2::scale_fill_manual(values = dist_colors)
    
    print(dist_plot)
    ggplot2::ggsave("model_visualizations/04.class_distribution.pdf", dist_plot, width = 8, height = 4, dpi = 300)  # 4:3
    cat("Class distribution plot saved to: model_visualizations/class_distribution.pdf\n")
  } else {
    cat("No training data information in model, cannot generate class distribution visualization\n")
  }
  
  # 5. TME adjustment information visualization (4:4)
  if (!is.null(model$TME_adjust) && !is.null(model$original_tme_ratio)) {
    tme_data <- data.frame(
      Type = c("Original Data", ifelse(is.numeric(model$TME_adjust), "Adjusted Data", "Unadjusted Data")),
      TME_Ratio = c(model$original_tme_ratio, 
                    ifelse(is.numeric(model$TME_adjust), model$target_tme_ratio, model$original_tme_ratio)) * 100
    )
    
    # 为TME比例对比图指定特定颜色: #F8C3CD、#9ACD32
    tme_colors <- rep(c("#F8C3CD", "#9ACD32"), length.out = nrow(tme_data))
    
    tme_plot <- ggplot2::ggplot(tme_data, ggplot2::aes(x = Type, y = TME_Ratio, fill = Type)) +
      ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", TME_Ratio)), vjust = -0.3, size = 4) +
      ggplot2::ylim(0, max(tme_data$TME_Ratio) * 1.2) +
      ggplot2::labs(title = "TME Cell Proportion Comparison", y = "TME Cell Proportion (%)") +
      ggplot2::scale_fill_manual(values = tme_colors) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_line(color = "grey95"),
        axis.line = ggplot2::element_line(color = "black")
      )
    
    print(tme_plot)
    ggplot2::ggsave("model_visualizations/05.tme_ratio_comparison.pdf", tme_plot, width = 3, height = 6, dpi = 300)
    cat("TME ratio comparison saved to: model_visualizations/tme_ratio_comparison.pdf\n")
  } else {
    cat("No TME adjustment information in model, skipping this visualization\n")
  }
  
  cat("\nAll visualization results have been saved to 'model_visualizations' folder\n")
  invisible(model)  # Return model object invisibly
}