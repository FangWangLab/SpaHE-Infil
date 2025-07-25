#Set data directory and filename (modify according to your actual situation)
#data_directory <- "/path/to/your/spatial/data" # Path to the folder where data is stored
#data_filename <- "filtered_feature_bc_matrix.h5" # For H5 format; or "filtered_feature_bc_matrix" for MTX format
#Load the function (first copy the SpatialAnno_auto function code into the R environment and execute it)
#source("your_function_script.R") # Or directly paste the function code in the R console
#Run the analysis
#result <- SpatialAnno_auto( data.dir = data_directory, filename = data_filename)

SpatialAnno_auto <- function(data.dir, filename, assay = "Spatial", slice = "slice1") {
  # Check if required packages are installed and load them
  required_packages <- c(
    "Seurat", "ggplot2", "patchwork", "dplyr", "hdf5r", 
    "SingleR", "celldex", "scRNAseq", "scuttle", "scran",
    "AUCell", "GSEABase", "reshape2", "scales"
  )
  
  # Install missing packages
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing_packages, update = FALSE)
  }
  
  # Load all required packages
  invisible(lapply(required_packages, library, character.only = TRUE))
  
  # Set working directory to data directory
  setwd(data.dir)
  message("Working directory set to: ", data.dir)
  
  #################### 1. Data Loading and Preprocessing ####################
  # Check file type and choose appropriate loading method
  if (grepl("h5$", filename, ignore.case = TRUE)) {
    # H5 format file (contains all data)
    spatial_data <- Load10X_Spatial(
      data.dir = data.dir,
      filename = filename,
      assay = assay,
      slice = slice
    )
  } else {
    # MTX format files (matrix.mtx, barcodes.tsv and features.tsv)
    spatial_data <- Load10X_Spatial(
      data.dir = data.dir,
      assay = assay,
      slice = slice,
      filter.matrix = TRUE,  # Filter low-quality data
      to.upper = FALSE       # Whether to convert gene names to uppercase
    )
  }
  
  # Save raw data
  saveRDS(spatial_data, "01.raw_spatial_data.rds")
  message("Raw data saved as: 01.raw_spatial_data.rds")
  
  #################### 2. Data Quality Control and Normalization ####################
  message("Starting data quality control...")
  spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^MT-")
  
  # Generate QC visualization
  pdf("02.QC_metrics.pdf", width = 10)
  print(VlnPlot(spatial_data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), 
                pt.size = 0.1, ncol = 3) & 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
  message("QC results saved as: 02.QC_metrics.pdf")
  
  #################### 3. Data Normalization and Feature Selection ####################
  message("Starting data normalization and feature selection...")
  spatial_data <- SCTransform(spatial_data, assay = assay, verbose = FALSE)
  spatial_data <- FindVariableFeatures(spatial_data, selection.method = "vst", nfeatures = 3000)
  spatial_data <- ScaleData(spatial_data)
  
  #################### 4. Dimensionality Reduction and Clustering ####################
  message("Starting dimensionality reduction and clustering analysis...")
  spatial_data <- RunPCA(spatial_data, features = VariableFeatures(object = spatial_data))
  
  set.seed(123)
  spatial_data <- RunTSNE(
    spatial_data,
    dims = 1:30,
    reduction = "pca",
    perplexity = 30,
    check_duplicates = FALSE
  )
  
  spatial_data <- FindNeighbors(spatial_data, dims = 1:30)
  spatial_data <- FindClusters(spatial_data, resolution = 1.2)
  
  # Save tSNE clustering results
  pdf("03.tSNE_clusters.pdf", width = 5, height = 5)
  print(DimPlot(spatial_data, reduction = "tsne", label = TRUE) + NoLegend())
  dev.off()
  message("tSNE clustering results saved as: 03.tSNE_clusters.pdf")
  
  #################### 5. Spatial Visualization and Image Optimization ####################
  message("Starting spatial clustering visualization...")
  if (!is.null(spatial_data@images[[slice]])) {
    message("High-resolution image loaded successfully, dimensions:", dim(spatial_data@images[[slice]]@image))
  } else {
    stop("Image loading failed, please check path and files")
  }
  
  pdf("04.spatial_clusters.pdf", width = 6, height = 5)
  print(SpatialDimPlot(
    spatial_data,
    images = slice,
    stroke = NA,
    pt.size.factor = 1.6,
    label = TRUE,
    label.size = 4
  ) + 
    theme(legend.position = "right") +
    ggtitle("Spatial Clustering"))
  dev.off()
  message("Spatial clustering results saved as: 04.spatial_clusters.pdf")
  
  #################### 6. AUCell Automatic Cell Type Annotation ####################
  message("Starting automatic cell type annotation...")
  
  # Define cell type marker gene sets
  markers <- list(
    "Tumor cells" = c("EPCAM", "KRT19", "MUC1"),
    "Normal epithelial cells" = c("FABP1", "CFTR", "KRT7","PAX8"),
    "CD8 T cells" = c("CD8A", "CD8B", "GZMB", "IFNG"),
    "CD4 T cells" = c("CD4", "IL7R"),
    "Macrophages" = c("CD68", "CD163", "CSF1R"),
    "Neutrophils" = c("S100A8", "S100A9", "FCGR3B"),
    "NK cells" = c("NKG7", "GNLY", "NCR1"),
    "Tregs" = c("FOXP3", "IL2RA"),
    "Fibroblasts" = c("COL1A1", "DCN", "PDGFRA"),
    "Endothelial cells" = c("PECAM1", "VWF", "CDH5"),
    "B cells" = c("CD19", "CD79A", "MS4A1", "CXCL13", "CD27", "CCL19"),
    "Dendritic cells" = c("CD1C", "CD141", "CD209", "HLA-DRA", "IDO1")
  )
  
  # Check and filter genes not present in data
  valid_markers <- lapply(markers, function(genes) {
    genes[genes %in% rownames(spatial_data)]
  })
  valid_markers <- valid_markers[sapply(valid_markers, length) > 0]
  
  if (length(valid_markers) == 0) {
    stop("All marker genes are missing from the data, please check gene names")
  }
  
  # Prepare input data
  spatial_counts <- GetAssayData(spatial_data, assay = assay, slot = "counts")
  
  # Run AUCell analysis pipeline
  set.seed(123)
  cells_rankings <- AUCell_buildRankings(spatial_counts, plotStats = FALSE)
  
  gene_sets <- GeneSetCollection(
    lapply(names(valid_markers), function(ct) {
      GeneSet(valid_markers[[ct]], setName = ct)
    })
  )
  
  cells_auc <- AUCell_calcAUC(gene_sets, cells_rankings)
  auc_matrix <- t(getAUC(cells_auc))
  
  # Boost AUC scores for CD8 T cells
  cd8_label <- "CD8 T cells"
  if (cd8_label %in% colnames(auc_matrix)) {
    auc_matrix[, cd8_label] <- auc_matrix[, cd8_label] * 1.3
    message("Boosted AUC scores for CD8 T cells to increase sensitivity")
  } else {
    warning("CD8 T cell marker gene set not found, please check naming")
  }
  
  # Assign cell type with highest AUC score to each cell
  spatial_data$auc_celltype <- colnames(auc_matrix)[apply(auc_matrix, 1, which.max)]
  
  #################### 7. Visualize AUC Scores and Annotation Results ####################
  message("Generating annotation visualizations...")
  
  # Define cell type color mapping (predefined for consistency)
  celltype_colors <- c(
    "Tumor cells" = "#FF6F61",                 # Pink (primary visual color, most prominent)
    "Normal epithelial cells" = "#FFD800",     # Gold (high saturation)
    "CD8 T cells" = "#7FDBFF",                 # Mint blue (cool tone)
    "CD4 T cells" = "#1E90FF",                 # Sea green (deep cool green, contrasts with tumor pink)
    "Macrophages" = "#9370DB",                 # Medium violet (medium saturation)
    "Neutrophils" = "#00CED1",                 # Dark turquoise (high saturation cool color)
    "NK cells" = "#DA70D6",                    # Orchid purple (reddish purple, distinct from macrophages)
    "Tregs" = "#CD853F",                       # Peru brown (warm brown)
    "Fibroblasts" = "#4CAF50",                 # Tomato red (high saturation warm color)
    "Endothelial cells" = "#9ACD32",           # Yellow-green (transition color)
    "B cells" = "#F8C3CD",                     # Light coral (high saturation pink)
    "Dendritic cells" = "#20B2AA"              # Light sea green (cool color)
  )
  
  # Add missing color definitions
  if(!all(unique(spatial_data$auc_celltype) %in% names(celltype_colors))) {
    missing_types <- setdiff(unique(spatial_data$auc_celltype), names(celltype_colors))
    warning(paste("The following cell types have no color definition, using defaults:", paste(missing_types, collapse = ", ")))
    celltype_colors <- c(celltype_colors, scales::hue_pal()(length(missing_types)))
    names(celltype_colors)[(length(celltype_colors)-length(missing_types)+1):length(celltype_colors)] <- missing_types
  }
  
  # AUC score distribution boxplot (using consistent colors)
  pdf("05.AUC_scores_boxplot.pdf", width = 6, height = 5)
  # Convert to dataframe for ggplot2
  auc_df <- reshape2::melt(auc_matrix)
  colnames(auc_df) <- c("Cell", "CellType", "AUC")
  # Ensure plotting order matches color order
  auc_df$CellType <- factor(auc_df$CellType, levels = names(celltype_colors))
  print(ggplot(auc_df, aes(x = CellType, y = AUC, fill = CellType)) +
          geom_boxplot() +
          scale_fill_manual(values = celltype_colors) +
          labs(title = "AUC Scores for Each Cell Type", y = "AUC Score") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # Annotation results on tSNE (using consistent colors)
  pdf("06.AUCell_tsne.pdf", width = 6, height = 5)
  print(DimPlot(spatial_data, reduction = "tsne", group.by = "auc_celltype", label = TRUE) +
          ggtitle("AUCell Annotation (tSNE)") +
          theme(legend.position = "right") +
          scale_color_manual(values = celltype_colors))
  dev.off()
  
  # Spatial distribution of annotations (using consistent colors)
  pdf("06.AUCell_spatial.pdf", width = 6, height = 5)
  print(SpatialDimPlot(
    spatial_data,
    images = slice,
    group.by = "auc_celltype",
    pt.size.factor = 1.5,
    stroke = NA
  ) + ggtitle("AUCell Annotation (Spatial)") +
    theme(legend.position = "right") +
    scale_fill_manual(values = celltype_colors))
  dev.off()
  
  #################### 8. Cell Type Proportion Visualization ####################
  message("Generating cell type proportion statistics...")
  celltype_counts <- table(spatial_data$auc_celltype)
  celltype_ratio <- prop.table(celltype_counts) * 100
  
  ratio_df <- data.frame(
    CellType = names(celltype_ratio),
    Proportion = as.numeric(celltype_ratio)
  )
  # Ensure plotting order matches color order
  ratio_df$CellType <- factor(ratio_df$CellType, levels = names(celltype_colors))
  
  # Cell type proportion plot (using consistent colors)
  pdf("08.celltype_proportion.pdf", width = 6, height = 4)
  print(ggplot(ratio_df, aes(x = CellType, y = Proportion, fill = CellType)) +
          geom_col(width = 0.7) +
          geom_text(aes(label = sprintf("%.1f%%", Proportion)), vjust = -0.3, size = 3.5) +
          scale_fill_manual(values = celltype_colors) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                plot.title = element_text(hjust = 0.5, size = 14)) +
          labs(title = "Cell type distribution", x = "Cell type", y = "Proportion (%)") +
          ylim(0, max(ratio_df$Proportion) + 5))
  dev.off()
  
  #################### 9. Save Final Results ####################
  message("Saving final results...")
  
  # Save final annotation results (PDF version)
  pdf("11.final_spatial_annotation.pdf", width = 6, height = 5)
  print(SpatialDimPlot(
    spatial_data,
    images = slice,
    group.by = "auc_celltype",
    pt.size.factor = 2,
    stroke = NA
  ) + ggtitle("Spatial distribution of cell types") +
    theme(legend.position = "right") +
    scale_fill_manual(values = celltype_colors))
  dev.off()
  
  # Save high-resolution TIFF version
  width_inches <- 3072 / 300
  height_inches <- 2048 / 300
  tiff("11.final_spatial_annotation.tif", 
       width = width_inches, height = height_inches, 
       units = "in", res = 300, compression = "lzw")
  print(SpatialDimPlot(
    spatial_data,
    images = slice,
    group.by = "auc_celltype",
    pt.size.factor = 2,
    stroke = NA
  ) + ggtitle("Spatial distribution of cell types") +
    theme(legend.position = "right") +
    scale_fill_manual(values = celltype_colors))
  dev.off()
  
  # Save final annotated data
  saveRDS(spatial_data, "10.final_annotated_data.rds")
  message("Analysis complete! Final annotated data saved as: 10.final_annotated_data.rds")
  
  # Return annotated object
  return(spatial_data)
}