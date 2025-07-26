# SpaHE-Infil

SpaHE-Infil is a spatial analysis tool for the tumor microenvironment, enabling the reconstruction of spatial distributions of key cell types from routine HE-stained sections. It trains models using spatial transcriptomics data and constructs a multimodal predictive model by integrating nuclear-cytoplasmic morphology, texture features, and cell density from HE images. Its innovation lies in the introduction of a dynamic TME proportion calibration mechanism, capable of outputting parameter-adjustable visualization results, thus providing an effective tool for pathological diagnosis and tumor immunology research.

<img width="727" height="727" alt="image" src="https://github.com/user-attachments/assets/eb2f58ca-7eb7-46c3-8ec9-51fb4d60c465" />

## Core Workflow

This pipeline consists of 4 main scripts and 2 supplementary scripts, designed to be run in sequence:

### Main Scripts
 **01.SpatialAnno_auto.R**  
   - Automated spatial domain annotation  
   - Initial cell type identification
   - This pipeline supports other annotation methods (including manual annotation) for integration into subsequent processes.

 **02.train_model-TMEadj.R**  
   - Train TME-adjusted classification model  
   - Integrate spatial context features

 **03.predict_celltypes-TMEadj.R**  
   - Predict cell types using trained model  
   - Generate spatial distribution maps

 **04.spa-HE evaluate Y-correction.R**  
   - Evaluate prediction accuracy  
   - Perform Y-correction analysis  
   - Integrate with HE staining data

### Supplementary Scripts
- **Sup.train_model-TMEadj_HE feature get.R**  
  Extract features from HE images for model training

- **Sup.train_model-TMEadj_model features.R**  
  Additional feature engineering for model optimization

## Usage Steps
1. Run spatial annotation  
   `Rscript 01.SpatialAnno_auto.R`

2. Train the model  
   `Rscript 02.train_model-TMEadj.R`
   
   Supplementary Script. Prepare features using supplementary scripts  
      `Rscript Sup.train_model-TMEadj_HE feature get.R`  
      `Rscript Sup.train_model-TMEadj_model features.R`
   
4. Predict cell types  
   `Rscript 03.predict_celltypes-TMEadj.R`

5. Evaluate and correct results  
   `Rscript 04.spa-HE evaluate Y-correction.R`

## Usage method
**1.Spatial annotation**   
   source("01.SpatialAnno_auto.R")   
   result <- SpatialAnno_auto( data.dir = data_directory, filename = data_filename)

**2. SpaHE-Infil train model**   
   source("02.train_model-TMEadj.R")   
   model <- train_celltype_model(rds_files = sample_files, output_model = "celltype_model.rds",
   min_samples_per_type = 5, allow_unknown = TRUE, density_radius = 15,
   TME_adjust = FALSE/Calibration value)

**3. SpaHE-Infil predict model**   
   source("03.predict_celltypes-TMEadj.R")   
   predictions <- generate_predictions(he_image_path = "path/to/your/he_image.tif", model_path = "celltype_model.rds",
   output_prefix = "he_prediction_original_tme",
   TME_adjust = FALSE/Calibration value, step_size = 10)
   
**4. Evaluate the analysis results of SpaHE-Infil using ST**   
   source("04.spa-HE evaluate Y-correction.R")   
   evaluation_results <- evaluate_spahe_accuracy(st_rds_path = "/path/to/st_data.rds", 
   spahe_pred_path = "/path/to/spahe_predictions.RData", output_prefix = "path/to/your/spahe_predictions_directory",
   grid_size = 150, min_cells_per_grid = 3, scale_spahe_coords = TRUE)


## Requirements
- R (v4.3.0+)
- dplyr (1.1.0+)
- lubridate (1.8.0+)
- Seurat (5.0.0+)
- readr (2.1.0+)
- Required packages specified in individual scripts
- Requirements for input H&E image quality: no less than 900×600 pixels, with a resolution of at least 300 dpi.

## Notes
Refer to comments within each script for parameter configuration and detailed usage instructions.

## Author contact
wangfang_lukas@jiangnan.edu.cn or wangfang_lukas@zju.edu.cn
