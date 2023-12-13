# Pre-processing plasma proteomics for SPRINT
This project performed a number of pre-processing pipelines for plasma data generated for Sputum and plasma proteomes in response to CFTR therapy (SPRINT).

## Running analysis
- ```main.R``` performs all analysis, initially:
- All proteins > 20% missing are filtered out.
- Abundance values are log2 transformed.
- Analysis pipelines are run on this data, these steps include: performing batch correction using [Combat](https://rdrr.io/bioc/sva/man/ComBat.html) or "housekeeping" proteins (proteins with the lowest coefficient of variation within plates), performing normalisation of data - either within batches (prior to batch correction) or across all batches (post-batch correction), and performing these steps with and without removing depleted proteins.

## Results
- Mean IQR of all proteins in QC samples was calculated following running each pipeline, the results of this are in ```sorted_IQR.csv```.
- Each folder in ```figures``` contains the plots and ```.csv```'s for an analysis pipeline.
- Each folder name indicates the pipeline steps and order performed. e.g.
	- folder ```depl_withinvsn_combat``` contains data and plots following the steps of: depleted proteins being removed, then [VSN](https://bioconductor.org/packages/devel/bioc/vignettes/vsn/inst/doc/A-vsn.html) being performed within batches, before combat is run across all batches.
	- folder ```hk_vsn``` contains data and plots following the steps of: depleted proteins are NOT removed, "housekeeping" proteins are identified across plates and used to batch correct, before VSN is performed on all data to normalise.
- ```log_filter``` contains data/ figures from filtered and log2 transformed data only, prior to any other processing steps.
- ```*_df.csv``` is the processed data, including missing values at a stage. 
- ```*_mf_df.csv``` is the processed data, imputed using [MissForest](https://academic.oup.com/bioinformatics/article/28/1/112/219101?).
- Within each folder the prefix to each ```.csv``` or figure refers to the data or figure produced following that specific step during the pipeline. e.g. in the folder ```depl_withinvsn_combat```, ```depl_df.csv``` is the processed data following the step of removing depleted proteins, however ```combat_df.csv``` is the processed data following the final step in the pipeline of performing combat.

## Figures
- ```*_boxplot_cohort.png``` / ```*_boxplot_plate.png``` are boxplots of all proteins in each sample coloured by cohort/ plate respectively.
- ```*_heatmap_cohort.png``` / ```*_heatmap_plate.png``` are heatmaps of all proteins in each sample coloured by cohort/ plate respectively. Missing values are imputed using MissForest.
- ```*_boxplot_selected*.png``` are boxplots of selected proteins defined in ```line 52``` onwards of ```main.R```.
- ```*_deplete.png``` are boxplots of depleted proteins in each sample (if they have not been removed).
- ```*_pca_cohort.png``` / ```*_pca_plate.png``` are PCA plots of each sample coloured by cohort/ plate respectively. Missing values are imputed using MissForest.
- ```*_pca_var_imp.png``` are barplots of top contributing proteins to the PCA, coloured by whether these proteins were depleted or not.
- ```p(1-3)_low_cv.png``` are plot if "housekeeping" proteins were used to batch correct. These plots are the log2 filtered values of the proteins used for this purpose within plate 1, 2 and 3 respectively.

## Recommendations
- Combat appears to have corrected batch sufficiently, using "housekeeping" proteins does not.
- Removing depleted proteins (```deplr```), or not, does not seem to have substantially impacted preprocessing. 
- Performing VSN prior to combat does not impact batch correction, but results in slightly reduced QC protein IQR and therefore may be preferred.
- Therefore, either the ```depl_withinvsn_combat``` or ```withinvsn_combat``` pipelines would be recommended for further processing, depending on if you preferred to process with or without depleted proteins. Within these folders, the ```combat_df.csv``` contains the processed data, including missing values, and ```combat_mf_df.csv``` includes the processed data with missing values imputed using MissForest.