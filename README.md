# Processing plasma proteomics for SPRINT
This project contains pre-processing and processing for plasma proteomics for SPRINT. It involved performing a number of pre-processing pipelines followed by DE and GSEA for plasma data generated for Sputum and plasma proteomes in response to CFTR therapy (SPRINT).

# Pre-processing:

## Running analysis
- ```main.R``` performs all analysis, initially:
- All proteins > 20% missing are filtered out.
- Abundance values are log2 transformed.
- Analysis pipelines are run on this data, these steps include: performing batch correction using [Combat](https://rdrr.io/bioc/sva/man/ComBat.html) or "housekeeping" proteins (the most abundant proteins with the lowest coefficient of variation within plates), performing normalisation of data - either within batches (prior to batch correction) or across all batches (post-batch correction), and performing these steps with and without removing depleted proteins.
- ```model_metadata.R``` performs variance partition analysis and saves figures.

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
- For pipelines where depleted proteins were not removed during pre-processing, ```"*_rm_depl_df.csv"``` is the same data as ```*_df.csv``` processed data, without depleted proteins. ```"*_rm_depl_mf_df.csv"``` is the equivalent with the data imputed using MissForest.

## Figures
- ```*_boxplot_cohort.png``` / ```*_boxplot_plate.png``` are boxplots of all proteins in each sample coloured by cohort/ plate respectively.
- ```*_heatmap_cohort.png``` / ```*_heatmap_plate.png``` are heatmaps of all proteins in each sample coloured by cohort/ plate respectively. Missing values are imputed using MissForest.
- ```*_boxplot_selected*.png``` are boxplots of selected proteins defined in ```line 52``` onwards of ```main.R```.
- NB. in cases where the proteins shown on ```*_boxplot_selected_least.png``` are the same as ```*_boxplot_selected_most.png```. Then the proteins in ```*_boxplot_selected_least.png``` have been filtered from the dataset, due to missingness, and R is saving the most recent figure.
- ```*_deplete.png``` are boxplots of depleted proteins in each sample (if they have not been removed).
- ```*_pca_cohort.png``` / ```*_pca_plate.png``` are PCA plots of each sample coloured by cohort/ plate respectively. Missing values are imputed using MissForest.
- ```*_pca_var_imp.png``` are barplots of top contributing proteins to the PCA, coloured by whether these proteins were depleted or not.
- ```p(1-3)_low_cv.png``` are plot if "housekeeping" proteins were used to batch correct. These plots are the log2 filtered values of the proteins used for this purpose within plate 1, 2 and 3 respectively.
 - A linear mixed model was fit for each protein to estimate the contribution of Cohort and Plate. In ```*_var_part.png```, A) shows the variance explained by Cohort and Plate of the top 30 proteins ordered by Cohort. B) shows violin plots representing the distribution of variance explained of all proteins.

## QC-RSC for Signal Drift and Batch Correction
- ```check_signal_drift.R``` runs testing of using the [QC-RSC algorithm](https://bioconductor.org/packages/release/bioc/vignettes/pmp/inst/doc/pmp_vignette_signal_batch_correction_mass_spectrometry.html) to correct batch effect and signal drift. 
- Folder ```QCRSC``` contains the results.
- The algorithm involves fitting a curve to QC samples to check for, and then attempt to remove, signal drift and batch effect at the same time.
- There were not enough QC samples in Plate 3 (4 are required) to run the algorithm as is intended. Instead, fitting the curve to healthy samples only and then to all samples was tested.
- Using healthy samples generally over-corrected for signal drift, worsening data quality.
- Using all samples, often did not improve signal drift and potentially removed biological variation.

## Recommendations
- Attempting to correct for signal drift using QC-RSC seemed to result in poorer data quality.
- Combat appears to have corrected batch sufficiently, using "housekeeping" proteins does not.
- Removing depleted proteins (```deplr```), or not, does not seem to have substantially impacted pre-processing, though it may be best to include them during pre-processing to reduce the risk of biasing results. 
- Performing VSN prior to combat, compared to afterwards, does not impact batch correction, but results in slightly reduced QC protein IQR and therefore may be preferred.
- Therefore, the ```withinvsn_combat``` pipeline would be recommended for further processing. Within these folders, the ```combat_rm_depl_df.csv``` contains the processed data, including missing values, and ```combat_rm_depl_mf_df.csv``` includes the processed data with missing values imputed using MissForest.

# Processing:
- ```differential_expression.R```, ```funcs_DE.R``` and ```DE_figs.R``` contain differential expression analysis and figures for pre vs post ETI and pre ETI vs healthy and post ETI vs healthy.
- ```fpsea.R``` runs functional gene set enrichment analysis and produces figures.