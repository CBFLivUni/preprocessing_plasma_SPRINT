# functions for assisting with `differential_expression.R`

library(missForest)

check_data_order <- function(meta_df_, pro_df_) {

	# check that column order of meta data and protein data aligns.
	# if it doesn't throw error.

	if (all(as.character(meta_df_$Sample) == colnames(pro_df_)) == FALSE) {
			# if all values aren't true, then throw error

			stop("meta_data sample column and data protein column names don't align")
		} else {
			print("Samples in meta data align with columns in protein data")
		}
}