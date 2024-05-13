library(S4Vectors)
library(reshape2)
library(gridExtra)
library(DEP)
library(stringr)

source("funcs_preprocess.R")

setwd(wd)
dir_path <- paste(wd, "/figures", sep = "")

# load data
df <- read.csv("FINAL_sprint_plasma_spectronaut_protein_output.csv")
meta_df <- read.csv("FINAL_sprint_plasma_spectronaut_output_grouping.csv")

# clean
meta_df$Cohort[meta_df$Cohort == "QC "] <- "QC"
meta_df$Cohort[meta_df$Cohort == "healthy "] <- "healthy"
colnames(df)[colnames(df) == "ï..PG.Genes"] <- "PG.Genes"
colnames(meta_df)[colnames(meta_df) == "ï..Order"] <- "Order"
meta_cols <- c("PG.Genes", "PG.ProteinDescriptions", "PG.UniProtIds")
data_cols <- setdiff(names(df), meta_cols)

plot_missingness(df)
ggsave(paste(dir_path, "/pre_filter_missing.png", sep = ""))

# filter - drop proteins than are more than 20% missing
filter_missing_df <-
  df[rowSums(is.na(df[, data_cols])) / length(data_cols) <= 0.20, ]

plot_missingness(filter_missing_df)
ggsave(paste(dir_path, "/post_filter_missing.png", sep = ""))

# log2
log2_df <- filter_missing_df %>%
  mutate(across(where(is.numeric), ~log2(. + 1)))  # add constant

# table of missing values
na_count <- rowSums(is.na(log2_df))
na_count_df <- data.frame(
        PG.Genes = log2_df$PG.Genes,
        PG.ProteinDescriptions = log2_df$PG.ProteinDescriptions,
        PG.UniProtIds = log2_df$PG.UniProtIds,
        missing_count = na_count
        )

write.csv(x = na_count_df,
        file = paste(dir_path, "/missing_count.csv", sep = ""))

# selected proteins for boxplots
# suspected to differ
selected_p <- c("Protein S100-A8", "Protein S100-A9",
        "C-reactive protein", "Fibrinogen-like protein 1")

# most and least differences in previous paper
# https://www.cysticfibrosisjournal.com/article/S1569-1993(23)01669-7/fulltext
most_diff <- c("CD59 glycoprotein", "Cystatin-B", "Cathepsin D",
        "Polymeric immunoglobulin receptor", "Antileukoproteinase",
        "Prothrombin", "Gamma-glutamyl hydrolase", "Protein S100-A8")

least_diff <- c("Cathepsin G", "Myeloblastin",
        "Peptidoglycan recognition protein 1", "Glycogen phosphorylase, liver form",
        "Matrix metalloproteinase-9")

# proteins with lowest coeff variation across batches
low_coeff_res <- low_coeff_var_p(log2_df, pc = 0.05, max_p = 5, data_cols)
low_coeff_var <- low_coeff_res$names

# plot selected proteins in boxplots
selected_p <- list("selected_p" = selected_p,
                "most_diff" = most_diff,
                "least_diff" = least_diff,
                "low_coeff_var" = low_coeff_var)

# process different pipelines, each string corresponds to the order
# that each process is run
pipelines <- c(
                "depl_combat_vsn",
                "combat_vsn",
                "depl_hk_vsn",
                "hk_vsn",
                "depl_withinvsn_combat",
                "withinvsn_combat",
                "depl_withinvsn_hk",
                "withinvsn_hk"

                ## arsyn is not running reliably due to incompatible experimental design
                #"depl_arsyn_vsn",
                #"arsyn_vsn"
                #"depl_withinvsn_arsyn",
                #"withinvsn_arsyn"
                )

iqr_df <- run_pipelines(pipelines,
                log2_df,
                df,
                meta_df,
                data_cols,
                meta_cols,
                dir_path,
                selected_p_ = selected_p)

# sort IQR
iqr_df <- iqr_df[order(iqr_df$IQR), ]
print(iqr_df)

write.csv(x = iqr_df, file = paste(dir_path, "/sorted_IQR.csv", sep = ""))
