library(pmp)
library(tidyverse)

setwd(wd)

# load data
df <- read.csv("FINAL_sprint_plasma_spectronaut_protein_output.csv")
meta_df <- read.csv("FINAL_sprint_plasma_spectronaut_output_grouping.csv")

# creat folder
dir.create(file.path(paste(wd, "/figures/QCRSC", sep = "")), showWarnings = FALSE)

# clean
meta_df$Cohort[meta_df$Cohort == "QC "] <- "QC"
meta_df$Cohort[meta_df$Cohort == "healthy "] <- "healthy"
colnames(df)[colnames(df) == "ï..PG.Genes"] <- "PG.Genes"
colnames(meta_df)[colnames(meta_df) == "ï..Order"] <- "Order"
meta_cols <- c("PG.Genes", "PG.ProteinDescriptions", "PG.UniProtIds")
data_cols <- setdiff(names(df), meta_cols)

# filter - drop proteins than are more than 20% missing
filter_missing_df <-
  df[rowSums(is.na(df[, data_cols])) / length(data_cols) <= 0.20, ]

# log2
log2_df <- filter_missing_df %>%
  mutate(across(where(is.numeric), ~log2(. + 1)))  # add constant

# signal drift correction PMP
# https://bioconductor.org/packages/release/bioc/vignettes/pmp/inst/doc/pmp_vignette_signal_batch_correction_mass_spectrometry.html
# JUST USING QC'S AS RECOMMENDED
qc_rm_depl_df <- QCRSC(df = t(log2_df[, data_cols]),
                        order = meta_df$Run.order,
                        batch = meta_df$Plate,
                        classes = meta_df$Cohort,
                        log = FALSE,
                        spar = 0,
                        minQC = 4)

# requires 4 QC data points per plate, only 2 in P3
plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = qc_rm_depl_df,
                  classes = meta_df$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/qc.pdf")


# QC AND HEALTHY as QC
healthy_as_qc_meta <- meta_df
healthy_as_qc_meta$Cohort[healthy_as_qc_meta$Cohort == "healthy"] <- "QC"

qc_healthy_depl_df <- QCRSC(df = t(log2_df[, data_cols]),
                        order = meta_df$Run.order,
                        batch = meta_df$Plate,
                        classes = healthy_as_qc_meta$Cohort,
                        log = FALSE,
                        spar = 0,
                        minQC = 4)

# requires 4 QC data points per plate, only 2 in P3
plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = qc_healthy_depl_df,
                  classes = healthy_as_qc_meta$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/qc_healthy.pdf")

# USING HEALTHY AS QC, BUT PROCESS BATCHES SEPARATELY, SO CAN DO BATCH CORRECTION AFTER
healthy_as_qc_meta <- meta_df
healthy_as_qc_meta$condition[healthy_as_qc_meta$condition == "healthy"] <- "QC"

p1_df <- log2_df[, grep("^P1", names(log2_df))]
p1_meta <- meta_df[healthy_as_qc_meta$Plate == "P1", ]
p2_df <- log2_df[, grep("^P2", names(log2_df))]
p2_meta <- meta_df[healthy_as_qc_meta$Plate == "P2", ]
p3_df <- log2_df[, grep("^P3", names(log2_df))]
p3_meta <- meta_df[healthy_as_qc_meta$Plate == "P3", ]

corrected_data_all_p1 <- QCRSC(df = t(p1_df),
                        order = p1_meta$Run.order,
                        batch = p1_meta$Plate,
                        classes = p1_meta$Cohort,
                        spar = 0,
                        log = FALSE,
                        minQC = 4)

corrected_data_all_p2 <- QCRSC(df = t(p2_df),
                        order = p2_meta$Run.order,
                        batch = p2_meta$Plate,
                        classes = p2_meta$Cohort,
                        spar = 0,
                        log = FALSE,
                        minQC = 4)

corrected_data_all_p3 <- QCRSC(df = t(p3_df),
                        order = p3_meta$Run.order,
                        batch = p3_meta$Plate,
                        classes = p3_meta$Cohort,
                        spar = 0,
                        log = FALSE,
                        minQC = 4)

# merge df's
merged_corrected_data_all_df <- cbind(corrected_data_all_p1, corrected_data_all_p2, corrected_data_all_p3)

plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = merged_corrected_data_all_df,
                  classes = healthy_as_qc_meta$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/healthyasQC_diff_batch.pdf")


# USE ALL AS QC
corrected_data_all <- QCRSC(df = t(log2_df[, data_cols]),
                        order = meta_df$Run.order,
                        batch = meta_df$Plate,
                        classes = rep("QC", length(meta_df$Cohort)),
                        log = FALSE,
                        spar = 0,
                        minQC = 4)

plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = corrected_data_all,
                  classes = meta_df$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/allasQC.pdf")

# ALL AS QC, AND SAME BATCH. TO SEE IF ANY OTHER BATCH CORRECTION STEP BEING PERFORMED
corrected_data_all_same_batch <- QCRSC(df = t(log2_df[, data_cols]),
                        order = meta_df$Run.order,
                        batch = rep("P1", length(meta_df$Plate)),
                        classes = rep("QC", length(meta_df$Cohort)),
                        spar = 0,
                        minQC = 4)

plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = corrected_data_all_same_batch,
                  classes = meta_df$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/allasQC_same_batch.pdf")


# ALL AS QC, BUT SEPARATE BATCHES
p1_df <- log2_df[, grep("^P1", names(log2_df))]
p1_meta <- meta_df[meta_df$Plate == "P1", ]
p2_df <- log2_df[, grep("^P2", names(log2_df))]
p2_meta <- meta_df[meta_df$Plate == "P2", ]
p3_df <- log2_df[, grep("^P3", names(log2_df))]
p3_meta <- meta_df[meta_df$Plate == "P3", ]

corrected_data_all_p1 <- QCRSC(df = t(p1_df),
                        order = p1_meta$Run.order,
                        batch = p1_meta$Plate,
                        classes = rep("QC", length(p1_meta$Cohort)),
                        spar = 0,
                        log = FALSE,
                        minQC = 4)

corrected_data_all_p2 <- QCRSC(df = t(p2_df),
                        order = p2_meta$Run.order,
                        batch = p2_meta$Plate,
                        classes = rep("QC", length(p2_meta$Cohort)),
                        spar = 0,
                        log = FALSE,
                        minQC = 4)

corrected_data_all_p3 <- QCRSC(df = t(p3_df),
                        order = p3_meta$Run.order,
                        batch = p3_meta$Plate,
                        classes = rep("QC", length(p3_meta$Cohort)),
                        spar = 0,
                        log = FALSE,
                        minQC = 4)

# merge df's
merged_corrected_data_all_df <- cbind(corrected_data_all_p1, corrected_data_all_p2, corrected_data_all_p3)

plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = merged_corrected_data_all_df,
                  classes = meta_df$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/allasQC_diff_batch.pdf")

# USING HEALTHY ONLY AS QC
corrected_data_healthy <- QCRSC(df = t(log2_df[, data_cols]),
                        order = meta_df$Run.order,
                        batch = meta_df$Plate,
                        qc_label = "healthy",
                        log = FALSE,
                        classes = meta_df$Cohort,
                        spar = 0,
                        minQC = 4)

plots <- sbc_plot(df = t(log2_df[, data_cols]),
                  corrected_df = corrected_data_healthy,
                  classes = meta_df$Cohort,
                  batch = meta_df$Plate,
                  output = "figures/QCRSC/healthyasQC.pdf")
                  