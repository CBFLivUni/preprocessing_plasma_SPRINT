library(limma)
library(tidyverse)
library(Biobase)

source("funcs_DE.R")

# load data
pro_df <- read.csv("figures/withinvsn_combat/combat_rm_depl_df.csv")
meta_df <- read.csv("cleaned_meta.csv", stringsAsFactors=TRUE)  # written during `model_metadata.R`

# drop ivacaftor samps
iva_samps <- read.csv("ivacaftor_samps.csv")[, 1]
pro_df <- pro_df[, !(names(pro_df) %in% iva_samps)]

# both 6021 labelled as pre, drop for now.
pro_df <- pro_df %>% select(-contains("6021"))
meta_df <- meta_df[which(meta_df["ID"] != 6021), ]

# check samples align in data and meta_data
# sort meta data sample alphabetically
meta_df <- meta_df[order(meta_df$Sample), ]
# sort colnames alphabetically
pro_df <- pro_df[, order(colnames(pro_df))]

# RUN PRE- POST ANALYSIS
pro_pre_post_df <- pro_df %>% select((contains("ETI_post") | contains("ETI_pre")))
meta_pre_post_df <- meta_df[which(meta_df["Cohort"] == "ETI post" | meta_df["Cohort"] == "ETI pre"), ]

pro_pre_post_df <- pro_pre_post_df[, order(colnames(pro_pre_post_df))]  # select changes order

# drop unused levels
meta_pre_post_df <- droplevels(meta_pre_post_df)

table(meta_pre_post_df$Cohort)

# only keep cols need for imputation
meta_for_imp_df <- subset(meta_pre_post_df, select = -c(NEW, ID, Sample))

# impute meta
imp_res <- missForest(meta_for_imp_df)
imp_df <- as.data.frame(imp_res$ximp)

# add id back as factor
# https://support.bioconductor.org/p/9153666/
imp_df$ID <- as.factor(meta_pre_post_df$ID)

# pre vs post

# check AIC of combinations of covariates/ factors for pre vs post
# needs to use imputed protein data, as can't have NA's for AIC
pro_imp_df <- read.csv("figures/withinvsn_combat/combat_rm_depl_mf_df.csv")

# drop ivacaftor samps
pro_imp_df <- pro_imp_df[, !(names(pro_imp_df) %in% iva_samps)]
# keep only pre and post
pro_imp_df <- pro_imp_df %>% select((contains("ETI_post") | contains("ETI_pre")))

# both 6021 labelled as pre, drop for now.
pro_imp_df <- pro_imp_df %>% select(-contains("6021"))

# sort colnames alphabetically to match meta
pro_imp_df <- pro_imp_df[, order(colnames(pro_imp_df))]

check_data_order(meta_pre_post_df, pro_imp_df)

designs <- list(
	# where A=age, F=FEV, W=Weight, T=V1-V2

	# AFWT
	design_AFWT = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit + imp_df$Baseline.FEV1. +
	imp_df$Baseline.weight + imp_df$V1.to.V2),

	# AFW, AWT, FWT, AFT
	design_FWT = model.matrix(~0 + imp_df$Cohort + imp_df$Baseline.FEV1. +
	imp_df$Baseline.weight + imp_df$V1.to.V2),
	design_AWT = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit +
	imp_df$Baseline.weight + imp_df$V1.to.V2),
	design_AFT = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit + imp_df$Baseline.FEV1. +
	imp_df$V1.to.V2),
	design_AFW = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit + imp_df$Baseline.FEV1. +
	imp_df$Baseline.weight),

	# A, F, W, T
	design_A = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit),
	design_F = model.matrix(~0 + imp_df$Cohort + imp_df$Baseline.FEV1.),
	design_W = model.matrix(~0 + imp_df$Cohort + imp_df$Baseline.weight),
	design_T = model.matrix(~0 + imp_df$Cohort + imp_df$V1.to.V2),

	# AF, AW, AT, FW, FT, WT
	design_AF = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit + imp_df$Baseline.FEV1.),
	design_AW = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit + imp_df$Baseline.weight),
	design_AT = model.matrix(~0 + imp_df$Cohort + imp_df$Age.at.baseline.visit + imp_df$V1.to.V2),
	design_FW = model.matrix(~0 + imp_df$Cohort + imp_df$Baseline.FEV1. + imp_df$Baseline.weight),
	design_FT = model.matrix(~0 + imp_df$Cohort + imp_df$Baseline.FEV1. + imp_df$V1.to.V2),
	design_WT = model.matrix(~0 + imp_df$Cohort + imp_df$Baseline.weight + imp_df$V1.to.V2),

	# None
	design_none = model.matrix(~0 + imp_df$Cohort)
)

aics <- selectModel(pro_imp_df, designlist = designs, criterion = "aic")

# get best model for each protein
pref_aic <- data.frame(table(aics$pref))
colnames(pref_aic)[1] <- "Model"
pref_aic$Model <- factor(pref_aic$Model, levels = pref_aic$Model[order(-pref_aic$Freq)])
pref_aic$Model <- gsub("design_", "", pref_aic$Model)

# plot times each combo of vars/ factors preferable

p_aic <- ggplot(data = pref_aic, aes(x = reorder(Model, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  theme_light() +
  labs(x = "Covariates included in model",
       y = "Frequency of model having lowest AIC for each protein") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylim(0, NA)  # set lower lim as 0

ggsave(plot = p_aic, paste("figures/DE/freq_aic_per_protein_pre_post.png", sep = ""))

# calculate to difference of not having covars/ factors
aic_df <- as.data.frame(aics) |>
	mutate(none_minus_AFWT = aics$IC[, "design_none"] - aics$IC[, "design_AFWT"],
		none_minus_FWT = aics$IC[, "design_none"] - aics$IC[, "design_FWT"],
		none_minus_AWT = aics$IC[, "design_none"] - aics$IC[, "design_AWT"],
		none_minus_AFT = aics$IC[, "design_none"] - aics$IC[, "design_AFT"],
		none_minus_AFW = aics$IC[, "design_none"] - aics$IC[, "design_AFW"],
		none_minus_A = aics$IC[, "design_none"] - aics$IC[, "design_A"],
		none_minus_F = aics$IC[, "design_none"] - aics$IC[, "design_F"],
		none_minus_W = aics$IC[, "design_none"] - aics$IC[, "design_W"],
		none_minus_T = aics$IC[, "design_none"] - aics$IC[, "design_T"],
		none_minus_AF = aics$IC[, "design_none"] - aics$IC[, "design_AF"],
		none_minus_AW = aics$IC[, "design_none"] - aics$IC[, "design_AW"],
		none_minus_AT = aics$IC[, "design_none"] - aics$IC[, "design_AT"],
		none_minus_FW = aics$IC[, "design_none"] - aics$IC[, "design_FW"],
		none_minus_FT = aics$IC[, "design_none"] - aics$IC[, "design_FT"],
		none_minus_WT = aics$IC[, "design_none"] - aics$IC[, "design_WT"],
	)

# calc diff in AIC vs no covars
aic_df <- bind_rows(aic_df) |>
  pivot_longer(
    cols = contains("_minus_"),
    names_to = "model_comparison",
    values_to = "aic_diff"
  )

aic_df$model_comparison <- gsub("none_minus_", "", aic_df$model_comparison)

# Plot AIC differences for each model comparison
bxplt_aic_diff <- aic_df |>
  ggplot(aes(temp, x = reorder(model_comparison, -aic_diff, median), y = aic_diff)) +
  geom_boxplot(alpha = 0, size = 0.5) +
  geom_violin(alpha = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0)) +
  labs(x = "Covariates included in model",
       y = "AIC difference vs no covariates (> 0 is better fit model than no covariates)")

ggsave(plot = bxplt_aic_diff, paste("figures/DE/aic_model_diff_pre_post.png", sep = ""))

# all models < 0 diff to no covariates. Therefore, best fit model with AIC is using no covariates

# plot AIC overall
raw_aics_col_df <- subset(as.data.frame(aics), select = -c(pref, criterion))
raw_aics_df <- stack(raw_aics_col_df)
raw_aics_df$ind <- gsub("IC.design_", "", raw_aics_df$ind)
names(raw_aics_df)[names(raw_aics_df) == "ind"] <- "Model"

bxplt_aic_raw <- raw_aics_df |>
  ggplot(aes(temp, x = reorder(Model, -values, median), y = values)) +
  geom_boxplot(alpha = 0, size = 0.5) +
  geom_violin(alpha = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0)) +
  labs(x = "Model",
       y = "AIC")

ggsave(plot = bxplt_aic_raw, paste("figures/DE/aic_model_pre_post.png", sep = ""))

meds_aic <- t(as.data.frame(lapply(raw_aics_col_df, median)))

# too difficult to see differences in raw AIC, better to use difference in AIC between models.

# run selected model for pre vs post. i.e. Using no covariates
design <- model.matrix(~0 + imp_df$Cohort)

# Rename the columns of the design matrix
colnames(design)[1:2] <- c("ETI_post", "ETI_pre")
colnames(design) <- gsub("imp_df\\$", "", colnames(design))

contrast <- makeContrasts(ETI_post - ETI_pre, levels = design)

check_data_order(meta_pre_post_df, pro_pre_post_df)

# ID as random effect.
# low correlation, but safer to include in model
corfit <- duplicateCorrelation(pro_pre_post_df, design, block = imp_df$ID)

# Fit model
fit <- lmFit(pro_pre_post_df, design, block = imp_df$ID,
			correlation = corfit$consensus.correlation)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Number of differentially expressed proteins
summary(decideTests(fit2))

tt_pre_post <- topTable(fit2, coef = NULL, n = Inf)
tt_pre_post$protein <- pro_df[match(rownames(tt_pre_post), rownames(pro_df)), "PG.Genes"]
write.csv(tt_pre_post, "figures/DE/topranked_pre_post.csv")

# RUN HEALTHY VS CF ANALYSIS (INC CF CONTROL)

# keep only pre, healthy and control
pro_cf_healthy_df <- pro_df %>% select((contains("healthy") | contains("ETI_pre") | contains("control_1st")))
meta_cf_healthy_df <- meta_df[which(meta_df["Cohort"] == "healthy" | meta_df["Cohort"] == "ETI pre" |
			meta_df["Cohort"] == "control 1st"), ]

# missing meta_data, and already have far more CF than healthy.
pro_cf_healthy_df <- pro_cf_healthy_df %>% select(-contains("P1_A6_28uL_2586_V1_control_1st"))

levels(meta_cf_healthy_df$Cohort)[match("ETI pre", levels(meta_df$Cohort))] <- "CF"
levels(meta_cf_healthy_df$Cohort)[match("control 1st", levels(meta_df$Cohort))] <- "CF"

# drop unused levels
meta_cf_healthy_df <- droplevels(meta_cf_healthy_df)

table(meta_cf_healthy_df$Cohort)

# check samples align in data and meta_data
meta_cf_healthy_df <- meta_cf_healthy_df[order(meta_cf_healthy_df$Sample), ]
pro_cf_healthy_df <- pro_cf_healthy_df[, order(colnames(pro_cf_healthy_df))]

check_data_order(meta_cf_healthy_df, pro_cf_healthy_df)

# don't do any imputation as no meta data for healthy and hasn't made a difference for prev model
# run selected model. Using no covariates
design <- model.matrix(~0 + meta_cf_healthy_df$Cohort)

# Rename the columns of the design matrix
colnames(design)[1:2] <- c("CF", "healthy")

contrast <- makeContrasts(CF - healthy, levels = design)

# Fit model
fit <- lmFit(pro_cf_healthy_df, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Number of differentially expressed proteins
summary(decideTests(fit2))

tt_healthy_cf <- topTable(fit2, coef = NULL, n = Inf)
tt_healthy_cf$protein <- pro_df[match(rownames(tt_healthy_cf), rownames(pro_df)), "PG.Genes"]
write.csv(tt_pre_post, "figures/DE/topranked_healthy_cf.csv")

# RUN HEALTHY VS CF ANALYSIS (NOT INC CF CONTROL)

# keep only pre, healthy and control
pro_eti_healthy_df <- pro_df %>% select((contains("healthy") | contains("ETI_pre")))
meta_eti_healthy_df <- meta_df[which(meta_df["Cohort"] == "healthy" | meta_df["Cohort"] == "ETI pre"), ]

levels(meta_eti_healthy_df$Cohort)[match("ETI pre", levels(meta_df$Cohort))] <- "CF"

# drop unused levels
meta_eti_healthy_df <- droplevels(meta_eti_healthy_df)

table(meta_eti_healthy_df$Cohort)

# check samples align in data and meta_data
meta_eti_healthy_df <- meta_eti_healthy_df[order(meta_eti_healthy_df$Sample), ]
pro_eti_healthy_df <- pro_eti_healthy_df[, order(colnames(pro_eti_healthy_df))]

check_data_order(meta_eti_healthy_df, pro_eti_healthy_df)

# don't do any imputation as no meta data for healthy and hasn't made a difference for prev model
# run selected model. Using no covariates
design <- model.matrix(~0 + meta_eti_healthy_df$Cohort)

# Rename the columns of the design matrix
colnames(design)[1:2] <- c("CF", "healthy")

contrast <- makeContrasts(CF - healthy, levels = design)

# Fit model
fit <- lmFit(pro_eti_healthy_df, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Number of differentially expressed proteins
summary(decideTests(fit2))
tt_healthy_pre <- topTable(fit2, coef = NULL, n = Inf)
tt_healthy_pre$protein <- pro_df[match(rownames(tt_healthy_pre), rownames(pro_df)), "PG.Genes"]
write.csv(tt_healthy_pre, "figures/DE/topranked_healthy_eti.csv")

# healthy 8 sig diff. row = 23, 64, 319, 202, 281, 57, 136, 91

# GET SIG DIFF P's ON PRE VS HEALTHY AND ADJUST ONLY THOSE ON ETI PRE VS POST

# p adjust on sig diff between healthy vs pre
sig_healthy_pre_rows <- rownames(tt_healthy_pre[tt_healthy_pre$adj.P.Val < 0.05, ])

# sig diff of pre vs healthy, in pre vs post
tt_pre_post_selected <- tt_pre_post[match(sig_healthy_pre_rows, rownames(tt_pre_post)),]

# adjust selected only
tt_pre_post_selected$sel_adjust <-
	p.adjust(tt_pre_post_selected$P.Value, method = "BH", n = nrow(tt_pre_post_selected))

# still nothing significant post selection

# RUN POST vs HEALTHY

# keep only post healthy
pro_post_healthy_df <- pro_df %>% select((contains("healthy") | contains("ETI_post")))
meta_post_healthy_df <- meta_df[which(meta_df["Cohort"] == "healthy" | meta_df["Cohort"] == "ETI post"), ]

levels(meta_post_healthy_df$Cohort)[match("ETI post", levels(meta_df$Cohort))] <- "post_ETI"

# drop unused levels
meta_post_healthy_df <- droplevels(meta_post_healthy_df)

table(meta_post_healthy_df$Cohort)

# check samples align in data and meta_data
meta_post_healthy_df <- meta_post_healthy_df[order(meta_post_healthy_df$Sample), ]
pro_post_healthy_df <- pro_post_healthy_df[, order(colnames(pro_post_healthy_df))]

check_data_order(meta_post_healthy_df, pro_post_healthy_df)

# don't do any imputation as no meta data for healthy and hasn't made a difference for prev model
# run selected model. Using no covariates
design <- model.matrix(~0 + meta_post_healthy_df$Cohort)

# Rename the columns of the design matrix
colnames(design)[1:2] <- c("post_ETI", "healthy")

contrast <- makeContrasts(post_ETI - healthy, levels = design)

# Fit model
fit <- lmFit(pro_post_healthy_df, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Number of differentially expressed proteins
summary(decideTests(fit2))
tt_healthy_post <- topTable(fit2, coef = NULL, n = Inf)
tt_healthy_post$protein <- pro_df[match(rownames(tt_healthy_post), rownames(pro_df)), "PG.Genes"]
write.csv(tt_healthy_post, "figures/DE/topranked_healthy_post.csv")
