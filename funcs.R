library(ggplot2)
library(tidyverse)
library(pheatmap)
library(missForest)
library(pvca)
library(corrr)
library(sva)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(stringr)
library(vsn)
library(MultiBaC)
library(variancePartition)
library(ggpubr)

plot_missingness <- function(df) {

  # plot barplot of % missingness of each protein
  missing_percentage <- rowSums(is.na(df)) / ncol(df) * 100

  plot_data <- data.frame(
    row_index = factor(seq_len(nrow(df)), levels = seq_len(nrow(df))),
    missing_percentage = missing_percentage
  )

  # Order rows by missing percentage
  plot_data <- plot_data[order(missing_percentage, decreasing = TRUE), ]

  # Reorder levels of the factor variable
  plot_data$row_index <- factor(plot_data$row_index, levels = plot_data$row_index)

  ggplot(plot_data, aes(x = row_index, y = missing_percentage)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
    labs(
      title = "Rows with the Most Missing Data",
      x = "Proteins",
      y = "Percentage of Missing Data"
    ) +
    theme_classic() +
    theme(
      panel.grid.minor = element_line(color = "gray", linetype = "dotted"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}


gg_boxplot <- function(df, x_lab, y_lab) {

  m_df <- melt(df)

  ggplot(m_df, aes(x = variable, y = value)) +
    geom_boxplot() +
    labs(x = x_lab, y = y_lab) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


missingness_heatmap <- function(df, labels, f_name) {

  # Replace NaN values with 0
  heatmap_df <- t(df)

  heatmap_df[heatmap_df == "NaN"] <- 0
  heatmap_df[is.na(heatmap_df)] <- 0

  pheatmap(
    heatmap_df,
    labels_row = labels,
    width = 20,
    height = 20,
    filename = f_name
  )
}


missing_per_sample_bar <- function(df) {

  missing_counts <- colSums(is.na(df))

  plot_data <- data.frame(
    Column = names(missing_counts),
    MissingCount = missing_counts
  )

  # Create a bar plot
  ggplot(plot_data, aes(x = Column, y = MissingCount)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    labs(title = "Number of Missing Values per Sample", x = "Column", y = "Missing Count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

}


CBF_PCA <- function(data, groups, useLabels = F, labels = "",
      pcs = c(1, 2), type = "scores", scale = T, legendName = "Treatment") {

  # CBF PCA function taken from R for data science course material.

  colores <- rainbow(length(unique(groups)))

  if (scale) { 
    pc <- prcomp(data, scale = T)
  } else { 
    pc <- prcomp(data) 
  }

  if (type == "scores") {
    if (useLabels & length(labels) != nrow(data)) {
      print("Warning: The labels not given or given incorrectly. Using rownames.") 
      labels <- rownames(data) 
    }

    pcdf <- data.frame(pc1 = pc$x[, pcs[1]], pc2 = pc$x[, pcs[2]])

    if (useLabels) 
      pcdf$labels <- labels

    perc_accounted <- summary(pc)$importance[2, pcs] * 100

    .e <- environment()
    p <- ggplot(data = pcdf, aes(x = pc1, y = pc2), environment = .e) +
      geom_point(aes(fill = groups), colour = "black", size = 5.5,
                 pch = 21) + scale_fill_manual(values = colores, name = legendName)

    if (useLabels)
      p <- p + geom_text_repel(aes(label = labels))

    p <- p + xlab(paste("PC", pcs[1], " (", round(perc_accounted[1],
                  2), "%)", sep = "")) + ylab(paste("PC", pcs[2], " (", round(perc_accounted[2],
                  2), "%)", sep = "")) + theme_bw(base_size = 20) + theme(legend.position = "bottom")

    p

  } else if (type == "loadings") {
 
    if (useLabels & length(labels) != nrow(pc$rotation)) { 
      print("Warning: loadings labels not given or given incorrectly. Using the column names.") 
      labels <- colnames(data) 

    }

    pcdf <- data.frame(load1 = pc$rotation[, pcs[1]], load2 = pc$rotation[,
            pcs[2]], var = labels)

    label_offset_x <- 0.035 * (range(pcdf$load1)[2] - range(pcdf$load1)[1]) 
    label_offset_y <- 0.035 * (range(pcdf$load2)[2] - range(pcdf$load2)[1])

    .e <- environment()

    p <- ggplot(data = pcdf, aes(x = load1, y = load2), environment = .e) + geom_point()

    if (useLabels)

      p <- p + geom_text_repel(aes(x = load1, y = load2), label = labels)

    p <- p + xlab(paste("Loadings for PC", pcs[1], sep = "")) + ylab(paste("Loadings for PC",
        pcs[2], sep = "")) + ggtitle("PCA loadings plot") + theme_bw(base_size = 20) 
    p

  } else if (type == "varAcc") {
    perc_accounted <- (pc$sdev/sum(pc$sdev) * 100) 
    perc_with_cumsum <- data.frame(pc = as.factor(1:length(perc_accounted)),
          perc_acc = perc_accounted, perc_cumsum = cumsum(perc_accounted))
    p <- ggplot(data = perc_with_cumsum, aes(x = pc, y = perc_cumsum)) +
      geom_bar(stat = "identity", col = "black", fill = "white") +
      geom_hline(yintercept = 95, col = "red") + geom_hline(yintercept = 0,
      col = "black") + xlab("PC") + ylab("% Variance") + ggtitle("% Variance accounted for by principle components") 

    print(p)

  } else {
    cat(sprintf("\nError: no type %s", type))

  }

  return(p)
}

CBFpvcaFunction <- function(data,phenotypedata){

  # taken from R for data science course material
  #The only thing the function now assumes is that your phenotype table 
  #is in the same order as your columns of your data.

  #packages <- c("pvca", "ggplot2", "Biobase") 
  #lapply(packages, library, character.only = TRUE, quietly = TRUE)

  if(length(seq_len(nrow(phenotypedata))) != length(seq_len(nrow(data)))) {
    warning(paste0("Data needs to be a wide", "\n", "Transposing data", "\n"))
    data <- t(data)
  }

  if(!identical(rownames(phenotypedata), rownames(data))) {
    warning(paste0("Rownames in phenotypedata do not match the rownames of data", "\n"))
    warning(paste0("Matching rownames in phenotypedata to data", "\n"))

      if(nrow(phenotypedata) != nrow(data)) {
        stop("Phenotypedata and data have different amounts of rows")
      } else {
    rownames(phenotypedata) <- rownames(data)
      }
  }

  data1_pheno_formatted<-new("AnnotatedDataFrame",data=phenotypedata)

  data_forPVCA<-ExpressionSet(assayData=t(data),phenoData=data1_pheno_formatted)

  pvcaObj_data <- pvcaBatchAssess(abatch = data_forPVCA, batch.factors = colnames(data1_pheno_formatted@data), threshold = 0.7)

  data_pvca_res<-data.frame(as.data.frame(pvcaObj_data$label),t(as.data.frame(pvcaObj_data$dat)))
  
  colnames(data_pvca_res)<-c("effect","variance")


  p<-ggplot((data_pvca_res[-nrow(data_pvca_res),]), aes(x= effect, y = variance)) + 
    geom_bar(stat = 'identity',  position = 'dodge', col ='transparent')+ 
    scale_fill_discrete(guide = 'none') + 
    theme_bw(base_size = 16)+
    theme(plot.title =  element_text(hjust=0.5),axis.text.x=element_text(angle=45,hjust=1))+
    labs(x = 'Effects', y = 'Weighted average proportion variance')+
    ggtitle("PVCA estimation bar chart corrected data")

  return(p)

}


pca_var_imp_plot <- function(df, og_df, colour_deplete = TRUE, top = 40, selected_p = FALSE) {
  # df - current df, data only
  # og_df - original data with protein description column
  # selected_p - FALSE or vector of selected proteins to plot

  # adapted from
  # https://github.com/kassambara/factoextra/blob/master/R/fviz_contrib.R

  # use protein names
  replace_name_list <- setNames(og_df$PG.ProteinDescriptions, rownames(og_df))

  existing_columns <- intersect(colnames(df), names(replace_name_list))

  df <- df %>%
    rename_at(vars(existing_columns), ~ unlist(replace_name_list)[.])

  # defaults
  fill = "steelblue"
  color = "steelblue"
  sort.val = "desc"
  top = top
  ggtheme = theme_bw()
  title = "Contribution of variables to Dim-1-2"
  xtickslab.rt = 45

  res <- pca_var_imp_df(df, choice = "var",
                         axes=1:2, fill=fill, color = color, top = top,
                         xtickslab.rt = 45, ggtheme = ggtheme)

  df_ <- res$df
  theo_contrib <- res$theo_contrib

  if (selected_p != FALSE){
    # then keep only selected proteins

    df_ <- df_[df_$name %in% selected_p, ]

  }

  # colour depleted proteins different colour, by adding new column
  if (colour_deplete == TRUE) {
    depl_proteins = c("immunoglobulin",
                  "albumin",
                  "transferrin",
                  "haptoglobin",
                  "fibrinogen",
                  "apolipoprotein",
                  "Alpha-2-macroglobulin",
                  "Alpha-1-antitrypsin",
                  "Alpha-1-acid glycoprotein"
    )

    df_ <- df_ %>%
      rowwise %>%
      mutate(., Depleted = ifelse(any(grepl(paste(depl_proteins, collapse='|'), name, ignore.case = TRUE)), "Depleted", "Not Depleted"))

    fill <- "Depleted"
  }
  
  p <- ggpubr::ggbarplot(df_, x = "name", y = "contrib", fill = fill, color = color,
                         sort.val = sort.val, top = top,
                         main = title, xlab = FALSE, ylab ="Contributions (%)",
                         xtickslab.rt = xtickslab.rt, ggtheme = ggtheme,
                         sort.by.groups = FALSE
                         )+
    geom_hline(yintercept = theo_contrib, linetype = 2, color = "red")

  return(p)
}

pca_var_imp_df <- function(df, choice = c("row", "col", "var", "ind", "quanti.var", "quali.var", "group", "partial.axes"),
                         axes=1, fill="steelblue", color = "steelblue", 
                         sort.val = c("desc", "asc", "none"), top = Inf,
                         xtickslab.rt = 45, ggtheme = theme_minimal(), ...) {

  # adapted from
  # https://github.com/kassambara/factoextra/blob/master/R/fviz_contrib.R

  X <- PCA(df)

  dd <- facto_summarize(X, element = choice, result = "contrib", axes = axes)
  contrib <- dd$contrib
  names(contrib) <-rownames(dd)

  # expected Average contribution
  theo_contrib <- 100/length(contrib)
  if(length(axes) > 1) {
    # Adjust variable contributions by the Dimension eigenvalues
    eig <- get_eigenvalue(X)[axes,1]
    theo_contrib <- sum(theo_contrib*eig)/sum(eig)
  }
  df <- data.frame(name = factor(names(contrib), levels = names(contrib)), contrib = contrib, stringsAsFactors = TRUE)

  return(list("df" = df, "theo_contrib" = theo_contrib))
}


plot_selected_boxplot <- function(df_, og_df, meta_df, selected_p){

  # requires df_ to be matrix with rownames and data cols only
  #ie.
  #mat <- as.matrix(df[, data_cols])
  #rownames(mat) <- rownames(df)

  if (length(rownames(df_)) == 0){
    stop("requires df_ to be matrix with rownames and data cols only")
  }

  if ("PG.ProteinDescriptions" %in% colnames(df_)){
    stop("requires df_ to be matrix with rownames and data cols only")
  }

  # change rownames to protein description
  replace_name_list <- setNames(og_df$PG.ProteinDescriptions, rownames(og_df))
  existing_rows <- intersect(rownames(df_), names(replace_name_list))
  rownames(df_) <- unname(replace_name_list[rownames(df_)])

  # keep only selected rows
  selected_df <- df_[rownames(df_) %in% selected_p, ]

  if (nrow(selected_df) == 0){
    print("None of the selected proteins are in the dataframe")
    # then leave function
    return()
  }

  # add meta data to df
  sel_m_df <- melt(selected_df)
  sel_m_df$Sample <- sel_m_df$Var2
  result_df <- sel_m_df %>%
    left_join(meta_df, by = "Sample")

  p <- ggplot(result_df, aes(x = Cohort, y = value, fill = Cohort)) +
  geom_boxplot() +
  facet_wrap(~Var1, scales = "free_y") +
  theme_bw() +
  labs(title = "Boxplots of Selected Proteins Across Conditions",
       x = "Condition",
       y = "Expression")
  
  return(p)
}


remove_depleted <- function(df){
  # remove depleted proteins

  depl_proteins = c("immunoglobulin",
                  "albumin",
                  "transferrin",
                  "haptoglobin",
                  "fibrinogen",
                  "apolipoprotein",
                  "Alpha-2-macroglobulin",
                  "Alpha-1-antitrypsin",
                  "Alpha-1-acid glycoprotein"
    )

    df_ <- df %>%
      rownames_to_column(var="row") %>%  # force tidyverse to preserve row names
      filter(!grepl(paste(depl_proteins, collapse='|'), PG.ProteinDescriptions, ignore.case = TRUE)) %>%
      column_to_rownames(var="row")  # force tidyverse to preserve row names

    return(df_)
}

df_to_matrix <- function(df_, d_cols){
  # keep rownames and remove data cols, for easier boxplotting

  df_mat <- as.matrix(df_[, d_cols])
  rownames(df_mat) <- rownames(df_)

  return(df_mat)
}


plot_deplete <- function(df_){

  # plot depleted proteins only
  # https://www.thermofisher.com/document-connect/document-connect.html?url=https://assets.thermofisher.com/TFS-Assets%2FLSG%2Fmanuals%2FMAN0017285_HighSlctTop_14AbundProteinDepletResin_PI.pdf
  depl_proteins = c("immunoglobulin",
                    "albumin",
                    "transferrin",
                    "haptoglobin",
                    "fibrinogen",
                    "apolipoprotein",
                    "Alpha-2-macroglobulin",
                    "Alpha-1-antitrypsin",
                    "Alpha-1-acid glycoprotein"
  )

  df_ig <- df_ %>%
    filter(grepl(paste(depl_proteins, collapse='|'), PG.ProteinDescriptions, ignore.case = TRUE))

  df_ig_m <- melt(df_ig)

  colnames(df_ig_m)[colnames(df_ig_m) == "variable"] <- "Sample"

  df_ig_m <- merge(df_ig_m, meta_df, by = "Sample", all.x = TRUE)

  df_ig_m$cond_rank <- dense_rank(df_ig_m$Cohort)

  p <- ggplot(df_ig_m, aes(x = reorder(Sample, cond_rank), y = value, fill = Cohort))  +
    geom_boxplot() +
    labs(x = "Samples", y = "log2") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_brewer(palette = "Accent")

  return(p)
}


plot_all_proteins <- function(df_, meta_df_, grouped){

  # plots all proteins in df, grouped by cohort or plate
  # either grouped and coloured by "cohort" or "plate"

  df_m <- melt(df_)

  colnames(df_m)[colnames(df_m) == "variable"] <- "Sample"

  # add meta
  df_m <- merge(df_m, meta_df_, by = "Sample", all.x = TRUE)

  if (grouped == "cohort"){

    # force columns to be ordered by condition
    df_m$Cohort <- factor(df_m$Cohort, levels = unique(df_m$Cohort))

    df_m$cond_rank <- dense_rank(df_m$Cohort)

    fill_ <- df_m$Cohort

  } else if (grouped == "plate"){

    # force columns to be ordered by plate
    df_m$Plate <- factor(df_m$Plate, levels = unique(df_m$Plate))

    df_m$cond_rank <- dense_rank(df_m$Plate)

    fill_ <- df_m$Plate

  } else {
    stop("'grouped' must be 'cohort' or 'plate'")
  }

  p <- ggplot(df_m, aes(x = reorder(Sample, cond_rank), y = value, fill = fill_))  +
    geom_boxplot() +
    labs(x = "Samples", y = "log2") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_brewer(palette = "Accent")

  return(p)
}

vsn_one_batch <- function(df, batch){

  # df - all data
  # perform vsn on one batches only. return normalised df.

  # df of that batch only
  p_df <- df[, grep(paste("^", batch, sep = ""), names(df))]

  p_mat <- df_to_matrix(p_df)

  # replace nan's so vsn is happy.
  p_mat[is.nan(p_mat) | is.infinite(p_mat)] <- NA

  # run vsn
  inv_log <- 2 ^ p_mat  # inverse log, as vsn log2 automatically
  vsn_res <- vsn2(inv_log)
  norm_mat <- exprs(vsn_res)

  return(norm_mat)
}


vsn_within_batch <- function(df){

  # normalise per batch
  norm_p1 <- vsn_one_batch(df, "P1")
  norm_p2 <- vsn_one_batch(df, "P2")
  norm_p3 <- vsn_one_batch(df, "P3")

  all_norm <- cbind(norm_p1, norm_p2, norm_p3)

  return(all_norm)
}


idx_low_cv_p <- function(df_, pc = 0.05){
  # return idx of the pc% lowest co var proteins in the df

  cv <- apply(df_, 1, function(x) sd(x) / mean(x) * 100)

  # keep idx of protein and sort in asc order
  cv_df <- as.data.frame(x = cv)
  cv_df$protein_idx <- rownames(cv_df)
  cv_df <- cv_df[order(cv), ]

  # Take pc% smallest value
  top_cv_df <- head(cv_df, n = pc * nrow(cv_df))

  protein_idxs <- top_cv_df$protein_idx

  return(protein_idxs)
}

plot_figures <- function(df_, og_df, meta_df_, data_cols, f_name_prefix, selected_p_){
  # generate and save all figures

  # check folder created
  subdir_path <- sub("/[^/]*$", "/", f_name_prefix)  # keep only path
  dir.create(file.path(subdir_path), showWarnings = FALSE)

  # impute using missforest for PCA
  imp_res <- missForest(t(as.matrix(df_[, data_cols])))
  imp_pca <- as.data.frame(imp_res$ximp)

  # check pca
  CBF_PCA(data = imp_pca,
          groups = meta_df_$Plate,
          useLabels = F,
          labels = "",
          pcs = c(1, 2),
          type = "scores",
          scale = T,
          legendName = "Plate")
  ggsave(paste(f_name_prefix, "pca_plate.png", sep = ""))

  CBF_PCA(data = imp_pca,
          groups = meta_df_$Cohort,
          useLabels = F,
          labels = "",
          pcs = c(1, 2),
          type = "scores",
          scale = T,
          legendName = "Plate")
  ggsave(paste(f_name_prefix, "pca_cohort.png", sep = ""))

  pheatmap(
      imp_pca,
      labels_row = meta_df_$Cohort,
      width = 20,
      height = 20,
      filename = paste(f_name_prefix, "heatmap_cohort.png", sep = "")
    )
  
  pheatmap(
      imp_pca,
      labels_row = meta_df_$Plate,
      width = 20,
      height = 20,
      filename = paste(f_name_prefix, "heatmap_plate.png", sep = "")
    )

  pca_var_imp_plot(imp_pca, og_df, colour_deplete = TRUE, top = 40)
  ggsave(paste(f_name_prefix, "pca_var_imp.png", sep = ""), width = 20)

  # boxplots
  # selected proteins to boxplot
  selected_p <- selected_p_$selected_p
  most_diff <- selected_p_$most_diff
  least_diff <- selected_p_$least_diff
  low_coeff_var <- selected_p_$low_coeff_var

  m_df <- df_to_matrix(df_[, data_cols])

  # plot selected proteins from main.R
  plot_selected_boxplot(m_df, og_df, meta_df_, selected_p)
  ggsave(paste(f_name_prefix, "boxplot_selected.png", sep = ""))

  plot_selected_boxplot(m_df, og_df, meta_df_, most_diff)
  ggsave(paste(f_name_prefix, "boxplot_selected_most.png", sep = ""))

  plot_selected_boxplot(m_df, og_df, meta_df_, least_diff)
  ggsave(paste(f_name_prefix, "boxplot_selected_least.png", sep = ""))

  plot_selected_boxplot(m_df, og_df, meta_df_, low_coeff_var)
  ggsave(paste(f_name_prefix, "boxplot_selected_low_coeff_var.png", sep = ""))

  plot_all_proteins(df_, meta_df_, grouped = "cohort")
  ggsave(paste(f_name_prefix, "boxplot_cohort.png", sep = ""), width = 20)

  plot_all_proteins(df_, meta_df_, grouped = "plate")
  ggsave(paste(f_name_prefix, "boxplot_plate.png", sep = ""), width = 20)

  plot_deplete(df_)
  ggsave(paste(f_name_prefix, "deplete.png", sep = ""), width = 20)

  # Plot variance importance.
  # get row names for plotting
  df_for_varpart <- as.data.frame(t(imp_pca))  # transpose again after missforest

  # add non-data cols
  sel_idx_only <- og_df[rownames(df_for_varpart), ]
  df_for_varpart$PG.ProteinDescriptions <- sel_idx_only$PG.ProteinDescriptions
  df_for_varpart$PG.Genes <- sel_idx_only$PG.Genes

  # cut rownames so fit in plot
  rownames(df_for_varpart) <-  paste0(substr(df_for_varpart$PG.ProteinDescriptions, 1, 20),
                          " - ",
                          df_for_varpart$PG.Genes)

  # recommended cat vars are modelled as random effects
  form <- ~(1|Plate) + (1|Cohort)
  var_part <- fitExtractVarPartModel(exprObj = df_to_matrix(df_for_varpart, data_cols),
                          formula = form,
                          data = meta_df_)

  vp <- sortCols(var_part, last = "Residuals")

  # top 30 proteins sorted by cohort variance %
  vp <- vp[order(-vp$Cohort), ]

  # top 30 proteins
  vp_bars <- plotPercentBars(vp[1:30, ]) + theme(legend.position = "bottom")  #top 30 so we can see some detail
  vp_violin <- plotVarPart(vp) 

  ggarrange(vp_bars, vp_violin, 
            labels = c("A", "B"),
            ncol = 2, nrow = 1)

  ggsave(paste(f_name_prefix, "var_part.png", sep = ""), width = 20)

  return(imp_pca)
}

run_pipelines <- function(pipelines, df_, og_df, meta_df_, data_cols, meta_cols, dir_path, selected_p_){
  # run the pipelines
  # df is filtered and log2 transformed.
  # og_df the original dataframe to plot against

  # df to record the iqr of QC's across all methods
  iqr_df <- data.frame(Method = character(), IQR = double())

  # run processing pipelines one at a time
  for (i in 1:length(pipelines)){

    p <- pipelines[[i]]

    # reset df at the start of each pipeline
    df_to_process <- df_

    # plot and write csv's of filtered and log transformed df to overall dir
    mf_data <- plot_figures(df_to_process, og_df, meta_df_, data_cols, f_name_prefix = paste(dir_path, "/log_filter/", sep=""), selected_p_)  
    write.csv(x = df_to_process, file = paste(dir_path, "/log_filter/", "df.csv", sep = ""))

    # if not depleted during pipeline, then also write csv of dfs with depleted removed.
    if (grepl("depl", p, fixed = TRUE) == FALSE){
      depl_df_ <- remove_depleted(df_to_process)
      write.csv(x = depl_df_, file = paste(dir_path, "/log_filter/", "rm_depl_df.csv", sep = ""))
    }

    # write misforest imputed data to csv
    mf_df <- cbind(df_to_process[, meta_cols], t(as.matrix(mf_data)))
    write.csv(x = mf_df, file = paste(dir_path, "/log_filter/", "mf_df.csv", sep = ""))

    # if not depleted during pipeline, then also write csv of dfs with depleted removed.
    if (grepl("depl", p, fixed = TRUE) == FALSE){
      depl_mf_df_ <- remove_depleted(mf_df)
      write.csv(x = depl_mf_df_, file = paste(dir_path, "/log_filter/", "rm_depl_mf_df.csv", sep = ""))
    }

    # get one pipeline
    all_processes <- str_split(p, "_")

    for (j in 1:length(all_processes[[1]])){

      process <- all_processes[[1]][[j]]

      # put plot in correct folder, formatted as "figures/(pipeline)/(process)_(file)"
      f_name_prefix <- paste(dir_path, "/", p, "/", process, "_", sep="")

      # run each process
      if (process == "combat"){

        # run combat
        combat_df <- ComBat(dat = df_to_process[, data_cols],
                  batch = meta_df_$Plate,
                  mod = NULL,
                  par.prior = TRUE,
                  prior.plots = TRUE)

        # add back non data cols
        df_to_process <- cbind(combat_df, df_to_process[, meta_cols])

      } else if (process == "withinvsn"){
        # vsn within batches ie. if performed before batch correction
        vsn_w_df <- vsn_within_batch(df_to_process[, data_cols])

        # add back non data cols
        df_to_process <- cbind(vsn_w_df, df_to_process[, meta_cols])

      } else if (process == "vsn"){
        # run vsn but over all batches ie. if performed after batch correction

        vsn_mat <- df_to_matrix(df_to_process[, data_cols])

        vsn_mat[is.nan(vsn_mat) | is.infinite(vsn_mat)] <- NA
        inv_log <- 2 ^ vsn_mat  # inverse log, as vsn log2 automatically
        vsn_normalized_data <- vsn2(inv_log)

        vsn_df <- exprs(vsn_normalized_data)

        # add back non data cols
        df_to_process <- cbind(vsn_df, df_to_process[, meta_cols])

      } else if (process == "depl"){
        # remove depleted proteins
        df_to_process <- remove_depleted(df_to_process)

      } else if (process == "hk"){
        # correct batches based on stable proteins with low coeff var.
        norm_df <- norm_to_low_coeff(df_to_process, og_df, pc = 0.05, max_p = 5, meta_df_, paste(dir_path, "/", p, sep=""), data_cols)
        # add back non data cols
        df_to_process <- cbind(norm_df, df_to_process[, meta_cols])

      } else if (process == "arsyn"){
        # batch correction using arsyn

        # divide df into an experiment per plate
        p1_df <- df_to_process[, grep("^P1", names(df_to_process))]
        p1_meta <- meta_df_[meta_df_$Plate == "P1", ]
        p2_df <- df_to_process[, grep("^P2", names(df_to_process))]
        p2_meta <- meta_df_[meta_df_$Plate == "P2", ]
        p3_df <- df_to_process[, grep("^P3", names(df_to_process))]
        p3_meta <- meta_df_[meta_df_$Plate == "P3", ]

        # requires all conditions to be in each plate, which they aren't, so fudge it.
        # model doesn't converge error as data input is incorrect
        # https://www.bioconductor.org/packages/devel/bioc/manuals/pcaMethods/man/pcaMethods.pdf
        ed_random <- list("p1" = rep(c("H", "C", "E"), length(p1_meta$Cohort))[1:length(p1_meta$Cohort)],
                          "p2" = rep(c("H", "C", "E"), length(p2_meta$Cohort))[1:length(p2_meta$Cohort)],
                          "p3" = rep(c("H", "C", "E"), length(p3_meta$Cohort))[1:length(p3_meta$Cohort)])

        # group control
        p1_r_c <- replace_control(p1_meta$Cohort)
        p2_r_c <- replace_control(p2_meta$Cohort)
        p3_r_c <- replace_control(p3_meta$Cohort)

        ed_group_controls <- list("p1" = p1_r_c,
                                  "p2" = p2_r_c,
                                  "p3" = p3_r_c)

        mbac <- createMbac(inputOmics = list(p1_df, p2_df, p3_df),
                batchFactor = c("p1", "p2", "p3"),
                experimentalDesign = ed_group_controls,
                omicNames = c("Protein"))

        arsyn_res <- ARSyNbac(mbac, batchEstimation = TRUE, filterNoise = FALSE, Interaction = FALSE,
                  Variability = 0.90, beta = 2, modelName = "Model 1", showplot = TRUE)

        df_to_process <- arsyn_res$CorrectedData

      }

      # plots
      mf_data <- plot_figures(df_to_process, og_df, meta_df_, data_cols, f_name_prefix, selected_p_)

      # write csv data
      write.csv(x = df_to_process, file = paste(f_name_prefix, "df.csv", sep = ""))

      # if not depleted during pipeline, then also write csv of dfs with depleted removed.
      if (grepl("depl", p, fixed = TRUE) == FALSE){
        depl_df_ <- remove_depleted(df_to_process)
        write.csv(x = depl_df_, file = paste(f_name_prefix, "rm_depl_df.csv", sep = ""))
      }

      # write misforest imputed data to csv
      mf_df <- cbind(df_to_process[, meta_cols], t(as.matrix(mf_data)))
      write.csv(x = mf_df, file = paste(f_name_prefix, "mf_df.csv", sep = ""))

      # if not depleted during pipeline, then also write csv of dfs with depleted removed.
      if (grepl("depl", p, fixed = TRUE) == FALSE){
        depl_mf_df_ <- remove_depleted(mf_df)
        write.csv(x = depl_mf_df_, file = paste(f_name_prefix, "rm_depl_mf_df.csv", sep = ""))
      }

    }

    # at the end of the pipeline, calculate the IQR of all proteins in QC samples
    all_qc_df <- df_to_process[, grep("QC", names(df_to_process))]
    iqr_values <- apply(all_qc_df, 1, calculate_iqr)
    mean_iqr <- mean(iqr_values)

    iqr_df <- iqr_df %>% add_row(Method = p, IQR = mean_iqr)

    # Plot variance importance.
    # get row names for plotting
    # df_for_varpart <- mf_df
    # rownames(df_for_varpart) <- df_for_varpart$PG.ProteinDescriptions

    # # cut rownames so fit in plot
    # rownames(df_for_varpart) <-  paste0(substr(rownames(df_for_varpart), 1, 20),
    #                         " - ",
    #                         df_for_varpart$PG.Genes)

    # # recommended cat vars are modelled as random effects
    # form <- ~(1|Plate) + (1|Cohort)
    # var_part <- fitExtractVarPartModel(exprObj = df_to_matrix(df_for_varpart, data_cols),
    #                         formula = form,
    #                         data = meta_df_)

    # vp <- sortCols(var_part, last = "Residuals")

    # # top 30 proteins sorted by cohort variance %
    # vp <- vp[order(-vp$Cohort), ]

    # # top 30 proteins
    # vp_bars <- plotPercentBars(vp[1:30, ]) + theme(legend.position = "bottom")  #top 30 so we can see some detail
    # vp_violin <- plotVarPart(vp) 

    # ggarrange(vp_bars, vp_violin, 
    #           labels = c("A", "B"),
    #           ncol = 2, nrow = 1)

    # f_name_prefix <- paste(dir_path, "/", p, "/", process, "_", sep="")
    # ggsave(paste(f_name_prefix, "var_part.png", sep = ""), width = 20)

  }
  return(iqr_df)
}

low_coeff_var_p <- function(df_, pc = 0.05, max_p = 5, data_cols) {
  # find the lowest coefficient of variation proteins across plates
  # Then select the most abundant "max_p" proteins.
  # pc is the top x % of variation to select

  p1_df <- df_[, grep("^P1", names(df_))]
  p2_df <- df_[, grep("^P2", names(df_))]
  p3_df <- df_[, grep("^P3", names(df_))]

  # Rank proteins for each plate in decreasing order of coeff var
  p1_pro_idx <- idx_low_cv_p(p1_df, pc = pc)
  p2_pro_idx <- idx_low_cv_p(p2_df, pc = pc)
  p3_pro_idx <- idx_low_cv_p(p3_df, pc = pc)

  # intersect across all plates
  intersect_idx <- intersect(intersect(p1_pro_idx, p2_pro_idx), p3_pro_idx)

  intersect_df <- df_[intersect_idx, data_cols]

  intersect_df$mean <- rowMeans(intersect_df)

  # rank by mean
  intersect_df <- intersect_df[order(intersect_df$mean, decreasing = TRUE), ]

  # drop proteins with mean < 10, slightly arbitrary threshold to ensure well above LOD
  intersect_df <- intersect_df[intersect_df$mean >= 10, ]

  # take only top 'max_p' number of proteins.
  low_idx_pro_idx <- rownames(intersect_df)[1 : max_p]

  # get protein name from idx
  low_cv_pro_names <- df_[low_idx_pro_idx, "PG.ProteinDescriptions"]

  return(list("names" = low_cv_pro_names,
            "idx" = low_idx_pro_idx))

}


norm_to_low_coeff <- function(df_, og_df, pc = 0.05, max_p = 5, meta_df_, dir_path, data_cols){

  # calculate proteins with lowest coefficient of variation
  low_coeff_var_res <- low_coeff_var_p(df_, pc = pc, max_p = 5, data_cols)
  low_cv_pro_names <- low_coeff_var_res$names
  low_idx_pro_idx <- low_coeff_var_res$idx

  p1_df <- df_[, grep("^P1", names(df_))]
  p2_df <- df_[, grep("^P2", names(df_))]
  p3_df <- df_[, grep("^P3", names(df_))]

  # normalise to low_cv proteins
  # https://nanostring.com/wp-content/uploads/Gene_Expression_Data_Analysis_Guidelines.pdf

  # geometric mean of low_cv_proteins across each plate
  p1_geo_mean <- exp(mean(unlist(log(p1_df[low_idx_pro_idx, ]))))
  p2_geo_mean <- exp(mean(unlist(log(p2_df[low_idx_pro_idx, ]))))
  p3_geo_mean <- exp(mean(unlist(log(p3_df[low_idx_pro_idx, ]))))

  # arithmetic mean of these means across all plates
  mean_all_plates <- mean(p1_geo_mean, p2_geo_mean, p3_geo_mean)

  # plate specific norm factor
  p1_norm_fact <- mean_all_plates / p1_geo_mean
  p2_norm_fact <- mean_all_plates / p2_geo_mean
  p3_norm_fact <- mean_all_plates / p3_geo_mean

  # multiply counts by norm factor
  p1_norm <- p1_df * p1_norm_fact
  p2_norm <- p2_df * p2_norm_fact
  p3_norm <- p3_df * p3_norm_fact

  all_norm_df <- cbind(p1_norm, p2_norm, p3_norm)

  p1_mat <- df_to_matrix(p1_df)
  p2_mat <- df_to_matrix(p2_df)
  p3_mat <- df_to_matrix(p3_df)

  plot_selected_boxplot(p1_mat, og_df, meta_df_, low_cv_pro_names)
  ggsave(paste(dir_path, "/p1_low_cv.png", sep = ""), width = 20)

  plot_selected_boxplot(p2_mat, og_df, meta_df_, low_cv_pro_names)
  ggsave(paste(dir_path, "/p2_low_cv.png", sep = ""), width = 20)

  plot_selected_boxplot(p3_mat, og_df, meta_df_, low_cv_pro_names)
  ggsave(paste(dir_path, "/p3_low_cv.png", sep = ""), width = 20)

  return(all_norm_df)
}

replace_control <- function(my_list) {
  # Replace elements containing "control" with "control"
  control_elements <- grep("control", my_list)
  my_list[control_elements] <- "control"
  return(my_list)
}

calculate_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr_value <- q[2] - q[1]
  return(iqr_value)
}