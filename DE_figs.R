library(pheatmap)
library(tidyverse)
library(reshape2)
library(EnhancedVolcano)
library(ggh4x)
library(fgsea)

source("funcs_DE.R")
source("funcs_preprocess.R")

# load data
pro_with_gene_df <- read.csv("figures/withinvsn_combat/combat_rm_depl_df.csv", row.names = 1)  # inc NA
pro_with_gene_imp_df <- read.csv("figures/withinvsn_combat/combat_rm_depl_mf_df.csv", row.names = 1)  # NA imp
meta_df <- read.csv("cleaned_meta.csv", stringsAsFactors=TRUE)  # written during `model_metadata.R`
DE_pro <- read.csv("figures/DE/topranked_healthy_eti.csv", row.names = 1)  # healthy_ETI_DE
idx_DE_pro <- rownames(DE_pro[DE_pro$adj.P.Val < 0.05, ])

# pre process

# reset rownames on df's, as the idx that matches selected DE proteins has been reset by limma
rownames(pro_with_gene_df) <- NULL
rownames(pro_with_gene_imp_df) <- NULL

# drop ivacaftor samps
iva_samps <- read.csv("ivacaftor_samps.csv")[, 1]
pro_df <- pro_with_gene_df[, !(names(pro_with_gene_df) %in% iva_samps)]
pro_imp_df <- pro_with_gene_imp_df[, !(names(pro_with_gene_imp_df) %in% iva_samps)]

# both 6021 labelled as pre, drop for now.
pro_df <- pro_df %>% select(-contains("6021"))
pro_imp_df <- pro_imp_df %>% select(-contains("6021"))
meta_df <- meta_df[which(meta_df["ID"] != 6021), ]

# keep only pre, post, healthy
pro_df <- pro_df %>% select((contains("ETI_post") | contains("ETI_pre") | contains("healthy")))
pro_imp_df <- pro_imp_df %>% select((contains("ETI_post") | contains("ETI_pre") | contains("healthy")))
meta_df <- meta_df[which(meta_df["Cohort"] == "ETI post" | meta_df["Cohort"] == "ETI pre" | meta_df["Cohort"] == "healthy"), ]

pro_df <- droplevels(pro_df)
pro_imp_df <- droplevels(pro_imp_df)
meta_df <- droplevels(meta_df)

# check samples align in data and meta_data
# sort meta data sample alphabetically
meta_df <- meta_df[order(meta_df$Sample), ]
# sort colnames alphabetically
pro_df <- pro_df[, order(colnames(pro_df))]
pro_imp_df <- pro_imp_df[, order(colnames(pro_imp_df))]

# DE selected proteins only
pro_df_sel <- pro_df[match(idx_DE_pro, rownames(pro_df)), ]
pro_imp_df_sel <- pro_imp_df[match(idx_DE_pro, rownames(pro_imp_df)), ]

# check orders
check_data_order(meta_df, pro_df)
check_data_order(meta_df, pro_imp_df)
check_data_order(meta_df, pro_df_sel)
check_data_order(meta_df, pro_imp_df_sel)

# heatmaps, pre vs post vs healthy

# all

# sort order force, ETI pre, post, healthy
ord_meta_df <- meta_df[order(meta_df$Cohort, decreasing = TRUE), ]
ord_meta_df <- rbind(ord_meta_df[!ord_meta_df$Cohort == "healthy" , ], ord_meta_df[ord_meta_df$Cohort == "healthy", ])

pheatmap(
      pro_imp_df,
	    labels_row = pro_with_gene_df[match(rownames(pro_imp_df), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_col = meta_df$Cohort,
      scale = "row",
      width = 20,
      height = 50,
      filename = "figures/DE/heatmap_pre_post_healthy.png"
    )

# force order
pheatmap(
      pro_imp_df[, match(ord_meta_df$Sample, colnames(pro_imp_df))],
	    labels_row = pro_with_gene_df[match(rownames(pro_imp_df), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_col = ord_meta_df$Cohort,
      cluster_cols = FALSE,
      scale = "row",
      width = 20,
      height = 50,
      filename = "figures/DE/heatmap_pre_post_healthy_no_cluster_coh.png"
    )

pheatmap(
      t(pro_imp_df),
	    labels_col = pro_with_gene_df[match(rownames(pro_imp_df), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_row = meta_df$Cohort,
      scale = "column",
      width = 50,
      height = 20,
      filename = "figures/DE/heatmap_pre_post_healthy_t.png"
    )

# force order
pheatmap(
      t(pro_imp_df[, match(ord_meta_df$Sample, colnames(pro_imp_df))]),
	    labels_col = pro_with_gene_df[match(rownames(pro_imp_df), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_row = ord_meta_df$Cohort,
      cluster_rows = FALSE,
      scale = "column",
      width = 50,
      height = 20,
      filename = "figures/DE/heatmap_pre_post_healthy_t_no_cluster_coh.png"
    )

# DE -  plot only protein's that were DE healthy vs pre.
pheatmap(
      pro_imp_df_sel,
	    labels_row = pro_with_gene_df[match(rownames(pro_imp_df_sel), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_col = meta_df$Cohort,
      scale = "row",
      width = 20,
      height = 10,
      filename = "figures/DE/heatmap_pre_post_healthy_DE.png"
    )

# force order
pheatmap(
      pro_imp_df_sel[, match(ord_meta_df$Sample, colnames(pro_imp_df))],
	    labels_row = pro_with_gene_df[match(rownames(pro_imp_df_sel), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_col = ord_meta_df$Cohort,
      cluster_cols = FALSE,
      scale = "row",
      width = 20,
      height = 10,
      filename = "figures/DE/heatmap_pre_post_healthy_DE_no_cluster_coh.png"
    )

pheatmap(
      t(pro_imp_df_sel),
	  labels_col = pro_with_gene_df[match(rownames(pro_imp_df_sel), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_row = meta_df$Cohort,
      scale = "column",
      width = 10,
      height = 20,
      filename = "figures/DE/heatmap_pre_post_healthy_DE_t.png"
    )

pheatmap(
      t(pro_imp_df_sel[, match(ord_meta_df$Sample, colnames(pro_imp_df))]),
	  labels_col = pro_with_gene_df[match(rownames(pro_imp_df_sel), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_row = ord_meta_df$Cohort,
      cluster_rows = FALSE,
      scale = "column",
      width = 10,
      height = 20,
      filename = "figures/DE/heatmap_pre_post_healthy_DE_t_no_cluster_coh.png"
    )

# mean
agg_pro_imp_df <- t(aggregate(t(pro_imp_df_sel), by = list(meta_df$Cohort), mean))

# tidy
colnames(agg_pro_imp_df) <- agg_pro_imp_df[1, ]
agg_pro_imp_df <- as.data.frame(agg_pro_imp_df[-1, ])
agg_pro_imp_df <- mutate_all(agg_pro_imp_df, function(x) as.numeric(as.character(x)))

pheatmap(
      agg_pro_imp_df,
	    labels_row = pro_with_gene_df[match(rownames(agg_pro_imp_df), rownames(pro_with_gene_df)), "PG.Genes"],
      labels_col = colnames(agg_pro_imp_df),
      scale = "row",
      width = 5,
      height = 5,
      filename = "figures/DE/heatmap_pre_post_healthy_DE_agg.png"
    )


# boxplots selected DE proteins, pre vs post vs healthy

df_for_bxplt <- pro_df_sel
df_for_bxplt$Protein <- DE_pro[match(rownames(df_for_bxplt), rownames(DE_pro)), "protein"]
df_for_bxplt <- melt(df_for_bxplt, by = Protein)

# change col names
df_for_bxplt$variable <- gsub(".*ETI_pre.*", "ETI_pre", df_for_bxplt$variable)
df_for_bxplt$variable <- gsub(".*ETI_post.*", "ETI_post", df_for_bxplt$variable)
df_for_bxplt$variable <- gsub(".*healthy.*", "Healthy", df_for_bxplt$variable)

names(df_for_bxplt)[names(df_for_bxplt) == "variable"] <- "Cohort"

df_for_bxplt$Cohort <- factor(df_for_bxplt$Cohort,
    levels = c("ETI_pre", "ETI_post", "Healthy"), ordered = TRUE)

bxp <- ggplot(df_for_bxplt, aes(x = Cohort, y = value, fill = Cohort)) +
  geom_boxplot() +
  facet_wrap(~Protein, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  labs(x = "Cohort",
       y = "Expression")

#ggsave("figures/DE/bxplt_DE_healthy_pre_post.png")

# bit more publication ready
df_for_bxplt_publ <- df_for_bxplt
levels(df_for_bxplt_publ$Cohort) <- list(`Pre` = "ETI_pre", `Post` = "ETI_post", Healthy = "Healthy")
# set order same as other figures
df_for_bxplt_publ$Cohort <- factor(df_for_bxplt_publ$Cohort, levels=c("Healthy", "Pre", "Post"))
df_for_bxplt_publ$Cohort <- droplevels(df_for_bxplt_publ$Cohort)

# calc which pathway each DE protein is part of
reactome_pway <- gmtPathways("../msigdb/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
de_pros <- sort(unique(df_for_bxplt_publ$Protein))
react_pway_stack <- stack(reactome_pway)

# if complement background, set background as red, else white
pway_to_pro <- list()
backgrounds <- list()

for (p in 1:length(de_pros)){
  protein <- de_pros[[p]]
  # get pathway protein is in
  pway_p_in <- as.character(react_pway_stack[which(react_pway_stack["values"] == protein), "ind"])

  print(protein)
  print(pway_p_in)
  print("----")

  # most are in multiple pathways, so groupo complement pathways together and check if complement in name at all
    if (any(sapply(pway_p_in, function(x) grepl("IMMUNE|COMPLEMENT", x)))) {
      pway_to_pro[protein] <- "COMPLEMENT"
      backgrounds <- append(backgrounds, list(element_rect(fill = "#777777y")))
    } else {
      pway_to_pro[protein] <- "NO COMPLEMENT"
      backgrounds <- append(backgrounds, list(element_rect(fill = "white")))
    }
}

# generate ind plot for each boxplot, so can display geom bracket correctly
plot_list <- list()

for (p in 1:length(de_pros)){
  protein <- de_pros[[p]]

  df_for_protein <- df_for_bxplt_publ[which(df_for_bxplt_publ["Protein"] == protein), ]
  df_for_protein$title <- protein  # force facet for easy strip colour

  pro_bxplt <- ggplot(df_for_protein, aes(x = Cohort, y = value, color = Cohort)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, shape = 1) +
    #facet_wrap(~Protein, scales = "fixed", nrow = 1) +  # fixed scales as otherwise geom_bracket doens't work
    # set colour of background on facet wrap
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank()
            ) +
    # should work as vector, but isn't, so plot brackets individually.
    # can't set geom_bracket y position for each facet!?
    geom_bracket(
      data = df_for_protein,
      xmin = "Pre", xmax = "Post",
      y.position = max(na.omit(df_for_protein$value)) * 1.01,
      label = "ns",
      tip.length = 0.01,
      color = "black",
    ) +
    geom_bracket(
      data = df_for_protein,
      xmin = "Post", xmax = "Healthy",
      y.position = max(na.omit(df_for_protein$value)) * 1.04,
      label = "ns",
      tip.length = 0.01,
      color = "black",
    ) +
    geom_bracket(
      data = df_for_protein,
      xmin = "Pre", xmax = "Healthy",
      y.position = max(na.omit(df_for_protein$value)) * 1.07,
      label = "*",
      tip.length = 0.01,
      color = "black",
    ) +
    #scale_color_manual(values=c(Healthy = "#000000", Pre = "#6c76ad", Post = "#e4d551")) +
    scale_color_manual(values=c(Healthy = "#c679ff", Pre = "#ffa91d", Post = "#7bad00")) +  # colour
    labs(x = "Cohort",
        y = "Expression") +
    facet_grid(. ~title)
      
    # only show legend if last plot
    if (protein != tail(de_pros, n=1)) {
      pro_bxplt$theme$legend.position <- "none"
    }

    # only show x axis title if first plot
    if (protein != head(de_pros, n=1)) {
      pro_bxplt$theme$axis.title.y <- element_blank()
    }

    # change facet strip colour to match protein
    if (pway_to_pro[protein] == "COMPLEMENT") {
      pro_bxplt$theme$strip.background <- element_rect(fill="#838383")
    }
    
  plot_list <- append(plot_list, list(pro_bxplt))
}

# give final fig with legend more space
widths_ <- c()
n_ <- 0.45
for (i in 1:length(de_pros)){
  widths_ <- append(widths_, 1 - (n_/length(de_pros)))
}
widths_[length(widths_)] <- (widths_[length(widths_)] + n_)

pub_bxp <- ggarrange(plotlist = plot_list, nrow = 1, ncol=length(de_pros), widths = widths_)
#ggsave("figures/DE/test_bxp.png", pub_bxp, width = 20)

# volcano plot
vp <- EnhancedVolcano::EnhancedVolcano(DE_pro,
	lab = DE_pro$protein,
	x = "logFC",
	y = "adj.P.Val",
	pCutoff = 0.05,
  FCcutoff = 0.5,
  ylim = c(0,2),
  xlim = c(-2.2, 2.2),
  selectLab = DE_pro[DE_pro$adj.P.Val < 0.05, "protein"],
	subtitle = "",  #Differential Expression
  legendLabels = c("Not Sig.", expression(Log[2] ~ FC), "Adj. p-value", expression(~ Adj. ~ p - value ~ and ~ log[2] ~ FC)),
  legendLabSize = 10,
	title = "Pre vs. Healthy",  #pre ETI vs Healthy
  drawConnectors = TRUE,
  caption="",
  col=c('black', '#c4c410', 'blue', 'red3'),
	legendPosition = "bottom"
)

#ggsave("figures/DE/volcano_healthy_pre.png")

# PCA selected proteins
pca_de <- CBF_PCA(data = t(pro_imp_df_sel),
        groups = meta_df$Cohort,
        useLabels = F,
        labels = "",
        pcs = c(1, 2),
        type = "scores",
        scale = T,
        ellipse = TRUE,
        legendName = "Cohort")
#ggsave("figures/DE/pca_de.png", pca_de)

pca_de_uns <- CBF_PCA(data = t(pro_imp_df_sel),
        groups = meta_df$Cohort,
        useLabels = F,
        labels = "",
        pcs = c(1, 2),
        type = "scores",
        scale = FALSE,
        ellipse = TRUE,
        legendName = "Cohort")
#ggsave("figures/DE/pca_de_unscaled.png", pca_de_uns)

# panel fig for paper.
vp$theme$plot.title$face <- "plain"

top_plot <- ggarrange(vp, pub_bxp,
        labels = c("A", "B"),
        ncol = 2, nrow = 1,
        widths = c(0.65, 1.35),
        heights = c(0.5, 0.5))

# fsea
# plot generated in fpsea.R
pre_post_plt <- readRDS("figures/GSEA/reactome_pre_post.rds")
pre_healthy_plt <- readRDS("figures/GSEA/reactome_pre_healthy.rds")

pre_post_plt$labels$title <- "Pre vs. Post"
pre_healthy_plt$labels$title <- "Pre vs. Healthy"

bottom_plot <- ggarrange(pre_post_plt, pre_healthy_plt,
        labels = c("C", "D"),
        ncol = 2, nrow = 1,
        widths = c(1.25, 1.25),
        heights = c(0.5, 0.5))

all_plot <- ggarrange(top_plot, bottom_plot,
            ncol = 1, nrow = 2)

ggsave("figures/pub_fig.png", all_plot, width=20, height=15)
