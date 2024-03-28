library(fgsea)
library(tibble)
library(tidyverse)
library(data.table)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(ggpubr)

# post in relation to pre.
pre_post_DE_res <- read.csv("figures/DE/topranked_pre_post.csv", row.names = 1)  # pre vs post DE results
pre_healthy_DE_res <- read.csv("figures/DE/topranked_healthy_eti.csv", row.names = 1)  # pre vs post DE results

hm_all_pway <- gmtPathways("../msigdb/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.symbols.gmt")
reactome_pway <- gmtPathways("..msigdb/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
kegg_pway <- gmtPathways("../msigdb/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
wp_pway <- gmtPathways("../msigdb/msigdb_v2023.2.Hs_GMTs/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt")

all_pways <- list("Hallmark all" = hm_all_pway, "Reactome" = reactome_pway, "KEGG" = kegg_pway, "Wikipathways" = wp_pway)

# store results
results_gsea <- list()

set.seed(1)

# run GSEA for pre vs post, and pre vs healthy.
analysis <- list("pre_vs_post" = pre_post_DE_res, "pre_vs_healthy" = pre_healthy_DE_res)


# get dataframe of proteins in dataset that are present in each sig pathway
reactome_leading_edge_df <- data.frame(matrix(ncol = 0, nrow = 0))

for (j in 1:length(analysis)){

	analysis_name <- names(analysis)[[j]]
	DE_res <- analysis[[j]]

	names(DE_res)[names(DE_res) == 'protein'] <- 'SYMBOL'

	# named list gene name: rank
	ranks <- DE_res$t
	names(ranks) <- DE_res$SYMBOL

	# check against all pathways.
	for (i in 1:length(all_pways)){

		# check against each pathway.
		pway <- all_pways[[i]]
		n_pway <- names(all_pways[i])

		fgseaRes <- fgsea(pathways = pway,
					stats = ranks,
					minSize = 15,
					maxSize = 500)
		
		# store results
		results_gsea[[paste(n_pway, "_", analysis_name, sep="")]] <- fgseaRes
		
		ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
			geom_col(aes(fill = padj < 0.05)) +
			coord_flip() +
			labs(x = "Pathway", y = "Normalized Enrichment Score",
				title = paste(n_pway, ": NES from GSEA", sep = "")) +
			theme_bw()

		if (n_pway == 'Reactome'){
			#ggsave(paste("figures/GSEA/", n_pway, analysis_name, "_pathway.png", sep = ""), width = 20)
		} else {
			#ggsave(paste("figures/GSEA/", n_pway, analysis_name, "_pathway.png", sep = ""))
		}

		topPathways <- fgseaRes[head(order(padj), n = 15)][order(NES), pathway]

		if (n_pway == 'Reactome'){
			gsea_all_tab <- plotGseaTable(pway[topPathways], ranks, fgseaRes, gseaParam = 0.5, render = FALSE,
							colwidths=c(0.75, 0.15, 0.06, 0.06, 0.06))
			#ggsave(paste("figures/GSEA/gsea_table_", n_pway, analysis_name, ".png", sep = ""), gsea_all_tab, width = 21)
		} else {
			gsea_all_tab <- plotGseaTable(pway[topPathways], ranks, fgseaRes, gseaParam = 0.5, render = FALSE,
							colwidths=c(0.6, 0.2, 0.1, 0.1, 0.1))
			#ggsave(paste("figures/GSEA/gsea_table_", n_pway, analysis_name, ".png", sep = ""), gsea_all_tab, width = 15)
		}

		# if sig, then plot enrichment
		sig <- fgseaRes[padj < 0.05]

		if (nrow(sig) > 0){
			for (j in 1:nrow(sig)){

			pathway_to_plot <- sig[j, ]$pathway

			plotEnrichment(pway[[pathway_to_plot]], ranks) +
					labs(title = pathway_to_plot,
						x = "Rank",
						y = "Enrichment Score") +
					theme_bw()
			
			#ggsave(paste("figures/GSEA/enrichment_", pathway_to_plot, analysis_name, ".png", sep=""), height = 5)

			if (n_pway == "Reactome") {
				# record leading edge genes, that are driving differences
				row_to_add <- sig[j, c("pathway", "leadingEdge")]
				row_to_add$analysis <- analysis_name
				row_to_add$leadingEdge <- paste(row_to_add$leadingEdge[[1]], collapse=', ')
				reactome_leading_edge_df <- rbind(reactome_leading_edge_df, row_to_add)
			}
			}
		}

		write.csv(reactome_leading_edge_df, "figures/GSEA/leading_edge_reactome.csv", row.names = FALSE)



		# dotplot
		both_plt <- ggplot(fgseaRes[!is.na(fgseaRes[["pathway"]]),], aes(x=.data[["NES"]], 
																		y=reorder(.data[["pathway"]], .data[["NES"]]), size=log10(.data[["size"]]))) + 
			geom_point(shape = 21, aes(fill = .data[["padj"]])) +
			xlab("NES") + ylab("") + ggtitle(n_pway) +
			geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
			scale_fill_continuous(high = "blue", low = "red") +
			scale_y_discrete(
				labels = ~ ifelse(.x %in% fgseaRes[padj < 0.05]$pathway, paste(.x, "*"), .x)
			) +
			theme_bw(base_size = 16) + theme(axis.title = element_text(size = 28),
										axis.text.y = element_text(size = 10))

		if (n_pway == 'Reactome'){
			#ggsave(paste0("figures/GSEA/gsea_dotplot_", n_pway, analysis_name, ".png"),	both_plt, width = 20)
		} else {
			#ggsave(paste0("figures/GSEA/gsea_dotplot_", n_pway, analysis_name, ".png"),	both_plt)
		}

		# sig dotplot only
		if (nrow(sig) > 0){
			sig_plt <- ggplot(sig[!is.na(sig[["pathway"]]),], aes(x=.data[["NES"]], 
																			y=reorder(.data[["pathway"]], .data[["NES"]]), size=log10(.data[["size"]]))) + 
				geom_point(shape = 21, aes(fill = .data[["padj"]])) +
				xlab("NES") + ylab("") + ggtitle(n_pway) +
				geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
				scale_fill_continuous(high = "#aaaaf3", low = "#ff0000") +  #, limits=c(0.001, 0.005), breaks=seq(0.001, 0.005,by=0.001)
				scale_size_continuous(range = c(10^log(1.3), 10^log(2.0))) +  # inverse log limits  #breaks arg kills the legend for some reason
				theme_bw(base_size = 16) + theme(axis.title = element_text(size = 28),
											axis.text.y = element_text(size = 10)) +
				guides(size = guide_legend(order = 1))  # force legend order

			#ggsave(paste0("figures/GSEA/gsea_dotplot_sig", n_pway, analysis_name, ".png"), sig_plt)

			# Save pre vs post hallmarks for publication figure
			if (n_pway == "Hallmark all" & analysis_name == "pre_vs_post"){
				#saveRDS(sig_plt, "figures/GSEA/hallmark_pre_post.rds")
			}
			if (n_pway == "Hallmark all" & analysis_name == "pre_vs_healthy"){
				#saveRDS(sig_plt, "figures/GSEA/hallmark_pre_healthy.rds")
			}
			if (n_pway == "Wikipathways" & analysis_name == "pre_vs_post"){
				#saveRDS(sig_plt, "figures/GSEA/wiki_pre_post.rds")
			}
			if (n_pway == "Wikipathways" & analysis_name == "pre_vs_healthy"){
				#saveRDS(sig_plt, "figures/GSEA/wiki_pre_healthy.rds")
			}
			if (n_pway == "Reactome" & analysis_name == "pre_vs_post"){
				#saveRDS(sig_plt, "figures/GSEA/reactome_pre_post.rds")
			}
			if (n_pway == "Reactome" & analysis_name == "pre_vs_healthy"){
				#saveRDS(sig_plt, "figures/GSEA/reactome_pre_healthy.rds")
			}
		}

		up_plt <- ggplot(fgseaRes[!is.na(fgseaRes[["pathway"]]) & fgseaRes$NES > 0 ,], aes(x=.data[["NES"]], 
																		y=reorder(.data[["pathway"]], .data[["NES"]]), size=log10(.data[["size"]]))) + 
			geom_point(shape = 21, aes(fill = .data[["padj"]])) +
			xlab("NES") + ylab("") +
			geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
			scale_fill_continuous(high = "blue", low = "red") +
			scale_y_discrete(
				labels = ~ ifelse(.x %in% fgseaRes[padj < 0.05 & NES > 0]$pathway, paste(.x, "*"), .x)
			) +
			theme_bw(base_size = 16) + theme(axis.title = element_text(size = 28),
										axis.text.y = element_text(size = 12))

		down_plt <- ggplot(fgseaRes[!is.na(fgseaRes[["pathway"]]) & fgseaRes$NES < 0 ,], aes(x=.data[["NES"]], 
																		y=reorder(.data[["pathway"]], .data[["NES"]]), size=log10(.data[["size"]]))) + 
			geom_point(shape = 21, aes(fill = .data[["padj"]])) +
			xlab("NES") + ylab("") +
			geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
			scale_fill_continuous(high = "blue", low = "red") +
			scale_y_discrete(
				labels = ~ ifelse(.x %in% fgseaRes[padj < 0.05 & NES < 0]$pathway, paste(.x, "*"), .x)
			) +
			theme_bw(base_size = 16) + theme(axis.title = element_text(size = 28),
										axis.text.y = element_text(size = 12))

		comb_plt <- ggarrange(up_plt, down_plt, ncol = 2)

		if (n_pway == 'Reactome'){
			#ggsave(paste0("figures/GSEA/gsea_dotplot_window_", n_pway, analysis_name, ".png"), comb_plt, width = 40)
		} else {
			#ggsave(paste0("figures/GSEA/gsea_dotplot_window_", n_pway, analysis_name, ".png"), comb_plt, width = 20)
		}

	}

	# coagulation and inflammation pathways significantly downregulated.

}
