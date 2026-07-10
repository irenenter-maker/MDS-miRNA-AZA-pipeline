suppressWarnings({
  suppressMessages({
    library(DESeq2)
    library(tidyverse)
    library(fgsea)
    library(ggrepel)
    library(ggraph)
    library(tidygraph)
    library(scales)
    library(writexl)
    library(patchwork)
    library(yaml)
  })
})

cfg <- yaml::read_yaml("config.yaml")

setwd(cfg$working_dir)

network_purple <- "#742881"
network_green  <- "#1b7939"
bar_fill_low   <- "#C3C3E5"
cluster_palette <- c("#E69F00", "#56B4E9", "#009E73",
                     "#0072B2", "#F0E442", "#D55E00", "#CC79A7")

consistent_theme <- theme_classic(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle = element_text(size = 14, color = "grey30"),
    axis.text     = element_text(size = 14, color = "black"),
    axis.title    = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.title  = element_text(size = 14, face = "bold"),
    legend.text   = element_text(size = 12),
    plot.margin   = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3)
  )

save_plot <- function(plot, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height)
  invisible(plot)
}

volcano_plot <- function(res, title, subtitle, x_limits, max_overlaps) {
  res %>%
    mutate(
      log10padj   = -log10(padj),
      significant = padj < 0.05,
      label       = ifelse(padj < 0.05 & abs(NES) > 1.5, pathway, NA)
    ) %>%
    ggplot(aes(x = NES, y = log10padj)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 2.5) +
    geom_text_repel(aes(label = label), size = 4,
                    max.overlaps = max_overlaps, box.padding = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = network_purple)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.5, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.5, color = "black") +
    labs(x = "Normalized Enrichment Score (NES)",
         y = expression(-log[10](padj)),
         color = "Sig.", title = title, subtitle = subtitle) +
    scale_y_continuous(limits = c(0, 2.2)) +
    scale_x_continuous(limits = x_limits) +
    consistent_theme
}

bar_plot <- function(plot_data, title, subtitle) {
  ggplot(plot_data, aes(x = ES, y = pathway, fill = neg_log10_padj)) +
    geom_col(width = 0.75) +
    scale_fill_gradient(low = bar_fill_low, high = network_purple,
                        name = expression(-log[10](padj))) +
    labs(x = "Enrichment Score (ES)", y = NULL,
         title = title, subtitle = subtitle) +
    consistent_theme +
    theme(panel.grid.major.y = element_blank())
}

make_enrich_plot <- function(pathway_name, pathways_list, stats, title_text,
                             y_lab = "Enrichment Score (ES)", compact = FALSE) {
  if (!pathway_name %in% names(pathways_list)) {
    warning("Pathway not found: ", pathway_name)
    return(ggplot() + theme_void())
  }
  p <- plotEnrichment(pathways_list[[pathway_name]], stats) +
    labs(title = title_text, x = "Ranked miRNAs", y = y_lab) +
    consistent_theme +
    theme(panel.grid = element_blank())

  if (compact) {
    p <- p + theme(plot.title  = element_text(face = "bold", size = 12),
                   axis.title  = element_text(size = 16),
                   axis.text   = element_text(size = 8),
                   plot.margin = margin(5, 5, 5, 5))
  }
  p
}

counts   <- read.csv(cfg$inputs$counts, row.names = 1, check.names = FALSE)
metadata <- read.csv(cfg$inputs$metadata, stringsAsFactors = FALSE)

metadata <- metadata %>%
  mutate(time     = case_when(treatment == "pre"  ~ "before",
                              treatment == "post" ~ "after"),
         response = case_when(treat_respondance == "responding"     ~ "R",
                              treat_respondance == "non_responding" ~ "NR"))

rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(counts), ]

metadata$time     <- factor(metadata$time,     levels = c("before", "after"))
metadata$response <- factor(metadata$response, levels = c("R", "NR"))
metadata$patient  <- factor(metadata$patient)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                              design = ~ response + time + response:time)

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds$response <- relevel(dds$response, ref = "NR")
dds$time     <- relevel(dds$time,     ref = "before")
dds <- DESeq(dds)

dds_res_df <- as.data.frame(DESeq2::results(dds)) %>%
  mutate(negative_log10_pvalue = -log10(pvalue),
         negative_log10_pvalue = ifelse(is.infinite(negative_log10_pvalue), 1000,
                                        negative_log10_pvalue),
         signed_rank_stats = sign(log2FoldChange) * negative_log10_pvalue +
           runif(n(), -1e-10, 1e-10),
         miRNA = rownames(.))

miRNA_list <- setNames(dds_res_df$signed_rank_stats, dds_res_df$miRNA) %>%
  sort(decreasing = TRUE)

gsea_specs <- list(
  miRWalk_Pathways_mature = list(
    min_size         = 15, max_size = 500,
    volcano_title    = "Signaling Pathways Modulated by AZA Response",
    volcano_subtitle = "miRWalk Pathways",
    volcano_xlim     = c(-2.2, 2.2),
    volcano_overlaps = Inf,
    bar_title        = "Top Signaling Pathways",
    bar_subtitle     = "Ranked by Enrichment Score",
    bar_size         = c(12, 8),
    bar_top_n        = NULL,
    enrich_title     = "Prostagladin Synthesis & Regulation"
  ),
  miRTarBase_mature = list(
    min_size         = 10, max_size = 500,
    volcano_title    = "Enrichment of miRNA Target Genes",
    volcano_subtitle = "miRTarBase",
    volcano_xlim     = c(-2.2, 2.2),
    volcano_overlaps = 10,
    bar_title        = "Top Targeted Genes Associated with Response",
    bar_subtitle     = "Top 20 Ranked by Enrichment Score",
    bar_size         = c(12, 8),
    bar_top_n        = 20,
    enrich_title     = "KCNB1"
  ),
  miRWalk_GO_mature = list(
    min_size         = 15, max_size = 300,
    volcano_title    = "Significant enriched GO terms",
    volcano_subtitle = "miRWalk GO",
    volcano_xlim     = c(-2.2, 2.2),
    volcano_overlaps = 15,
    bar_title        = "Significant enriched GO terms",
    bar_subtitle     = "miRWalk GO",
    bar_size         = c(14, 8),
    bar_top_n        = NULL,
    enrich_title     = "Prostaglandin Biosynthetic Process"
  ),
  miRPathDB_GO_Biological_process_mature = list(
    min_size         = 20, max_size = 500,
    volcano_title    = "Biological Processes Enriched in AZA Response",
    volcano_subtitle = "miRPathDB GO Bio Process",
    volcano_xlim     = c(-2.3, 2.3),
    volcano_overlaps = 15,
    bar_title        = "Top Biological Processes (GO)",
    bar_subtitle     = "Distinguishing Responders",
    bar_size         = c(14, 8),
    bar_top_n        = NULL,
    enrich_title     = "Response to starvation"
  )
)

gsea_specs <- imap(gsea_specs, ~ c(.x, cfg$gsea[[.y]]))

run_gsea_analysis <- function(spec, name, stats) {
  pathways <- gmtPathways(spec$gmt)

  set.seed(42)
  res <- fgsea(pathways = pathways, stats = stats,
               minSize = spec$min_size, maxSize = spec$max_size,
               eps = 1e-300, nPerm = 10000, nproc = 1)

  sig <- res %>% filter(padj < 0.05) %>% arrange(padj)

  save_plot(volcano_plot(res, spec$volcano_title, spec$volcano_subtitle,
                         spec$volcano_xlim, spec$volcano_overlaps),
            spec$volcano_stem, width = 12, height = 8)

  plot_data <- sig %>%
    mutate(neg_log10_padj = -log10(padj),
           pathway = reorder(pathway, -ES))

  plot_data %>%
    mutate(leadingEdge = sapply(leadingEdge, toString)) %>%
    write_xlsx(spec$xlsx)

  bar_data <- if (is.null(spec$bar_top_n)) plot_data else
    plot_data %>% arrange(padj) %>% slice(seq_len(spec$bar_top_n))

  save_plot(bar_plot(bar_data, spec$bar_title, spec$bar_subtitle),
            spec$bar_stem,
            width = spec$bar_size[1], height = spec$bar_size[2])

  top_hit <- sig %>% slice(1) %>% pull(pathway)
  save_plot(make_enrich_plot(top_hit, pathways, stats, spec$enrich_title),
            spec$enrich_stem, width = 8, height = 6)

  list(res = res, sig = sig, plot_data = plot_data,
       pathways = pathways, top_hit = top_hit)
}

gsea_results <- imap(gsea_specs, ~ run_gsea_analysis(.x, .y, miRNA_list))

gsea_data1 <- gsea_results$miRWalk_Pathways_mature$sig
gsea_data2 <- gsea_results$miRTarBase_mature$sig
gsea_data3 <- gsea_results$miRWalk_GO_mature$sig
gsea_data4 <- gsea_results$miRPathDB_GO_Biological_process_mature$sig
pathways   <- gsea_results$miRWalk_Pathways_mature$pathways
pathways4  <- gsea_results$miRPathDB_GO_Biological_process_mature$pathways

network_input <- bind_rows(
  gsea_data1 %>% arrange(padj) %>% slice(1:7)  %>% mutate(database = "Signaling"),
  gsea_data4 %>% mutate(database = "Biological Process"),
  gsea_data3 %>% arrange(padj) %>% slice(1:3)  %>% mutate(database = "GO"),
  gsea_data2 %>% arrange(padj) %>% slice(1:10) %>% mutate(database = "Target Gene")
)

graph_data <- network_input %>%
  select(pathway, leadingEdge, database) %>%
  unnest(leadingEdge) %>%
  rename(from = pathway, to = leadingEdge) %>%
  mutate(from = str_replace_all(from, "_", " "),
         from = str_wrap(from, width = 20))

miRNA_info <- dds_res_df %>%
  select(miRNA, log2FoldChange) %>%
  rename(name = miRNA, val = log2FoldChange)

build_base_nodes <- function(graph_data, miRNA_info) {
  data.frame(name = unique(c(graph_data$from, graph_data$to))) %>%
    left_join(miRNA_info, by = "name") %>%
    mutate(type = ifelse(name %in% graph_data$from, "Biological Theme", "miRNA"),
           val  = ifelse(type == "Biological Theme", 0, val))
}

network_common_layers <- function(mirna_size) {
  list(
    geom_edge_link(aes(alpha = after_stat(index)), color = "grey70",
                   show.legend = FALSE),
    geom_node_point(aes(filter = type == "Biological Theme",
                        size = type, shape = type), color = "grey20"),
    geom_node_text(
      aes(label = name,
          size          = I(ifelse(type == "Biological Theme", 3, 2)),
          fontface      = I(ifelse(type == "Biological Theme", "bold", "plain")),
          segment.color = I(ifelse(type == "Biological Theme", "grey20", "grey95"))),
      color = "black", repel = TRUE, check_overlap = TRUE,
      bg.color = "white", bg.r = 0.15, box.padding = 0.5, force = 11,
      min.segment.length = 0, max.overlaps = Inf),
    scale_size_manual(values = c("Biological Theme" = 4, "miRNA" = mirna_size)),
    scale_shape_manual(values = c("Biological Theme" = 15, "miRNA" = 19)),
    theme_graph(base_family = "sans", background = "white")
  )
}

graph_lfc <- tbl_graph(nodes = build_base_nodes(graph_data, miRNA_info),
                       edges = graph_data)

set.seed(42)
p <- ggraph(graph_lfc, layout = "nicely") +
  network_common_layers(mirna_size = 2.5) +
  geom_node_point(aes(filter = type == "miRNA", color = val,
                      size = type, shape = type)) +
  scale_color_gradient2(low = network_green, high = network_purple, mid = "white",
                        midpoint = 0, limits = c(-3, 3), oob = scales::squish,
                        name = "Log2 Fold Change")
print(p)
save_plot(p, cfg$outputs$network_lfc, width = 18, height = 10)

plot_single_mir <- function(dds_object, mir_name, title_text) {
  if (!(mir_name %in% rownames(dds_object))) {
    return(ggplot() + theme_void() +
             geom_text(aes(0, 0, label = paste("Not found:", mir_name))))
  }

  plotCounts(dds_object, gene = mir_name,
             intgroup = c("response", "time", "patient"), returnData = TRUE) %>%
    ggplot(aes(x = time, y = count, color = response, group = patient)) +
    geom_line(aes(linetype = response), linewidth = 0.8, alpha = 0.6) +
    geom_point(size = 3) +
    scale_y_log10() +
    scale_color_manual(values = c("R" = "#5cae63", "NR" = "#986eac")) +
    labs(title = title_text, subtitle = mir_name,
         y = "Norm. Counts (log10)", x = NULL) +
    consistent_theme +
    theme(legend.position = "none",
          plot.title    = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12, color = "grey40"),
          axis.title.y  = element_text(size = 12))
}

mir_groups <- list(
  "Negative log2 FC" = c("hsa-miR-18a-3p", "hsa-miR-423-5p", "hsa-miR-99b-5p",
                         "hsa-let-7e-5p", "hsa-miR-92a-3p"),
  "Positive log2 FC" = c("hsa-miR-301a-3p", "hsa-miR-30b-5p", "hsa-miR-20a-5p",
                         "hsa-miR-331-3p", "hsa-miR-106a-5p")
)

final_panel <- imap(mir_groups, function(mirs, title) {
  wrap_plots(map(mirs, ~ plot_single_mir(dds, .x, title)), nrow = 1)
}) %>%
  wrap_plots(ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

save_plot(final_panel, cfg$outputs$key_drivers_panel, width = 25, height = 10)

top_n_hits <- 8

dot_data <- bind_rows(
  gsea_data1 %>% arrange(padj) %>% slice(1:top_n_hits) %>% mutate(Category = "Pathways (miRWalk)"),
  gsea_data3 %>% arrange(padj) %>% slice(1:top_n_hits) %>% mutate(Category = "GO (miRWalk)"),
  gsea_data4 %>% arrange(padj) %>% slice(1:top_n_hits) %>% mutate(Category = "Biological Processes (miRPathDB)")
) %>%
  mutate(neg_log10_padj = -log10(padj),
         pathway_clean  = str_wrap(str_replace_all(pathway, "_", " "), width = 50),
         Category = factor(Category, levels = c("Pathways (miRWalk)",
                                                "GO (miRWalk)",
                                                "Biological Processes (miRPathDB)")))

plot_A <- ggplot(dot_data, aes(x = NES, y = reorder(pathway_clean, NES))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_col(aes(fill = neg_log10_padj), color = "black", width = 0.7) +
  facet_wrap(~ Category, ncol = 1, strip.position = "top", scales = "free_y") +
  scale_fill_gradient(low = "#b8d5e6", high = "#0072b2", name = "-log10(padj)") +
  labs(x = "NES", y = NULL) +
  consistent_theme +
  theme(strip.background   = element_rect(fill = "grey95", color = NA),
        strip.text         = element_text(face = "bold", size = 13),
        axis.text.y        = element_text(size = 11, lineheight = 0.8),
        plot.margin        = margin(10, 10, 10, 10),
        legend.position    = "right")

panel_b_specs <- list(
  list(pathway = "WP98_Prostaglandin_Synthesis_and_Regulation",
       set = pathways,  title = "Prostaglandin Synthesis & Regulation"),
  list(pathway = "response_to_starvation",
       set = pathways4, title = "Response to Starvation"),
  list(pathway = "hsa03010_Ribosome",
       set = pathways,  title = "Ribosome")
)

plot_B_combined <- map(panel_b_specs, ~ make_enrich_plot(
  .x$pathway, .x$set, miRNA_list, .x$title, y_lab = "ES", compact = TRUE
)) %>% wrap_plots(ncol = 1)

plot_C <- p + theme(plot.margin = margin(20, 20, 20, 20))

final_figure_v3 <- ((plot_A | wrap_elements(plot_B_combined)) +
                      plot_layout(widths = c(1.2, 0.8))) / plot_C +
  plot_layout(heights = c(1.5, 1))

save_plot(final_figure_v3, cfg$outputs$composite_figure,
          width = 18, height = 18)

get_mir_stem <- function(x) {
  x %>%
    tolower() %>%
    str_remove_all("hsa-") %>%
    str_remove_all("mir-|mir") %>%
    str_remove_all("-3p|-5p") %>%
    str_remove_all("-[0-9]+$")
}

read_cluster_gmt <- function(path) {
  readLines(path) %>%
    str_split("\t") %>%
    map_dfr(function(x) {
      if (length(x) < 3) return(NULL)
      data.frame(Cluster_ID = x[1], precursor = x[3:length(x)],
                 stringsAsFactors = FALSE)
    }) %>%
    mutate(match_key = get_mir_stem(precursor))
}

plot_cluster_network <- function(gmt_path, legend_name, out_stem) {
  cluster_df <- read_cluster_gmt(gmt_path)

  nodes <- build_base_nodes(graph_data, miRNA_info) %>%
    mutate(match_key = get_mir_stem(name)) %>%
    left_join(cluster_df %>% select(Cluster_ID, match_key), by = "match_key") %>%
    group_by(name) %>% slice(1) %>% ungroup()

  large_clusters <- nodes %>%
    filter(type == "miRNA", !is.na(Cluster_ID)) %>%
    count(Cluster_ID) %>%
    filter(n > 2) %>%
    pull(Cluster_ID)

  graph <- tbl_graph(
    nodes = nodes %>%
      mutate(Cluster_Color_Group = ifelse(Cluster_ID %in% large_clusters,
                                          Cluster_ID, NA)),
    edges = graph_data
  )

  set.seed(42)
  plt <- ggraph(graph, layout = "nicely") +
    network_common_layers(mirna_size = 3) +
    geom_node_point(aes(filter = type == "miRNA", color = Cluster_Color_Group,
                        size = type, shape = type)) +
    scale_color_manual(name = legend_name, values = cluster_palette,
                       na.value = "grey50") +
    theme(legend.position = "right")

  ggsave(paste0(out_stem, ".png"), plot = plt, width = 18, height = 10, dpi = 300)
  print(plt)
  invisible(plt)
}

plot_cluster_network(cfg$inputs$cluster_gmt,
                     "Precursor Clusters (>2 miRs)", cfg$outputs$network_clusters)

plot_cluster_network(cfg$inputs$family_gmt,
                     "Precursor Families (>2 miRs)", cfg$outputs$network_families)