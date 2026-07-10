suppressWarnings({
  suppressMessages({
    library("DESeq2")
    library("tidyverse")
    library("ggrepel")
    library("tidygraph")
    library("kableExtra")
    library("yaml")
  })
})

config <- read_yaml("config.yaml")
setwd(config$working_dir)  

counts   <- read.csv(config$inputs$counts, row.names = 1, check.names = FALSE)
metadata <- read.csv(config$inputs$metadata, stringsAsFactors = FALSE)

metadata <- metadata %>%
  mutate(time = case_when(treatment == "pre" ~ "before",
                          treatment == "post" ~ "after"),
         response = case_when(treat_respondance == "responding" ~ "R",
                              treat_respondance == "non_responding" ~ "NR"))

rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(counts), ]

metadata$time <- factor(metadata$time, levels = c("before", "after"))
metadata$response <- factor(metadata$response, levels = c("R", "NR"))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ response + time + response:time)

keep <- rowSums(counts(dds) >= config$params$min_count) >= config$params$min_samples
dds <- dds[keep, ]
dds$response <- relevel(dds$response, ref = "NR")
dds$time <- relevel(dds$time, ref = "before")
dds <- DESeq(dds)
dds_res <- DESeq2::results(dds)

dds_res_df <- as.data.frame(dds_res) %>%
  mutate(negative_log10_pvalue = -log10(pvalue)) %>%
  mutate(negative_log10_pvalue = ifelse(
    is.infinite(negative_log10_pvalue),
    1000,
    negative_log10_pvalue)
  ) %>%
  mutate(signed_rank_stats = sign(log2FoldChange) * negative_log10_pvalue + 
           runif(n(), -1e-10, 1e-10))

dds_res_df$miRNA <- rownames(dds_res_df)
miRNA_list <- dds_res_df$signed_rank_stats
names(miRNA_list) <- dds_res_df$miRNA
miRNA_list <- sort(miRNA_list, decreasing = TRUE)

dds_res_df <- dds_res_df %>%
  mutate(
    Significance = "Not Significant",
    Label = if_else(miRNA %in% top20_miRNAs$miRNA, miRNA, ""),
    category = case_when(
      Label != "" ~ "Highlighted",
      # Added $params$ here:
      !is.na(pvalue) & pvalue < config$params$sig_cut & log2FoldChange > config$params$lfc_cut ~ "Upregulated",
      !is.na(pvalue) & pvalue < config$params$sig_cut & log2FoldChange < -config$params$lfc_cut ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    category = factor(category,
      levels = c("Not significant", "Downregulated", "Upregulated", "Highlighted"))
  )

top20_miRNAs <- dds_res_df %>%
  filter(!is.na(pvalue)) %>%
  arrange(pvalue) %>%
  head(config$params$n_labelled)

max_y_val <- max(dds_res_df$negative_log10_pvalue, na.rm = TRUE)
y_limit_upper <- if(max_y_val > 0) max_y_val * 1.3 else 0.05

dds_res_df <- dds_res_df %>%
  mutate(
    Significance = "Not Significant",
    Label = if_else(miRNA %in% top20_miRNAs$miRNA, miRNA, ""),
    category = case_when(
      Label != "" ~ "Highlighted",
      # Added $params$ here:
      !is.na(pvalue) & pvalue < config$params$sig_cut & log2FoldChange > config$params$lfc_cut ~ "Upregulated",
      !is.na(pvalue) & pvalue < config$params$sig_cut & log2FoldChange < -config$params$lfc_cut ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    category = factor(category,
                      levels = c("Not significant", "Downregulated", "Upregulated", "Highlighted"))
  )

  network_purple <- "#742881"
network_green  <- "#1b7939"

consistent_theme <- theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle = element_text(size = 14, color = "grey30"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3)
  )

volcano_plot <- ggplot(
  dds_res_df %>% arrange(category),
  aes(x = log2FoldChange, y = negative_log10_pvalue, colour = category)
) +
  # Fixed below: config$params$lfc_cut
  geom_vline(xintercept = c(-config$params$lfc_cut, config$params$lfc_cut), linetype = "dashed",
             colour = "grey50", alpha = 0.5) +
  # Fixed below: config$params$sig_cut
  geom_hline(yintercept = -log10(config$params$sig_cut), linetype = "dashed",
             colour = "grey50", alpha = 0.5) +
  geom_point(aes(size = category, alpha = category)) +
  scale_colour_manual(values = c(
    "Not significant" = "grey70",
    "Downregulated"   = "#2C7BB6",
    "Upregulated"     = "#D7191C",
    "Highlighted"     = network_purple
  )) +
  scale_size_manual(values = c("Not significant" = 2, "Downregulated" = 2.5,
                               "Upregulated" = 2.5, "Highlighted" = 3.5),
                    guide = "none") +
  scale_alpha_manual(values = c("Not significant" = 0.4, "Downregulated" = 0.7,
                                "Upregulated" = 0.7, "Highlighted" = 1),
                     guide = "none") +
  geom_text_repel(
    data = filter(dds_res_df, Label != ""),
    aes(label = Label),
    colour = "black", size = 3.5, fontface = "bold",
    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2,
    segment.color = "grey60", segment.size = 0.4
  ) +
  coord_cartesian(ylim = c(0, y_limit_upper)) +
  labs(
    x = expression(log[2] ~ "(fold change)"),
    y = expression(-log[10] ~ "(p-value)"),
    colour = NULL
  ) +
  consistent_theme


ggsave(
  config$outputs$suppl_volcano, 
  plot = volcano_plot,
  width = 9,
  height = 8,
  dpi = 300
)

supp_kbl <- top20_miRNAs %>%
  dplyr::select(
    baseMean, stat, log2FoldChange, pvalue, padj
  ) %>%
  kbl(
    caption = "Top 20 miRNAs ranked from miRNA set enrichment analysis.",
    align = "c",
    col.names = c(
      "miRNA", "Base mean", "Wald statistic",
      "log2(Fold change)", "p-value", "Adjusted p-value"
    )
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive")
  )

save_kable(supp_kbl, file = config$outputs$suppl_top20_table)
