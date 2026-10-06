# Visualize final PAML results and the ASTRAL species tree.

required <- c(
  "readr", "dplyr", "tidyr", "ggplot2",
  "ape", "tibble", "stringr", "ggtree"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ape)
  library(tibble)
  library(stringr)
  library(ggtree)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  SCRIPT_DIR <- dirname(
    normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
  )
} else {
  SCRIPT_DIR <- getwd()
}

ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)

RESULTS <- file.path(ROOT, "output", "paml", "paml_master_results.tsv")
BEB <- file.path(ROOT, "output", "paml", "paml_beb_sites.tsv")
TREE_FILE <- file.path(ROOT, "output", "species_tree", "astral_species_tree.tre")
OUTDIR <- file.path(ROOT, "output", "figures")

for (path in c(RESULTS, BEB, TREE_FILE)) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.3),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      plot.title = element_blank()
    )
}

save_figure <- function(plot, name, width_mm, height_mm) {
  base <- file.path(OUTDIR, name)

  ggsave(
    paste0(base, ".pdf"),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    paste0(base, ".png"),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      paste0(base, ".svg"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      device = svglite::svglite
    )
  }
}

# Read final PAML summary tables.
results <- read_tsv(RESULTS, show_col_types = FALSE) %>%
  mutate(gene = recode(gene, "SDHA1" = "SDHA"))

beb <- read_tsv(BEB, show_col_types = FALSE) %>%
  mutate(gene = recode(gene, "SDHA1" = "SDHA"))

hyp_order <- c(
  "nilssonii_terminal",
  "ept_stem",
  "ept_terminals",
  "myotis_stem",
  "myotis_all7_terminals",
  "myotis_wholeclade",
  "vespertilionidae_stem"
)

hyp_labels <- paste0("H", 1:7)

site_order <- c(
  "M1a_vs_M2a",
  "M7_vs_M8",
  "M8a_vs_M8"
)

site_labels <- c(
  "M1a vs M2a",
  "M7 vs M8",
  "M8a vs M8"
)

# Match the gene order used in the manuscript heatmap.
gene_order <- c(
  "PPARGC1A", "PPARG", "PLIN1", "PRDM16", "CPT1B", "ACSL1",
  "SLC25A20", "CPT1A", "CPT2", "SIRT6", "PRKAA1", "SDHB",
  "ATP5F1B", "ATP5F1A", "SDHA", "COX4I1", "NDUFV2", "UQCRB",
  "NDUFS2", "NDUFB8", "ATP5MC1", "NDUFA1", "COX7A1", "CREB5",
  "CREB3L1", "FGFR1", "SOS2", "SOS1", "ADCY3", "ADCY7",
  "ADCY8", "CREB3L4", "UCP1"
)

missing_genes <- setdiff(unique(results$gene), gene_order)

if (length(missing_genes) > 0) {
  stop(
    "Genes present in PAML results but absent from manuscript heatmap order: ",
    paste(missing_genes, collapse = ", ")
  )
}

site_beb <- bind_rows(
  beb %>%
    filter(family == "site", test == "M2a") %>%
    mutate(comparison = "M1a_vs_M2a"),
  beb %>%
    filter(family == "site", test == "M8") %>%
    mutate(comparison = "M7_vs_M8"),
  beb %>%
    filter(family == "site", test == "M8") %>%
    mutate(comparison = "M8a_vs_M8")
) %>%
  group_by(gene, comparison) %>%
  summarise(
    n95 = sum(probability >= 0.95, na.rm = TRUE),
    n99 = sum(probability >= 0.99, na.rm = TRUE),
    .groups = "drop"
  )

branchsite_beb <- beb %>%
  filter(
    family == "branchsite",
    str_ends(test, "_BS_alt")
  ) %>%
  mutate(
    comparison = str_remove(test, "_BS_alt$")
  ) %>%
  group_by(gene, comparison) %>%
  summarise(
    n95 = sum(probability >= 0.95, na.rm = TRUE),
    n99 = sum(probability >= 0.99, na.rm = TRUE),
    .groups = "drop"
  )

classify <- function(data, use_beb = TRUE) {
  data %>%
    mutate(
      n95 = coalesce(n95, 0L),
      n99 = coalesce(n99, 0L),
      raw_sig = !is.na(p_value) & p_value < 0.05,
      fdr_sig = !is.na(q_value) & q_value < 0.05,
      fill = case_when(
        use_beb & fdr_sig & n99 > 0 ~ "beb99",
        use_beb & fdr_sig & n95 > 0 ~ "beb95",
        fdr_sig ~ "fdr",
        raw_sig ~ "raw",
        TRUE ~ "none"
      ),
      star = case_when(
        use_beb & fdr_sig & n99 > 0 ~ "**",
        use_beb & fdr_sig & n95 > 0 ~ "*",
        TRUE ~ ""
      )
    )
}

site_data <- expand_grid(
  gene = gene_order,
  comparison = site_order
) %>%
  left_join(
    results %>%
      filter(test_family == "site") %>%
      select(gene, comparison, p_value, q_value),
    by = c("gene", "comparison")
  ) %>%
  left_join(
    site_beb,
    by = c("gene", "comparison")
  ) %>%
  classify(TRUE) %>%
  mutate(
    panel = "Site models",
    xcol = factor(comparison, levels = site_order, labels = site_labels)
  )

branchsite_data <- expand_grid(
  gene = gene_order,
  comparison = hyp_order
) %>%
  left_join(
    results %>%
      filter(test_family == "branchsite") %>%
      select(gene, comparison, p_value, q_value),
    by = c("gene", "comparison")
  ) %>%
  left_join(
    branchsite_beb,
    by = c("gene", "comparison")
  ) %>%
  classify(TRUE) %>%
  mutate(
    panel = "Branch-site models",
    xcol = factor(comparison, levels = hyp_order, labels = hyp_labels)
  )

branch_data <- expand_grid(
  gene = gene_order,
  comparison = hyp_order
) %>%
  left_join(
    results %>%
      filter(test_family == "branch") %>%
      select(gene, comparison, p_value, q_value),
    by = c("gene", "comparison")
  ) %>%
  mutate(
    n95 = 0L,
    n99 = 0L
  ) %>%
  classify(FALSE) %>%
  mutate(
    panel = "Branch models",
    xcol = factor(comparison, levels = hyp_order, labels = hyp_labels)
  )

plot_data <- bind_rows(
  site_data %>% select(gene, xcol, panel, fill, star),
  branchsite_data %>% select(gene, xcol, panel, fill, star),
  branch_data %>% select(gene, xcol, panel, fill, star)
) %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    panel = factor(
      panel,
      levels = c("Site models", "Branch-site models", "Branch models")
    ),
    fill = factor(
      fill,
      levels = c("none", "raw", "fdr", "beb95", "beb99")
    ),
    text_colour = if_else(fill %in% c("beb95", "beb99"), "white", "black")
  )

fill_colours <- c(
  none = "#FFFFFF",
  raw = "#F8DFA8",
  fdr = "#8EC5E8",
  beb95 = "#6E3B8B",
  beb99 = "#183A5A"
)

fill_labels <- c(
  none = "No evidence",
  raw = "Raw p < 0.05",
  fdr = "FDR q < 0.05",
  beb95 = "FDR + BEB >= 0.95",
  beb99 = "FDR + BEB >= 0.99"
)

# Reproduce the manuscript PAML evidence heatmap.
heatmap_plot <- ggplot(
  plot_data,
  aes(x = xcol, y = gene, fill = fill)
) +
  geom_tile(colour = "#D9D9D9", linewidth = 0.25) +
  geom_text(
    aes(label = star, colour = text_colour),
    fontface = "bold",
    size = 3.1,
    show.legend = FALSE
  ) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = fill_colours,
    labels = fill_labels,
    drop = FALSE,
    name = "Selection evidence"
  ) +
  scale_colour_identity() +
  labs(x = NULL, y = NULL) +
  theme_publication(base_size = 9) +
  theme(
    axis.text.x = element_text(
      size = 7.5,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 7.5, colour = "black"),
    strip.text = element_text(size = 9.5),
    panel.spacing = unit(4, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(colour = "#D9D9D9")
    )
  )

save_figure(
  heatmap_plot,
  "Figure_PAML_selection_heatmap",
  235,
  150
)

# Reproduce the manuscript species-tree panel.
tree <- read.tree(TREE_FILE)

yinpterochiroptera <- c(
  "Pteropus_vampyrus",
  "Pteropus_medius",
  "Pteropus_alecto",
  "Rousettus_aegyptiacus",
  "Rhinolophus_sinicus",
  "Rhinolophus_hipposideros",
  "Rhinolophus_ferrumequinum",
  "Hipposideros_armiger"
)

outgroup <- intersect(yinpterochiroptera, tree$tip.label)

if (
  length(outgroup) >= 2 &&
  isTRUE(tryCatch(is.monophyletic(tree, outgroup), error = function(e) FALSE))
) {
  tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
}

tree <- ladderize(tree)

myotis7 <- c(
  "Myotis_lucifugus",
  "Myotis_yumanensis",
  "Myotis_brandtii",
  "Myotis_davidii",
  "Myotis_daubentonii",
  "Myotis_myotis",
  "Myotis_nattereri"
)

vespertilionidae <- c(
  myotis7,
  "Cnephaeus_nilssonii",
  "Eptesicus_fuscus",
  "Pipistrellus_kuhlii",
  "Plecotus_auritus"
)

display_name <- function(x) {
  str_replace_all(x, "_", " ")
}

tree_plot <- ggtree(
  tree,
  layout = "rectangular",
  branch.length = "none",
  size = 0.45,
  colour = "black"
) +
  geom_tiplab(
    aes(label = display_name(label)),
    fontface = "italic",
    size = 3.2,
    align = TRUE,
    linetype = 0,
    offset = 0.12,
    colour = "black"
  ) +
  theme_tree() +
  theme(
    plot.margin = margin(5, 115, 5, 5),
    plot.background = element_rect(fill = "white", colour = NA)
  )

tree_data <- tree_plot$data %>%
  mutate(
    support = suppressWarnings(as.numeric(label))
  )

support_data <- tree_data %>%
  filter(
    !isTip,
    !is.na(support),
    support < 0.95
  )

if (nrow(support_data) > 0) {
  tree_plot <- tree_plot +
    geom_text(
      data = support_data,
      aes(
        x = x,
        y = y,
        label = sprintf("%.2f", support)
      ),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = -0.55,
      fontface = "bold",
      size = 3
    )
}

tip_positions <- tree_data %>%
  filter(isTip) %>%
  select(label, y)

range_for_taxa <- function(taxa) {
  positions <- tip_positions %>%
    filter(label %in% taxa)

  if (nrow(positions) == 0) {
    return(c(NA_real_, NA_real_))
  }

  range(positions$y)
}

vesper_range <- range_for_taxa(vespertilionidae)
yinptero_range <- range_for_taxa(yinpterochiroptera)
yangochiroptera <- setdiff(tree$tip.label, yinpterochiroptera)
yango_range <- range_for_taxa(yangochiroptera)

x_max <- max(tree_data$x, na.rm = TRUE)
x_vesper <- x_max * 1.58
x_suborder <- x_max * 1.83

add_vertical_clade_label <- function(plot, x, y_range, label) {
  if (any(is.na(y_range))) {
    return(plot)
  }

  plot +
    annotate(
      "segment",
      x = x,
      xend = x,
      y = y_range[1] - 0.25,
      yend = y_range[2] + 0.25,
      linewidth = 0.65,
      colour = "black"
    ) +
    annotate(
      "text",
      x = x + 0.12,
      y = mean(y_range),
      label = label,
      angle = 90,
      fontface = "bold",
      size = 3.6
    )
}

tree_plot <- add_vertical_clade_label(
  tree_plot,
  x_vesper,
  vesper_range,
  "Vespertilionidae"
)

tree_plot <- add_vertical_clade_label(
  tree_plot,
  x_suborder,
  yango_range,
  "Yangochiroptera"
)

tree_plot <- add_vertical_clade_label(
  tree_plot,
  x_suborder,
  yinptero_range,
  "Yinpterochiroptera"
)

tree_plot <- tree_plot +
  xlim(0, x_max * 2.12) +
  coord_cartesian(clip = "off")

save_figure(
  tree_plot,
  "Figure_species_tree_foregrounds",
  195,
  125
)