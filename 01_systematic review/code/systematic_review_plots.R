#!/usr/bin/env Rscript

# Analyse the filtered 137-study systematic-review dataset and reproduce review figures/tables.

# Load packages required for review summaries and figures.
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(scales)
  library(maps)
  library(writexl)
  library(patchwork)

# Resolve paths relative to this script so the workflow is portable.
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

# Read the locked 137-study curated dataset as df_final.
INPUT_FILE <- Sys.getenv(
  "REVIEW_DATASET",
  file.path(ROOT, "input", "Table_S1_curated_review_dataset.xlsx")
)

OUT_DIR <- file.path(ROOT, "output", "review")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")
FIGURE_INPUT_DIR <- file.path(OUT_DIR, "figure_inputs")

for (directory in c(OUT_DIR, TABLE_DIR, FIGURE_DIR, FIGURE_INPUT_DIR)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(INPUT_FILE)) {
  stop("Review dataset not found: ", INPUT_FILE)
}

extension <- tolower(tools::file_ext(INPUT_FILE))

df_final <- switch(
  extension,
  xlsx = read_excel(INPUT_FILE, sheet = 1),
  xls = read_excel(INPUT_FILE, sheet = 1),
  csv = read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE),
  tsv = read.delim(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE),
  stop("Unsupported review dataset format: ", extension)
)

# Confirm that the supplied dataset is the final 137-study set.
if (nrow(df_final) != 137L) {
  stop("Expected exactly 137 studies, found ", nrow(df_final), ".")
}

# Check that columns required by the reported analyses are present.
required_columns <- c(
  "Study ID",
  "Publication Year",
  "Continent",
  "Country",
  "Biogeographical region",
  "species",
  "genus",
  "family",
  "order",
  "genetic/physiological",
  "physiological trait",
  "Gene(s) / molecular marker(s) studied",
  "study_system"
)

missing_columns <- setdiff(required_columns, names(df_final))

if (length(missing_columns) > 0) {
  stop(
    "Final review dataset is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Define delimiters and fixed analysis categories.
split_comma <- "\\s*[,;]\\s*"
split_region <- "\\s*[,;/|]\\s*"

valid_regions <- c(
  "Nearctic", "Neotropical", "Palearctic",
  "Afrotropical", "Indomalayan", "Australasian"
)

# Nuclear candidate-gene panel used in the bat selection analysis.
candidate_genes_nuclear <- c(
  "ACSL1", "ADCY3", "ADCY7", "ADCY8",
  "ATP5F1A", "ATP5F1B", "ATP5MC1",
  "COX4I1", "COX7A1",
  "CPT1A", "CPT1B", "CPT2",
  "CREB3L1", "CREB3L4", "CREB5",
  "FGFR1",
  "NDUFA1", "NDUFB8", "NDUFS2", "NDUFV2",
  "PLIN1", "PPARG", "PPARGC1A", "PRDM16",
  "PRKAA1", "SDHA", "SDHB",
  "SIRT6", "SLC25A20", "SOS1", "SOS2",
  "UCP1", "UQCRB"
)

# Labels retained in Table S1 but excluded from gene-frequency summaries.
non_specific_gene_labels <- c(
  "CAMP-PKA", "CRTCS-CREB", "A-COA", "PLT",
  "PDH", "PK", "PFK", "ALPHA-GPO", "COX", "CS",
  "SNRK2", "PP2C", "Mitochondrial D-loop"
)

# Figure palettes.
study_type_palette <- c(
  "Physiological" = "#3E7CB1",
  "Genetic" = "#E9B872",
  "Both" = "#6E3B8B"
)

gene_status_palette <- c(
  "In panel" = "#2A9D8F",
  "Not in panel" = "#F0A7A7"
)

study_system_palette <- c(
  "Wild / non-model" = "#3E7CB1",
  "Model organism" = "#E9B872",
  "Mixed" = "#6E3B8B",
  "Unspecified" = "#C9C9C9"
)

order_palette_map <- c(
  "Chiroptera" = "#F002E9",
  "Rodentia" = "#6E3B8B",
  "Artiodactyla" = "#E9B872",
  "Lagomorpha" = "#3F8DB3",
  "Primates" = "#74D9B8",
  "Carnivora" = "#DDA0DD",
  "Other orders" = "#BDBDBD"
)

# Shared plotting theme for review figures.
theme_publication <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(
        fill = "grey90", colour = "black", linewidth = 0.3
      ),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.title = element_blank()
    )
}

# Standardize trait labels used in the trait-group figure.
clean_trait_label <- function(x) {
  x <- str_to_lower(as.character(x))
  x <- str_replace_all(x, "-", " ")
  x <- str_squish(x)
  
  case_when(
    x %in% c(
      "nst", "non shivering thermogenesis",
      "nonshivering thermogenesis", "non shiveringthermogenesis"
    ) ~ "NST",
    x == "bat" ~ "BAT",
    x == "metabolic rate" ~ "Metabolic Rate",
    x == "metabolism" ~ "Metabolism",
    x == "body mass" ~ "Body Mass",
    x == "thermogenesis" ~ "Thermogenesis",
    x == "thermogenic capacity" ~ "Thermogenic Capacity",
    x == "torpor" ~ "Torpor",
    x == "hibernation" ~ "Hibernation",
    x == "acclimatization" ~ "Acclimatization",
    x == "acclimation" ~ "Acclimation",
    x == "bat acclimatization" ~ "Bat Acclimatization",
    x == "body temperature regulation" ~ "Body Temperature Regulation",
    x == "heat loss reduction" ~ "Heat Loss Reduction",
    x == "thermoregulation" ~ "Thermoregulation",
    x == "insulation" ~ "Insulation",
    x == "fur" ~ "Fur",
    x == "blubber" ~ "Blubber",
    x == "vasoconstriction" ~ "Vasoconstriction",
    x == "shivering" ~ "Shivering",
    x == "white adipose tissue" ~ "White Adipose Tissue",
    x == "respiratory rate" ~ "Respiratory Rate",
    TRUE ~ str_to_title(x)
  )
}

# Collapse detailed traits into the five reported functional groups.
assign_trait_group <- function(trait) {
  case_when(
    trait %in% c(
      "Metabolic Rate", "Metabolism", "Body Mass"
    ) ~ "Metabolic Regulation",
    trait %in% c(
      "NST", "BAT", "White Adipose Tissue",
      "Thermogenesis", "Thermogenic Capacity", "Shivering"
    ) ~ "Heat Production",
    trait %in% c(
      "Body Temperature Regulation", "Thermoregulation",
      "Insulation", "Fur", "Heat Loss Reduction", "Vasoconstriction"
    ) ~ "Thermoregulation",
    trait %in% c("Torpor", "Hibernation") ~ "Torpor / Hibernation",
    trait %in% c(
      "Acclimation", "Acclimatization", "Bat Acclimatization"
    ) ~ "Acclimation",
    TRUE ~ "Other"
  )
}

# Identify sequence-accession labels that should not count as gene symbols.
is_accession_label <- function(x) {
  str_detect(x, "^[A-Z]{2}[0-9]{6}$")
}

# Build pie-slice polygons used on the regional map.
draw_pie_df <- function(x, y, values, radius = 8, start = 0) {
  values <- as.numeric(values)
  
  if (all(is.na(values)) || sum(values, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  values[is.na(values)] <- 0
  props <- values / sum(values)
  ends <- cumsum(props) * 2 * pi + start
  starts <- c(start, head(ends, -1))
  
  map_dfr(seq_along(starts), function(index) {
    theta <- seq(starts[index], ends[index], length.out = 100)
    
    tibble(
      x = c(x, x + radius * cos(theta)),
      y = c(y, y + radius * sin(theta)),
      slice = index
    )
  })
}

# Export each figure in journal and editable formats.
save_figure <- function(plot, filename, width_mm, height_mm) {
  base <- file.path(FIGURE_DIR, filename)
  
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

# Expand the final gene/marker column to one token per study.
gene_tokens_all <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    gene_field = `Gene(s) / molecular marker(s) studied`
  ) %>%
  filter(!is.na(gene_field), str_squish(gene_field) != "") %>%
  separate_rows(gene_field, sep = split_comma) %>%
  transmute(
    paper_id,
    gene = str_squish(gene_field)
  ) %>%
  filter(gene != "") %>%
  distinct(paper_id, gene)

# Record labels intentionally omitted from gene-level synthesis.
gene_synthesis_exclusions <- gene_tokens_all %>%
  filter(
    gene %in% non_specific_gene_labels |
      is_accession_label(gene)
  ) %>%
  mutate(
    exclusion_reason = case_when(
      is_accession_label(gene) ~ "Sequence accession; not a gene symbol",
      gene %in% non_specific_gene_labels ~
        "Non-specific pathway/enzyme/marker label",
      TRUE ~ "Other"
    )
  )

# Retain one valid gene symbol per study for frequency calculations.
gene_by_paper <- gene_tokens_all %>%
  filter(
    !gene %in% non_specific_gene_labels,
    !is_accession_label(gene)
  ) %>%
  distinct(paper_id, gene)

# Count the number of studies reporting each gene.
gene_counts <- gene_by_paper %>%
  count(gene, sort = TRUE, name = "n")

# Count studies containing gene or molecular-marker information.
n_gene_reporting_studies <- df_final %>%
  filter(
    !is.na(`Gene(s) / molecular marker(s) studied`),
    `Gene(s) / molecular marker(s) studied` != ""
  ) %>%
  nrow()

# Build paper-level long tables used by summaries and figures.
region_by_paper <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    region_field = `Biogeographical region`
  ) %>%
  filter(
    !is.na(region_field),
    region_field != "",
    region_field != "Global"
  ) %>%
  separate_rows(region_field, sep = split_region) %>%
  mutate(region = str_squish(region_field)) %>%
  filter(region %in% valid_regions) %>%
  distinct(paper_id, region)

continent_by_paper <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    continent_field = Continent
  ) %>%
  filter(!is.na(continent_field), continent_field != "") %>%
  separate_rows(continent_field, sep = split_comma) %>%
  mutate(continent = str_squish(continent_field)) %>%
  filter(continent != "") %>%
  distinct(paper_id, continent)

country_by_paper <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    country_field = Country
  ) %>%
  filter(!is.na(country_field), country_field != "") %>%
  separate_rows(country_field, sep = split_comma) %>%
  mutate(Country = str_squish(country_field)) %>%
  filter(Country != "") %>%
  distinct(paper_id, Country)

order_by_paper <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    order_field = order
  ) %>%
  filter(!is.na(order_field), order_field != "") %>%
  separate_rows(order_field, sep = split_comma) %>%
  mutate(order = str_to_title(str_squish(order_field))) %>%
  filter(order != "") %>%
  distinct(paper_id, order)

trait_by_paper <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    trait_field = `physiological trait`
  ) %>%
  filter(!is.na(trait_field), trait_field != "") %>%
  separate_rows(trait_field, sep = split_comma) %>%
  mutate(trait = clean_trait_label(trait_field)) %>%
  filter(!is.na(trait), trait != "") %>%
  distinct(paper_id, trait)

trait_group_by_paper <- trait_by_paper %>%
  mutate(trait_group = assign_trait_group(trait)) %>%
  distinct(paper_id, trait_group)

study_system_by_paper <- df_final %>%
  transmute(
    paper_id = `Study ID`,
    study_system
  ) %>%
  filter(!is.na(study_system), study_system != "") %>%
  distinct(paper_id, study_system)

# Validate the locked review counts used in the manuscript.
expected_region_counts <- tibble(
  region = c(
    "Nearctic", "Neotropical", "Palearctic",
    "Afrotropical", "Indomalayan", "Australasian"
  ),
  expected_n = c(36L, 5L, 72L, 7L, 10L, 4L)
)

region_check <- expected_region_counts %>%
  left_join(
    region_by_paper %>% count(region, name = "observed_n"),
    by = "region"
  ) %>%
  mutate(match = expected_n == observed_n)

if (any(is.na(region_check$match)) || any(!region_check$match)) {
  stop("Regional study counts differ from the locked final dataset.")
}

n_global <- sum(
  df_final$`Biogeographical region` == "Global",
  na.rm = TRUE
)

if (n_global != 2L) {
  stop("Expected exactly two intentional Global studies.")
}

if (n_gene_reporting_studies != 52L) {
  stop(
    "Expected 52 studies with gene/marker information, found ",
    n_gene_reporting_studies, "."
  )
}

n_unique_gene_labels <- n_distinct(gene_by_paper$gene)
singleton_gene_n <- gene_counts %>%
  filter(n == 1L) %>%
  nrow()

singleton_gene_pct <- 100 * singleton_gene_n / n_unique_gene_labels

if (n_unique_gene_labels != 478L) {
  stop(
    "Expected 478 unique gene symbols, found ",
    n_unique_gene_labels, "."
  )
}

if (
  singleton_gene_n != 426L ||
  round(singleton_gene_pct, 1) != 89.1
) {
  stop(
    "Expected 426/478 singleton genes (89.1%), found ",
    singleton_gene_n, "/", n_unique_gene_labels, "."
  )
}

ppargc1a_n <- gene_counts %>%
  filter(gene == "PPARGC1A") %>%
  pull(n)

if (length(ppargc1a_n) != 1L || ppargc1a_n != 7L) {
  stop("Expected PPARGC1A in 7 studies.")
}

# Derive manuscript-level review summary statistics.
family_count <- df_final %>%
  filter(!is.na(family), family != "") %>%
  separate_rows(family, sep = split_comma) %>%
  transmute(family = str_squish(family)) %>%
  filter(family != "") %>%
  summarise(n = n_distinct(family)) %>%
  pull(n)

genus_count <- df_final %>%
  filter(!is.na(genus), genus != "") %>%
  separate_rows(genus, sep = split_comma) %>%
  transmute(genus = str_squish(genus)) %>%
  filter(genus != "") %>%
  summarise(n = n_distinct(genus)) %>%
  pull(n)

species_entry_count <- df_final %>%
  filter(!is.na(species), species != "") %>%
  separate_rows(species, sep = split_comma) %>%
  transmute(species = str_squish(species)) %>%
  filter(species != "") %>%
  nrow()

summary_stats <- tibble(
  Metric = c(
    "Total included studies",
    "Continents represented",
    "Countries represented",
    "Biogeographical regions represented in regional plots",
    "Global studies excluded from regional plots",
    "Mammalian orders represented",
    "Families represented",
    "Genera represented",
    "Species entries reported",
    "Genetic studies",
    "Physiological studies",
    "Both-type studies",
    "Wild / non-model studies",
    "Model organism studies",
    "Mixed studies",
    "Unspecified study system",
    "Studies reporting genes/markers",
    "Unique gene symbols used in gene-level synthesis",
    "Genes reported in only one study",
    "Percentage of genes reported in only one study",
    "Earliest publication year",
    "Latest publication year"
  ),
  Value = c(
    nrow(df_final),
    n_distinct(continent_by_paper$continent),
    n_distinct(country_by_paper$Country),
    n_distinct(region_by_paper$region),
    n_global,
    n_distinct(order_by_paper$order),
    family_count,
    genus_count,
    species_entry_count,
    sum(df_final$`genetic/physiological` == "Genetic", na.rm = TRUE),
    sum(df_final$`genetic/physiological` == "Physiological", na.rm = TRUE),
    sum(df_final$`genetic/physiological` == "Both", na.rm = TRUE),
    sum(df_final$study_system == "Wild / non-model", na.rm = TRUE),
    sum(df_final$study_system == "Model organism", na.rm = TRUE),
    sum(df_final$study_system == "Mixed", na.rm = TRUE),
    sum(df_final$study_system == "Unspecified", na.rm = TRUE),
    n_gene_reporting_studies,
    n_unique_gene_labels,
    singleton_gene_n,
    round(singleton_gene_pct, 1),
    min(df_final$`Publication Year`, na.rm = TRUE),
    max(df_final$`Publication Year`, na.rm = TRUE)
  )
)

# Prepare figure input tables from the final study-level dataset.
year_type_counts <- df_final %>%
  filter(
    !is.na(`Publication Year`),
    !is.na(`genetic/physiological`),
    `genetic/physiological` != ""
  ) %>%
  count(`Publication Year`, `genetic/physiological`, name = "n")

gene_counts_main <- gene_counts %>%
  slice_head(n = 12) %>%
  mutate(
    candidate_status = if_else(
      gene %in% candidate_genes_nuclear,
      "In panel",
      "Not in panel"
    ),
    candidate_status = factor(
      candidate_status,
      levels = c("In panel", "Not in panel")
    ),
    gene = factor(gene, levels = gene)
  )

order_region_long <- region_by_paper %>%
  inner_join(order_by_paper, by = "paper_id") %>%
  distinct(paper_id, region, order)

main_map_orders <- c(
  "Chiroptera", "Rodentia", "Artiodactyla",
  "Lagomorpha", "Primates", "Carnivora"
)

order_region_counts <- order_region_long %>%
  mutate(
    order_group = if_else(
      order %in% main_map_orders,
      order,
      "Other orders"
    ),
    order_group = factor(
      order_group,
      levels = c(main_map_orders, "Other orders")
    )
  ) %>%
  count(region, order_group, name = "n")

region_totals_order <- order_region_long %>%
  distinct(paper_id, region) %>%
  count(region, name = "total_studies")

trait_group_levels <- c(
  "Heat Production",
  "Metabolic Regulation",
  "Thermoregulation",
  "Torpor / Hibernation",
  "Acclimation"
)

traitgroup_region_system <- trait_group_by_paper %>%
  inner_join(region_by_paper, by = "paper_id") %>%
  inner_join(study_system_by_paper, by = "paper_id") %>%
  distinct(paper_id, trait_group, region, study_system) %>%
  count(region, trait_group, study_system, name = "n") %>%
  filter(
    trait_group %in% trait_group_levels,
    study_system != "Unspecified"
  ) %>%
  mutate(
    region = factor(region, levels = valid_regions),
    trait_group = factor(
      trait_group,
      levels = rev(trait_group_levels)
    ),
    study_system = factor(
      study_system,
      levels = c(
        "Wild / non-model",
        "Model organism",
        "Mixed"
      )
    )
  )

# Plot publication trends by evidence type.
p_year <- ggplot(
  year_type_counts,
  aes(
    x = `Publication Year`,
    y = n,
    fill = `genetic/physiological`
  )
) +
  geom_col(width = 0.85, colour = NA) +
  scale_fill_manual(
    values = study_type_palette,
    name = "Evidence type"
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(year_type_counts$`Publication Year`, na.rm = TRUE) / 5) * 5,
      ceiling(max(year_type_counts$`Publication Year`, na.rm = TRUE) / 5) * 5,
      by = 5
    )
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Year", y = "Number of studies") +
  theme_publication(base_size = 9) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot the twelve most frequently reported genes.
p_gene <- ggplot(
  gene_counts_main,
  aes(x = gene, y = n, fill = candidate_status)
) +
  geom_col(width = 0.75, colour = NA) +
  scale_fill_manual(
    values = gene_status_palette,
    name = "Candidate panel",
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Gene", y = "Number of studies") +
  theme_publication(base_size = 9) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Define plotting coordinates for regional pies and labels.
region_coords <- tribble(
  ~region, ~lon, ~lat,
  "Nearctic", -105, 50,
  "Neotropical", -60, -15,
  "Palearctic", 35, 55,
  "Afrotropical", 20, 0,
  "Indomalayan", 95, 20,
  "Australasian", 135, -25
)

region_label_coords <- tribble(
  ~region, ~label_lon, ~label_lat,
  "Nearctic", -130, 63,
  "Neotropical", -88, -14,
  "Palearctic", 12, 78,
  "Afrotropical", -8, 12,
  "Indomalayan", 122, 29,
  "Australasian", 152, -17
)

order_region_plot <- order_region_counts %>%
  left_join(region_totals_order, by = "region") %>%
  left_join(region_coords, by = "region") %>%
  mutate(
    radius = rescale(
      total_studies,
      to = c(11, 28),
      from = range(
        region_totals_order$total_studies,
        na.rm = TRUE
      )
    )
  )

pie_polys <- list()

for (region_name in unique(order_region_plot$region)) {
  region_data <- order_region_plot %>%
    filter(region == region_name) %>%
    arrange(order_group)
  
  pie_data <- draw_pie_df(
    x = region_data$lon[1],
    y = region_data$lat[1],
    values = region_data$n,
    radius = region_data$radius[1]
  )
  
  if (!is.null(pie_data)) {
    slice_sizes <- sapply(
      split(pie_data, pie_data$slice),
      nrow
    )
    
    pie_data$order_group <- rep(
      region_data$order_group,
      times = slice_sizes
    )[seq_len(nrow(pie_data))]
    
    pie_data$group_id <- paste(
      region_name,
      pie_data$slice,
      sep = "_"
    )
    
    pie_polys[[region_name]] <- pie_data
  }
}

pie_polys <- bind_rows(pie_polys)

region_labels <- region_totals_order %>%
  left_join(region_label_coords, by = "region") %>%
  mutate(
    label = paste0(
      region,
      "\n(",
      total_studies,
      ")"
    )
  )

# Draw the regional mammalian-order map.
world_map <- map_data("world")

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    colour = "white",
    linewidth = 0.15
  ) +
  geom_polygon(
    data = pie_polys,
    aes(
      x = x,
      y = y,
      group = group_id,
      fill = order_group
    ),
    colour = "grey10",
    linewidth = 0.15
  ) +
  geom_text(
    data = region_labels,
    aes(
      x = label_lon,
      y = label_lat,
      label = label
    ),
    size = 3.2,
    lineheight = 0.9
  ) +
  coord_quickmap(
    xlim = c(-180, 180),
    ylim = c(-58, 84)
  ) +
  scale_fill_manual(
    values = order_palette_map,
    breaks = c(main_map_orders, "Other orders"),
    name = "Mammalian order",
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9, colour = "black")
  )

# Combine the map, publication trend and top-gene panels.
figure_1_combined <- (
  p_map + labs(tag = "(a)")
) / (
  (p_year + labs(tag = "(b)")) |
    (p_gene + labs(tag = "(c)"))
) +
  patchwork::plot_layout(
    heights = c(1.12, 1.00)
  )

# Plot trait groups across regions and study systems.
p_traits <- ggplot(
  traitgroup_region_system,
  aes(
    x = n,
    y = trait_group,
    fill = study_system
  )
) +
  geom_col(width = 0.75, colour = NA) +
  facet_wrap(~ region, ncol = 3, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(
    values = study_system_palette[
      c(
        "Wild / non-model",
        "Model organism",
        "Mixed"
      )
    ],
    name = "Study system",
    drop = TRUE
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Number of studies",
    y = "Trait group"
  ) +
  theme_publication(base_size = 9) +
  theme(legend.position = "bottom")

# Export analysis tables used for reporting and supplementary material.
write_xlsx(
  list(
    summary_statistics = summary_stats,
    region_check = region_check,
    gene_counts = gene_counts,
    gene_synthesis_exclusions = gene_synthesis_exclusions
  ),
  file.path(TABLE_DIR, "review_analysis_tables.xlsx")
)

write.csv(
  summary_stats,
  file.path(TABLE_DIR, "summary_statistics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  gene_counts,
  file.path(TABLE_DIR, "gene_counts_for_synthesis.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  gene_synthesis_exclusions,
  file.path(TABLE_DIR, "gene_synthesis_exclusions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Export figure input tables for transparent figure reproduction.
write.csv(
  year_type_counts,
  file.path(
    FIGURE_INPUT_DIR,
    "Figure_1B_year_evidence_type_counts.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  gene_counts_main,
  file.path(
    FIGURE_INPUT_DIR,
    "Figure_1C_top_gene_counts.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  order_region_plot,
  file.path(
    FIGURE_INPUT_DIR,
    "Figure_1A_order_region_counts.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  traitgroup_region_system,
  file.path(
    FIGURE_INPUT_DIR,
    "Figure_2_traitgroup_region_studysystem_counts.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Export the final review figures.
save_figure(
  figure_1_combined,
  "Figure_1_combined_ABC",
  width_mm = 180,
  height_mm = 150
)

save_figure(
  p_map,
  "Figure_1A_region_order_map",
  width_mm = 180,
  height_mm = 100
)

save_figure(
  p_year,
  "Figure_1B_publication_trend",
  width_mm = 95,
  height_mm = 78
)

save_figure(
  p_gene,
  "Figure_1C_top_reported_genes",
  width_mm = 85,
  height_mm = 78
)

save_figure(
  p_traits,
  "Figure_2_traitgroup_region_studysystem",
  width_mm = 180,
  height_mm = 105
)