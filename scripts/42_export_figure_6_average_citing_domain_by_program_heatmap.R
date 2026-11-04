# Create a three-panel heatmap of average citing papers per original publication.
# This is the average (rather than total-count) companion to Figure 12.
# Run from the project root for Figure 13 (log10 color scale):
#   Rscript scripts/42_export_figure_13_average_citing_domain_by_program_heatmap.R
# Add --linear for Figure 6 disciplinary distribution (linear color scale):
#   Rscript scripts/42_export_figure_13_average_citing_domain_by_program_heatmap.R --linear

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)])
project_root <- if (length(script_file) == 1L) dirname(dirname(normalizePath(script_file))) else normalizePath(".")
local_r_library <- file.path(project_root, ".Rlib")
if (dir.exists(local_r_library)) .libPaths(c(local_r_library, .libPaths()))
use_log10_scale <- !("--linear" %in% commandArgs(trailingOnly = TRUE))
output_stem <- if (use_log10_scale) "figure_13" else "figure_6_disciplinary_distribution"
y_axis_label <- if (use_log10_scale) "Original Publication Group" else "Collaboration Type"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

program_order <- c("HTAN", "HuBMAP", "SenNet")
group_order <- c("Intradivision Publications", "Interdivision Publications", "Other Publications")
domain_order <- c("Life Sciences", "Health Sciences", "Physical Sciences", "Social Sciences")
domain_colours <- c(
  "Life Sciences" = "#008042",
  "Health Sciences" = "#6D99D5",
  "Physical Sciences" = "#9964A8",
  "Social Sciences" = "#D55E00"
)
x_axis_label_colour <- if (use_log10_scale) unname(domain_colours[domain_order]) else "black"
x_axis_label_face <- if (use_log10_scale) "bold" else "plain"
common_font_size <- if (use_log10_scale) 11 else 22
common_font_face <- "plain"
strip_text_face <- if (use_log10_scale) "bold" else common_font_face
na_label_size <- if (use_log10_scale) 3 else common_font_size / ggplot2::.pt
colourbar_height <- if (use_log10_scale) grid::unit(4, "cm") else grid::unit(7, "cm")
output_width <- if (use_log10_scale) 10 else 20
output_height <- if (use_log10_scale) 8.5 else 12

read_selected_csv <- function(path, columns) {
  header <- names(read.csv(path, nrows = 0, check.names = FALSE, encoding = "UTF-8"))
  missing_columns <- setdiff(columns, header)
  if (length(missing_columns) > 0) stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
  classes <- rep("NULL", length(header)); names(classes) <- header; classes[columns] <- "character"
  read.csv(path, check.names = FALSE, colClasses = classes, encoding = "UTF-8", na.strings = c("", "NA"))
}

publication_path <- file.path(project_root, "data", "processed", "NIHReporter-map_consortia-publications-validated-citations.csv")
citation_path <- file.path(project_root, "outputs", "tables", "openalex-all-publication-citing-works-with-collab-type.csv")
for (path in c(publication_path, citation_path)) if (!file.exists(path)) stop("Required input does not exist: ", path)

group_from_collab_type <- function(x) {
  case_when(
    x == "Intradivision Collaboration" ~ "Intradivision Publications",
    x == "Interdivision Collaboration" ~ "Interdivision Publications",
    x == "Single Project" ~ "Other Publications",
    TRUE ~ NA_character_
  )
}

# The denominator is every retained original publication in the specified
# program and collaboration group, including those with no citing work in a domain.
original_publications <- read_selected_csv(
  publication_path, c("consortia", "pmid", "preprint_matched", "collab_type")
) %>%
  mutate(
    consortia = factor(consortia, levels = program_order),
    publication_group = group_from_collab_type(collab_type)
  ) %>%
  filter(
    consortia %in% program_order,
    preprint_matched == "FALSE",
    !is.na(publication_group)
  ) %>%
  distinct(consortia, pmid, publication_group, .keep_all = TRUE)

publication_denominators <- original_publications %>%
  count(consortia, publication_group, name = "number_original_publications")

# Numerator: same citation eligibility and OpenAlex-work de-duplication as Figure 10.
citing_edges <- read_selected_csv(
  citation_path,
  c("original_pmid", "citing_openalex_id", "journal", "field", "domain",
    "tag_consortia", "tag_collab_type", "collab_type_match_status")
) %>%
  mutate(
    consortia = factor(tag_consortia, levels = program_order),
    publication_group = group_from_collab_type(tag_collab_type)
  ) %>%
  filter(
    collab_type_match_status == "matched",
    consortia %in% program_order,
    !is.na(citing_openalex_id),
    !is.na(field), nzchar(trimws(field)), tolower(trimws(field)) != "nan",
    !is.na(domain), nzchar(trimws(domain)), domain %in% domain_order,
    !is.na(journal), nzchar(trimws(journal)), tolower(trimws(journal)) != "nan",
    !is.na(publication_group)
  ) %>%
  transmute(consortia, pmid = original_pmid, citing_openalex_id, domain, publication_group) %>%
  distinct(consortia, pmid, citing_openalex_id, .keep_all = TRUE) %>%
  inner_join(
    original_publications %>% select(consortia, pmid, publication_group),
    by = c("consortia", "pmid", "publication_group"),
    relationship = "many-to-one"
  )

domain_counts <- citing_edges %>%
  group_by(consortia, publication_group, domain) %>%
  summarise(number_citing_papers = n_distinct(citing_openalex_id), .groups = "drop")

all_cells <- expand.grid(
  consortia = program_order,
  publication_group = group_order,
  domain = domain_order,
  stringsAsFactors = FALSE
)
heatmap_data <- all_cells %>%
  left_join(domain_counts, by = c("consortia", "publication_group", "domain")) %>%
  left_join(publication_denominators, by = c("consortia", "publication_group")) %>%
  mutate(
    number_citing_papers = coalesce(number_citing_papers, 0L),
    # A program/group with no retained original publications has no defined mean.
    average_citing_papers_per_original_publication = if_else(
      is.na(number_original_publications) | number_original_publications == 0,
      NA_real_,
      number_citing_papers / number_original_publications
    ),
    consortia = factor(consortia, levels = program_order),
    publication_group = factor(publication_group, levels = group_order),
    domain = factor(domain, levels = domain_order)
  )

table_path <- file.path(project_root, "outputs", "tables", "average-citing-paper-domain-by-program-and-collaboration-group-heatmap.csv")
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
write.csv(heatmap_data, table_path, row.names = FALSE, na = "")

max_average <- max(heatmap_data$average_citing_papers_per_original_publication, na.rm = TRUE)
fill_breaks <- c(0, 1, 10, 100, 1000)
fill_breaks <- fill_breaks[fill_breaks <= max_average]
if (!(0 %in% fill_breaks)) fill_breaks <- c(0, fill_breaks)
linear_fill_breaks <- pretty(c(0, max_average), n = 5)
linear_fill_breaks <- linear_fill_breaks[linear_fill_breaks >= 0 & linear_fill_breaks <= max_average]
if (!(0 %in% linear_fill_breaks)) linear_fill_breaks <- c(0, linear_fill_breaks)

domain_key <- data.frame(
  consortia = factor(program_order[1], levels = program_order),
  publication_group = factor(group_order[1], levels = group_order),
  domain = factor(domain_order, levels = domain_order)
)

# Figure 14 omits the HTAN intradivision row because no retained HTAN original
# publication belongs to that group; the row remains in the source-data table.
plot_heatmap_data <- if (use_log10_scale) {
  heatmap_data
} else {
  heatmap_data %>%
    filter(!(consortia == "HTAN" & publication_group == "Intradivision Publications")) %>%
    mutate(publication_group = factor(publication_group, levels = rev(group_order)))
}
plot_domain_key <- if (use_log10_scale) domain_key else plot_heatmap_data[0, c("consortia", "publication_group", "domain")]

fill_scale <- if (use_log10_scale) {
  scale_fill_gradient(
    low = "#F7FBFF", high = "#084081",
    trans = "log10",
    na.value = "#D9D9D9",
    breaks = fill_breaks,
    labels = format(fill_breaks, big.mark = ",", trim = TRUE),
    name = "Average Number of Citing Papers\nper Original Publication\n(log10 scale)"
  )
} else {
  scale_fill_gradient(
    low = "#F7FBFF", high = "#084081",
    na.value = "#D9D9D9",
    breaks = linear_fill_breaks,
    labels = format(linear_fill_breaks, big.mark = ",", trim = TRUE),
    name = "Average Citations"
  )
}

domain_legend_guide <- if (use_log10_scale) {
  guide_legend(order = 1, override.aes = list(alpha = 1, shape = 15, size = 4))
} else {
  "none"
}

figure_13 <- ggplot(
  plot_heatmap_data,
  aes(x = domain, y = publication_group, fill = average_citing_papers_per_original_publication)
) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    data = filter(plot_heatmap_data, is.na(average_citing_papers_per_original_publication)),
    aes(x = domain, y = publication_group, label = "N/A"),
    inherit.aes = FALSE, colour = "#555555", size = na_label_size
  ) +
  geom_point(
    data = plot_domain_key,
    aes(x = domain, y = publication_group, colour = domain),
    inherit.aes = FALSE, alpha = 0, show.legend = TRUE
  ) +
  facet_wrap(~consortia, ncol = 1, scales = if (use_log10_scale) "fixed" else "free_y") +
  scale_x_discrete(labels = domain_order) +
  scale_y_discrete(limits = if (use_log10_scale) rev(group_order) else NULL, drop = TRUE) +
  fill_scale +
  scale_colour_manual(values = domain_colours, name = "Citing-paper Domain") +
  labs(x = "Domain", y = y_axis_label) +
  theme_minimal(base_size = common_font_size) +
  theme(
    axis.text.x = element_text(face = x_axis_label_face, size = common_font_size, colour = x_axis_label_colour),
    axis.text.y = element_text(face = common_font_face, size = common_font_size),
    axis.title = element_text(face = common_font_face, size = common_font_size),
    strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text = element_text(face = strip_text_face, size = common_font_size),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(face = common_font_face, size = common_font_size),
    legend.title = element_text(face = common_font_face, size = common_font_size),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)
  ) +
  guides(
    colour = domain_legend_guide,
    fill = guide_colorbar(order = 2, barheight = colourbar_height)
  )

figure_dir <- file.path(project_root, "outputs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, paste0(output_stem, ".svg")), figure_13, device = "svg", width = output_width, height = output_height, units = "in")
ggsave(file.path(figure_dir, paste0(output_stem, ".jpg")), figure_13, device = "jpeg", width = output_width, height = output_height, units = "in", dpi = 300)
message("Created ", output_stem, " and ", basename(table_path), ".")
