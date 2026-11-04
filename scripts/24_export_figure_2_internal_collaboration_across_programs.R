# Export Figure 2: internal collaboration across Mapping Programs.
#
# Run from the project root:
#   Rscript scripts/29_export_figure_2_internal_collaboration_across_programs.R
#
# Inputs are the three publication/project tables and the six coauthorship-node
# tables listed in outputs/sharing/figure_6-reproducibility/README.md.

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)])
project_root <- if (length(script_file) == 1L) {
  dirname(dirname(normalizePath(script_file)))
} else {
  normalizePath(".")
}

local_r_library <- file.path(project_root, ".Rlib")
if (dir.exists(local_r_library)) {
  .libPaths(c(local_r_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

data_path <- function(...) file.path(project_root, "data", ...)
figure_path <- function(...) file.path(project_root, "outputs", "figures", ...)

program_order <- c("HTAN", "HuBMAP", "SenNet")
program_colours <- c(HTAN = "#008042", HuBMAP = "#6D99D5", SenNet = "#9964A8")

theme_atlas <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(colour = "#555555", size = base_size),
      axis.text = element_text(colour = "black"),
      legend.position = "right",
      strip.background = element_rect(fill = "#f0f0f0", colour = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "#eeeeee"),
      plot.caption = element_text(colour = "#444444", size = 9, hjust = 0)
    )
}

theme_set(theme_atlas())

publications <- read.csv(
  data_path("processed", "NIHReporter-map_consortia-publications-validated-citations.csv"),
  encoding = "UTF-8", colClasses = c(pmid = "character", pmc_id = "character")
) %>%
  mutate(
    consortia = factor(consortia, levels = program_order),
    division_collaboration = ifelse(
      collab_type == "Single Project", "Other Publications",
      "Internal Collaborative Publications"
    )
  ) %>%
  filter(preprint_matched == FALSE)

publication_metadata <- read.csv(
  data_path("metadata", "NIH-map_consortia-publication-metadata.csv"),
  encoding = "UTF-8", colClasses = c(pmid = "character")
) %>%
  left_join(
    publications %>% select(consortia, pmid, division_collaboration) %>% distinct(),
    by = c("consortia", "pmid")
  )

awards <- read.csv(
  data_path("processed", "NIHReporter-map_consortia-awards-clean.csv"),
  encoding = "UTF-8"
) %>%
  select(
    project_id, opportunity_number, division, division_type, ic,
    project_serial_number, core_project_number
  ) %>%
  distinct() %>%
  transmute(consortia = factor(project_id, levels = program_order))

# Historical Figure 6 denominator: unique award records after the selected
# award metadata fields are de-duplicated. This is intentionally not the later
# Table 6 core-project-number denominator.
projects_overall <- awards %>%
  count(consortia, name = "Total_Projects")

project_participation <- publication_metadata %>%
  filter(!is.na(division_collaboration)) %>%
  select(consortia, division_collaboration, core_project_number) %>%
  distinct() %>%
  group_by(consortia, division_collaboration) %>%
  summarise(Count_of_Projects = n(), .groups = "drop") %>%
  left_join(projects_overall, by = "consortia") %>%
  mutate(Percent_of_Projects = round(100 * Count_of_Projects / Total_Projects, 2))

publication_participation <- publications %>%
  count(consortia, division_collaboration, name = "Count_of_Publications") %>%
  group_by(consortia) %>%
  mutate(Percent_of_Publications = round(100 * Count_of_Publications / sum(Count_of_Publications), 2)) %>%
  ungroup()

read_nodes <- function(program, file_suffix) {
  read.csv(data_path("networks", paste0(program, file_suffix)), encoding = "UTF-8") %>%
    mutate(consortia = factor(program, levels = program_order))
}

baseline_nodes <- bind_rows(lapply(
  program_order,
  read_nodes,
  file_suffix = "-coauthorship-baseline-nodes.csv"
))
collaboration_nodes <- bind_rows(lapply(
  program_order,
  read_nodes,
  file_suffix = "-coauthorship-collaboration_partition-nodes.csv"
))

pi_participation <- baseline_nodes %>%
  select(consortia, id, role) %>%
  left_join(
    collaboration_nodes %>% select(consortia, id, ic_publications),
    by = c("consortia", "id")
  ) %>%
  mutate(ic_publications = coalesce(ic_publications, 0)) %>%
  filter(role == "Principal Investigator") %>%
  group_by(consortia) %>%
  summarise(
    Percent_of_PIs = round(100 * sum(ic_publications > 0) / n(), 2),
    .groups = "drop"
  )

format_percentage <- function(x, digits) {
  formatted <- formatC(round(x, digits), format = "f", digits = digits)
  sub("[.]0+$", "", formatted)
}

make_bar_plot <- function(data, value_column, tag, y_label, label_digits) {
  ggplot(data, aes(x = factor(consortia, levels = program_order), y = .data[[value_column]], fill = consortia)) +
    geom_col(alpha = 0.95) +
    geom_text(aes(label = paste0(format_percentage(.data[[value_column]], label_digits), "%")), vjust = -0.5, fontface = "bold") +
    scale_fill_manual(values = program_colours) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(tag = tag, y = y_label, fill = "Mapping Program") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

panel_a <- publication_participation %>%
  filter(division_collaboration == "Internal Collaborative Publications") %>%
  make_bar_plot("Percent_of_Publications", "A", "Percentage of Publications (%)", 1)
panel_b <- project_participation %>%
  filter(division_collaboration == "Internal Collaborative Publications") %>%
  make_bar_plot("Percent_of_Projects", "B", "Percentage of Projects (%)", 1)
panel_c <- pi_participation %>%
  make_bar_plot("Percent_of_PIs", "C", "Percentage of PIs (%)", 2)

figure_2 <- panel_a + panel_b + panel_c
dir.create(dirname(figure_path()), recursive = TRUE, showWarnings = FALSE)

ggsave(
  figure_path("figure_2_internal_collaboration_across_programs.svg"), figure_2,
  width = 6.5, height = 2.75, units = "in", device = "svg"
)
ggsave(
  figure_path("figure_2_internal_collaboration_across_programs.jpg"), figure_2,
  width = 6.5, height = 2.75, units = "in", dpi = 300
)
