# Re-create Figure 3: temporal publication trends by Mapping Program.

.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))

library(dplyr)
library(ggplot2)

project_root <- normalizePath(getwd(), winslash = "/")
program_order <- c("HTAN", "HuBMAP", "SenNet")
consortium_colours <- c(HuBMAP = "#6D99D5", HTAN = "#008042", SenNet = "#9964A8")

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
  file.path(project_root, "data", "processed",
            "NIHReporter-map_consortia-publications-validated-citations.csv"),
  colClasses = c(pmid = "character"), check.names = FALSE
) %>%
  mutate(
    consortia = factor(consortia, levels = program_order),
    division_collaboration = if_else(
      collab_type == "Single Project",
      "Other Publications",
      "Internal Collaborative Publications"
    )
  ) %>%
  filter(preprint_matched == FALSE)

# Program year is the ordinal publication year within each mapping program.
project_year <- publications %>%
  group_by(consortia) %>%
  select(consortia, pub_year) %>%
  distinct() %>%
  mutate(project_year = row_number())

publications_by_year_overall <- publications %>%
  filter(pub_year != max(pub_year)) %>%
  left_join(project_year, by = c("consortia", "pub_year")) %>%
  group_by(consortia, project_year) %>%
  summarise(total_publications = n(), .groups = "drop") %>%
  mutate(division_collaboration = "All Publications")

publications_by_year_ic <- publications %>%
  filter(pub_year != max(pub_year)) %>%
  left_join(project_year, by = c("consortia", "pub_year")) %>%
  group_by(consortia, project_year, division_collaboration) %>%
  summarise(total_publications = n(), .groups = "drop") %>%
  filter(division_collaboration != "Other Publications")

publications_by_year <- bind_rows(
  publications_by_year_overall,
  publications_by_year_ic
) %>%
  mutate(
    division_collaboration = factor(
      division_collaboration,
      levels = c("All Publications", "Internal Collaborative Publications")
    )
  ) %>%
  arrange(consortia, division_collaboration)

figure_3 <- ggplot(
  publications_by_year,
  aes(
    x = project_year, y = total_publications, colour = consortia,
    linetype = division_collaboration
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  facet_grid(~ consortia) +
  scale_linetype_manual(
    values = c(
      "All Publications" = "solid",
      "Internal Collaborative Publications" = "2121"
    )
  ) +
  scale_x_continuous(breaks = seq(1, max(project_year$project_year), 1)) +
  scale_colour_manual(values = consortium_colours) +
  labs(
    x = "Publication Year",
    y = "Number of Publications",
    colour = "Atlas Program",
    linetype = "Publication Group"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = grid::unit(0.5, "pt"),
    panel.grid.minor.x = element_blank()
  )

output_dir <- file.path(project_root, "outputs", "figures")
ggsave(file.path(output_dir, "figure_3_temporal_trends.svg"), figure_3,
       device = "svg", width = 6.5, height = 2.75, units = "in")
ggsave(file.path(output_dir, "figure_3_temporal_trends.jpg"), figure_3,
       device = "jpeg", width = 6.5, height = 2.75, units = "in", dpi = 300)
message("Wrote Figure 3 temporal trends to: ", output_dir)
