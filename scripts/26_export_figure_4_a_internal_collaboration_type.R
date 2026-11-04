#!/usr/bin/env Rscript

# Export Figure 4A: internal collaboration type scores by consortium component.
# This is the standalone version of the corresponding figure code in
# notebooks/data_analysis/NIH_map-publication_analysis.Rmd.

source(here::here("scripts/00_setup.R"))

library(tidyverse)

options(scipen = 999)

output_dir <- here::here("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

consortium_colours <- c(
  HuBMAP = "#6D99D5",
  HTAN   = "#008042",
  SenNet = "#9964A8"
)

publications <- read.csv(
  file.path(project_root, dirs[4],
            "NIHReporter-map_consortia-publications-validated-citations.csv"),
  header = TRUE,
  encoding = "UTF8",
  colClasses = c(pmid = "character", pmc_id = "character")
)

publication_metadata <- read.csv(
  file.path(project_root, dirs[2],
            "NIH-map_consortia-publication-metadata.csv"),
  header = TRUE,
  encoding = "UTF8",
  colClasses = c(pmid = "character")
)

awards <- read.csv(
  file.path(project_root, dirs[4],
            "NIHReporter-map_consortia-awards-clean.csv"),
  header = TRUE,
  encoding = "UTF8"
)

# Reproduce the publication metadata preparation used by the Rmd analysis.
publication_metadata <- publication_metadata %>%
  left_join(
    publications %>%
      select(consortia, pmid, is_preprint, preprint_matched, collab_type) %>%
      distinct(),
    by = c("consortia", "pmid")
  ) %>%
  mutate(
    division_collaboration = ifelse(collab_type == "Single Project", 0, 1),
    division_collaboration = factor(
      division_collaboration,
      levels = c(1, 0),
      labels = c("Internal Collaborative Publications", "Other Publications")
    ),
    consortia = factor(consortia, levels = c("HTAN", "HuBMAP", "SenNet")),
    division = factor(
      division,
      levels = c("DCC", "HTA", "PCA", "HIVE", "DEMO", "TMC_H", "TTD",
                 "RIT", "CODCC", "TMC_S", "TDA")
    )
  ) %>%
  filter(preprint_matched == FALSE)

awards <- awards %>%
  select(project_id, division, project_serial_number, core_project_number) %>%
  distinct() %>%
  rename(consortia = project_id) %>%
  mutate(
    consortia = factor(consortia, levels = c("HTAN", "HuBMAP", "SenNet")),
    division = factor(
      division,
      levels = c("DCC", "HTA", "PCA", "HIVE", "DEMO", "TMC_H", "TTD",
                 "RIT", "CODCC", "TMC_S", "TDA")
    )
  )

projects_overall <- awards %>%
  group_by(consortia, division) %>%
  summarise(Total_Projects = n(), .groups = "drop")

division_projects_overall <- publication_metadata %>%
  select(consortia, division, core_project_number) %>%
  distinct() %>%
  group_by(consortia, division) %>%
  summarise(Count_of_All_Projects = n(), .groups = "drop") %>%
  left_join(projects_overall, by = c("consortia", "division")) %>%
  mutate(Percent_of_All_Projects = round(Count_of_All_Projects / Total_Projects * 100, 2)) %>%
  select(consortia, division, Count_of_All_Projects, Percent_of_All_Projects)

division_projects_grouped <- publication_metadata %>%
  select(consortia, division, division_collaboration, core_project_number) %>%
  distinct() %>%
  group_by(consortia, division, division_collaboration) %>%
  summarise(Count_of_Projects_with_Internal_Collaboration = n(), .groups = "drop") %>%
  left_join(projects_overall, by = c("consortia", "division")) %>%
  mutate(Percent_of_Projects_with_Internal_Collaboration =
           round(Count_of_Projects_with_Internal_Collaboration / Total_Projects * 100, 2)) %>%
  filter(division_collaboration == "Internal Collaborative Publications") %>%
  select(consortia, division, Count_of_Projects_with_Internal_Collaboration,
         Percent_of_Projects_with_Internal_Collaboration)

division_projects <- left_join(
  division_projects_overall,
  division_projects_grouped,
  by = c("consortia", "division")
)

division_overall_publications <- publication_metadata %>%
  select(consortia, division, pmid) %>%
  group_by(consortia, division) %>%
  summarise(Count_of_Publications = n_distinct(pmid), .groups = "drop")

division_grouped_publications <- publication_metadata %>%
  select(consortia, division, division_collaboration, pmid) %>%
  mutate(division_collaboration = str_replace_all(as.character(division_collaboration), " ", "_")) %>%
  group_by(consortia, division, division_collaboration) %>%
  summarise(Count_of_Publications = n_distinct(pmid), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(consortia, division),
    names_from = division_collaboration,
    names_prefix = "Count_of_",
    values_from = Count_of_Publications,
    values_fill = 0
  ) %>%
  select(-any_of("Count_of_Other_Publications"))

division_collab_type_publications <- publication_metadata %>%
  select(consortia, division, collab_type, pmid) %>%
  mutate(collab_type = str_replace_all(as.character(collab_type), " ", "_")) %>%
  group_by(consortia, division, collab_type) %>%
  summarise(Count_of_Publications = n_distinct(pmid), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(consortia, division),
    names_from = collab_type,
    names_prefix = "Count_of_",
    values_from = Count_of_Publications,
    values_fill = 0
  ) %>%
  select(-any_of("Count_of_Single_Project"))

division_publication_statistics <- division_projects %>%
  left_join(division_overall_publications, by = c("consortia", "division")) %>%
  left_join(division_grouped_publications, by = c("consortia", "division")) %>%
  left_join(division_collab_type_publications, by = c("consortia", "division")) %>%
  mutate(
    Percent_of_Internal_Collaborations = round(
      Count_of_Internal_Collaborative_Publications / Count_of_Publications * 100, 2
    ),
    Percent_of_Interdivision_Collaborations = round(
      Count_of_Interdivision_Collaboration /
        Count_of_Internal_Collaborative_Publications * 100, 2
    ),
    Percent_of_Intradivision_Collaborations = round(
      Count_of_Intradivision_Collaboration /
        Count_of_Internal_Collaborative_Publications * 100, 2
    ),
    Internal_Collaboration_Type_Score = round(
      (Percent_of_Interdivision_Collaborations -
         Percent_of_Intradivision_Collaborations) / 100, 2
    ),
    Division = recode(as.character(division), RIT = "RTI")
  )

gradient_data <- data.frame(
  collaboration_type_concentration = seq(-1.07, 1.07, 0.001),
  color = seq(-1.07, 1.07, 0.001)
)

figure_4_a <- division_publication_statistics %>%
  ggplot(aes(x = Internal_Collaboration_Type_Score,
             y = forcats::fct_rev(Division))) +
  geom_segment(
    data = gradient_data,
    alpha = 0.15,
    aes(x = color, xend = color, y = -Inf, yend = Inf, color = color),
    show.legend = FALSE
  ) +
  scale_color_gradient2(low = "#ffeda0", mid = "white", high = "#87C980", midpoint = 0) +
  annotate("text", x = -0.85, y = "DCC", size = 1.65, color = "#6B6443",
           label = "Intradivision\nCollaborations") +
  annotate("text", x = 0.85, y = "DCC", size = 1.65, color = "#3F6B40",
           label = "Interdivision\nCollaborations") +
  annotate("text", x = 0, y = "DCC", size = 1.75, color = "#444444",
           label = "Mixed Collaboration Types") +
  ggnewscale::new_scale_color() +
  geom_point(aes(colour = consortia, size = Count_of_Internal_Collaborative_Publications)) +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = seq(-1, 1, 0.1), expand = c(0, 0)) +
  scale_color_manual(values = consortium_colours) +
  scale_size_area(max_size = 3) +
  labs(
    x = "Internal Collaboration Type Score",
    y = "Consortium Component",
    size = "Collaborative Publications"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = unit(8, "pt")),
    legend.text = element_text(size = unit(6, "pt")),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = unit(6, "pt")),
    axis.text.y = element_text(size = unit(6, "pt")),
    axis.title = element_text(size = unit(8, "pt"))
  )

ggsave(file.path(output_dir, "figure_4_a.svg"), figure_4_a,
       device = "svg", width = 5, height = 3, units = "in")
ggsave(file.path(output_dir, "figure_4_a.jpg"), figure_4_a,
       device = "jpeg", width = 5, height = 3, units = "in", dpi = 300)

message("Saved Figure 4A SVG and JPG to ", output_dir)
