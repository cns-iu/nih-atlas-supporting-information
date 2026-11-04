#!/usr/bin/env Rscript

# Generate table_3_atlas_program_authors_and participation.csv.
# This is the standalone version of the author-participation calculation
# previously embedded in notebooks/data_analysis/NIH_map-publication_analysis.Rmd.

source(here::here("scripts/00_setup.R"))

library(tidyverse)

programs <- c("HTAN", "HuBMAP", "SenNet")
network_dir <- file.path(project_root, "data", "networks")
output_file <- file.path(
  project_root,
  "outputs", "tables",
  "table_3_atlas_program_authors_and participation.csv"
)

read_nodes <- function(program, suffix) {
  path <- file.path(network_dir, paste0(program, suffix))
  read_csv(path, show_col_types = FALSE) %>%
    mutate(consortia = program)
}

baseline_nodes <- map_dfr(
  programs,
  read_nodes,
  suffix = "-coauthorship-baseline-nodes.csv"
)

collaboration_nodes <- map_dfr(
  programs,
  read_nodes,
  suffix = "-coauthorship-collaboration_partition-nodes.csv"
)

# Total authors by program and role.
author_role_count <- baseline_nodes %>%
  group_by(consortia, role) %>%
  summarise(n_authors = n(), .groups = "drop") %>%
  mutate(role = str_replace_all(role, " ", "_")) %>%
  pivot_wider(
    id_cols = consortia,
    names_from = role,
    names_prefix = "Count_of_",
    values_from = n_authors,
    values_fill = 0
  ) %>%
  mutate(Count_of_Authors = Count_of_Principal_Investigator + Count_of_Researcher) %>%
  rename(
    Count_of_PIs = Count_of_Principal_Investigator,
    Count_of_Researchers = Count_of_Researcher
  )

# Authors with at least one internal-collaboration publication, by role.
author_role_counts_ic <- collaboration_nodes %>%
  filter(ic_publications > 0) %>%
  group_by(consortia, role) %>%
  summarise(Count_of_Authors = n(), .groups = "drop") %>%
  mutate(role = str_replace_all(role, " ", "_") %>% tolower()) %>%
  pivot_wider(
    id_cols = consortia,
    names_from = role,
    names_prefix = "Count_of_",
    values_from = Count_of_Authors,
    values_fill = 0
  ) %>%
  mutate(
    Count_of_Authors_involved_in_Internal_Collaborative_Research =
      Count_of_principal_investigator + Count_of_researcher
  ) %>%
  rename(
    Count_PIs_involved_in_Internal_Collaborative_Research = Count_of_principal_investigator,
    Count_Researchers_involved_in_Internal_Collaborative_Research = Count_of_researcher
  )

table_3 <- author_role_count %>%
  left_join(author_role_counts_ic, by = "consortia") %>%
  mutate(
    Percent_of_PIs_involved_in_Collaborative_Research = round(
      Count_PIs_involved_in_Internal_Collaborative_Research / Count_of_PIs * 100, 2
    ),
    Percent_of_Researchers_involved_in_Collaborative_Research = round(
      Count_Researchers_involved_in_Internal_Collaborative_Research /
        Count_of_Researchers * 100, 2
    ),
    Percent_of_Authors_involved_in_Collaborative_Research = round(
      Count_of_Authors_involved_in_Internal_Collaborative_Research /
        Count_of_Authors * 100, 2
    )
  ) %>%
  rename(
    `Atlas Program` = consortia,
    `Count of Authors` = Count_of_Authors,
    `Count of IC Authors` = Count_of_Authors_involved_in_Internal_Collaborative_Research,
    `Percent of Authors in IC publications` = Percent_of_Authors_involved_in_Collaborative_Research,
    `Count of PI Authors` = Count_of_PIs,
    `Count of PI IC Authors` = Count_PIs_involved_in_Internal_Collaborative_Research,
    `Percent of PI in IC publications` = Percent_of_PIs_involved_in_Collaborative_Research,
    `Percent of Other researcher in IC publications` = Percent_of_Researchers_involved_in_Collaborative_Research
  ) %>%
  select(
    `Atlas Program`, `Count of Authors`, `Count of IC Authors`,
    `Percent of Authors in IC publications`, `Count of PI Authors`,
    `Count of PI IC Authors`, `Percent of PI in IC publications`,
    `Percent of Other researcher in IC publications`
  ) %>%
  arrange(match(`Atlas Program`, programs))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(table_3, output_file, row.names = FALSE, fileEncoding = "UTF8")
print(table_3)
message("Saved: ", output_file)
