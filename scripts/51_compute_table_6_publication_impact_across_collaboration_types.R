#!/usr/bin/env Rscript

# Generate Table 6: publication impact across collaboration types.
#
# The publication cohort and RCR statistics come from the classified publication
# records. Citation-count summaries use the versioned OpenAlex refresh saved in
# results/table7_openalex_updated_20261005_citation_metrics.csv. This keeps the
# standalone Table 6 generator aligned with the current citation snapshot.

source(here::here("scripts/00_setup.R"))

library(tidyverse)

programs <- c("HTAN", "HuBMAP", "SenNet")
groups <- c("All Publications", "Internal Collaborative Publications", "Other Publications")
output_file <- file.path(
  project_root,
  "outputs", "tables",
  "table_6_publication_impact_across_collaboration_types.csv"
)

classified <- read.csv(
  file.path(project_root, "results", "publication_functional_classification.csv"),
  header = TRUE,
  encoding = "UTF8",
  colClasses = c(pmid = "character")
)

# Match the legacy Table 7/Table 6 eligible cohort: journal publications only,
# excluding preprints and matched preprints.
classified <- classified %>%
  mutate(
    is_preprint_bool = as.logical(is_preprint_bool),
    preprint_matched = as.logical(preprint_matched)
  ) %>%
  filter(
    !is_preprint_bool,
    !preprint_matched,
    !is.na(journal),
    str_trim(journal) != ""
  ) %>%
  mutate(
    Publication_Group = if_else(
      collab_type == "Single Project",
      "Other Publications",
      "Internal Collaborative Publications"
    )
  )

# Citation-count metrics were refreshed from OpenAlex and are retained as a
# separate versioned input so rerunning this script does not silently revert to
# the older WoS/Google Scholar citation field.
openalex_metrics <- read.csv(
  file.path(project_root, "results", "table7_openalex_updated_20261005_citation_metrics.csv"),
  header = TRUE,
  encoding = "UTF8"
)

# Calculate RCR summaries for the same eligible cohort.
rcr_metrics <- classified %>%
  mutate(
    Publication_Group = factor(Publication_Group, levels = groups),
    Mapping_Program = factor(consortia, levels = programs)
  ) %>%
  group_by(Mapping_Program, Publication_Group) %>%
  summarise(
    Mean_Relative_Citation_Ratio = round(mean(relative_citation_ratio, na.rm = TRUE), 2),
    SD_Relative_Citation_Ratio = round(sd(relative_citation_ratio, na.rm = TRUE), 3),
    Var_Relative_Citation_Ratio = round(var(relative_citation_ratio, na.rm = TRUE), 3),
    .groups = "drop"
  )

table_6 <- openalex_metrics %>%
  mutate(
    Mapping_Program = factor(Mapping_Program, levels = programs),
    Publication_Group = factor(Publication_Group, levels = groups)
  ) %>%
  left_join(rcr_metrics, by = c("Mapping_Program", "Publication_Group")) %>%
  arrange(Mapping_Program, Publication_Group) %>%
  mutate(
    Mapping_Program = as.character(Mapping_Program),
    Publication_Group = as.character(Publication_Group)
  ) %>%
  select(
    Mapping_Program, Publication_Group,
    Total_Citations, Mean_Citations, SD_Citations, Var_Citations,
    Mean_Relative_Citation_Ratio, SD_Relative_Citation_Ratio,
    Var_Relative_Citation_Ratio
  )

stopifnot(nrow(table_6) == 9)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(table_6, output_file, row.names = FALSE, fileEncoding = "UTF8")
print(table_6)
message("Saved: ", output_file)
