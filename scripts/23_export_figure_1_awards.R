# Re-export the current five-panel Figure 1 from the award data. The code
# preserves the R Markdown theme and fixes the E-panel Mapping Program order.

.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

project_root <- normalizePath(getwd(), winslash = "/")
program_order <- c("HTAN", "HuBMAP", "SenNet")
consortium_colours <- c(HuBMAP = "#6D99D5", HTAN = "#008042", SenNet = "#9964A8")

# This is the original global theme defined in NIH_map-award_funding_analysis.Rmd.
theme_atlas <- function(base_size = 24) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(colour = "#555555", size = base_size),
      axis.text = element_text(colour = "black"),
      legend.position = "right",
      strip.background = element_rect(fill = "#f0f0f0", colour = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "#eeeeee"),
      plot.caption = element_text(colour = "#444444", size = 18, hjust = 0)
    )
}
theme_set(theme_atlas())

awards_files <- file.path(
  project_root, "data", "processed",
  paste0("NIHReporter-", program_order, "-awards-clean.csv")
)
awards <- bind_rows(lapply(awards_files, function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
           colClasses = c(duns_number = "character"))
})) %>%
  rename(consortium = project_id) %>%
  mutate(
    consortium = factor(consortium, levels = program_order),
    is_subproject = !is.na(subproject_number) & subproject_number != ""
  )

inventory <- awards %>%
  group_by(consortium) %>%
  summarise(
    core_projects = n_distinct(core_project_number),
    total_awards = n_distinct(application_id),
    parent_awards = n_distinct(project_number),
    sub_awards = n_distinct(project_number[is_subproject]),
    .groups = "drop"
  )

inv_total <- inventory %>%
  select(consortium, total_awards) %>%
  pivot_longer(-consortium, names_to = "type", values_to = "count")
inv_long <- inventory %>%
  select(consortium, core_projects, parent_awards, sub_awards) %>%
  pivot_longer(-consortium, names_to = "type", values_to = "count")

bar_panel <- function(data, panel_tag, y_label) {
  ggplot(data, aes(x = consortium, y = count, fill = consortium)) +
    geom_col(position = "dodge", width = 0.55, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = consortium_colours) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
    labs(tag = panel_tag, x = NULL, y = y_label, fill = "Mapping Program") +
    guides(fill = guide_legend(order = 1)) +
    theme(legend.position = "bottom", legend.direction = "vertical")
}

p1 <- bar_panel(inv_total, "A", "RePORTER Project Records")
p2 <- bar_panel(filter(inv_long, type == "core_projects"), "B", "Core Project Numbers")
p_parent <- bar_panel(filter(inv_long, type == "parent_awards"), "C", "Full Project Numbers")
p_sub <- bar_panel(filter(inv_long, type == "sub_awards"), "D", "Sub Project Numbers")

funding_real <- awards %>%
  filter(!is_subproject) %>%
  group_by(consortium) %>%
  summarise(Amount = sum(total_cost_ic, na.rm = TRUE), .groups = "drop")

p_total <-
  ggplot(funding_real,
         aes(x = factor(consortium, levels = program_order), y = Amount, fill = consortium)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = consortium_colours) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1),
                     limits = c(0, 225000000)) +
  labs(tag = "E", x = NULL, y = "Total Cost (USD, million)", fill = "Mapping Program") +
  guides(fill = guide_legend(order = 1)) +
  theme(legend.position = "bottom", legend.direction = "vertical")

fig1 <-
  p1 + p2 + p_parent + p_sub + p_total +
  plot_layout(ncol = 5) &
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.tag = element_text(size = 24)
  )

output_dir <- file.path(project_root, "outputs", "figures")
ggsave(file.path(output_dir, "figure_1.svg"),
       plot = fig1, device = "svg", width = 20, height = 5, units = "in")
ggsave(file.path(output_dir, "figure_1.jpg"),
       plot = fig1, device = "jpeg", width = 20, height = 5, units = "in", dpi = 300)
message("Updated Figure 1. E-panel order: ", paste(program_order, collapse = ", "))
