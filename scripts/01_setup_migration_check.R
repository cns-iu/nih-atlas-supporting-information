# Test that files are detected
source("scripts/00_setup.R")

# Should find your files
awards <- load_nih_data("awards", project_id = "HuBMAP")
pubs <- load_nih_data("pubs", project_id = "SenNet")

# Should load all projects at once
all_awards <- load_all_projects("awards")
all_pubs <- load_all_projects("pubs")

View(all_awards$HTAN)
