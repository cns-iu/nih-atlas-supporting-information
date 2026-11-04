"""Count NIH opportunity numbers, project serial numbers, and award records.

The opportunity-number count follows the project metadata definition of
primary funding (Core == 1). Project serial numbers and award records are
counted from the cleaned combined NIH RePORTER award table.
"""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
AWARDS = ROOT / "data/processed/NIHReporter-map_consortia-awards-clean.csv"
SEARCH_VARIABLES = ROOT / "data/metadata/NIHReporter_Search_Variables.xlsx"
OUTPUT = ROOT / "outputs/tables/table_2_program_award_counts.csv"
PROGRAMS = ["HTAN", "HuBMAP", "SenNet"]


def main() -> None:
    awards = pd.read_csv(AWARDS, low_memory=False)
    awards["project_id"] = awards["project_id"].astype(str).str.strip()

    # Count only primary-funding opportunity numbers, matching the table
    # definition in NIHReporter_Search_Variables.xlsx (Core == 1).
    primary = pd.read_excel(SEARCH_VARIABLES, sheet_name="NIH_Opportunity_Numbers")
    primary = primary.loc[primary["Core"].eq(1)].copy()
    primary["Project"] = primary["Project"].astype(str).str.strip()
    primary["Opportunity Number"] = primary["Opportunity Number"].astype(str).str.strip()
    primary_counts = (
        primary.groupby("Project")["Opportunity Number"]
        .nunique()
        .reindex(PROGRAMS, fill_value=0)
    )

    counts_by_program = {
        "Opportunity Numbers": [int(primary_counts[p]) for p in PROGRAMS],
        "Project Serial Numbers": [
            int(awards.loc[awards["project_id"].eq(p), "project_serial_number"].nunique())
            for p in PROGRAMS
        ],
        "Award Records": [int(awards["project_id"].eq(p).sum()) for p in PROGRAMS],
    }
    # Match the screenshot layout: metrics in rows and NIH programs in columns.
    counts = pd.DataFrame(counts_by_program, index=PROGRAMS).T.reset_index()
    counts = counts.rename(columns={"index": "NIH Program"})

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(OUTPUT, index=False)
    print(counts.to_string(index=False))
    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
