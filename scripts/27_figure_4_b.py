#!/usr/bin/env python3
"""Render a component-level alternative to Figure 1 Panel A.

Each Panel-A point is one consortium component.  Its ternary coordinates are
the shares of that component's associated within-atlas internal-collaboration
publications classified as A (within component), B (cross component within a
functional type), or C (cross functional type).  A publication associated with
multiple components contributes to every relevant component profile.
"""
from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "figure1_candidates"
FIGURE_OUT = ROOT / "outputs" / "figures"
RES = ROOT / "results" / "figure1_candidates"
OUT.mkdir(parents=True, exist_ok=True)
FIGURE_OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

PANEL_B = runpy.run_path(str(ROOT / "scripts" / "96_figure1_panelC_publication_compositions.py"))
PROGRAMS = ["HTAN", "HuBMAP", "SenNet"]
PROGRAM_COLORS = PANEL_B["PROGRAM_COLORS"]
TRIANGLE_HEIGHT = np.sqrt(3) / 2
TYPE_LABEL = {"technology": "Technology", "experimental": "Experimental", "coordination": "Coordination"}
TYPE_MARKERS = {"technology": "o", "experimental": "s", "coordination": "^"}
PATTERN_MAP = {"within_cc": "A", "cross_cc_same_type": "B", "cross_type": "C"}


def build_component_profiles():
    """Calculate one A/B/C composition for every consortium component."""
    classified = pd.read_csv(ROOT / "data" / "results" / "publication_functional_classification.csv", dtype={"pmid": str})
    metadata = pd.read_csv(
        ROOT / "data" / "metadata" / "NIH-map_consortia-publication-metadata.csv",
        dtype={"pmid": str, "core_project_number": str},
    )
    component_types = pd.read_csv(ROOT / "results" / "component_type_mapping.csv")
    z = metadata.merge(
        classified[["consortia", "pmid", "functional_collaboration_class"]],
        on=["consortia", "pmid"], how="inner",
    )
    z = z[z.functional_collaboration_class.isin(PATTERN_MAP)].copy()
    z["pattern"] = z.functional_collaboration_class.map(PATTERN_MAP)
    z = z.merge(component_types[["consortia", "division", "component_type"]], on=["consortia", "division"], how="left")
    assert z.component_type.notna().all()
    # One program-component-PMID association is the unit for this panel.
    z = z.drop_duplicates(["consortia", "division", "pmid"])
    counts = z.groupby(["consortia", "division", "component_type", "pattern"]).pmid.nunique().unstack(fill_value=0)
    counts = counts.reindex(columns=list("ABC"), fill_value=0).reset_index()
    counts["n_IC"] = counts[["A", "B", "C"]].sum(axis=1)
    for pattern in "ABC":
        counts[f"pct_{pattern}"] = 100 * counts[pattern] / counts.n_IC
    counts = counts.rename(columns={"consortia": "atlas_program", "division": "component"})
    counts["structural_constraint_flag"] = counts[["atlas_program", "component"]].apply(tuple, axis=1).isin({("HTAN", "DCC"), ("SenNet", "CODCC")})
    counts = counts.sort_values(["atlas_program", "component"]).reset_index(drop=True)
    assert len(counts) == 11
    assert counts.n_IC.ge(2).all()
    assert np.allclose(counts[["pct_A", "pct_B", "pct_C"]].sum(axis=1), 100)
    return counts


def ternary_xy(a, b, c):
    return b + .5 * c, TRIANGLE_HEIGHT * c


def draw_component_panel(ax, profiles):
    # DCC and CODCC are singleton coordination structures.  They remain in
    # the exported source data but are deliberately omitted from Panel A.
    profiles = profiles[~profiles.component.isin(["DCC", "CODCC"])].copy()
    triangle = np.array([[0, 0], [1, 0], [.5, TRIANGLE_HEIGHT]])
    ax.add_patch(Polygon(triangle, fill=False, lw=1.2, color="0.25"))
    for fraction in [.25, .50, .75]:
        ax.plot([.5 * fraction, 1 - .5 * fraction], [TRIANGLE_HEIGHT * fraction, TRIANGLE_HEIGHT * fraction], color="0.9", lw=.7)
        ax.plot([1 - fraction, .5 * (1 - fraction)], [0, TRIANGLE_HEIGHT * (1 - fraction)], color="0.9", lw=.7)
        ax.plot([fraction, .5 + .5 * fraction], [0, TRIANGLE_HEIGHT * (1 - fraction)], color="0.9", lw=.7)
    # Absolute label positions place every label close to its point while
    # keeping nearby labels separate; the observations are not moved.
    label_positions = {
        ("HTAN", "DCC"): (.430, .860, "right"),
        ("HTAN", "HTA"): (.980, .110, "left"),
        ("HTAN", "PCA"): (.980, .035, "left"),
        ("HuBMAP", "DEMO"): (.895, .190, "left"),
        ("HuBMAP", "HIVE"): (.190, .370, "right"),
        ("HuBMAP", "RIT"): (.490, .810, "right"),
        ("HuBMAP", "TMC_H"): (.435, .625, "left"),
        ("HuBMAP", "TTD"): (.070, .190, "right"),
        ("SenNet", "CODCC"): (.535, .835, "left"),
        ("SenNet", "TDA"): (.130, .085, "right"),
        ("SenNet", "TMC_S"): (.205, .440, "right"),
    }
    for _, row in profiles.iterrows():
        x, y = ternary_xy(row.pct_A / 100, row.pct_B / 100, row.pct_C / 100)
        color = PROGRAM_COLORS[row.atlas_program]
        face = "white" if row.structural_constraint_flag else color
        ax.scatter(x, y, s=67, marker=TYPE_MARKERS[row.component_type], facecolor=face,
                   edgecolor=color, linewidth=1.35, zorder=3)
        label_x, label_y, alignment = label_positions[(row.atlas_program, row.component)]
        marker = " †" if row.structural_constraint_flag else ""
        label = {"RIT": "RTI"}.get(row.component, row.component.replace("_", "-"))
        ax.text(label_x, label_y, f"{label}{marker}", fontsize=14.4, ha=alignment, va="center", zorder=4)
    ax.text(-.04, -.06, "Pattern A", ha="center", va="top", fontsize=16.4, fontweight="bold")
    ax.text(1.04, -.06, "Pattern B", ha="center", va="top", fontsize=16.4, fontweight="bold")
    ax.text(.5, TRIANGLE_HEIGHT + .04, "Pattern C", ha="center", va="bottom", fontsize=16.4, fontweight="bold")
    color_handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=PROGRAM_COLORS[p], markeredgecolor=PROGRAM_COLORS[p], label=p, markersize=7)
        for p in PROGRAMS
    ]
    shape_handles = [
        Line2D([], [], marker=TYPE_MARKERS[t], color="0.25", markerfacecolor="white", label=TYPE_LABEL[t], linestyle="None", markersize=7)
        for t in ["technology", "experimental", "coordination"]
    ]
    color_legend = ax.legend(handles=color_handles, title="Atlas program", title_fontsize=14.0,
                             ncol=1, fontsize=14.0, frameon=False,
                             loc="upper left", bbox_to_anchor=(1.12, .78), borderaxespad=0,
                             labelspacing=.45, handletextpad=.45)
    ax.add_artist(color_legend)
    ax.legend(handles=shape_handles, title="CC type", title_fontsize=14.0,
              ncol=1, fontsize=14.0, frameon=False,
              loc="upper left", bbox_to_anchor=(1.12, .44), borderaxespad=0,
              labelspacing=.45, handletextpad=.45)
    ax.set(xlim=(-.18, 1.14), ylim=(-.17, 1.00), aspect="equal", xticks=[], yticks=[])
    ax.set_frame_on(False)


def save(fig, stem, out_dir=OUT):
    for extension in ("png", "svg"):
        fig.savefig(out_dir / f"{stem}.{extension}", dpi=360, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax, label):
    ax.text(-.03, 1.04, label, transform=ax.transAxes, fontsize=26, fontweight="bold", va="top")


def main():
    profiles = build_component_profiles()
    profiles.to_csv(RES / "Fig1_component_ABC_profiles.csv", index=False)

    # Standalone, corrected component-composition panel for Figure 1.
    panel_a_figure = plt.figure(figsize=(7.2, 6.9))
    panel_a = panel_a_figure.add_subplot(1, 1, 1)
    draw_component_panel(panel_a, profiles)
    panel_label(panel_a, "A")
    panel_a_figure.subplots_adjust(left=.05, right=.78, bottom=.08, top=.92)
    save(panel_a_figure, "figure_4_b", out_dir=FIGURE_OUT)

    publication_data = PANEL_B["build_plot_data"]()
    fig = plt.figure(figsize=(15.0, 6.9))
    grid = fig.add_gridspec(1, 2, wspace=.52)
    left = fig.add_subplot(grid[0, 0])
    draw_component_panel(left, profiles)
    panel_label(left, "A")
    right = fig.add_subplot(grid[0, 1])
    PANEL_B["draw_triangle"](right, publication_data, "project_share", "", "", minimal=True,
                              show_corner_labels=True, legend_right=True, font_scale=2,
                              aggregate_coincident=True)
    panel_label(right, "B")
    fig.subplots_adjust(left=.04, right=.90, bottom=.08, top=.91)
    save(fig, "Fig1_component_profiles_panelB_within_atlas_composition")
    print("Component-level Figure 1 alternative complete")


if __name__ == "__main__":
    main()
