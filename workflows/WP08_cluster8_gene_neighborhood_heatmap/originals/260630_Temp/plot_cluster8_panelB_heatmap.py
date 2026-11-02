#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch


MODULES = {
    "Hydroxymethyl-base\nprecursor": [
        "C5_MT",
        "TET_5mYOX",
        "dCMP_hydroxymethylase",
        "hmC_hmdU_related",
    ],
    "Pyrimidine /\nnucleotide metabolism": [
        "ThyX",
        "ThyS",
        "dCMP_dTMP_enzyme",
        "nucleotidyltransferase",
    ],
    "DNA metabolism": [
        "DNA_polymerase",
        "helicase",
        "ligase",
        "nuclease",
        "primase",
        "recombinase",
    ],
    "Glycan pathway": [
        "cluster8_GT_like",
        "TET_associated_GT_like",
        "additional_GT",
        "NDP_sugar_enzyme",
        "epimerase",
        "dehydratase",
        "reductase",
        "aminotransferase",
        "methyltransferase",
        "acetyltransferase",
        "sulfotransferase",
        "APS_PAPS_kinase",
    ],
    "Virion structural": [
        "capsid",
        "portal",
        "tail",
        "tape_measure",
        "terminase",
    ],
    "Unknown /\ndark ORFs": [
        "conserved_hypothetical",
        "DUF",
    ],
}


GROUP_ORDER = {
    "C5-MT + 5mYOX-associated": 1,
    "dCMP hydroxymethylase-associated": 2,
    "ThyX/ThyS / nucleotide-metabolism-associated": 3,
    "additional-GT/glycan-associated": 4,
    "ambiguous/no DNA-context": 5,
}


GROUP_LABELS = {
    "C5-MT + 5mYOX-associated": "C5-MT + 5mYOX",
    "dCMP hydroxymethylase-associated": "dCMP-HMase",
    "ThyX/ThyS / nucleotide-metabolism-associated": "ThyX/ThyS",
    "additional-GT/glycan-associated": "Glycan-associated",
    "ambiguous/no DNA-context": "Ambiguous",
}


GROUP_COLORS = {
    "C5-MT + 5mYOX-associated": "#5B2A86",
    "dCMP hydroxymethylase-associated": "#1F77B4",
    "ThyX/ThyS / nucleotide-metabolism-associated": "#2CA02C",
    "additional-GT/glycan-associated": "#FF7F0E",
    "ambiguous/no DNA-context": "#8C8C8C",
}


def short_col_label(col):
    label_map = {
        "C5_MT": "C5-MT",
        "TET_5mYOX": "5mYOX/\nTET",
        "dCMP_hydroxymethylase": "dCMP\nHMase",
        "hmC_hmdU_related": "5hmC/\n5hmdU",
        "ThyX": "ThyX",
        "ThyS": "ThyS",
        "dCMP_dTMP_enzyme": "dCMP/\ndTMP",
        "nucleotidyltransferase": "NTase/\nPPase",
        "DNA_polymerase": "Pol",
        "helicase": "Hel",
        "ligase": "Lig",
        "nuclease": "Nuc",
        "primase": "Pri",
        "recombinase": "Rec",
        "cluster8_GT_like": "Cluster 8\nGT-like",
        "TET_associated_GT_like": "TET-assoc.\nGT-like",
        "additional_GT": "Other\nGT",
        "NDP_sugar_enzyme": "NDP-sugar",
        "epimerase": "Epi",
        "dehydratase": "Dehyd",
        "reductase": "Red",
        "aminotransferase": "AminoT",
        "methyltransferase": "MT",
        "acetyltransferase": "AcT",
        "sulfotransferase": "ST",
        "APS_PAPS_kinase": "APS/PAPS",
        "capsid": "Capsid",
        "portal": "Portal",
        "tail": "Tail",
        "tape_measure": "Tape\nmeasure",
        "terminase": "Terminase",
        "conserved_hypothetical": "Cons.\nhyp.",
        "DUF": "DUF",
    }
    return label_map.get(col, col)


def main():
    parser = argparse.ArgumentParser(
        description="Plot Panel B global neighborhood heatmap from cluster8_neighborhood_matrix.tsv"
    )
    parser.add_argument(
        "--matrix",
        required=True,
        help="cluster8_neighborhood_matrix.tsv"
    )
    parser.add_argument(
        "--outdir",
        default="panelB_visualization",
        help="Output directory"
    )
    parser.add_argument(
        "--label-col",
        default="candidate_id",
        help="Column to use as row label"
    )
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(matrix_path, sep="\t")

    # Add missing columns as zeros
    heatmap_cols = []
    module_boundaries = []
    module_centers = []
    start = 0

    for module_name, cols in MODULES.items():
        present_cols = []
        for col in cols:
            if col not in df.columns:
                df[col] = 0
            present_cols.append(col)

        heatmap_cols.extend(present_cols)
        end = start + len(present_cols)
        module_boundaries.append(end)
        module_centers.append((start + end - 1) / 2)
        start = end

    # Sort rows biologically
    if "primary_group" not in df.columns:
        df["primary_group"] = "ambiguous/no DNA-context"

    df["group_rank_for_plot"] = df["primary_group"].map(GROUP_ORDER).fillna(99)

    for score_col in ["DNA_context_score", "glycan_score", "virion_score"]:
        if score_col not in df.columns:
            df[score_col] = 0

    df = df.sort_values(
        by=[
            "group_rank_for_plot",
            "DNA_context_score",
            "glycan_score",
            "virion_score",
            args.label_col,
        ],
        ascending=[True, False, False, False, True],
    ).reset_index(drop=True)

    data = df[heatmap_cols].fillna(0).astype(float).clip(0, 3).values

    row_labels = df[args.label_col].astype(str).tolist()

    # Figure size scales with row number
    n_rows, n_cols = data.shape
    fig_width = max(16, 0.34 * n_cols)
    fig_height = max(6, 0.28 * n_rows + 2.0)

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Layout:
    # left group bar, main heatmap, right model text
    gs = fig.add_gridspec(
        nrows=1,
        ncols=3,
        width_ratios=[0.35, 12, 3.2],
        wspace=0.03,
    )

    ax_group = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_model = fig.add_subplot(gs[0, 2])

    # Discrete colormap for scores 0,1,2,3
    cmap = ListedColormap(["#F2F2F2", "#C6DBEF", "#6BAED6", "#08519C"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)

    # X/Y labels
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([short_col_label(c) for c in heatmap_cols], rotation=90, fontsize=8)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=8)

    # Minor grid
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Module separators
    for boundary in module_boundaries[:-1]:
        ax.axvline(boundary - 0.5, color="black", linewidth=1.2)

    # Row group separators
    groups = df["primary_group"].tolist()
    for i in range(1, n_rows):
        if groups[i] != groups[i - 1]:
            ax.axhline(i - 0.5, color="black", linewidth=1.2)
            ax_group.axhline(i - 0.5, color="black", linewidth=1.2)
            ax_model.axhline(i - 0.5, color="black", linewidth=1.2)

    # Module labels on top
    module_names = list(MODULES.keys())
    for center, module_name in zip(module_centers, module_names):
        ax.text(
            center,
            -2.0,
            module_name,
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            clip_on=False,
        )

    # Left group annotation bar
    group_rgb = []
    for g in df["primary_group"]:
        group_rgb.append(GROUP_COLORS.get(g, "#8C8C8C"))

    group_index = np.arange(n_rows).reshape(-1, 1)
    group_cmap = ListedColormap(group_rgb)

    ax_group.imshow(group_index, aspect="auto", cmap=group_cmap)
    ax_group.set_xticks([])
    ax_group.set_yticks(np.arange(n_rows))
    ax_group.set_yticklabels([])
    ax_group.set_title("Class", fontsize=8)

    # Right inferred model text
    ax_model.set_xlim(0, 1)
    ax_model.set_ylim(n_rows - 0.5, -0.5)
    ax_model.set_xticks([])
    ax_model.set_yticks([])

    if "inferred_model" in df.columns:
        for y, model in enumerate(df["inferred_model"].astype(str)):
            ax_model.text(
                0.02,
                y,
                model,
                va="center",
                ha="left",
                fontsize=7,
            )

    ax_model.set_title("Inferred model", fontsize=8)

    # Titles
    ax.set_title(
        "Panel B. Global gene-neighborhood context of cluster 8 members",
        fontsize=13,
        pad=45,
    )

    # Colorbar
    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.018,
        pad=0.01,
        ticks=[0, 1, 2, 3],
    )
    cbar.ax.set_yticklabels(["0\nAbsent", "1\nWeak", "2\nMedium", "3\nHigh"])
    cbar.set_label("Annotation evidence", fontsize=8)

    # Legend for group bar
    legend_handles = [
        Patch(facecolor=color, edgecolor="none", label=GROUP_LABELS.get(group, group))
        for group, color in GROUP_COLORS.items()
    ]

    ax.legend(
        handles=legend_handles,
        title="Neighborhood class",
        bbox_to_anchor=(0.0, -0.32),
        loc="upper left",
        fontsize=8,
        title_fontsize=8,
        frameon=False,
        ncol=2,
    )

    # Clean spines
    for axis in [ax, ax_group, ax_model]:
        for spine in axis.spines.values():
            spine.set_visible(False)

    plt.tight_layout()

    pdf_out = outdir / "PanelB_cluster8_global_neighborhood_heatmap.pdf"
    png_out = outdir / "PanelB_cluster8_global_neighborhood_heatmap.png"

    fig.savefig(pdf_out, dpi=300, bbox_inches="tight")
    fig.savefig(png_out, dpi=300, bbox_inches="tight")

    print(f"[OK] Wrote {pdf_out}")
    print(f"[OK] Wrote {png_out}")


if __name__ == "__main__":
    main()
