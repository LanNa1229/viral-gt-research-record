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
    "Pyrimidine /\nnucleotide\nmetabolism": [
        "ThyX",
        "ThyS",
        "dCMP_dTMP_enzyme",
        "nucleotidyltransferase",
    ],
    "DNA\nmetabolism": [
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
    "Virion\nstructural": [
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


GROUP_SHORT = {
    "C5-MT + 5mYOX-associated": "C5-MT + 5mYOX",
    "dCMP hydroxymethylase-associated": "dCMP-HMase",
    "ThyX/ThyS / nucleotide-metabolism-associated": "ThyX/ThyS",
    "additional-GT/glycan-associated": "Glycan-associated",
    "ambiguous/no DNA-context": "Ambiguous",
}


GROUP_COLORS = {
    "C5-MT + 5mYOX-associated": "#6A3D9A",
    "dCMP hydroxymethylase-associated": "#1F78B4",
    "ThyX/ThyS / nucleotide-metabolism-associated": "#33A02C",
    "additional-GT/glycan-associated": "#FF7F00",
    "ambiguous/no DNA-context": "#999999",
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
        "NDP_sugar_enzyme": "NDP-\nsugar",
        "epimerase": "Epi",
        "dehydratase": "Dehyd",
        "reductase": "Red",
        "aminotransferase": "AminoT",
        "methyltransferase": "MT",
        "acetyltransferase": "AcT",
        "sulfotransferase": "ST",
        "APS_PAPS_kinase": "APS/\nPAPS",

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--outdir", default="panelB_visualization_clean")
    parser.add_argument("--label-col", default="candidate_id")
    parser.add_argument("--no-title", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.matrix, sep="\t")

    # Prepare heatmap columns
    heatmap_cols = []
    module_spans = []
    start = 0

    for module_name, cols in MODULES.items():
        for col in cols:
            if col not in df.columns:
                df[col] = 0
        heatmap_cols.extend(cols)
        end = start + len(cols)
        module_spans.append((module_name, start, end))
        start = end

    # Sort rows
    if "primary_group" not in df.columns:
        df["primary_group"] = "ambiguous/no DNA-context"

    df["group_rank_for_plot"] = df["primary_group"].map(GROUP_ORDER).fillna(99)

    for col in ["DNA_context_score", "glycan_score", "virion_score"]:
        if col not in df.columns:
            df[col] = 0

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
    n_rows, n_cols = data.shape
    row_labels = df[args.label_col].astype(str).tolist()

    # Figure size
    fig_width = max(13, 0.36 * n_cols + 3.5)
    fig_height = max(5.5, 0.32 * n_rows + 2.8)

    fig = plt.figure(figsize=(fig_width, fig_height))

    # More top room, no right text panel
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        height_ratios=[0.55, 8],
        width_ratios=[0.35, 12, 0.45],
        hspace=0.02,
        wspace=0.03,
    )

    ax_module = fig.add_subplot(gs[0, 1])
    ax_group = fig.add_subplot(gs[1, 0])
    ax = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])

    # Main heatmap
    cmap = ListedColormap(["#F2F2F2", "#C6DBEF", "#6BAED6", "#08519C"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)

    # Axes
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([short_col_label(c) for c in heatmap_cols], rotation=90, fontsize=7)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=8)

    # Grid
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Module separators
    for module_name, start, end in module_spans[:-1]:
        ax.axvline(end - 0.5, color="black", linewidth=1.1)

    # Row group separators
    groups = df["primary_group"].tolist()
    for i in range(1, n_rows):
        if groups[i] != groups[i - 1]:
            ax.axhline(i - 0.5, color="black", linewidth=1.2)
            ax_group.axhline(i - 0.5, color="black", linewidth=1.2)

    # Module header axis
    ax_module.set_xlim(-0.5, n_cols - 0.5)
    ax_module.set_ylim(0, 1)
    ax_module.axis("off")

    for module_name, start, end in module_spans:
        center = (start + end - 1) / 2
        ax_module.text(
            center,
            0.45,
            module_name,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
        ax_module.plot([start - 0.5, end - 0.5], [0.05, 0.05], color="black", linewidth=1.0)

    # Left group color bar
    group_codes = []
    group_palette = []
    group_to_code = {}
    for g in df["primary_group"]:
        if g not in group_to_code:
            group_to_code[g] = len(group_to_code)
            group_palette.append(GROUP_COLORS.get(g, "#999999"))
        group_codes.append(group_to_code[g])

    ax_group.imshow(np.array(group_codes).reshape(-1, 1), aspect="auto",
                    cmap=ListedColormap(group_palette))
    ax_group.set_xticks([])
    ax_group.set_yticks([])
    ax_group.set_title("Class", fontsize=8, pad=4)

    # Colorbar
    cbar = fig.colorbar(im, cax=ax_cbar, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["0", "1", "2", "3"])
    cbar.set_label("Evidence\nscore", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    # Clean spines
    for axis in [ax, ax_group, ax_module]:
        for spine in axis.spines.values():
            spine.set_visible(False)

    # Legend at bottom
    handles = [
        Patch(facecolor=GROUP_COLORS[g], edgecolor="none", label=GROUP_SHORT.get(g, g))
        for g in GROUP_ORDER
        if g in set(df["primary_group"])
    ]

    fig.legend(
        handles=handles,
        title="Neighborhood class",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=min(5, len(handles)),
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )

    if not args.no_title:
        fig.suptitle(
            "Global gene-neighborhood context of cluster 8 members",
            fontsize=12,
            y=0.99,
        )

    fig.subplots_adjust(
        top=0.86,
        bottom=0.24,
        left=0.10,
        right=0.93,
    )

    pdf_out = outdir / "PanelB_cluster8_global_neighborhood_heatmap_clean.pdf"
    png_out = outdir / "PanelB_cluster8_global_neighborhood_heatmap_clean.png"

    fig.savefig(pdf_out, dpi=300, bbox_inches="tight")
    fig.savefig(png_out, dpi=300, bbox_inches="tight")

    print(f"[OK] Wrote {pdf_out}")
    print(f"[OK] Wrote {png_out}")


if __name__ == "__main__":
    main()
