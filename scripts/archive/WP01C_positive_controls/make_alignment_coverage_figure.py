#!/usr/bin/env python3
"""Render a primary horizontal alignment-coverage schematic."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


CONTROL_ORDER = [
    ("2C0N", "A197", "2C0N", "PDB"),
    ("7SP6_A", "CvHAS", "7SP6_A", "PDB"),
    ("Che8_gp110_model", "gp110", "Che8_gp110", "AF3"),
    ("Q84432_model", "A111/114R", "Q84432", "AF3"),
    ("pbcv1_a064r_1_model", "A064R", "pbcv1_a064r", "AF3"),
    ("Q5UQC3_apo_model", "L230 (apo)", "Q5UQC3_apo", "AF3"),
]

REP_COLOR = "#4F86B5"
BEST_COLOR = "#E5A900"
TRACK_COLOR = "#D8DEE5"
TEXT = "#17212B"
MUTED = "#5F6C7B"
GRID = "#E5E9ED"


def short_id(value: str, max_len: int = 20) -> str:
    return value if len(value) <= max_len else value[: max_len - 1] + "…"


def configure_style():
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_formats(fig, stem: Path):
    fig.canvas.draw()
    for suffix, kwargs in [("pdf", {}), ("svg", {}), ("png", {"dpi": 600})]:
        tmp = stem.parent / f"{stem.name}.tmp.{suffix}"
        final = stem.with_suffix(f".{suffix}")
        fig.savefig(tmp, bbox_inches="tight", facecolor="white", **kwargs)
        tmp.replace(final)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()

    data = pd.read_csv(args.input, sep="\t")
    required = {
        "comparison_type", "status", "query", "target", "alnlen", "qstart",
        "qend", "qlen", "qcov", "tcov", "qtmscore", "ttmscore",
    }
    missing_columns = required - set(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")
    if (data["status"] != "found").any():
        missing = data.loc[data["status"] != "found", ["query", "comparison_type"]]
        raise ValueError(f"Incomplete comparison table:\n{missing.to_string(index=False)}")

    comparison_order = [
        "representative",
        "best_nonrepresentative_giant_member",
    ]
    indexed = data.set_index(["query", "comparison_type"], drop=False)
    for control_id, _, _, _ in CONTROL_ORDER:
        for comparison in comparison_order:
            if (control_id, comparison) not in indexed.index:
                raise ValueError(f"Missing {comparison} row for {control_id}")

    fig, ax = plt.subplots(figsize=(11.2, 8.2))
    ax.set_xlim(-0.76, 1.70)
    ax.set_ylim(-0.6, len(CONTROL_ORDER) * 2 + 2.25)
    ax.axis("off")

    title_y = len(CONTROL_ORDER) * 2 + 1.95
    legend_y = len(CONTROL_ORDER) * 2 + 1.43
    header_y = len(CONTROL_ORDER) * 2 + 0.82
    n_c_y = len(CONTROL_ORDER) * 2 + 0.39
    track_end = 0.82
    ax.text(-0.76, title_y, "Direct alignment extent and structural support", fontsize=12, fontweight="bold", color=TEXT)
    ax.plot([-0.10, -0.02], [legend_y, legend_y], color=REP_COLOR, lw=6, solid_capstyle="butt")
    ax.text(0.00, legend_y, "versus A0A3Q8Q3U2", fontsize=7, color=TEXT, va="center")
    ax.plot([0.45, 0.53], [legend_y, legend_y], color=BEST_COLOR, lw=6, solid_capstyle="butt")
    ax.text(0.55, legend_y, "versus best non-representative baseline giant member", fontsize=7, color=TEXT, va="center")
    ax.text(-0.76, header_y, "Control", fontsize=7, fontweight="bold", color=MUTED)
    ax.text(-0.21, header_y, "Direct comparison", fontsize=7, fontweight="bold", color=MUTED, ha="center")
    ax.text(0.00, n_c_y, "N", fontsize=6.5, color=MUTED, ha="center")
    ax.text(track_end, n_c_y, "C", fontsize=6.5, color=MUTED, ha="center")
    ax.text(0.91, header_y, "Target", fontsize=7, fontweight="bold", color=MUTED)
    headers = [(1.28, "aligned/qlen"), (1.41, "qcov"), (1.52, "tcov"), (1.64, "minTM")]
    for x, label in headers:
        ax.text(x, header_y, label, fontsize=6.5, fontweight="bold", color=MUTED, ha="center")
    ax.plot([-0.76, 1.69], [len(CONTROL_ORDER) * 2 + 0.62] * 2, color=GRID, lw=0.8)

    start_y = len(CONTROL_ORDER) * 2 - 0.05
    for group_index, (control_id, protein, display_id, source) in enumerate(CONTROL_ORDER):
        top_y = start_y - group_index * 2
        ax.text(-0.76, top_y - 0.35, protein, fontsize=8.2, fontweight="bold", color=TEXT, va="center")
        ax.text(-0.76, top_y - 0.71, f"{display_id} · {source}", fontsize=6.4, color=MUTED, va="center")

        for offset, comparison in enumerate(comparison_order):
            y = top_y - offset * 0.78
            row = indexed.loc[(control_id, comparison)]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            color = REP_COLOR if comparison == "representative" else BEST_COLOR
            comparison_label = "Representative" if comparison == "representative" else "Best giant member"
            ax.text(-0.05, y, comparison_label, fontsize=6.5, color=MUTED, ha="right", va="center")
            ax.plot([0, track_end], [y, y], color=TRACK_COLOR, lw=6.0, solid_capstyle="butt", zorder=1)
            qlen = float(row["qlen"])
            qstart = float(row["qstart"])
            qend = float(row["qend"])
            start = max(0.0, (qstart - 1.0) / qlen)
            end = min(1.0, qend / qlen)
            start_x = start * track_end
            end_x = end * track_end
            ax.plot([start_x, end_x], [y, y], color=color, lw=7.2, solid_capstyle="butt", zorder=2)
            ax.plot([start_x, start_x], [y - 0.10, y + 0.10], color=color, lw=0.8)
            ax.plot([end_x, end_x], [y - 0.10, y + 0.10], color=color, lw=0.8)
            ax.text(0.91, y, short_id(str(row["target"]), max_len=17), fontsize=6.0, color=TEXT, va="center")
            ax.text(1.28, y, f"{int(row['alnlen'])}/{int(qlen)}", fontsize=6.2, color=TEXT, ha="center", va="center")
            ax.text(1.41, y, f"{float(row['qcov']):.2f}", fontsize=6.2, color=TEXT, ha="center", va="center")
            ax.text(1.52, y, f"{float(row['tcov']):.2f}", fontsize=6.2, color=TEXT, ha="center", va="center")
            min_tm = min(float(row["qtmscore"]), float(row["ttmscore"]))
            ax.text(1.64, y, f"{min_tm:.2f}", fontsize=6.2, color=TEXT, ha="center", va="center")

        if group_index < len(CONTROL_ORDER) - 1:
            ax.plot([-0.76, 1.69], [top_y - 1.35] * 2, color=GRID, lw=0.55)

    ax.text(
        -0.76,
        -0.38,
        "Bars show query-normalized positions of the direct Foldseek alignment. Best-member hits are ranked by min(qTM, tTM).",
        fontsize=6.6,
        color=MUTED,
    )

    stem = args.output_dir / "direct_alignment_coverage"
    save_formats(fig, stem)
    plt.close(fig)

    legend = (
        "Direct alignment extent and structural support for the corrected viral-GT positive controls. "
        "For each control, the upper track shows its direct Foldseek alignment to the A0A3Q8Q3U2 "
        "cluster representative and the lower track shows the highest-scoring non-representative member "
        "of the original baseline giant cluster, ranked by the smaller of the query- and target-normalized "
        "TM-scores. Grey tracks denote full control length and colored segments denote the aligned query "
        "interval. Values report aligned residues relative to query length, query coverage (qcov), target "
        "coverage (tcov), and min(qTM, tTM). E-values, bit scores, sequence identity, LDDT, probability, "
        "alignment TM-score and both directional TM-scores are provided in Source Data."
    )
    (args.output_dir / "direct_alignment_coverage_legend.txt").write_text(legend + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
