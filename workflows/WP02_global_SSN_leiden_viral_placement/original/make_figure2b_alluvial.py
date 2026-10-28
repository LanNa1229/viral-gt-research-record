#!/usr/bin/env python3
"""Create Figure 2b: CAZy annotation -> structural placement -> evidence outcome."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle
import pandas as pd


TOTAL_EXPECTED = 687

ANNOTATION_ORDER = ["CAZy assigned", "CAZy unassigned"]

PLACEMENT_ORDER = [
    "Established-family context",
    "Family-peripheral context",
    "Mixed-family community",
    "Intercommunity connector",
    "Sparsely annotated community",
    "Annotation-depleted community",
    "Virus-enriched dark community",
    "Singleton or isolate",
]

OUTCOME_ORDER = [
    "Concordant known-family placement",
    "Discordant assigned-family context",
    "Directly supported unassigned protein",
    "Limited, ambiguous, or alternative support",
    "No direct annotated-neighbour support",
]

PLACEMENT_LABELS = {
    "established_family_context": "Established-family context",
    "peripheral_established_family_context": "Family-peripheral context",
    "genuine_mixed_family_community": "Mixed-family community",
    "intercommunity_bridge": "Intercommunity connector",
    "sparsely_annotated_community": "Sparsely annotated community",
    "annotation_depleted_community": "Annotation-depleted community",
    "virus_enriched_dark_community": "Virus-enriched dark community",
    "singleton_or_isolate": "Singleton or isolate",
}

COLORS = {
    "CAZy assigned": "#365F91",
    "CAZy unassigned": "#A7ADB4",
    "Established-family context": "#365F91",
    "Family-peripheral context": "#87A9CF",
    "Mixed-family community": "#C99A24",
    "Intercommunity connector": "#D87532",
    "Sparsely annotated community": "#879B50",
    "Annotation-depleted community": "#A7ADB4",
    "Virus-enriched dark community": "#B04B78",
    "Singleton or isolate": "#D2D5D8",
    "Concordant known-family placement": "#365F91",
    "Discordant assigned-family context": "#D87532",
    "Directly supported unassigned protein": "#879B50",
    "Limited, ambiguous, or alternative support": "#C99A24",
    "No direct annotated-neighbour support": "#BFC4C8",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("Figure2B_alluvial_v3"))
    return parser.parse_args()


def derive_outcome(row: pd.Series) -> str:
    if row["cazy_status"] == "CAZy_assigned":
        relation = row["cazy_relation_to_community"]
        if relation == "concordant":
            return "Concordant known-family placement"
        if relation == "discordant":
            return "Discordant assigned-family context"
        raise ValueError(f"Unexpected assigned-family relation: {relation!r}")

    call = row["unassigned_family_support_call"]
    if call == "direct_support_concordant_with_community":
        return "Directly supported unassigned protein"
    if call in {
        "direct_support_discordant_with_community",
        "ambiguous_multi_family_direct_support",
        "insufficient_independent_annotated_neighbors",
    }:
        return "Limited, ambiguous, or alternative support"
    if call == "no_direct_annotated_neighbor_support":
        return "No direct annotated-neighbour support"
    raise ValueError(f"Unexpected unassigned support call: {call!r}")


def prepare_data(path: Path) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t")
    required = {
        "name",
        "cazy_status",
        "placement_major",
        "cazy_relation_to_community",
        "unassigned_family_support_call",
    }
    missing = required.difference(x.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if len(x) != TOTAL_EXPECTED:
        raise ValueError(f"Expected {TOTAL_EXPECTED} proteins, found {len(x)}")
    if x["name"].duplicated().any():
        duplicates = x.loc[x["name"].duplicated(), "name"].tolist()[:10]
        raise ValueError(f"Duplicate protein names detected: {duplicates}")

    x = x.copy()
    x["annotation"] = x["cazy_status"].map(
        {"CAZy_assigned": "CAZy assigned", "CAZy_unassigned": "CAZy unassigned"}
    )
    x["placement"] = x["placement_major"].map(PLACEMENT_LABELS)
    x["outcome"] = x.apply(derive_outcome, axis=1)

    for column in ["annotation", "placement", "outcome"]:
        if x[column].isna().any():
            raise ValueError(f"Unmapped values in {column}")
    return x


def node_intervals(
    counts: pd.Series,
    order: list[str],
    *,
    height: float = 1.0,
    gap: float = 0.018,
) -> dict[str, tuple[float, float]]:
    usable = height - gap * (len(order) - 1)
    if usable <= 0:
        raise ValueError("Node gaps exceed the available height")
    total = counts.reindex(order, fill_value=0).sum()
    intervals: dict[str, tuple[float, float]] = {}
    cursor = height
    for label in order:
        node_height = usable * counts.get(label, 0) / total
        intervals[label] = (cursor - node_height, cursor)
        cursor -= node_height + gap
    return intervals


def flow_offsets(
    flow: pd.DataFrame,
    source_col: str,
    target_col: str,
    source_order: list[str],
    target_order: list[str],
    source_intervals: dict[str, tuple[float, float]],
    target_intervals: dict[str, tuple[float, float]],
) -> list[dict[str, float | str | int]]:
    lookup = flow.set_index([source_col, target_col])["n"].to_dict()
    source_totals = flow.groupby(source_col)["n"].sum().to_dict()
    target_totals = flow.groupby(target_col)["n"].sum().to_dict()
    source_cursor = {label: source_intervals[label][1] for label in source_order}
    target_cursor = {label: target_intervals[label][1] for label in target_order}
    records = []

    for source in source_order:
        source_height = source_intervals[source][1] - source_intervals[source][0]
        for target in target_order:
            n = int(lookup.get((source, target), 0))
            if n == 0:
                continue
            target_height = target_intervals[target][1] - target_intervals[target][0]
            sh = source_height * n / source_totals[source]
            th = target_height * n / target_totals[target]
            sy1 = source_cursor[source]
            sy0 = sy1 - sh
            ty1 = target_cursor[target]
            ty0 = ty1 - th
            source_cursor[source] = sy0
            target_cursor[target] = ty0
            records.append(
                {
                    "source": source,
                    "target": target,
                    "n": n,
                    "sy0": sy0,
                    "sy1": sy1,
                    "ty0": ty0,
                    "ty1": ty1,
                }
            )
    return records


def ribbon_path(
    x0: float,
    x1: float,
    sy0: float,
    sy1: float,
    ty0: float,
    ty1: float,
) -> MplPath:
    dx = (x1 - x0) * 0.46
    verts = [
        (x0, sy1),
        (x0 + dx, sy1),
        (x1 - dx, ty1),
        (x1, ty1),
        (x1, ty0),
        (x1 - dx, ty0),
        (x0 + dx, sy0),
        (x0, sy0),
        (x0, sy1),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def draw_ribbons(
    ax: plt.Axes,
    flows: list[dict[str, float | str | int]],
    x0: float,
    x1: float,
    *,
    color_by: str,
) -> None:
    for record in sorted(flows, key=lambda r: int(r["n"]), reverse=True):
        color_label = str(record[color_by])
        patch = PathPatch(
            ribbon_path(
                x0,
                x1,
                float(record["sy0"]),
                float(record["sy1"]),
                float(record["ty0"]),
                float(record["ty1"]),
            ),
            facecolor=COLORS[color_label],
            edgecolor="none",
            alpha=0.30,
            zorder=1,
        )
        ax.add_patch(patch)


def display_label(label: str, n: int, total: int) -> str:
    return f"{label}\n{n} ({100*n/total:.1f}%)"


def draw_nodes(
    ax: plt.Axes,
    x: float,
    intervals: dict[str, tuple[float, float]],
    counts: pd.Series,
    order: list[str],
    *,
    width: float,
    label_side: str,
    label_y: dict[str, float] | None = None,
    fontsize: float = 7.6,
) -> None:
    label_y = label_y or {}
    for label in order:
        y0, y1 = intervals[label]
        ax.add_patch(
            Rectangle(
                (x - width / 2, y0),
                width,
                y1 - y0,
                facecolor=COLORS[label],
                edgecolor="#34383C",
                linewidth=0.55,
                zorder=3,
            )
        )
        n = int(counts.get(label, 0))
        if label_side == "left":
            tx, ha = x - width / 2 - 0.065, "right"
        else:
            tx, ha = x + width / 2 + 0.065, "left"
        node_mid = (y0 + y1) / 2
        text_y = label_y.get(label, node_mid)
        if abs(text_y - node_mid) > 0.008:
            line_start = x - width / 2 if label_side == "left" else x + width / 2
            line_end = tx + 0.025 if label_side == "left" else tx - 0.025
            ax.plot(
                [line_start, line_end],
                [node_mid, text_y],
                color="#5A6066",
                linewidth=0.45,
                zorder=3.5,
            )
        ax.text(
            tx,
            text_y,
            display_label(label, n, TOTAL_EXPECTED),
            ha=ha,
            va="center",
            fontsize=fontsize,
            color="#24282C",
            linespacing=1.25,
            zorder=4,
        )


def make_figure(x: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    annotation_counts = x["annotation"].value_counts()
    placement_counts = x["placement"].value_counts()
    outcome_counts = x["outcome"].value_counts()

    annotation_intervals = node_intervals(annotation_counts, ANNOTATION_ORDER, gap=0.035)
    placement_intervals = node_intervals(placement_counts, PLACEMENT_ORDER, gap=0.016)
    outcome_intervals = node_intervals(outcome_counts, OUTCOME_ORDER, gap=0.025)

    flow_ap = (
        x.groupby(["annotation", "placement"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    flow_po = (
        x.groupby(["placement", "outcome"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    full_flow = (
        x.groupby(["annotation", "placement", "outcome"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    full_flow["fraction_of_all_viral"] = full_flow["n"] / len(x)
    full_flow.to_csv(output_dir / "figure2b_alluvial_flow_counts.tsv", sep="\t", index=False)

    summary = pd.concat(
        [
            annotation_counts.rename_axis("category").rename("n").reset_index().assign(stage="Existing annotation"),
            placement_counts.rename_axis("category").rename("n").reset_index().assign(stage="Structural placement"),
            outcome_counts.rename_axis("category").rename("n").reset_index().assign(stage="Evidence interpretation"),
        ],
        ignore_index=True,
    )
    summary["fraction_of_all_viral"] = summary["n"] / len(x)
    summary[["stage", "category", "n", "fraction_of_all_viral"]].to_csv(
        output_dir / "figure2b_node_counts.tsv", sep="\t", index=False
    )

    ap_offsets = flow_offsets(
        flow_ap,
        "annotation",
        "placement",
        ANNOTATION_ORDER,
        PLACEMENT_ORDER,
        annotation_intervals,
        placement_intervals,
    )
    po_offsets = flow_offsets(
        flow_po,
        "placement",
        "outcome",
        PLACEMENT_ORDER,
        OUTCOME_ORDER,
        placement_intervals,
        outcome_intervals,
    )

    mpl.rcParams.update(
        {
            "font.family": "Nimbus Sans",
            "font.size": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.linewidth": 0.6,
        }
    )

    fig, ax = plt.subplots(figsize=(12.1, 7.0))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    xs = [0.0, 1.95, 4.05]
    node_width = 0.105
    draw_ribbons(
        ax,
        ap_offsets,
        xs[0] + node_width / 2,
        xs[1] - node_width / 2,
        color_by="target",
    )
    draw_ribbons(
        ax,
        po_offsets,
        xs[1] + node_width / 2,
        xs[2] - node_width / 2,
        color_by="source",
    )

    draw_nodes(
        ax,
        xs[0],
        annotation_intervals,
        annotation_counts,
        ANNOTATION_ORDER,
        width=node_width,
        label_side="left",
    )
    draw_nodes(
        ax,
        xs[1],
        placement_intervals,
        placement_counts,
        PLACEMENT_ORDER,
        width=node_width,
        label_side="left",
        label_y={
            "Sparsely annotated community": 0.287,
            "Annotation-depleted community": 0.242,
            "Virus-enriched dark community": 0.194,
            "Singleton or isolate": 0.078,
        },
        fontsize=7.2,
    )
    draw_nodes(
        ax,
        xs[2],
        outcome_intervals,
        outcome_counts,
        OUTCOME_ORDER,
        width=node_width,
        label_side="right",
    )

    headings = [
        "Existing annotation",
        "Structural placement",
        "Evidence interpretation",
    ]
    for xpos, heading in zip(xs, headings):
        ax.text(
            xpos,
            1.045,
            heading,
            ha="center",
            va="bottom",
            fontsize=9.2,
            fontweight="bold",
            color="#202428",
        )

    ax.text(
        -0.92,
        1.145,
        "b",
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        color="#111315",
    )
    ax.text(
        -0.72,
        1.145,
        "Structural placement of viral proteins",
        fontsize=12.0,
        fontweight="bold",
        ha="left",
        va="top",
        color="#111315",
    )
    ax.text(
        -0.72,
        1.095,
        "CAZy annotation → weighted-Leiden placement → direct family-evidence outcome (n = 687)",
        fontsize=8.5,
        ha="left",
        va="top",
        color="#555B61",
    )
    ax.set_xlim(-0.95, 5.0)
    ax.set_ylim(-0.01, 1.16)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.99, top=0.99, bottom=0.03)

    stem = output_dir / "Figure2B_viral_structural_placement_alluvial_v3"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data = prepare_data(args.input_tsv)
    make_figure(data, args.output_dir)
    print(f"Created Figure 2b outputs in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
