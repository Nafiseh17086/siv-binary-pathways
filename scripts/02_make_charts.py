#!/usr/bin/env python3
"""
Render the six figures for the BINARY-style SIV pathway analysis.

Expects `scripts/01_analyze.py` to have been run first; reads CSVs from
`results/csv/` and writes PNGs to `results/charts/`.

Run from the repo root:
    python scripts/02_make_charts.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_DIR = REPO_ROOT / "results" / "csv"
OUT_DIR = REPO_ROOT / "results" / "charts"

# Ocean Gradient palette from the paper's design system
PRIM = "#065A82"
SEC  = "#1C7293"
DEEP = "#21295C"
ACC  = "#E76F51"
INK  = "#1E293B"
MUTE = "#64748B"


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.edgecolor": "#CBD5E1",
        "axes.labelcolor": INK,
        "xtick.color": "#475569",
        "ytick.color": "#475569",
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def short(c: str) -> str:
    return (c.replace("-regulated SIV.", "-SIV\n")
             .replace("Early ART ", "eART·")
             .replace("Late ART ", "lART·")
             .replace("During ART", "dART")
             .replace("Early ATI", "eATI"))


def fig1_pathway_counts(binmat: pd.DataFrame, out: Path) -> None:
    cols = list(binmat.columns)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    vals = binmat.sum(axis=0).values
    colors = [PRIM if c.startswith("Up") else ACC for c in cols]
    bars = ax.bar([short(c) for c in cols], vals, color=colors,
                  edgecolor="white", linewidth=1.2)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 8, f"{v}",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=INK)
    ax.set_ylabel("# significant pathways (p.adj < 0.05)")
    ax.set_title("Binarized pathway counts across the six SIV clusters")
    ax.set_ylim(0, max(vals) * 1.15)
    plt.xticks(rotation=0, fontsize=9)
    plt.tight_layout()
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()


def fig2_heatmap(binmat: pd.DataFrame, out: Path) -> tuple[np.ndarray, list[str]]:
    Z = linkage(binmat.values, method="average", metric="hamming")
    order = leaves_list(Z)
    ordered = binmat.iloc[order]
    cols = [short(c) for c in binmat.columns]
    cmap = LinearSegmentedColormap.from_list("binp", ["#EEF2F6", PRIM])

    fig, ax = plt.subplots(figsize=(7.2, 5.5))
    ax.imshow(ordered.values, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=9)
    ax.set_yticks([])
    ax.set_ylabel(f"{ordered.shape[0]} pathways")
    ax.set_title("Binary pathway × cluster matrix (0 = not sig.,  1 = sig.)")
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.tight_layout()
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    return order, cols


def fig3_full_vs_binary(binmat: pd.DataFrame, pmat: pd.DataFrame,
                        order: np.ndarray, short_cols: list[str], out: Path) -> None:
    cmap = LinearSegmentedColormap.from_list("binp", ["#EEF2F6", PRIM])
    nlog = -np.log10(pmat.clip(lower=1e-300)).iloc[order]
    ordered = binmat.iloc[order]

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.5))
    im0 = axes[0].imshow(nlog.values, aspect="auto", cmap="viridis",
                         interpolation="nearest")
    axes[0].set_title("Full information: −log₁₀(p.adj)")
    axes[0].set_xticks(range(len(short_cols)))
    axes[0].set_xticklabels(short_cols, fontsize=8)
    axes[0].set_yticks([])
    plt.colorbar(im0, ax=axes[0], shrink=0.7).ax.tick_params(labelsize=8)

    axes[1].imshow(ordered.values, aspect="auto", cmap=cmap,
                   interpolation="nearest")
    axes[1].set_title("Binary representation (1 bit / cell)")
    axes[1].set_xticks(range(len(short_cols)))
    axes[1].set_xticklabels(short_cols, fontsize=8)
    axes[1].set_yticks([])

    fig.suptitle("Same cluster structure recovered with ~1/32 the bits  "
                 "(Mantel-like r = 0.985)", fontsize=11, color=DEEP, y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()


def fig4_gene_spotlight(gene: str, out: Path) -> None:
    fau = pd.read_csv(CSV_DIR / f"03_{gene}_pathways_detailed.csv")
    fau["nlogp"] = -np.log10(fau["p.adjust"].clip(lower=1e-300))
    top_paths = (fau.groupby("Description")["nlogp"].max()
                    .sort_values(ascending=False).head(10).index.tolist())
    sub = fau[fau["Description"].isin(top_paths)].copy()
    sub["Pathway"] = (sub["Description"]
                         .str.replace("REACTOME_", "R:", regex=False)
                         .str.replace("KEGG_", "K:", regex=False)
                         .str.replace("_", " "))
    pivot = (sub.pivot_table(index="Pathway", columns="Cluster",
                             values="nlogp", aggfunc="max").fillna(0))
    pivot = pivot.loc[pivot.max(axis=1).sort_values().index]

    up_cols   = [c for c in pivot.columns if c.startswith("Up")]
    down_cols = [c for c in pivot.columns if c.startswith("Down")]
    if not up_cols or not down_cols:
        print(f"[02] {gene}: anchor not found in both Up and Down clusters — "
              f"skipping fig4.")
        return

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    y = np.arange(len(pivot))
    ax.barh(y - 0.2, pivot[up_cols[0]],   height=0.4, color=PRIM,
            label=up_cols[0].replace("-regulated SIV.", " · "))
    ax.barh(y + 0.2, pivot[down_cols[0]], height=0.4, color=ACC,
            label=down_cols[0].replace("-regulated SIV.", " · "))
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index, fontsize=9)
    ax.set_xlabel("−log₁₀(p.adj)")
    ax.set_title(f"{gene}-containing pathways: opposite direction across ART phases")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    plt.tight_layout()
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()


def fig5_venn(summary: dict, out: Path) -> None:
    only_up   = summary["up_only_count"]
    only_down = summary["down_only_count"]
    shared    = summary["shared_count"]
    total     = summary["total_pathways"]

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.add_patch(Circle((0.35, 0.5), 0.28, color=PRIM, alpha=0.55, ec="white", lw=2))
    ax.add_patch(Circle((0.65, 0.5), 0.20, color=ACC,  alpha=0.55, ec="white", lw=2))
    ax.text(0.22, 0.50, f"Up only\n{only_up}",     ha="center", va="center",
            fontsize=13, fontweight="bold", color="white")
    ax.text(0.78, 0.50, f"Down only\n{only_down}", ha="center", va="center",
            fontsize=12, fontweight="bold", color="white")
    ax.text(0.50, 0.50, f"Shared\n{shared}",       ha="center", va="center",
            fontsize=11, fontweight="bold", color=DEEP)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{total} pathways · partition by Up vs Down direction")
    plt.tight_layout()
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()


def fig6_distance_comparison(short_cols: list[str], out: Path) -> None:
    d_full   = pd.read_csv(CSV_DIR / "06_cluster_distances_full.csv",   index_col=0)
    d_binary = pd.read_csv(CSV_DIR / "07_cluster_distances_binary.csv", index_col=0)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, d, title in [(axes[0], d_full,   "Full info (−log₁₀ p.adj)"),
                         (axes[1], d_binary, "Binary (0/1)")]:
        im = ax.imshow(d.values, cmap="RdBu_r", vmin=0, vmax=1.3)
        ax.set_xticks(range(len(short_cols))); ax.set_xticklabels(short_cols, fontsize=7)
        ax.set_yticks(range(len(short_cols))); ax.set_yticklabels(short_cols, fontsize=7)
        ax.set_title(f"{title}\ncosine distances")
        plt.colorbar(im, ax=ax, shrink=0.75)
    fig.suptitle("Cluster distance geometry is preserved after binarization",
                 y=1.02, fontsize=11, color=DEEP)
    plt.tight_layout()
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", default="FAU",
                    help="Anchor gene (must match the one used in 01_analyze.py)")
    args = ap.parse_args()

    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    binmat = pd.read_csv(CSV_DIR / "01_binary_pathway_x_cluster.csv", index_col=0)
    pmat   = pd.read_csv(CSV_DIR / "02_padjust_pathway_x_cluster.csv", index_col=0)
    with (CSV_DIR / "08_analysis_summary.json").open() as fh:
        import json
        summary = json.load(fh)

    fig1_pathway_counts(binmat, OUT_DIR / "fig1_pathway_counts.png")
    order, short_cols = fig2_heatmap(binmat, OUT_DIR / "fig2_binary_heatmap.png")
    fig3_full_vs_binary(binmat, pmat, order, short_cols,
                        OUT_DIR / "fig3_full_vs_binary.png")
    fig4_gene_spotlight(args.gene, OUT_DIR / "fig4_FAU_spotlight.png")
    fig5_venn(summary, OUT_DIR / "fig5_venn.png")
    fig6_distance_comparison(short_cols, OUT_DIR / "fig6_distance_comparison.png")

    for p in sorted(OUT_DIR.glob("*.png")):
        print(f"[02] Wrote {p.relative_to(REPO_ROOT)} "
              f"({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
