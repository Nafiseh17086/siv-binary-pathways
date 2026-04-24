#!/usr/bin/env python3
"""
BINARY-style analysis of SIV pathway enrichment results.

The core idea from Lin et al. 2024 (Cell Genomics):
    "Complete quantitative information is not always necessary."

We adapt it to pathway enrichment: binarize (pathway x cluster) significance
at p.adj < 0.05, then show the binary matrix preserves the cluster geometry of
the full -log10(p.adj) matrix (Mantel-like r).

An anchor gene (default: FAU) is traced across all pathways it belongs to to
produce a simple two-state signature.

Run from the repo root:
    python scripts/01_analyze.py

Optional:
    python scripts/01_analyze.py --gene RPL10
    python scripts/01_analyze.py --alpha 0.01
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_IN = REPO_ROOT / "data" / "AllSets_SIV_DEGs_Pathways_Enrich_read.csv"
OUT_DIR = REPO_ROOT / "results" / "csv"


def load(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df[df["Cluster"].notna() & (df["Cluster"].astype(str).str.strip() != "")].copy()
    for col in ("pvalue", "p.adjust", "qvalue", "Count"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["Description"] = df["Description"].astype(str)
    df["geneID"] = df["geneID"].astype(str)

    def split_cluster(s: str):
        m = re.match(r"(Up|Down)-regulated SIV\.(.+)$", s)
        return (m.group(1), m.group(2)) if m else (None, None)

    df[["Direction", "Context"]] = df["Cluster"].apply(
        lambda s: pd.Series(split_cluster(s))
    )
    return df


def binarize(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    sig = (df["p.adjust"] < alpha).astype(int)
    return (
        df.assign(_sig=sig)
          .pivot_table(index="Description", columns="Cluster",
                       values="_sig", aggfunc="max", fill_value=0)
          .astype(int)
    )


def padjust_matrix(df: pd.DataFrame) -> pd.DataFrame:
    return df.pivot_table(index="Description", columns="Cluster",
                          values="p.adjust", aggfunc="min",
                          fill_value=1.0)


def cosine_distance(mat: pd.DataFrame) -> np.ndarray:
    """Cosine distance between the columns of a numeric DataFrame."""
    M = mat.values.T.astype(float)
    n = np.linalg.norm(M, axis=1, keepdims=True)
    n[n == 0] = 1.0
    M = M / n
    return 1.0 - (M @ M.T)


def contains_gene(s: str, gene: str) -> bool:
    return bool(re.search(rf"\b{re.escape(gene)}\b", str(s)))


def analyse(df: pd.DataFrame, gene: str, alpha: float) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- binary matrix
    binmat = binarize(df, alpha)
    binmat.to_csv(OUT_DIR / "01_binary_pathway_x_cluster.csv")

    # ---- full-info matrix
    pmat = padjust_matrix(df)
    pmat.to_csv(OUT_DIR / "02_padjust_pathway_x_cluster.csv")

    # ---- Mantel-like comparison
    nlog = -np.log10(pmat.clip(lower=1e-300))
    d_full = cosine_distance(nlog)
    d_bin = cosine_distance(binmat.astype(float))
    cols = list(binmat.columns)
    pd.DataFrame(d_full, index=cols, columns=cols).to_csv(OUT_DIR / "06_cluster_distances_full.csv")
    pd.DataFrame(d_bin,  index=cols, columns=cols).to_csv(OUT_DIR / "07_cluster_distances_binary.csv")
    iu = np.triu_indices(len(cols), k=1)
    mantel_r = float(np.corrcoef(d_full[iu], d_bin[iu])[0, 1])

    # ---- Up/Down partition
    up_cols   = [c for c in cols if c.startswith("Up")]
    down_cols = [c for c in cols if c.startswith("Down")]
    up_any   = set(binmat.index[(binmat[up_cols].sum(axis=1)   > 0)])
    down_any = set(binmat.index[(binmat[down_cols].sum(axis=1) > 0)])

    # ---- anchor gene
    df_gene = df[df["geneID"].apply(lambda s: contains_gene(s, gene))].copy()
    (df_gene
     .sort_values(["Cluster", "p.adjust"])
     [["Cluster", "Direction", "Context", "Description",
       "GeneRatio", "BgRatio", "p.adjust", "Count", "geneID"]]
     .to_csv(OUT_DIR / f"03_{gene}_pathways_detailed.csv", index=False))

    (df_gene
     .assign(_sig=lambda d: (d["p.adjust"] < alpha).astype(int))
     .pivot_table(index="Description", columns="Cluster",
                  values="_sig", aggfunc="max", fill_value=0)
     .astype(int)
     .to_csv(OUT_DIR / f"04_{gene}_binary_matrix.csv"))

    gene_summary = (df_gene
                    .groupby(["Direction", "Context"])
                    .agg(n_pathways=("Description", "nunique"),
                         min_padj=("p.adjust", "min"),
                         median_padj=("p.adjust", "median"))
                    .reset_index())
    gene_summary.to_csv(OUT_DIR / f"05_{gene}_summary_by_cluster.csv", index=False)

    # ---- summary JSON
    summary = dict(
        alpha=alpha,
        gene=gene,
        total_rows=int(len(df)),
        total_pathways=int(binmat.shape[0]),
        total_clusters=int(binmat.shape[1]),
        clusters=cols,
        density=float(binmat.values.mean()),
        sig_per_cluster={c: int(binmat[c].sum()) for c in cols},
        up_only_count=len(up_any - down_any),
        down_only_count=len(down_any - up_any),
        shared_count=len(up_any & down_any),
        mantel_r=mantel_r,
        gene_total_rows=int(len(df_gene)),
        gene_unique_pathways=int(df_gene["Description"].nunique()),
        gene_per_cluster={c: int((df_gene["Cluster"] == c).sum()) for c in cols},
        gene_min_padj=(float(df_gene["p.adjust"].min())
                       if len(df_gene) else None),
    )
    with (OUT_DIR / "08_analysis_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default=str(DEFAULT_IN),
                    help=f"Input CSV (default: {DEFAULT_IN.relative_to(REPO_ROOT)})")
    ap.add_argument("--gene", default="FAU", help="Anchor gene symbol (default: FAU)")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold on p.adjust (default: 0.05)")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"[01] Input not found: {in_path}\n"
                         f"     Place the CSV at {DEFAULT_IN} or pass --input.")

    print(f"[01] Loading {in_path}")
    df = load(in_path)
    print(f"[01] Rows after cleaning: {len(df)}")
    print(f"[01] Anchor gene: {args.gene} · alpha: {args.alpha}")

    summary = analyse(df, gene=args.gene, alpha=args.alpha)

    print(f"[01] Binary matrix: {summary['total_pathways']} pathways × "
          f"{summary['total_clusters']} clusters · density {summary['density']:.3f}")
    print(f"[01] Mantel-like r(full, binary) = {summary['mantel_r']:.3f}")
    print(f"[01] Up-only / Down-only / Shared = "
          f"{summary['up_only_count']} / {summary['down_only_count']} / {summary['shared_count']}")
    print(f"[01] Gene {args.gene!r}: {summary['gene_total_rows']} rows · "
          f"{summary['gene_unique_pathways']} unique pathways")
    print(f"[01] Wrote CSVs and JSON to {OUT_DIR.relative_to(REPO_ROOT)}/")


if __name__ == "__main__":
    main()
