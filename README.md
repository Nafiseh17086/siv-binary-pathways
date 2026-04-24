# siv-binary-pathways

**A BINARY-style analysis of SIV pathway enrichment data, anchored on the FAU gene.**

This repository adapts the core idea of Lin et al. 2024 (*Cell Genomics* 4, 100565 —
"Complete spatially resolved gene expression is not necessary for identifying spatial
domains") to a different data type: **pathway enrichment tables** from a SIV
differential-expression study.

The claim we borrow from the original paper:

> Binarized data (0/1 presence-absence) can preserve the structure that matters,
> even when quantitative values are available.

The question we ask with our data:

> If we throw away `p.adjust`, `qvalue`, `GeneRatio`, and `Count`, and keep only a
> 0/1 flag for "significantly enriched in this cluster (`p.adj < 0.05`)", do we still
> recover the same cluster-level biology?

**Headline result: yes — r = 0.985** between the cluster-distance matrices computed
on the full `−log₁₀(p.adj)` values and on the binary representation.

---

## Dataset

`AllSets_SIV_DEGs_Pathways_Enrich_read` — a clusterProfiler-style enrichment table
with 1,499 rows across **six clusters**, laid out as a 2 × 3 design:

|                       | Early ART · During ART | Early ART · Early ATI | Late ART · Early ATI |
| --------------------- | ---------------------: | --------------------: | -------------------: |
| **Up-regulated SIV**   | 356 pathway rows       | 376                   | 363                  |
| **Down-regulated SIV** | 139                    | 135                   | 127                  |

After collapsing to unique pathway descriptions, the analysis matrix is
**641 pathways × 6 clusters**. Up-regulated clusters carry ~2.7× the
significant pathways of down-regulated clusters — a real experimental
asymmetry, not an artifact.

## Method

Four steps, no deep learning, no hyperparameters beyond α = 0.05.

1. **Ingest** — parse the Numbers workbook sheet into a tidy DataFrame.
2. **Binarize** — for each (pathway × cluster) cell, set `1` if `p.adjust < 0.05` else `0`.
3. **Compare** — compute pairwise cosine distance between clusters on
   both the full `−log₁₀(p.adj)` matrix and the binary matrix. Correlate the
   two distance matrices (Mantel-like).
4. **Anchor** — trace the **FAU** ribosomal-protein gene across pathways and clusters.

The entire analytical "model" is essentially:

```python
mat = (df.assign(sig = df['p.adjust'] < 0.05)
         .pivot_table(index='Description',
                      columns='Cluster',
                      values='sig',
                      aggfunc='max',
                      fill_value=0)
         .astype(int))
```

## Results

### 1. Binarization preserves cluster geometry (r = 0.985)

Cosine distances between the six clusters, computed on `−log₁₀(p.adj)` vs on the
binary matrix, correlate at Pearson r = 0.985 across all 15 off-diagonal pairs.
The binary matrix is roughly 1/32 the size (1 bit vs a 32-bit float per cell) yet
captures ~98.5% of the cluster-level geometry.

See `results/charts/fig3_full_vs_binary.png` and `fig6_distance_comparison.png`.

### 2. Up vs Down pathway partition

- **Up-regulated only:** 439 pathways
- **Down-regulated only:** 156 pathways
- **Shared (Up & Down):** 46 pathways

### 3. FAU anchor gene — two-state signature

FAU appears in **32 enrichment rows** covering **16 unique pathways**, and is
significantly enriched in **exactly 2 of the 6 clusters**:

| Cluster                                 | FAU pathways | Min p.adj      | Median p.adj  |
| --------------------------------------- | -----------: | -------------: | ------------: |
| Up · Early ART · During ART             | 16           | 7.3 × 10⁻⁸²   | 1.7 × 10⁻²²  |
| Down · Late ART · Early ATI             | 16           | 5.7 × 10⁻¹²   | 5.6 × 10⁻⁹   |

Same pathways, opposite direction across ART phases. Biologically, the FAU-linked
ribosome / translation program is activated during controlled ART and suppressed
when treatment is interrupted after prolonged therapy — consistent with a host
translational response tied to viremia control.

See `results/charts/fig4_FAU_spotlight.png`.

## Repo layout

```
.
├── README.md                        # this file
├── requirements.txt                 # Python dependencies
├── LICENSE                          # MIT
├── .gitignore
├── data/
│   ├── README.md                    # place your input CSV here
│   └── AllSets_SIV_DEGs_Pathways_Enrich_read.csv   # not checked in
├── scripts/
│   ├── 00_export_numbers_to_csv.py  # .numbers → .csv (Apple Numbers → tabular)
│   ├── 01_analyze.py                # binarize, compare, anchor on FAU
│   └── 02_make_charts.py            # all six figures
├── results/
│   ├── csv/                         # analysis outputs (written by 01_analyze.py)
│   └── charts/                      # PNGs (written by 02_make_charts.py)
└── notebooks/
    └── exploration.ipynb            # optional interactive exploration
```

## Reproduce

```bash
# 1. Clone and enter
git clone https://github.com/<you>/siv-binary-pathways.git
cd siv-binary-pathways

# 2. Install
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 3. Put the input data in place
#    (a) If you have the Apple Numbers file:
cp /path/to/AllSets_SIV_DEGs_Pathways_Enrich_read.numbers data/
python scripts/00_export_numbers_to_csv.py \
    data/AllSets_SIV_DEGs_Pathways_Enrich_read.numbers \
    data/AllSets_SIV_DEGs_Pathways_Enrich_read.csv
#    (b) Or drop the CSV directly into data/ with the same filename.

# 4. Run the pipeline
python scripts/01_analyze.py
python scripts/02_make_charts.py

# 5. Outputs land in results/csv/ and results/charts/
```

## Anchor gene: changing FAU → any other gene

The anchor is a one-line change in `scripts/01_analyze.py`:

```python
GENE = "FAU"        # ← change here
```

Useful anchors to try on this dataset: `RPL10`, `EEF1A1`, `EIF3G`, `NFKB1`, `STAT1`.
Each will produce its own per-cluster signature table and binary matrix in `results/csv/`.

## Limitations

- **Cluster-level only.** Binarization preserves cluster geometry (r = 0.985) but the
  individual p-values do still matter for ranking within a cluster. If you need the
  most significant pathway in a single cluster, use `p.adjust`, not the binary table.
- **α is a knob.** We chose `α = 0.05` to match standard enrichment reporting.
  Sensitivity to α (tried 0.01, 0.10 informally) is low for cluster-distance
  correlation but will change pathway counts.
- **Direction sign is in the cluster label, not the matrix.** The Up/Down
  information is encoded by which cluster a pathway is on in, not by signed values.

## Reference

Lin, S., Cui, Y., Zhao, F., Yang, Z., Song, J., Yao, J., Zhao, Y., Qian, B.-Z.,
Zhao, Y., & Yuan, Z. (2024). Complete spatially resolved gene expression is not
necessary for identifying spatial domains. *Cell Genomics*, 4(6), 100565.
<https://doi.org/10.1016/j.xgen.2024.100565>

## License

MIT — see `LICENSE`.
