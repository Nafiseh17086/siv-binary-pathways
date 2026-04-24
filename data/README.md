# Data

This folder is intentionally empty in git (see `../.gitignore`).

Place your raw input file here before running the scripts:

- `AllSets_SIV_DEGs_Pathways_Enrich_read.numbers`  **or**
- `AllSets_SIV_DEGs_Pathways_Enrich_read.csv`

If you have the `.numbers` version, run:

```bash
python ../scripts/00_export_numbers_to_csv.py \
    AllSets_SIV_DEGs_Pathways_Enrich_read.numbers \
    AllSets_SIV_DEGs_Pathways_Enrich_read.csv
```

The downstream scripts expect the CSV at
`data/AllSets_SIV_DEGs_Pathways_Enrich_read.csv`.

## Expected columns

```
Cluster, DE * Group, ID, Description, GeneRatio, BgRatio,
pvalue, p.adjust, qvalue, geneID, Count
```

(Two blank spacer columns in the original workbook are ignored by the parser.)
