# Foodborne bacteria example

These notebooks reproduce the C18:2n-6 workflow for *Listeria monocytogenes*
described in the study referenced at the top of each notebook.

## Setup

From the repository root:

```bash
python -m pip install -e ".[notebooks]"
cd "Examples/Foodborne bacteria"
jupyter lab
```

Run all cells in `prepare_data.ipynb`, then all cells in `Analysis.ipynb`, using
this directory as the working directory. The notebooks require internet access
for KEGG organism/gene lookups and for uncached database downloads.

After providing the source spreadsheet, execute and save both notebooks in order:

```bash
jupyter nbconvert --to notebook --execute --inplace prepare_data.ipynb
jupyter nbconvert --to notebook --execute --inplace Analysis.ipynb
```

## Required source data

Obtain the original study's supplementary spreadsheet and set `SOURCE_FILE` in
`prepare_data.ipynb` to its local path. Its default name is
`aem.03027-20-s0002.xlsx`. The spreadsheet is not bundled with this repository.
The selected worksheet must contain these columns:

- `Gene_id` and `gene_name`;
- `C18:2n-6_fold_change` (a ratio, not log2 fold change);
- `C18:2n-6_pValue` and `C18:2n-6_fdr`.

The notebook retains the study's thresholds: fold change >= 2.5 or <= 0.4,
with FDR <= 0.05. Check these assumptions before adapting it to another dataset.

## Generated files and caches

The preparation notebook builds `df_all.xlsx`, `df_up.xlsx`, `df_down.xlsx`,
`background.txt`, the two gene-ID lists, and the two gene-name lists. It also
writes `unmapped_genes.xlsx`; review excluded genes before interpreting results.
The historical filenames `downregulted.txt` and `downregultedGNames.txt` are
retained for compatibility.

The analysis notebook creates `exemple_fb.gmt` if absent, downloads `aro.obo`
if absent, and caches VFDB files in `vfdb_data/`. Reuse these files to avoid
repeated downloads. Record their download dates/versions for reproducible work.
Enrichment reports are saved under `enrichment_exampleup/` and
`enrichment_exampledown/`. The two expression groups use separate color maps.

Annotations rely on names rather than sequence or phenotype validation. Saved
outputs were cleared after the code corrections; rerun with the original data
to obtain current results. Automated regression tests exercise local fixtures
and notebook transformations without claiming to reproduce the study's findings.
