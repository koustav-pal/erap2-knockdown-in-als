# Manuscript figure workflow

Sequential Snakemake workflow for Figure 1, S2, S5C, 3B–D, S6 and S7 of
*Multi-omic profiling reveals that reduction of ERAP2 corrects protein and lipid alterations in VCP mutant ALS astrocytes*.

Source manuscript: [Google Doc](https://docs.google.com/document/d/17fpl8B5VMb7OsAXWDgXepkA6fX_5bwcGkHms8ncPET0/edit), inspected 20 October 2026; reported document modification time 2026-10-20 05:30:45 UTC.

This workflow starts from existing processed matrices, patient summaries, fitted
module assignments and curated enrichment tables. It regenerates plots and the
BODIPY/RSF models needed by those plots. It does not rerun read alignment, variant
calling, single-cell reference annotation or WGCNA discovery. **No BAM files are
inputs or archive contents.**

## Run

Run commands from repository root. With existing R packages:

```sh
Rscript workflow/scripts/check_syntax.R
snakemake --cores 1 --dry-run
snakemake --cores 1
```

For a different R executable:

```sh
snakemake --cores 1 --config rscript=/path/to/Rscript
```

Alternatively create the Snakemake environment from `workflow/envs/snakemake.yaml`
and let Snakemake install the analysis packages:

```sh
conda env create -f workflow/envs/snakemake.yaml -n erap2-workflow
conda activate erap2-workflow
snakemake --cores 1 --use-conda
```

Environment specifications are dependency declarations, not a reconstruction of
the historical package lockfile. Each stage saves `sessionInfo.txt`. Syntax and
selected real-data checks used R 4.5; fresh conda environment solves have not been
tested. S7 needs Seurat 5. Output files and logs go under `workflow/results/`.

Every stage depends on the preceding stage's outputs, so analysis stages remain
sequential even when more cores are supplied. `workflow/config/steps.json` maps
stages to panels; `workflow/provenance/panels.tsv` records corresponding original
plots and sources.

## Inputs

`workflow/config/inputs.tsv` lists the required files and their original project
paths. Data live in `workflow/data/`, which is excluded from Git. Most original
local RDS/plot files were zero-byte placeholders; real inputs were retrieved from
`koustav@patanilab01.lsi.nus.edu.sg:/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD`.
The annotated single-cell object is approximately 3.7 GB and bulk counts 185 MB.

For an archive recipient, put the listed processed input files in `workflow/data/`.
For an existing project tree, this optional helper copies only manifest files:

```sh
python workflow/scripts/materialize_inputs.py \
  --source-root /home/koustav/Projects/012-Ben-Clarke-ERAP1-KD \
  --host koustav@patanilab01.lsi.nus.edu.sg
```

Omit `--host` for local files. `--files filename ...` stages a subset. `genes.gtf`
contains only `gene` rows of the original Ensembl GRCh38.104 GTF, preserving its
gene IDs, names and biotypes. The manifest and helper never select BAMs.

STRING enrichment uses the frozen grouped table in
`workflow/resources/protein_string_grouped.csv`. The existing saved STRING table
had no `direction_group` column and could not reproduce the grouped panels.
The replacement was queried against STRING **11.5**, with the original four
protein direction groups, species 9606 and FDR < 0.05. This is a regenerated
annotation snapshot, not a recovered historical result; exact numerical equality
to the manuscript's old STRING plot has not been established. Query details are
in `workflow/provenance/string_snapshot.json`; workflow execution is offline.

## Analysis choices and small fixes

- PCA follows the existing DEP/VSN/MinProb/PCAtools code. Cached HGNC–UniProt
  mapping replaces its live Ensembl query; the unused full GTF import is removed.
- GO plots use the existing manually annotated tables; they do not silently rerun
  enrichment against changing databases. Figure 3C selects the module containing
  ERAP2, confirmed as M3, from the saved module assignments.
- S2B retains the original GO:BP 2023 gene-set parsing, minimum size 10,
  `nPermSimple=10000`, and pathways with p < 0.1 in both comparisons. The erroneous
  extra protein-LFC layer in its NES plot is removed. A fixed seed is added.
- S5C retains the original cell-level `lme4::lmer` model and emmeans contrasts.
  An unused, broken well-summary block was omitted.
- 3D/S6D–E use the SNP-adjusted source in `analysis/18`, not the older SNP-free
  source in `analysis/11`. The original 1000 seeds, 1000 trees, mutation selection,
  median risk and `>= median` high-risk threshold are retained for cervical and
  lumbar regions. Models are processed one at a time to avoid keeping 2000
  forests in RAM. This full fit is expensive and has not been rerun here.
- S7 uses the user's current notebook and processed annotated object. The workflow
  copy removes pasted console output and an unused duplicated bulk block, supplies
  the missing colour vector, and fixes `patient_scores` → `patient_scores_df` in
  the contribution calculation. Its statistical models and module scoring remain
  from the source. S7 was syntax checked only, as requested; it was not rerun.
- The original S7 notebook had a stray `> erap2_summary` console prompt inside an
  R chunk. That one line was commented out; extracted R chunks now parse.

Original source files and line ranges are recorded in
`workflow/provenance/source_ranges.json`. Plot layouts may combine multiple
manuscript panels; final Illustrator assembly, lettering and microscopy images
are outside this workflow. Figure 1B's axis label now describes the plotted
proportion correctly instead of calling it log10 frequency.

## Validation and archive status

All workflow R scripts and the original S7 notebook's extracted R chunks parse.
A full Snakemake dry run resolves all 11 stages with real input paths. Real-data
plot checks passed for protein frequency/rescue, bulk expression/module and
grouped STRING plots. GO plots also ran using existing remote R packages.
PCA, GSEA, full RSF and S7 analyses have not been rerun. The available remote
library lacks emmeans, so the BODIPY execution check stopped at package loading.

This prepares code and an input manifest for a future Zenodo archive. It does not
publish a deposition. Package environments, full numerical reproduction and
final figure matching still need verification before describing the archive as
fully reproduced. Other manuscript panels and remotely discovered autophagy
candidates are recorded in `workflow/provenance/panels.tsv`.
