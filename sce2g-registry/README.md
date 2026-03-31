# scE2G predictions registry

A searchable registry of scE2G enhancer-to-gene predictions generated in the Engreitz lab.

## Directory structure

```
sce2g-registry/
  runs.tsv              One row per scE2G run
  predictions.tsv       One row per biosample × model
  scripts/
    ingest.py           Generates completed output from partial templates + scE2G configs
    merge.py            Merges reviewed output into the registry
    validate.py         Validates OAK paths and stream URLs
    explore.py          Marimo web app for browsing predictions
  templates/
    template_partial_runs.tsv         Fill in to describe a new run
    template_partial_predictions.tsv  Fill in to describe new biosamples
    template_runs.tsv                 Full runs.tsv schema (reference)
    template_predictions.tsv          Full predictions.tsv schema (reference)
  temp/                 Staging area — ingest.py writes here; review before merging
  input/                Put your filled-in partial templates here
```

---

## How to add a new run

### Step 1 — fill in the partial templates

Copy the templates from `templates/` and fill them in.

**`template_partial_runs.tsv`** — one row per run:

| Column | Notes |
|--------|-------|
| `run_id` | e.g. `2025_0318_all_clusters` |
| `config_file` | Full path to scE2G config YAML — used to derive `cell_clusters_config` and `igv_dir` |
| `results_dir` | Full OAK path to results dir (overrides config if provided) |
| `igv_dir` | Full OAK path to IGV directory (overrides config if provided) |
| `gene_reference_file` | URL or path to gene TSS reference BED |
| `scE2G_version` | e.g. `scE2G_v1.2` |
| `date_produced` | YYYY-MM-DD |
| `contact` | |
| `genome_assembly` | e.g. `hg38` |
| `transcriptome_reference` | e.g. `GENCODE v43` |
| `gene_symbols` | e.g. `RefSeq/UCSC` |
| `notes` | |

**`template_partial_predictions.tsv`** — one row per biosample:

| Column | Example |
|--------|---------|
| `run_id` | `2025_0318_all_clusters` |
| `biosample_id` | `BMMC22_B1_B` |
| `biosample_description` | `B1 B cell` |
| `dataset` | `BMMC22` |
| `data_provenance` | `Luecken et al.` |
| `dataset_description` | `NeurIPS BMMC dataset, individual cell type resolution` |
| `data_preprocessing` | `See scE2G paper Methods` |
| `citation` | `scE2G paper` |

If the `run_id` already exists in `runs.tsv`, you can omit the runs template — run metadata will be pulled from the registry automatically.

### Step 2 — run `ingest.py`

```bash
python scripts/ingest.py \
    --predictions_template input/my_predictions.tsv \
    --runs_template input/my_runs.tsv   # omit if run already in runs.tsv
```

Output goes to `temp/` by default. Auto-fills per (biosample × model):

| Field | Source |
|-------|--------|
| `model_name`, `score_threshold` | Results directory scan |
| `n_cells`, `n_fragments`, `n_umis` | Per-model stats file |
| `qc_pass`, `is_preferred` | Computed (see QC rules below) |
| `oak_predictions_full`, `oak_predictions_thresholded` | Verified to exist on OAK |
| `stream_predictions_bedpe`, `atac_bigwig_url` | Derived from `igv_dir` using `/oak/stanford/groups/engreitz` → `https://mitra.stanford.edu/engreitz/oak`; only set if OAK file exists |
| `rna_matrix_file`, `atac_frag_file` | From `cell_clusters_config` |
| `category` | Mitra metadata TSV |

### Step 3 — review the output

Check `temp/{run_id}_predictions.tsv` and `temp/{run_id}_runs.tsv`. Verify:
- All expected biosamples and models are present
- Stream URLs look correct (or are empty where not yet uploaded)
- Any OAK path warnings printed by ingest.py are resolved
- Many any manual corrections or additions to the data (e.g., add predictions that don't follow the expected pattern)

### Step 4 — merge into the registry

```bash
# New run (both files):
python scripts/merge.py temp/{run_id}_predictions.tsv temp/{run_id}_runs.tsv

# Adding biosamples to an existing run (predictions only):
python scripts/merge.py temp/{run_id}_predictions.tsv

# Overwrite existing entries for this run_id:
python scripts/merge.py temp/{run_id}_predictions.tsv temp/{run_id}_runs.tsv --force
```

---

## QC rules

**`qc_pass`:**
- Multiome models: n_fragments ≥ 2,000,000 AND n_umis ≥ 1,000,000
- ATAC-only models: n_fragments ≥ 2,000,000

**`is_preferred`** (one model marked True per biosample):
- Default order: `multiome_powerlaw_v3` > `multiome_megamap_v3` > `scATAC_powerlaw_v3` > `scATAC_megamap_v3`
- Exception: if n_fragments ≥ 2M but n_umis < 1M, the ATAC equivalent is preferred over multiome

Always check `qc_pass` before using predictions. Also check `gene_symbols` — many predictions
use RefSeq/UCSC symbols from the ABC reference, which are not consistent with GENCODE.
Use Ensembl IDs when uncertain.

---

## Validating the registry

```bash
python scripts/validate.py               # check all OAK paths
python scripts/validate.py --urls        # also check stream URLs (slow)
python scripts/validate.py --new temp/{run_id}_predictions.tsv   # validate before merging
```

---

## Browsing predictions

```bash
marimo run scripts/explore.py
```

Or open directly in the browser (no setup required):

https://marimo.app/github.com/EngreitzLab/eutensils/blob/main/sce2g-registry/scripts/explore.py

---

## Schema

### runs.tsv

| Column | Description |
|--------|-------------|
| `run_id` | Unique ID, typically the results directory basename |
| `scE2G_version` | e.g. `scE2G_v1.2` |
| `date_produced` | YYYY-MM-DD |
| `contact` | Who ran it |
| `results_dir` | Absolute OAK path to results directory |
| `cell_clusters_config` | Path to cell clusters TSV used as input (derived from config) |
| `gene_reference_file` | URL or path to gene TSS reference BED |
| `igv_dir` | Absolute OAK path to IGV directory (for stream URL derivation) |
| `genome_assembly` | e.g. `hg38` |
| `transcriptome_reference` | e.g. `GENCODE v43` |
| `gene_symbols` | Gene symbol convention used |
| `notes` | |

### predictions.tsv

| Column | Description | Source |
|--------|-------------|--------|
| `run_id` | Links to runs.tsv | manual |
| `biosample_id` | Cluster/cell type ID as used in file paths | manual |
| `biosample_description` | Human-readable cell type name | manual |
| `dataset` | Dataset group (e.g. `BMMC22`, `Islets`) | manual |
| `data_provenance` | Source publication or dataset | manual |
| `dataset_description` | Free text description of the dataset | manual |
| `genome_assembly` | e.g. `hg38` | from runs.tsv |
| `transcriptome_reference` | e.g. `GENCODE v43` | from runs.tsv |
| `gene_symbols` | Gene symbol convention used | from runs.tsv |
| `data_preprocessing` | How input data was processed | manual |
| `citation` | Publication to cite | manual |
| `model_name` | e.g. `multiome_powerlaw_v3` | auto |
| `score_threshold` | Score cutoff for thresholded predictions | auto |
| `is_preferred` | Recommended model for this biosample | auto |
| `n_cells` | Number of cells | auto |
| `n_fragments` | Total ATAC fragments | auto |
| `n_umis` | Total RNA UMIs | auto |
| `qc_pass` | Whether QC thresholds are met | auto |
| `oak_predictions_full` | OAK path to full (unthresholded) predictions TSV | auto |
| `oak_predictions_thresholded` | OAK path to thresholded predictions TSV | auto |
| `stream_predictions_bedpe` | Mitra stream URL for thresholded BEDPE | auto (if OAK file exists) |
| `atac_bigwig_url` | Mitra stream URL for normalized ATAC bigwig | auto (if OAK file exists) |
| `rna_matrix_file` | OAK path to input RNA count matrix | auto |
| `atac_frag_file` | OAK path to input ATAC fragment file | auto |
| `notes` | Cell-type-specific caveats | manual |
