#!/usr/bin/env python3
"""
Merge reviewed ingest output into the main registry.

Usage:
    python scripts/merge.py temp/{run_id}_predictions.tsv temp/{run_id}_runs.tsv

    # Merge predictions only (run already in registry):
    python scripts/merge.py temp/{run_id}_predictions.tsv

Options:
    --force   Overwrite existing entries for the same run_id (default: error)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

REGISTRY_DIR = Path(__file__).parent.parent
RUNS_TSV = REGISTRY_DIR / "runs.tsv"
PREDICTIONS_TSV = REGISTRY_DIR / "predictions.tsv"

RUNS_COLS = [
    "run_id", "scE2G_version", "date_produced", "contact",
    "results_dir", "cell_clusters_config", "gene_reference_file",
    "igv_dir", "genome_assembly", "transcriptome_reference", "gene_symbols",
    "notes",
]

PREDICTIONS_COLS = [
    "run_id", "biosample_id", "biosample_description", "dataset",
    "data_provenance", "dataset_description", "genome_assembly",
    "transcriptome_reference", "gene_symbols", "data_preprocessing", "citation",
    "model_name", "score_threshold", "is_preferred",
    "n_cells", "n_fragments", "n_umis", "qc_pass",
    "oak_predictions_full", "oak_predictions_thresholded",
    "stream_predictions_bedpe", "atac_bigwig_url",
    "rna_matrix_file", "atac_frag_file",
    "notes",
]


def load(path, cols):
    if Path(path).exists():
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        for col in cols:
            if col not in df.columns:
                df[col] = ""
        return df[cols]
    return pd.DataFrame(columns=cols)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("predictions_tsv", help="Reviewed predictions TSV to merge")
    parser.add_argument("runs_tsv", nargs="?", help="Reviewed runs TSV to merge (optional)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing entries for the same run_id")
    args = parser.parse_args()

    new_preds = pd.read_csv(args.predictions_tsv, sep="\t", dtype=str).fillna("")
    run_ids = new_preds["run_id"].unique().tolist()

    new_runs = None
    if args.runs_tsv:
        new_runs = pd.read_csv(args.runs_tsv, sep="\t", dtype=str).fillna("")

    # Load registry
    registry_preds = load(PREDICTIONS_TSV, PREDICTIONS_COLS)
    registry_runs = load(RUNS_TSV, RUNS_COLS)

    # Check for conflicts
    existing_run_ids = set(registry_runs["run_id"].tolist())
    existing_biosample_keys = set(
        zip(registry_preds["run_id"], registry_preds["biosample_id"], registry_preds["model_name"])
    )

    conflicts_runs = [r for r in run_ids if r in existing_run_ids and new_runs is not None]
    conflicts_preds = [
        (r, b, m) for r, b, m in zip(new_preds["run_id"], new_preds["biosample_id"], new_preds["model_name"])
        if (r, b, m) in existing_biosample_keys
    ]

    if (conflicts_runs or conflicts_preds) and not args.force:
        if conflicts_runs:
            print(f"ERROR: run(s) already in runs.tsv: {conflicts_runs}")
        if conflicts_preds:
            print(f"ERROR: {len(conflicts_preds)} (run, biosample, model) row(s) already in predictions.tsv")
            for r, b, m in conflicts_preds[:5]:
                print(f"  {r} / {b} / {m}")
            if len(conflicts_preds) > 5:
                print(f"  ... and {len(conflicts_preds) - 5} more")
        print("Use --force to overwrite.")
        sys.exit(1)

    # Remove existing entries for these run_ids if --force
    if args.force:
        registry_runs = registry_runs[~registry_runs["run_id"].isin(run_ids)]
        registry_preds = registry_preds[~registry_preds["run_id"].isin(run_ids)]

    # Merge runs
    if new_runs is not None:
        for col in RUNS_COLS:
            if col not in new_runs.columns:
                new_runs[col] = ""
        registry_runs = pd.concat([registry_runs, new_runs[RUNS_COLS]], ignore_index=True)

    # Merge predictions
    for col in PREDICTIONS_COLS:
        if col not in new_preds.columns:
            new_preds[col] = ""
    registry_preds = pd.concat([registry_preds, new_preds[PREDICTIONS_COLS]], ignore_index=True)

    registry_runs.to_csv(RUNS_TSV, sep="\t", index=False)
    registry_preds.to_csv(PREDICTIONS_TSV, sep="\t", index=False)

    print(f"Merged {len(new_preds)} prediction rows for run(s): {run_ids}")
    if new_runs is not None:
        print(f"Merged {len(new_runs)} run row(s) into runs.tsv")
    print(f"Registry now: {len(registry_preds)} predictions, {len(registry_runs)} runs")


if __name__ == "__main__":
    main()
