#!/usr/bin/env python3
"""
Remove a run_id from the registry (predictions.tsv and runs.tsv).

Usage:
    python scripts/remove.py <run_id> [<run_id> ...]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

REGISTRY_DIR = Path(__file__).parent.parent
RUNS_TSV = REGISTRY_DIR / "runs.tsv"
PREDICTIONS_TSV = REGISTRY_DIR / "predictions.tsv"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_ids", nargs="+", help="run_id(s) to remove")
    args = parser.parse_args()

    run_ids = args.run_ids

    runs = pd.read_csv(RUNS_TSV, sep="\t", dtype=str).fillna("")
    preds = pd.read_csv(PREDICTIONS_TSV, sep="\t", dtype=str).fillna("")

    missing = [r for r in run_ids if r not in runs["run_id"].values and r not in preds["run_id"].values]
    if missing:
        print(f"ERROR: run_id(s) not found in registry: {missing}")
        sys.exit(1)

    for run_id in run_ids:
        n_preds = (preds["run_id"] == run_id).sum()
        n_runs = (runs["run_id"] == run_id).sum()
        print(f"Removing {run_id}: {n_preds} prediction rows, {n_runs} run row(s)")

    runs = runs[~runs["run_id"].isin(run_ids)]
    preds = preds[~preds["run_id"].isin(run_ids)]

    runs.to_csv(RUNS_TSV, sep="\t", index=False)
    preds.to_csv(PREDICTIONS_TSV, sep="\t", index=False)

    print(f"Registry now: {len(preds)} predictions, {len(runs)} runs")


if __name__ == "__main__":
    main()
