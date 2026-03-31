#!/usr/bin/env python3
"""
Validate OAK paths and stream URLs in the scE2G registry.

Usage:
    # Validate the full registry
    python validate.py

    # Also check stream URLs (slow HTTP HEAD requests)
    python validate.py --urls

    # Validate a new template before adding to the registry
    python validate.py --new new_predictions.tsv [--new_runs new_runs.tsv]
    python validate.py --new new_predictions.tsv --urls
"""

import argparse
import sys
from pathlib import Path
import pandas as pd

REGISTRY_DIR = Path(__file__).parent.parent
RUNS_TSV = REGISTRY_DIR / "runs.tsv"
PREDICTIONS_TSV = REGISTRY_DIR / "predictions.tsv"

PRED_OAK_COLS = ["oak_predictions_full", "oak_predictions_thresholded", "rna_matrix_file", "atac_frag_file"]
PRED_URL_COLS = ["stream_predictions_bedpe", "atac_bigwig_url"]
RUN_OAK_COLS = ["results_dir", "cell_clusters_config"]


def check_url(url, session):
    try:
        r = session.head(url, timeout=10, allow_redirects=True)
        return r.status_code < 400
    except Exception:
        return False


def check_oak_paths(preds_df, runs_df, label="registry"):
    issues = []
    for _, row in preds_df.iterrows():
        entry = f"{row.get('biosample_id','')} / {row.get('model_name','')} (run: {row.get('run_id','')})"
        for col in PRED_OAK_COLS:
            val = str(row.get(col, "")).strip()
            if val and val != "nan" and not Path(val).exists():
                issues.append(f"  [predictions | {col}] {val}\n    -- {entry}")
    for _, row in runs_df.iterrows():
        entry = f"run: {row.get('run_id','')}"
        for col in RUN_OAK_COLS:
            val = str(row.get(col, "")).strip()
            if val and val != "nan" and not Path(val).exists():
                issues.append(f"  [runs | {col}] {val}\n    -- {entry}")
    n_checked = sum(
        (preds_df[c].astype(str).str.strip().replace("nan", "") != "").sum()
        for c in PRED_OAK_COLS if c in preds_df.columns
    )
    if issues:
        print(f"\n=== {label}: {len(issues)} missing OAK path(s) ===")
        for i in issues:
            print(i)
    else:
        print(f"{label} OAK paths: all {n_checked} checked OK")
    return issues


def check_urls(preds_df, label="registry"):
    import requests
    session = requests.Session()
    issues = []
    total = 0
    for col in PRED_URL_COLS:
        if col not in preds_df.columns:
            continue
        urls = preds_df[col].astype(str).str.strip().replace("nan", "")
        non_empty = urls[urls != ""]
        print(f"  Checking {len(non_empty)} URLs in {col} ...", flush=True)
        for idx, url in non_empty.items():
            total += 1
            if not check_url(url, session):
                row = preds_df.loc[idx]
                issues.append(
                    f"  [predictions | {col}] {url}\n"
                    f"    -- {row.get('biosample_id','')} / {row.get('model_name','')} (run: {row.get('run_id','')})"
                )
    if issues:
        print(f"\n=== {label}: {len(issues)} unreachable URL(s) ===")
        for i in issues:
            print(i)
    else:
        print(f"{label} stream URLs: all {total} checked OK")
    return issues


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--urls", action="store_true", help="Also check stream URLs (slow)")
    parser.add_argument("--new", metavar="PREDICTIONS_TSV", help="Validate a new predictions template file")
    parser.add_argument("--new_runs", metavar="RUNS_TSV", help="New runs template file to validate alongside --new")
    args = parser.parse_args()

    all_issues = []

    if args.new:
        # Validate a new template file
        new_preds = pd.read_csv(args.new, sep="\t", dtype=str).fillna("")
        new_runs = (
            pd.read_csv(args.new_runs, sep="\t", dtype=str).fillna("")
            if args.new_runs else pd.DataFrame()
        )
        print(f"Validating new file: {args.new} ({len(new_preds)} rows)")
        all_issues += check_oak_paths(new_preds, new_runs, label="new file")
        if args.urls:
            all_issues += check_urls(new_preds, label="new file")
    else:
        # Validate the full registry
        preds = pd.read_csv(PREDICTIONS_TSV, sep="\t", dtype=str).fillna("")
        runs = pd.read_csv(RUNS_TSV, sep="\t", dtype=str).fillna("")
        print(f"Validating registry ({len(preds)} prediction rows, {len(runs)} runs)")
        all_issues += check_oak_paths(preds, runs)
        if args.urls:
            print("Checking stream URLs ...")
            all_issues += check_urls(preds)

    if all_issues:
        print(f"\n{len(all_issues)} total issue(s) found.")
        sys.exit(1)
    else:
        print("\nAll checks passed.")


if __name__ == "__main__":
    main()
