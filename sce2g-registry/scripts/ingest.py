#!/usr/bin/env python3
"""
Generate completed predictions and run entries from partial templates + scE2G configs.

If --runs_template is omitted and the run_id already exists in runs.tsv,
existing run metadata is used automatically.

Input templates (see templates/ directory):
  template_partial_runs.tsv        -- one row per run
  template_partial_predictions.tsv -- one row per biosample

Partial runs template columns:
  run_id, config_file, results_dir, gene_reference_file, igv_dir,
  scE2G_version, date_produced, contact,
  genome_assembly, transcriptome_reference, gene_symbols, notes

  - config_file: full path to scE2G config YAML; used to derive cell_clusters_config
    and igv_dir if igv_dir is not provided directly. Not stored in runs.tsv.
  - results_dir: full OAK path. If blank, derived from config_file.
  - igv_dir: full OAK path to IGV directory. If blank, derived from config_file.

Partial predictions template columns:
  run_id, biosample_id, biosample_description, dataset,
  data_provenance, dataset_description, data_preprocessing, citation

Auto-filled per (biosample x model):
  model_name, score_threshold, is_preferred, n_cells, n_fragments, n_umis, qc_pass,
  oak_predictions_full, oak_predictions_thresholded,
  stream_predictions_bedpe, atac_bigwig_url,   <- only if OAK file exists
  rna_matrix_file, atac_frag_file

Stream URL mapping: /oak/stanford/groups/engreitz -> https://mitra.stanford.edu/engreitz/oak

Output goes to ../temp/ by default. Review, then run merge.py to add to registry.

Usage:
    python scripts/ingest.py \\
        --predictions_template my_predictions.tsv \\
        [--runs_template my_runs.tsv] \\
        [--output_dir temp/]
"""

import argparse
import glob
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

REGISTRY_DIR = Path(__file__).parent.parent
RUNS_TSV = REGISTRY_DIR / "runs.tsv"

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

PREFERRED_MODEL_ORDER = [
    "multiome_powerlaw_v3",
    "multiome_megamap_v3",
    "scATAC_powerlaw_v3",
    "scATAC_megamap_v3",
]

OAK_PREFIX = "/oak/stanford/groups/engreitz"
MITRA_PREFIX = "https://mitra.stanford.edu/engreitz/oak"


def oak_to_mitra(path):
    return str(path).replace(OAK_PREFIX, MITRA_PREFIX)


def resolve(path_str, base_dir):
    p = Path(str(path_str).strip())
    return p if p.is_absolute() else Path(base_dir) / p


def infer_scE2G_dir(config_file):
    """Infer repo root from config file (assumes config lives 2+ levels deep in repo)."""
    return Path(config_file).resolve().parent.parent.parent


def get_threshold(model_dir):
    stats = glob.glob(str(model_dir / "*_stats.tsv"))
    if not stats:
        return None
    m = re.search(r"threshold([\d.]+)", Path(stats[0]).name)
    return m.group(1) if m else None


def get_pred_paths(model_dir, threshold):
    def _find(*names):
        for name in names:
            p = model_dir / name
            if p.exists():
                return str(p)
        return ""
    full = _find("scE2G_predictions.tsv.gz", "encode_e2g_predictions.tsv.gz")
    thresh = _find(
        f"scE2G_predictions_threshold{threshold}.tsv.gz",
        f"encode_e2g_predictions_threshold{threshold}.tsv.gz",
    )
    return full, thresh


def get_stream_urls(igv_dir, biosample_id, model_name, threshold):
    base = Path(igv_dir) / biosample_id
    bedpe_oak = None
    for name in [
        f"scE2G_predictions_threshold{threshold}.bedpe",
        f"encode_e2g_predictions_threshold{threshold}.bedpe",
    ]:
        candidate = base / model_name / name
        if candidate.exists():
            bedpe_oak = candidate
            break
    bw_oak = base / "ATAC_norm.bw"
    return (
        oak_to_mitra(bedpe_oak) if bedpe_oak else "",
        oak_to_mitra(bw_oak) if bw_oak.exists() else "",
    )


def compute_qc_pass(n_fragments, n_umis, model_name):
    try:
        frags = int(float(n_fragments)) if n_fragments else 0
        umis = int(float(n_umis)) if n_umis else 0
    except (ValueError, TypeError):
        return False
    atac_only = "multiome" not in model_name.lower()
    return frags >= 2_000_000 if atac_only else (frags >= 2_000_000 and umis >= 1_000_000)


def pick_preferred(biosample_rows):
    models = [r["model_name"] for r in biosample_rows]
    if not models:
        return None
    by_model = {r["model_name"]: r for r in biosample_rows}
    for mname in [m for m in models if "multiome" in m.lower()]:
        mr = by_model[mname]
        try:
            frags = int(float(mr["n_fragments"])) if mr["n_fragments"] else 0
            umis = int(float(mr["n_umis"])) if mr["n_umis"] else 0
        except (ValueError, TypeError):
            frags = umis = 0
        if frags >= 2_000_000 and umis < 1_000_000:
            atac_equiv = mname.replace("multiome", "scATAC")
            if atac_equiv in models:
                return atac_equiv
    return next((m for m in PREFERRED_MODEL_ORDER if m in models), models[0])


def scan_biosample(results_dir, igv_dir, biosample_id, clusters_map):
    biosample_path = Path(results_dir) / biosample_id
    if not biosample_path.is_dir():
        return []
    rows = []
    for model_dir in sorted(biosample_path.iterdir()):
        if not model_dir.is_dir():
            continue
        model_name = model_dir.name
        threshold = get_threshold(model_dir)
        if threshold is None:
            continue
        full_pred, thresh_pred = get_pred_paths(model_dir, threshold)
        n_cells = n_fragments = n_umis = ""
        stats_files = glob.glob(str(model_dir / "*_stats.tsv"))
        if stats_files:
            try:
                s = pd.read_csv(stats_files[0], sep="\t").iloc[0]
                n_cells = str(int(s["cell_count"]))
                n_fragments = str(int(s["fragments_total"]))
                n_umis = str(int(s["umi_count"]))
            except Exception:
                pass
        stream_bedpe, atac_bw = "", ""
        if igv_dir:
            stream_bedpe, atac_bw = get_stream_urls(igv_dir, biosample_id, model_name, threshold)
        cluster_info = clusters_map.get(biosample_id, {})
        rows.append({
            "model_name": model_name,
            "score_threshold": threshold,
            "n_cells": n_cells,
            "n_fragments": n_fragments,
            "n_umis": n_umis,
            "oak_predictions_full": full_pred,
            "oak_predictions_thresholded": thresh_pred,
            "stream_predictions_bedpe": stream_bedpe,
            "atac_bigwig_url": atac_bw,
            "rna_matrix_file": cluster_info.get("rna_matrix_file", ""),
            "atac_frag_file": cluster_info.get("atac_frag_file", ""),
        })
    pref = pick_preferred(rows)
    for r in rows:
        r["is_preferred"] = str(r["model_name"] == pref)
        r["qc_pass"] = str(compute_qc_pass(r["n_fragments"], r["n_umis"], r["model_name"]))
    return rows


def validate_oak_paths(preds_df, runs_df):
    issues = []
    pred_oak = ["oak_predictions_full", "oak_predictions_thresholded",
                "rna_matrix_file", "atac_frag_file"]
    run_oak = ["results_dir", "cell_clusters_config"]
    for _, row in preds_df.iterrows():
        label = f"{row['biosample_id']} / {row['model_name']}"
        for col in pred_oak:
            val = str(row.get(col, "")).strip()
            if val and val != "nan" and not Path(val).exists():
                issues.append(f"  [predictions | {col}] {val}  -- {label}")
    for _, row in runs_df.iterrows():
        for col in run_oak:
            val = str(row.get(col, "")).strip()
            if val and val != "nan" and not Path(val).exists():
                issues.append(f"  [runs | {col}] {val}  -- run: {row['run_id']}")
    return issues


def resolve_run_info(run_id, runs_template_row, existing_runs):
    """
    Build a dict of run metadata from the template row (if provided)
    or from the existing runs.tsv entry.
    Returns (run_dict, results_dir, igv_dir, cell_clusters_path) or raises.
    """
    # Start from existing registry entry if available
    existing = {}
    if existing_runs is not None and run_id in existing_runs.index:
        existing = existing_runs.loc[run_id].to_dict()

    row = runs_template_row if runs_template_row is not None else {}

    def pick(key, default=""):
        # Template row takes precedence over existing
        v = str(row.get(key, "")).strip()
        if v:
            return v
        v = str(existing.get(key, "")).strip()
        return v if v and v != "nan" else default

    config_file = str(row.get("config_file", "")).strip()
    config = {}
    scE2G_dir = None
    if config_file and Path(config_file).exists():
        with open(config_file) as f:
            config = yaml.safe_load(f)
        scE2G_dir = infer_scE2G_dir(config_file)
    elif config_file:
        print(f"  WARNING: config_file not found: {config_file}", file=sys.stderr)

    # results_dir: template > existing > config
    results_dir_str = pick("results_dir")
    if not results_dir_str and config.get("results_dir") and scE2G_dir:
        results_dir_str = str(resolve(config["results_dir"], scE2G_dir))
    if not results_dir_str:
        raise ValueError(f"No results_dir found for run '{run_id}'")
    results_dir = Path(results_dir_str)

    # igv_dir: template > existing > config
    igv_dir_str = pick("igv_dir")
    if not igv_dir_str and config.get("IGV_dir") and scE2G_dir:
        igv_dir_str = str(resolve(config["IGV_dir"], scE2G_dir))
    igv_dir = Path(igv_dir_str) if igv_dir_str else None

    # cell_clusters: existing > config (never manually provided in template)
    cell_clusters_str = str(existing.get("cell_clusters_config", "")).strip()
    if (not cell_clusters_str or cell_clusters_str == "nan") and config.get("cell_clusters") and scE2G_dir:
        cell_clusters_str = str(resolve(config["cell_clusters"], scE2G_dir))
    cell_clusters_path = Path(cell_clusters_str) if cell_clusters_str else None

    run_dict = {
        "run_id": run_id,
        "scE2G_version": pick("scE2G_version"),
        "date_produced": pick("date_produced"),
        "contact": pick("contact"),
        "results_dir": str(results_dir),
        "cell_clusters_config": str(cell_clusters_path) if cell_clusters_path else "",
        "gene_reference_file": pick("gene_reference_file"),
        "igv_dir": str(igv_dir) if igv_dir else "",
        "genome_assembly": pick("genome_assembly", "hg38"),
        "transcriptome_reference": pick("transcriptome_reference"),
        "gene_symbols": pick("gene_symbols"),
        "notes": pick("notes"),
    }

    return run_dict, results_dir, igv_dir, cell_clusters_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions_template", required=True,
                        help="Partial predictions TSV (one row per biosample)")
    parser.add_argument("--runs_template",
                        help="Partial runs TSV (one row per run). Optional if run already in runs.tsv.")
    parser.add_argument("--output_dir", default=str(REGISTRY_DIR / "temp"),
                        help="Directory for output files (default: temp/)")
    args = parser.parse_args()

    preds_partial = pd.read_csv(args.predictions_template, sep="\t", dtype=str).fillna("")

    runs_partial = None
    if args.runs_template:
        runs_partial = pd.read_csv(args.runs_template, sep="\t", dtype=str).fillna("")

    # Load existing runs for fallback
    existing_runs = None
    if RUNS_TSV.exists():
        existing_runs = pd.read_csv(RUNS_TSV, sep="\t", dtype=str).fillna("").set_index("run_id")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_preds = []
    all_runs = []

    for run_id in preds_partial["run_id"].unique():
        print(f"\n── Run: {run_id} ──")

        # Get template row for this run (if runs_template provided)
        template_row = None
        if runs_partial is not None:
            match = runs_partial[runs_partial["run_id"] == run_id]
            if not match.empty:
                template_row = match.iloc[0].to_dict()
            else:
                print(f"  Note: run '{run_id}' not in runs_template; using existing runs.tsv entry")

        if template_row is None and (existing_runs is None or run_id not in existing_runs.index):
            print(f"  ERROR: no run info for '{run_id}' — provide --runs_template or add to runs.tsv first",
                  file=sys.stderr)
            continue

        try:
            run_dict, results_dir, igv_dir, cell_clusters_path = resolve_run_info(
                run_id, template_row, existing_runs
            )
        except ValueError as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            continue

        print(f"  results_dir:   {results_dir}")
        print(f"  igv_dir:       {igv_dir}")
        print(f"  cell_clusters: {cell_clusters_path}")

        # Load cell clusters for input file paths
        clusters_map = {}
        if cell_clusters_path and Path(cell_clusters_path).exists():
            try:
                cdf = pd.read_csv(cell_clusters_path, sep="\t", dtype=str)[
                    ["cluster", "rna_matrix_file", "atac_frag_file"]
                ]
                clusters_map = cdf.set_index("cluster").to_dict("index")
            except Exception as e:
                print(f"  WARNING: could not read cell_clusters: {e}", file=sys.stderr)

        all_runs.append(run_dict)

        run_preds = preds_partial[preds_partial["run_id"] == run_id]
        missing = []

        for _, bp in run_preds.iterrows():
            biosample_id = bp["biosample_id"]
            model_rows = scan_biosample(results_dir, igv_dir, biosample_id, clusters_map)
            if not model_rows:
                missing.append(biosample_id)
                continue
            print(f"  {biosample_id}: {', '.join(r['model_name'] for r in model_rows)}")
                for mr in model_rows:
                all_preds.append({
                    "run_id": run_id,
                    "biosample_id": biosample_id,
                    "biosample_description": bp.get("biosample_description", ""),
                    "dataset": bp.get("dataset", ""),
                    "data_provenance": bp.get("data_provenance", ""),
                    "dataset_description": bp.get("dataset_description", ""),
                    "genome_assembly": run_dict["genome_assembly"],
                    "transcriptome_reference": run_dict["transcriptome_reference"],
                    "gene_symbols": run_dict["gene_symbols"],
                    "data_preprocessing": bp.get("data_preprocessing", ""),
                    "citation": bp.get("citation", ""),
                    "notes": "",
                    **mr,
                })

        if missing:
            print(f"  WARNING: no results found for: {', '.join(missing)}")

    if not all_preds:
        print("\nNo predictions generated.", file=sys.stderr)
        sys.exit(1)

    preds_df = pd.DataFrame(all_preds)[PREDICTIONS_COLS]
    runs_df = pd.DataFrame(all_runs)[RUNS_COLS]

    issues = validate_oak_paths(preds_df, runs_df)
    if issues:
        print(f"\nWARNING: {len(issues)} OAK path(s) do not exist:")
        for i in issues:
            print(i)
    else:
        print(f"\nAll OAK paths verified.")

    stream_n = (preds_df["stream_predictions_bedpe"].str.strip() != "").sum()
    bw_n = (preds_df["atac_bigwig_url"].str.strip() != "").sum()
    print(f"Stream bedpe URLs: {stream_n}/{len(preds_df)} populated")
    print(f"ATAC bigwig URLs:  {bw_n}/{len(preds_df)} populated")

    run_ids = "_".join(runs_df["run_id"].tolist())
    out_preds = output_dir / f"{run_ids}_predictions.tsv"
    out_runs = output_dir / f"{run_ids}_runs.tsv"
    preds_df.to_csv(out_preds, sep="\t", index=False)
    runs_df.to_csv(out_runs, sep="\t", index=False)

    print(f"\nOutput ({len(preds_df)} prediction rows, {len(runs_df)} run row(s)):")
    print(f"  {out_preds}")
    print(f"  {out_runs}")
    print(f"\nReview, then run: python scripts/merge.py {out_preds} {out_runs}")


if __name__ == "__main__":
    main()
