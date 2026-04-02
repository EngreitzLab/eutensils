#!/usr/bin/env python3
"""
Marimo web app for browsing scE2G predictions.

Usage:
    marimo run scripts/explore.py

Or open in browser (no setup):
    https://marimo.app/github.com/EngreitzLab/eutensils/blob/main/sce2g-registry/scripts/explore.py
"""

import marimo as mo

app = mo.App(width="full")


@app.cell
def __():
    import pandas as pd
    import json
    from pathlib import Path
    import marimo as mo
    return pd, json, Path, mo


@app.cell
def __(pd, Path):
    GITHUB_RAW = "https://raw.githubusercontent.com/EngreitzLab/eutensils/main/sce2g-registry"

    def load_tsv(filename):
        try:
            local = Path(__file__).parent.parent / filename
            if local.exists():
                return pd.read_csv(local, sep="\t", dtype=str).fillna("")
        except NameError:
            pass  # __file__ not defined in WASM
        return pd.read_csv(f"{GITHUB_RAW}/{filename}", sep="\t", dtype=str).fillna("")

    preds = load_tsv("predictions.tsv")
    runs = load_tsv("runs.tsv")
    return GITHUB_RAW, preds, runs


@app.cell
def __(mo):
    mo.md("# scE2G predictions registry")


# ── Filters ───────────────────────────────────────────────────────────────────

@app.cell
def __(mo, preds):
    search = mo.ui.text(placeholder="Search biosample, dataset...", label="Search")

    # Map dataset_description → dataset id for filter labels
    dataset_map = {}
    for _, _row in preds[["dataset", "dataset_description"]].drop_duplicates().iterrows():
        _label = _row["dataset_description"] or _row["dataset"] or "(no dataset)"
        dataset_map[_label] = _row["dataset"]
    dataset_options = sorted(dataset_map)
    dataset_filter = mo.ui.multiselect(dataset_options, value=dataset_options, label="Dataset")

    model_options = sorted(preds["model_name"].unique())
    model_filter = mo.ui.multiselect(model_options, value=model_options, label="Model")

    qc_only = mo.ui.checkbox(label="QC pass only")
    preferred_only = mo.ui.checkbox(label="Preferred only", value=True)
    show_all = mo.ui.checkbox(label="Show all")
    return search, dataset_filter, dataset_map, model_filter, qc_only, preferred_only, show_all


@app.cell
def __(mo, search, dataset_filter, model_filter, qc_only, preferred_only, show_all):
    mo.hstack(
        [search, dataset_filter, model_filter, qc_only, preferred_only, show_all],
        gap=2,
        wrap=True,
    )


@app.cell
def __(preds, search, dataset_filter, dataset_map, model_filter, qc_only, preferred_only, show_all):
    df = preds.copy()
    if not show_all.value:
        if search.value.strip():
            q = search.value.strip().lower()
            mask = (
                df["biosample_id"].str.lower().str.contains(q, na=False)
                | df["biosample_description"].str.lower().str.contains(q, na=False)
                | df["dataset"].str.lower().str.contains(q, na=False)
                | df["dataset_description"].str.lower().str.contains(q, na=False)
            )
            df = df[mask]
        if dataset_filter.value:
            selected_datasets = {dataset_map[v] for v in dataset_filter.value}
            df = df[df["dataset"].isin(selected_datasets)]
        if model_filter.value:
            df = df[df["model_name"].isin(model_filter.value)]
        if qc_only.value:
            df = df[df["qc_pass"] == "True"]
        if preferred_only.value:
            df = df[df["is_preferred"] == "True"]
    return (df,)


# ── Table ─────────────────────────────────────────────────────────────────────

@app.cell
def __(mo, df):
    DISPLAY_COLS = [
        "run_id", "biosample_id", "biosample_description", "dataset",
        "model_name", "is_preferred", "qc_pass",
        "n_cells", "n_fragments", "n_umis",
    ]
    table = mo.ui.table(
        df[[c for c in DISPLAY_COLS if c in df.columns]].reset_index(drop=True),
        selection="multi",
    )
    return (table,)


@app.cell
def __(mo, table, df):
    mo.vstack([
        mo.md(f"**{len(table.value)} selected** · {len(df)} predictions shown"),
        table,
    ])


# ── Resolve full rows for selected entries ────────────────────────────────────

@app.cell
def __(table, preds, runs):
    sel = table.value
    if len(sel) > 0:
        sel_keys = set(zip(sel["run_id"], sel["biosample_id"], sel["model_name"]))
        full_sel = preds[
            [(r, b, m) in sel_keys
             for r, b, m in zip(preds["run_id"], preds["biosample_id"], preds["model_name"])]
        ].copy()
        sel_runs = runs[runs["run_id"].isin(full_sel["run_id"])].copy()
    else:
        full_sel = preds.iloc[:0].copy()
        sel_runs = runs.iloc[:0].copy()
    return full_sel, sel_runs


# ── Download ──────────────────────────────────────────────────────────────────

@app.cell
def __(mo, full_sel, sel_runs):
    if len(full_sel) > 0:
        download_ui = mo.vstack([
            mo.md(f"Download both files to get full metadata for the {len(full_sel)} selected prediction(s):"),
            mo.hstack([
                mo.download(
                    full_sel.to_csv(sep="\t", index=False).encode(),
                    filename="selected_predictions.tsv",
                    label=f"predictions.tsv ({len(full_sel)} rows)",
                ),
                mo.download(
                    sel_runs.to_csv(sep="\t", index=False).encode(),
                    filename="selected_runs.tsv",
                    label=f"runs.tsv ({len(sel_runs)} run(s))",
                ),
            ]),
        ])
    else:
        download_ui = mo.md("*Select rows above to download or build an IGV session.*")
    download_ui


# ── IGV session ───────────────────────────────────────────────────────────────

@app.cell
def __(mo):
    mo.md("---\n## IGV session")


@app.cell
def __(mo):
    locus_input = mo.ui.text(value="chr10:79,017,034-79,273,289", label="Locus")
    include_crispr = mo.ui.checkbox(label="K562 CRISPR tracks")
    include_gwas = mo.ui.checkbox(label="UKBB GWAS track")
    return locus_input, include_crispr, include_gwas


@app.cell
def __(mo, locus_input, include_crispr, include_gwas):
    mo.hstack([locus_input, include_crispr, include_gwas], gap=2)


@app.cell
def __(mo, full_sel, locus_input, include_crispr, include_gwas, json):
    # Nature color palette — cycling per dataset
    ATAC_COLORS     = ["#002359", "#00488d", "#5496ce", "#9bcae9", "#c5e5fb"]
    MULTIOME_COLORS = ["#430b4e", "#792374", "#a64791", "#d3a9ce", "#e9d3ea"]
    SCATAC_COLORS   = ["#003648", "#006479", "#0096a0", "#96ced3", "#cae5ee"]

    CRISPR_TRACKS = [
        {
            "name": "K562 CRISPR positives",
            "url": "https://mitra.stanford.edu/engreitz/oak/igv/hg38/scE2G/crispr_positives.bedpe",
            "color": "#9b241c",
            "height": 100,
        },
        {
            "name": "K562 CRISPR tested elements",
            "url": "https://mitra.stanford.edu/engreitz/oak/igv/hg38/scE2G/crispr_elements.bed",
            "color": "#9b241c",
            "height": 35,
        },
    ]
    GWAS_TRACK = {
        "name": "UKBB GWAS variants",
        "url": "https://mitra.stanford.edu/engreitz/oak/igv/hg38/scE2G/UKBB_GWAS/merged_variants.bed",
        "color": "#9b241c",
        "height": 50,
    }

    if len(full_sel) == 0:
        igv_ui = mo.md("*Select biosamples above to generate an IGV session.*")
    else:
        datasets = sorted(full_sel["dataset"].unique().tolist())
        color_idx = {ds: i % len(ATAC_COLORS) for i, ds in enumerate(datasets)}

        tracks = []
        if include_crispr.value:
            tracks.extend(CRISPR_TRACKS)

        for _, row in full_sel.iterrows():
            model = row.get("model_name", "")
            ci = color_idx.get(row.get("dataset", ""), 0)
            is_multiome = "multiome" in model.lower()
            pred_color = MULTIOME_COLORS[ci] if is_multiome else SCATAC_COLORS[ci]
            atac_color = ATAC_COLORS[ci]
            label = row.get("biosample_description") or row.get("biosample_id", "")
            bedpe = str(row.get("stream_predictions_bedpe", "")).strip()
            bw = str(row.get("atac_bigwig_url", "")).strip()
            if bedpe:
                tracks.append({
                    "name": label,
                    "url": bedpe,
                    "color": pred_color,
                    "arcType": "proportional",
                    "height": 100,
                })
            if bw:
                tracks.append({
                    "name": f"ATAC — {label}",
                    "url": bw,
                    "color": atac_color,
                    "height": 35,
                })

        if include_gwas.value:
            tracks.append(GWAS_TRACK)

        session = {"genome": "hg38", "locus": locus_input.value, "tracks": tracks}
        session_json = json.dumps(session, indent=2)

        div_id = f"igv-div-{abs(hash(session_json)) % 10**8}"
        igv_html = f"""<div id="{div_id}" style="padding-top:10px;"></div>
<script>(function() {{
  var options = {session_json};
  function launch() {{
    var el = document.getElementById('{div_id}');
    if (!el) return;
    el.innerHTML = '';
    if (typeof igv !== 'undefined') {{
      igv.createBrowser(el, options);
    }} else {{
      var s = document.createElement('script');
      s.src = 'https://cdn.jsdelivr.net/npm/igv@3.0.5/dist/igv.min.js';
      s.onload = function() {{ igv.createBrowser(el, options); }};
      document.head.appendChild(s);
    }}
  }}
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', launch);
  else launch();
}})();</script>"""

        igv_ui = mo.vstack([
            mo.md(f"**{len(tracks)} tracks** for {len(full_sel)} selected biosample(s)"),
            mo.download(
                session_json.encode(),
                filename="igv_session.json",
                label=f"Download IGV session JSON ({len(tracks)} tracks)",
            ),
            mo.Html(igv_html),
        ])

    igv_ui


if __name__ == "__main__":
    app.run()
