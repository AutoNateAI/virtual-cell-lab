#!/usr/bin/env python3
"""Run x02: put two real human cell types side by side and see which genes give them away.

Takes the Tabula Sapiens lung dataset from a pinned CZ CELLxGENE Census release, samples the same number of cells of each
type, and reports which genes each type has on. No model is involved and nothing is scored. Only summary tables are kept
here; the per-cell data stays in the Census.

    pip install cellxgene-census numpy pandas
    python3 -I runs/x02-two-cell-types/compare.py
"""
import json
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CENSUS_VERSION = "2025-11-08"
DATASET = "0d2ee4ac-05ee-40b2-afb6-ebb584caa867"  # Tabula Sapiens - Lung
ASSAY = "10x 3' v3"
TYPES = {"lung": "pulmonary alveolar type 2 cell", "immune": "CD4-positive, alpha-beta T cell"}
CELLS_PER_TYPE = 1000
SEED = 13
ON = 0.10  # a gene counts as "on" in a cell type if at least this share of its cells have any RNA from it
WATCH = ["SFTPC", "CD3E", "ACTB", "CFTR"]


def main():
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        human = census["census_data"]["homo_sapiens"]
        obs = human.obs.read(value_filter=f"dataset_id == '{DATASET}' and assay == \"{ASSAY}\"", column_names=["soma_joinid", "cell_type", "donor_id"]).concat().to_pandas()
        rng = np.random.default_rng(SEED)
        chosen, donors = {}, {}
        for key, name in TYPES.items():
            pool = obs[obs.cell_type.astype(str) == name]
            ids = np.sort(rng.choice(pool.soma_joinid.to_numpy(), size=min(CELLS_PER_TYPE, len(pool)), replace=False))
            chosen[key] = ids
            donors[key] = pool[pool.soma_joinid.isin(ids)].donor_id.astype(str).value_counts().to_dict()
        coords = np.concatenate(list(chosen.values()))
        adata = cellxgene_census.get_anndata(census, organism="Homo sapiens", obs_coords=np.sort(coords), var_value_filter="feature_type == 'protein_coding'")

    names = adata.var["feature_name"].astype(str).to_numpy()
    matrix = adata.X.tocsr()
    joinids = adata.obs["soma_joinid"].to_numpy()
    out = {"census_version": CENSUS_VERSION, "dataset_id": DATASET, "assay": ASSAY, "protein_coding_genes": int(len(names)), "seed": SEED, "on_threshold": ON, "types": {}}
    summary = {}
    # Added after the first result, and labelled as such: how many RNA molecules sit behind each "detected"?
    out["post_hoc_counts"] = {}
    for key, name in TYPES.items():
        sub = matrix[np.isin(joinids, chosen[key])]
        totals = np.asarray(sub.sum(axis=1)).ravel()
        detected_per_cell = np.asarray((sub > 0).sum(axis=1)).ravel()
        fraction_cells = np.asarray((sub > 0).mean(axis=0)).ravel()
        normalized = sub.multiply(1 / np.maximum(totals, 1)[:, None] * 1e4).tocsr()
        normalized.data = np.log1p(normalized.data)
        mean_log = np.asarray(normalized.mean(axis=0)).ravel()
        summary[key] = (fraction_cells, mean_log)
        watched = {g: np.asarray(sub[:, np.where(names == g)[0][0]].todense()).ravel() for g in WATCH}
        out["post_hoc_counts"][key] = {g: {"mean_molecules_per_cell": round(float(v.mean()), 2), "median_molecules_in_detecting_cells": float(np.median(v[v > 0])) if (v > 0).any() else 0.0, "share_of_cell_rna_percent": round(float(100 * v.sum() / totals.sum()), 3)} for g, v in watched.items()}
        out["types"][key] = {
            "cell_type": name, "cells": int(sub.shape[0]), "donors": donors[key],
            "median_genes_detected_per_cell": float(np.median(detected_per_cell)),
            "median_share_of_protein_coding_genes_per_cell": float(np.median(detected_per_cell) / len(names)),
            "genes_on": int((fraction_cells >= ON).sum()),
        }
    on_a, on_b = summary["lung"][0] >= ON, summary["immune"][0] >= ON
    out["genes_on_in_both"] = int((on_a & on_b).sum())
    out["genes_on_in_lung_only"] = int((on_a & ~on_b).sum())
    out["genes_on_in_immune_only"] = int((on_b & ~on_a).sum())
    out["genes_on_in_neither"] = int((~on_a & ~on_b).sum())

    table = pd.DataFrame({"gene": names, "lung_fraction": summary["lung"][0], "immune_fraction": summary["immune"][0], "lung_mean_log": summary["lung"][1], "immune_mean_log": summary["immune"][1]})
    table["difference"] = table.lung_mean_log - table.immune_mean_log
    out["watched"] = table[table.gene.isin(WATCH)].set_index("gene").round(4).to_dict("index")
    top_lung = table.nlargest(10, "difference").round(4)
    top_immune = table.nsmallest(10, "difference").round(4)
    out["top_lung"] = top_lung.gene.tolist()
    out["top_immune"] = top_immune.gene.tolist()
    pd.concat([top_lung, top_immune]).to_csv(HERE / "top_genes.csv", index=False)
    (HERE / "results.json").write_text(json.dumps(out, indent=2) + "\n")
    (HERE / "results.md").write_text(report(out))
    print(report(out))


def report(o):
    w, a, b = o["watched"], o["types"]["lung"], o["types"]["immune"]
    pct = lambda v: f"{100 * v:.1f}%"
    lines = ["# Run x02: two cells, side by side", "", f"Generated by `compare.py` from CZ CELLxGENE Census {o['census_version']}, Tabula Sapiens - Lung ({o['assay']}).", "", "| | Type 2 alveolar cell | CD4 T cell |", "|---|---|---|",
             f"| Cells sampled | {a['cells']:,} | {b['cells']:,} |", f"| Median genes detected per cell | {a['median_genes_detected_per_cell']:,.0f} | {b['median_genes_detected_per_cell']:,.0f} |",
             f"| Share of the {o['protein_coding_genes']:,} protein-coding genes | {pct(a['median_share_of_protein_coding_genes_per_cell'])} | {pct(b['median_share_of_protein_coding_genes_per_cell'])} |", f"| Genes on in at least {o['on_threshold']:.0%} of cells | {a['genes_on']:,} | {b['genes_on']:,} |",
             "", f"Genes on in both: {o['genes_on_in_both']:,}. Only in the lung cell: {o['genes_on_in_lung_only']:,}. Only in the T cell: {o['genes_on_in_immune_only']:,}. In neither: {o['genes_on_in_neither']:,}.", "", "## Share of cells with any RNA from the gene", "", "| Gene | Type 2 alveolar cell | CD4 T cell |", "|---|---|---|"]
    lines += [f"| {g} | {pct(w[g]['lung_fraction'])} | {pct(w[g]['immune_fraction'])} |" for g in WATCH]
    lines += ["", "## Added after the first result: molecules behind each detection", "", "| Gene | Lung cell: mean per cell | T cell: mean per cell | Share of the lung cell's RNA | Share of the T cell's RNA |", "|---|---|---|---|---|"]
    pc = o["post_hoc_counts"]
    lines += [f"| {g} | {pc['lung'][g]['mean_molecules_per_cell']:,} | {pc['immune'][g]['mean_molecules_per_cell']:,} | {pc['lung'][g]['share_of_cell_rna_percent']}% | {pc['immune'][g]['share_of_cell_rna_percent']}% |" for g in WATCH]
    lines += ["", "## Genes that most separate them", "", f"Higher in the lung cell: {', '.join(o['top_lung'])}", "", f"Higher in the T cell: {', '.join(o['top_immune'])}", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
