#!/usr/bin/env python3
"""Run 002: the same three simple models, scored on a leaky split and a fair one. A labelled leakage demonstration.

    python3 -I runs/002-leakage-demo/leakage.py <path to adata_Training.h5ad> <holdout file>

Data: Virtual Cell Challenge 2025 training file (CC0 1.0).
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "001-open-the-h1-data"))
import h5chunks as h5  # noqa: E402

PERT, CTRL, SEED, N_GENES = "target_gene", "non-targeting", 13, 1000


def main(path, holdout_path):
    held = {l.strip() for l in open(holdout_path) if l.strip() and not l.startswith("#")}
    with h5py.File(path, "r") as f:
        pert = h5.obs_column(f, PERT)
        n, m = h5.shape(f)
        groups = sorted(set(pert) - {""})
        assert CTRL in groups and held <= set(groups), "holdout perturbations or control label missing from the file"
        gi = {g: i for i, g in enumerate(groups)}
        code = np.array([gi.get(p, -1) for p in pert])
        leaky_test = np.random.default_rng(SEED).random(n) < 0.2
        fair_test = np.isin(pert, sorted(held))
        # sums[split][side] : groups x genes ; counts likewise
        sums = {s: np.zeros((2, len(groups), m)) for s in ("leaky", "fair")}
        cnts = {s: np.zeros((2, len(groups))) for s in ("leaky", "fair")}
        for start, block in h5.row_chunks(f):
            rows = slice(start, start + block.shape[0])
            totals = np.asarray(block.sum(axis=1)).ravel()
            scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
            block = block.multiply(scale[:, None]).tocsr()
            block.data = np.log1p(block.data)
            c = code[rows]
            for split, test in (("leaky", leaky_test[rows]), ("fair", fair_test[rows])):
                for side in (0, 1):
                    sel = (test == bool(side)) & (c >= 0)
                    if not sel.any():
                        continue
                    onehot = np.zeros((len(groups), int(sel.sum())))
                    onehot[c[sel], np.arange(int(sel.sum()))] = 1
                    sums[split][side] += np.asarray(block[sel].T @ onehot.T).T
                    cnts[split][side] += onehot.sum(axis=1)

    out = {"n_cells": n, "n_genes": m, "groups": len(groups) - 1, "held_out": sorted(held), "splits": {}}
    for split in ("leaky", "fair"):
        mean = sums[split] / np.maximum(cnts[split], 1)[:, :, None]
        ci = gi[CTRL]
        ctrl = mean[0, ci]
        genes = np.argsort(-ctrl)[:N_GENES]
        train_perts = [g for g in groups if g != CTRL and cnts[split][0, gi[g]] > 0]
        test_perts = [g for g in groups if g != CTRL and cnts[split][1, gi[g]] > 0]
        pmean = np.mean([mean[0, gi[g]] for g in train_perts], axis=0)
        per = []
        for g in test_perts:
            obs = mean[1, gi[g]][genes]
            preds = {"control-mean": ctrl[genes], "perturbation-mean": pmean[genes],
                     "lookup": (mean[0, gi[g]] if g in train_perts else pmean)[genes]}
            row = {"perturbation": g, "test_cells": int(cnts[split][1, gi[g]]), "train_cells": int(cnts[split][0, gi[g]])}
            od = obs - ctrl[genes]
            for name, p in preds.items():
                pd_ = p - ctrl[genes]
                row[f"{name}_pearson_delta"] = None if np.allclose(pd_, 0) else float(np.corrcoef(pd_, od)[0, 1])
                row[f"{name}_mse"] = float(np.mean((p - obs) ** 2))
            per.append(row)
        summ = {}
        for name in ("control-mean", "perturbation-mean", "lookup"):
            pv = [r[f"{name}_pearson_delta"] for r in per if r[f"{name}_pearson_delta"] is not None]
            summ[name] = {"pearson_delta": round(float(np.mean(pv)), 4) if pv else None, "mse": round(float(np.mean([r[f"{name}_mse"] for r in per])), 6)}
        lookup_best = sum(1 for r in per if r["lookup_mse"] < min(r["control-mean_mse"], r["perturbation-mean_mse"]))
        out["splits"][split] = {"test_perturbations": len(test_perts), "train_perturbations": len(train_perts), "summary": summ,
                                "lookup_beats_both_on_mse": lookup_best, "per_perturbation": per}

    L, F = out["splits"]["leaky"]["summary"], out["splits"]["fair"]["summary"]
    out["predictions"] = [
        {"statement": "lookup leaky Pearson delta >= 0.8", "held": L["lookup"]["pearson_delta"] >= 0.8},
        {"statement": "lookup drops by >= 0.4 leaky to fair", "held": L["lookup"]["pearson_delta"] - F["lookup"]["pearson_delta"] >= 0.4},
        {"statement": "perturbation-mean fair Pearson delta > 0.2", "held": F["perturbation-mean"]["pearson_delta"] > 0.2},
        {"statement": "perturbation-mean fair MSE < control-mean fair MSE", "held": F["perturbation-mean"]["mse"] < F["control-mean"]["mse"]},
        {"statement": "lookup beats both baselines on MSE for >= 90% of perturbations (leaky)", "held": out["splits"]["leaky"]["lookup_beats_both_on_mse"] >= 135},
    ]
    (HERE / "results.json").write_text(json.dumps(out, indent=2) + "\n")
    lines = ["# Run 002: leakage demonstration", "", "Generated by `leakage.py`. Same cells, same three models, two splits.", "",
             "| Split | Model | Mean Pearson delta | Mean MSE |", "|---|---|---|---|"]
    for split in ("leaky", "fair"):
        for name, s in out["splits"][split]["summary"].items():
            lines.append(f"| {split} | {name} | {s['pearson_delta'] if s['pearson_delta'] is not None else 'undefined'} | {s['mse']} |")
    lines += ["", f"Leaky split: lookup beats both baselines on MSE for {out['splits']['leaky']['lookup_beats_both_on_mse']} of {out['splits']['leaky']['test_perturbations']} perturbations.",
              "", "| Prediction | Held |", "|---|---|"] + [f"| {p['statement']} | {'yes' if p['held'] else 'NO'} |" for p in out["predictions"]]
    (HERE / "results.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
