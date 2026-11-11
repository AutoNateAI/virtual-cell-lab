# Run 002: how to cheat at predicting cells, by accident

A labelled leakage demonstration on the 2025 Challenge training file. The same three simple models are scored on a leaky split,
where every perturbation has cells on both sides, and a fair one, where 30 perturbations are held out whole (`holdout.txt`).
A model that only remembers each perturbation's average scored 0.76 (mean Pearson delta) on the leaky split and 0.21 on the fair one.
Three of five predictions held. See `run-card.yaml` and `results.md`.

```bash
# needs the training file from run 001's README, and runs/001-open-the-h1-data/h5chunks.py next to this folder
python3 -I runs/002-leakage-demo/leakage.py vcc/train/adata_Training.h5ad runs/002-leakage-demo/holdout.txt
```

About 7 minutes on a MacBook Air M1, peak memory about 8 GB. Data CC0 1.0, not copied here.
