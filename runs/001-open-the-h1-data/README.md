# Run 001: opening the 2025 Virtual Cell Challenge H1 data

Opens Arc's three 2025 Challenge files (training, validation, test) and checks them against Arc's own description before any model
sees them: targets, splits, depth, whether each knockdown lowered its target, and the control cells. Nothing is modelled or scored.

All five predictions held. Every file turned out to hold the same 38,176 control cells. See `run-card.yaml` and `results.md`.

```bash
B=https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025
mkdir -p vcc/train vcc/validation vcc/test
curl -L -C - -o vcc/train/adata_Training.h5ad $B/train/adata_Training.h5ad        # 15.5 GB
curl -L -C - -o vcc/validation/adata_Validation.h5ad $B/validation/adata_Validation.h5ad   # 6.9 GB
curl -L -C - -o vcc/test/adata_Test.h5ad $B/test/adata_Test.h5ad                  # 12.0 GB
pip install h5py numpy scipy pandas "anndata>=0.11"
python3 -I runs/001-open-the-h1-data/open.py vcc
```

Reads the files in row chunks; peak memory was about 5 GB. The data is CC0 1.0 (Arc Virtual Cell Atlas) and is not copied here.
Cite Arc Institute's Virtual Cell Challenge (Roohani et al., Cell 2025).
