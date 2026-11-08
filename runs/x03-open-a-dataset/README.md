# Run x03: opening a real Perturb-seq file

Opens the pilot experiment from Adamson et al. 2016, as harmonized by scPerturb, and asks how big the matrix is, how much of it
is zero, and whether the targeted genes went down. Nothing is modelled or scored.

Three of four predictions held and one missed. See `run-card.yaml` and `results.md`.

```bash
pip install h5py numpy scipy
curl -L -o a.h5ad "https://zenodo.org/records/13350497/files/AdamsonWeissman2016_GSM2406675_10X001.h5ad?download=1"
python3 -I runs/x03-open-a-dataset/open.py a.h5ad
```

The data is CC BY 4.0 and is not copied here. Cite Adamson et al., Cell 2016, and Peidli et al., scPerturb.
