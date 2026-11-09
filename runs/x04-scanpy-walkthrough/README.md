# Run x04: the standard Scanpy workflow on 2,700 blood cells

Takes a public single-cell file through quality control, PCA, UMAP and Leiden clustering, then asks how much of the answer depends
on our choices and how you tell a real difference in a gene from luck. No model is scored.

Six of seven predictions held and one missed: seven genes passed multiple-testing correction on shuffled labels, which carry no signal.
See `run-card.yaml` and `results.md`. Section 10 of the notebook is post-hoc and labelled so.

[Open in Colab](https://colab.research.google.com/github/AutoNateAI/virtual-cell-lab/blob/main/runs/x04-scanpy-walkthrough/walkthrough.ipynb)
and choose Run all. Not yet tested in Colab; it has been tested from a clean local install.

To run it locally on Python 3.12:

```bash
cd runs/x04-scanpy-walkthrough
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt jupyterlab
jupyter lab walkthrough.ipynb
```

The notebook downloads about 7.6 MB from 10x Genomics into `data/`, checks its SHA-256, and overwrites `results.json`, `results.md` and
`umap.json`, which should come out unchanged.

The data is not copied here. It is *3k PBMCs from a Healthy Donor*, 10x Genomics; cite that and Wolf, Angerer & Theis, Genome Biology 2018.
