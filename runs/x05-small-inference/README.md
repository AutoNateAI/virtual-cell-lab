# Run x05: ask a pretrained virtual cell model what happens next

Gives Arc Institute's State model real K562 control cells and asks what five gene knockdowns it never trained on in K562 would
do, plus one gene name it has never heard of. It looks at what comes back. Nothing is scored against measured data here.

**Registered, not yet run.** The four predictions are in `run-card.yaml`, committed before the model was downloaded.

[Open in Colab](https://colab.research.google.com/github/AutoNateAI/virtual-cell-lab/blob/main/runs/x05-small-inference/inference.ipynb),
choose a T4 GPU (Runtime → Change runtime type), then Run all. It downloads about 1.6 GB and saves `results.json`.

## Licence

State's weights and outputs are under the [Arc Research Institute State Model Non-Commercial License](https://github.com/ArcInstitute/state/blob/main/MODEL_LICENSE.md).
Use is allowed only for non-commercial purposes, and research or education engaged in or sponsored by a commercial entity does
not count. The notebook refuses to download anything until you set `ACCEPT_STATE_NONCOMMERCIAL_LICENSE = True`. Only do that if
your use qualifies; if unsure, ask Arc first. This is a summary, not legal advice.

Any output shown must cite: Adduri, A. et al. (2025) *Predicting cellular responses to perturbation across diverse contexts
with State.* bioRxiv 10.1101/2025.06.26.661135; Cell 2026, doi 10.1016/j.cell.2026.07.052. Data: Replogle et al., Cell 2022.

Everything this run produces is model output, not a biological finding.
