# Run x02: two cells, side by side

Compares 1,000 type 2 alveolar cells with 1,000 CD4 T cells from the Tabula Sapiens lung dataset, through a pinned
CZ CELLxGENE Census release. Nothing is modelled or scored.

Three of the five predictions held and two missed. `run-card.yaml` has the predictions, written before any expression
value was read, along with what we would not trust. `results.md` is the readable result.

```bash
pip install cellxgene-census numpy pandas
python3 -I runs/x02-two-cell-types/compare.py
```

Only summary tables are kept here. Cite: The Tabula Sapiens Consortium, Science 2022, doi 10.1126/science.abl4896.
