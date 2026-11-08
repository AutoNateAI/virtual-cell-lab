# Run x01: one real gene

Looks CFTR up in HGNC, NCBI Gene, RefSeq, Ensembl, UniProt and ClinVar, translates its coding letters three at a time,
compares that with the published protein, and removes the three letters of the delta F508 variant.

`run-card.yaml` holds the predictions, written before any record was opened, along with what was observed and what we
wouldn't trust. `results.md` is the readable result.

To re-run it (Python 3, no dependencies):

```bash
python3 -I runs/x01-one-real-gene/lookup.py
```

It downloads the current records into `snapshots/`, so results can change as the databases are updated.
