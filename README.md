# Virtual Cell Lab

The open code and data behind [Inside the Virtual Cell](https://autonateai.com/publications), an independent video series
and research project on predicting how cells respond to change.

Each folder in `runs/` is one experiment: the card written before the result, the script that produced it, the output,
and the raw records it used. What went wrong is written down on the card too.

| Run | What it asks | Result |
|---|---|---|
| [x01](runs/x01-one-real-gene/) | What do the public records say about one gene, CFTR, and do they agree? | [results](runs/x01-one-real-gene/results.md) |
| [x02](runs/x02-two-cell-types/) | How different are the genes switched on in a lung cell and a T cell? | [results](runs/x02-two-cell-types/results.md) |
| [x03](runs/x03-open-a-dataset/) | What is actually in a real Perturb-seq file? | [results](runs/x03-open-a-dataset/results.md) |

## Licences

The code here is MIT licensed. Saved records keep their original terms: see each run's `snapshots/README.md`.

Independent project, not affiliated with or reviewed by any researcher or institution named. General information, not
medical advice.
