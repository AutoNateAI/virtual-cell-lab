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
| [x04](runs/x04-scanpy-walkthrough/) | Can a computer find cell types nobody labelled, and how do you tell a real difference from luck? | [results](runs/x04-scanpy-walkthrough/results.md) |
| [x05](runs/x05-small-inference/) | What does it look like to run a virtual cell model, and what can't its answer tell you? | registered, not yet run |
| [001](runs/001-open-the-h1-data/) | What did Arc actually measure for the 2025 Challenge, and does it match the description? | [results](runs/001-open-the-h1-data/results.md) |
| [002](runs/002-leakage-demo/) | How much of a great score comes from splitting the data wrong? | [results](runs/002-leakage-demo/results.md) |

## Licences

The code here is MIT licensed. Saved records keep their original terms: see each run's `snapshots/README.md`.

Independent project, not affiliated with or reviewed by any researcher or institution named. General information, not
medical advice.
