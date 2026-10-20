# Coding Exercise: Cleaning RNA-seq Data

This exercise takes about **40 minutes**. There is no single correct solution: we care about your reasoning and how you document it more than about a perfect result. 

## The situation

We treated human cells with a drug and want to know which genes it switches on or off. Eight samples were sequenced (4 control, 4 treated). For every gene and sample, we counted how many sequencing reads came from that gene: **more reads means the gene is more active**.

You get two files, put together by hand and not very carefully:

- `data/counts_raw.csv`: one row per gene, one column per sample, read counts in the cells
- `data/samples_raw.csv`: one row per sample with its condition, sequencing date and RNA quality score (RIN, 1–10, higher is better)

The data is simulated; the gene IDs are not real.

## Setup

1. Click **"Use this template"** on GitHub to create your own **private** copy and add **namsaraeva** as a collaborator.
2. Clone it and create a branch called `solution`.
3. Use Python 3.10+ with `pandas` and `numpy`. Script or notebook, your choice.

## Tasks

If you run out of time, stop and write in `NOTES.md` what you would do next. Prioritising is part of the task.

**1. Clean.** Produce `output/counts_clean.csv` (gene IDs as rows, integer counts, one column per sample named like the sample sheet) and `output/samples_clean.csv` (one row per sample, consistent values). The sample columns of the count table must match the rows of the sample sheet. List the problems you found and every decision you made in `NOTES.md`.

**2. Normalise with numpy.** Samples were sequenced to different depths, so raw counts are not comparable. Write a function using only NumPy, without loops:

```python
def log_cpm(counts: np.ndarray) -> np.ndarray:
    """genes x samples counts -> log2(counts per million + 1)"""
```

where counts per million = count / total count of that sample × 1,000,000.

**3. Analyse.**

1. For each gene, compute the log2 fold change = mean log-CPM of treated samples − mean log-CPM of control samples.
2. Report the 5 most strongly up- and the 5 most strongly down-regulated genes (by symbol).
3. Does any sample look suspicious? How did you check, and what did you do about it?

Put your results in `NOTES.md`.

## Git

We read your commit history as part of your solution: commit in small logical steps with clear messages, and keep the repo free of clutter. When done, push your branch, **open a Pull Request from `solution` into `main`** with a short summary, and send the link to **a.namsaraeva@tum.de** by **28.10.2026**.

## Rules

Use any documentation you like. If you use AI assistants, say so in `NOTES.md` and be ready to explain every line. Do not edit the files in `data/`.
