# Notes

## Data problems found

- **Sample IDs:** Inconsistent formats in `samples_raw.csv` (`S1`, `S02`, `s3`, `Sample 4`, `sample_7`).
- **Condition labels:** Inconsistent synonyms (`control`, `Control`, `ctrl`, `untreated` vs `treated`, `Treated`, `drug`).
- **Duplicate sample row:** Sample `s3` was listed twice with identical values.
- **Unsequenced sample:** `S9` was present in `samples_raw.csv` with RIN 2.4, but had no column in `counts_raw.csv`.
- **Low RNA quality:** `S06` (RIN 3.1) and `S9` (RIN 2.4) have low RIN scores.
- **Count column headers:** Full BAM file paths (`/scratch/run1/S2.sorted.bam`) instead of sample IDs.
- **Scrambled sample order:** Column order in `counts_raw.csv` did not match sample sheet rows.
- **Number formatting:** Thousand-separators (`"1,086"`) forced string types, and some integers were stored as floats (`243.0`).
- **Duplicate gene rows:** `CDH2` and `TUBB` appeared twice with identical counts across all samples.
- **Non-gene rows:** Five HTSeq summary statistics rows at the bottom (`__no_feature`, `__ambiguous`, etc.).

---

## Cleaning decisions (what, why, alternatives considered)

- **Standardize sample IDs to S1..S8:** Extracted digits using regex to create consistent IDs matching BAM files. 
Alternative: manual mapping, rejected as less scalable.
- **Standardize conditions to 'control' and 'treated':** Mapped labels with a dictionary. 
Alternative: 0/1 encoding, rejected because string labels are self-documenting.
- **Remove duplicate s3 row:** Dropped the second identical row. 
Alternative: rename to s4, rejected because only one BAM file exists for S3.
- **Exclude S9 from clean sample sheet:** S9 was never sequenced (no BAM file), likely discarded due to RIN 2.4. Dropped so sample sheet rows match count columns. 
Alternative: keep S9 with NaN counts, rejected because task requires integer counts.
- **Strip BAM paths from count headers:** Extracted `S\d+` via regex. 
Alternative: string splitting on slashes, rejected because regex is cleaner.
- **Reorder count columns to S1..S8:** Matched column order to sample sheet rows.
- **Drop HTSeq `__*` rows:** Filtered out summary stats so only true genes remain. 
Alternative: keep at bottom, rejected because task specifies gene IDs as rows.
- **Drop duplicate CDH2 and TUBB rows:** Kept first occurrence since counts were identical. 
Alternative: sum rows, rejected because it would double their expression.
- **Cast counts to int64:** Stripped commas and converted floats to int. Task specifies integer counts.
- **Index by gene_id and keep symbol separate:** Stored gene symbols in a dictionary to keep count matrix purely numeric for NumPy.

---

## Results

**Top 5 Up-regulated genes:**
- **TXNRD1** (`ENSG90000000062`) - +2.6935
- **ATF4** (`ENSG90000000054`) - +2.6856
- **DDIT3** (`ENSG90000000055`) - +2.6284
- **DNAJB1** (`ENSG90000000053`) - +2.5091
- **MT2A** (`ENSG90000000068`) - +2.2952

**Top 5 Down-regulated genes:**
- **BUB1** (`ENSG90000000047`) - -3.8366
- **CCNB1** (`ENSG90000000039`) - -3.4751
- **TOP2A** (`ENSG90000000038`) - -3.2158
- **MKI67** (`ENSG90000000036`) - -2.7491
- **AURKA** (`ENSG90000000045`) - -2.5770

---

## Suspicious samples

- **Sample S6 is suspicious:**
  - Low RIN of 3.1 (normal samples are 8.2–9.8).
  - Low library size: only 1,754 total reads (~50x lower than ~64k–117k in others).
  - Low correlation with other treated samples (r is approx. 0.84 vs r > 0.98 among normal replicates).
  - Shallow depth caused zero-count dropouts for moderate genes like `CD44` (0 reads in S6 vs ~45 in other treated), which falsely exaggerated `CD44` downregulation (-3.12).
- **What I did:** Excluded S6 when computing treated means so technical dropouts do not skew the results. Without S6, `CD44` log2FC is -0.94 and the top downregulated genes consistently reflect cell cycle arrest.
- **Sample S9:** RIN 2.4, failed QC and was not sequenced; excluded from the sample sheet.

---

## With more time I would...

- Make PCA and volcano plots to visualize sample clustering and significant hits.

---

## Tools used (docs, AI assistants, ...)
