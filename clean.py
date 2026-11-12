import os
import re
import pandas as pd


def clean_samples(raw_path: str, out_path: str) -> pd.DataFrame:
    df = pd.read_csv(raw_path)

    # Standardize sample IDs (S02, s3, Sample 4 -> S1..S8)
    df["sample_id"] = df["sample_id"].apply(
        lambda x: f"S{int(re.search(r'\d+', str(x)).group())}"
    )

    # Standardize condition labels
    cond_map = {
        "control": "control", "Control": "control", "ctrl": "control", "untreated": "control",
        "treated": "treated", "Treated": "treated", "drug": "treated"
    }
    df["condition"] = df["condition"].map(cond_map)

    # Drop duplicate row (s3) and unsequenced sample (S9, degraded RIN 2.4)
    df = df.drop_duplicates(subset=["sample_id"])
    df = df[df["sample_id"].isin([f"S{i}" for i in range(1, 9)])].copy()

    df["rin"] = df["rin"].astype(float)
    df = df.sort_values(by="sample_id", key=lambda col: col.str.extract(r'(\d+)')[0].astype(int))
    df = df.reset_index(drop=True)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)
    return df


def clean_counts(raw_path: str, out_path: str, sample_order: list[str]) -> pd.DataFrame:
    df = pd.read_csv(raw_path)

    # Drop HTSeq alignment summary rows
    df = df[~df["gene_id"].str.startswith("__")].copy()

    # Drop duplicate gene rows (CDH2, TUBB)
    df = df.drop_duplicates(subset=["gene_id"]).copy()

    # Rename BAM paths to sample IDs (/scratch/run1/S2.sorted.bam -> S2)
    col_map = {c: re.search(r"(S\d+)", c).group(1) for c in df.columns if re.search(r"(S\d+)", c)}
    df = df.rename(columns=col_map)

    # Strip commas and cast to integer counts
    for col in sample_order:
        df[col] = df[col].astype(str).str.replace(",", "").astype(float).astype(int)

    # Gene IDs as index, matching sample sheet columns
    counts = df.set_index("gene_id")[sample_order]

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    counts.to_csv(out_path)
    return counts


if __name__ == "__main__":
    samples = clean_samples("data/samples_raw.csv", "output/samples_clean.csv")
    counts = clean_counts("data/counts_raw.csv", "output/counts_clean.csv", samples["sample_id"].tolist())
    print(f"Done: {len(samples)} samples, {len(counts)} genes cleaned.")
