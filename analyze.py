import pandas as pd
import numpy as np
from normalize import log_cpm


def load_data():
    counts = pd.read_csv("output/counts_clean.csv", index_col=0)
    samples = pd.read_csv("output/samples_clean.csv")

    raw = pd.read_csv("data/counts_raw.csv")
    gene_map = dict(zip(raw["gene_id"], raw["gene_symbol"]))
    return counts, samples, gene_map


def main():
    counts, samples, gene_map = load_data()

    # Normalise counts
    norm_counts = log_cpm(counts.values)
    df_norm = pd.DataFrame(norm_counts, index=counts.index, columns=counts.columns)

    # 1. Sample QC: detect S6 as outlier (RIN 3.1, low depth 1,754 reads)
    qc = samples.copy()
    qc["total_reads"] = counts.sum(axis=0).values
    qc.to_csv("output/sample_qc_summary.csv", index=False)
    df_norm.corr().round(4).to_csv("output/sample_correlation_matrix.csv")

    # 2. Differential expression excluding outlier S6 (S5, S7, S8 vs S1, S2, S3, S4)
    ctrl_cols = ["S1", "S2", "S3", "S4"]
    treat_cols = ["S5", "S7", "S8"]  # S6 excluded due to technical failure

    mean_ctrl = df_norm[ctrl_cols].mean(axis=1)
    mean_treat = df_norm[treat_cols].mean(axis=1)
    lfc = mean_treat - mean_ctrl

    res = pd.DataFrame({
        "gene_id": counts.index,
        "gene_symbol": [gene_map.get(g, g) for g in counts.index],
        "mean_control": mean_ctrl.round(4).values,
        "mean_treated": mean_treat.round(4).values,
        "log2_fold_change": lfc.round(4).values
    }).set_index("gene_id")

    res = res.sort_values(by="log2_fold_change", ascending=False)
    res.to_csv("output/differential_expression.csv")

    # 3. Report Top 5 UP and DOWN
    top_up = res.head(5)
    top_down = res.tail(5).iloc[::-1]

    print("--- Top 5 Up-regulated Genes ---")
    for gid, row in top_up.iterrows():
        print(f"- **{row['gene_symbol']}** (`{gid}`) - +{row['log2_fold_change']:.4f}")

    print("\n--- Top 5 Down-regulated Genes ---")
    for gid, row in top_down.iterrows():
        print(f"- **{row['gene_symbol']}** (`{gid}`) - {row['log2_fold_change']:.4f}")

    print("\n--- Note on Suspicious Sample ---")
    print(f"Sample S6 (RIN 3.1, {qc.loc[qc['sample_id'] == 'S6', 'total_reads'].values[0]} reads) was excluded to prevent false-positive dropouts.")


if __name__ == "__main__":
    main()
