import numpy as np
import pandas as pd


def log_cpm(counts: np.ndarray) -> np.ndarray:
    """genes x samples counts -> log2(counts per million + 1)"""
    library_sizes = counts.sum(axis=0, keepdims=True)
    cpm = (counts / library_sizes) * 1e6
    return np.log2(cpm + 1.0)


if __name__ == "__main__":
    # Load clean counts and normalize
    df_counts = pd.read_csv("output/counts_clean.csv", index_col=0)
    log_cpm_mat = log_cpm(df_counts.values)

    # Save normalized table
    df_log_cpm = pd.DataFrame(
        np.round(log_cpm_mat, 4),
        index=df_counts.index,
        columns=df_counts.columns
    )
    df_log_cpm.to_csv("output/counts_log_cpm.csv")

    # Quick sanity check: gene with 1 count in 1M reads = log2(1 + 1) = 1.0
    test = log_cpm(np.array([[1], [999999]]))
    assert np.isclose(test[0, 0], 1.0)

    print(f"Done: normalized {df_log_cpm.shape[0]} genes x {df_log_cpm.shape[1]} samples -> output/counts_log_cpm.csv")
