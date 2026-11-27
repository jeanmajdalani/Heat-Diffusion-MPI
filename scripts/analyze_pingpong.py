#!/usr/bin/env python3
"""Analyze MPI ping-pong CSVs and create model summary plus plots."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT = RESULTS / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "Same VM - SHM": "pingpong_same_vm_shm.csv",
    "Same VM - TCP": "pingpong_same_vm_tcp.csv",
    "Node 1 - Node 2": "pingpong_node1_node2.csv",
    "Node 1 - Jean": "pingpong_node1_jean.csv",
    "Node 2 - Jean": "pingpong_node2_jean.csv",
}
REQUIRED = {"bytes", "reps", "t_oneway_s", "bandwidth_MBps"}

def load_data(label, filename):
    path = RESULTS / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    df = pd.read_csv(path)
    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"{filename} missing columns: {sorted(missing)}")
    for col in REQUIRED:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(REQUIRED))
    df = df[(df["bytes"] > 0) & (df["reps"] > 0) & (df["t_oneway_s"] > 0)]
    df = df.sort_values("bytes").reset_index(drop=True)
    if len(df) < 2:
        raise ValueError(f"{filename} has fewer than two valid rows")
    df["dataset"] = label
    df["message_size_MB"] = df["bytes"] / 1_000_000
    return df

def fit_model(df):
    x = df["bytes"].to_numpy(dtype=float)
    y = df["t_oneway_s"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    ss_res = np.sum((y - predicted) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    model_bw = (1 / slope) / 1_000_000 if slope > 0 else np.nan
    return {
        "dataset": df["dataset"].iloc[0],
        "measurements": len(df),
        "repetition_counts": ", ".join(map(str, sorted(df["reps"].astype(int).unique()))),
        "startup_latency_s": intercept,
        "startup_latency_ms": intercept * 1000,
        "transfer_time_s_per_byte": slope,
        "model_bandwidth_MBps": model_bw,
        "r_squared": r2,
    }

def main():
    frames, summaries = [], []
    for label, filename in DATASETS.items():
        df = load_data(label, filename)
        frames.append(df)
        summaries.append(fit_model(df))
        reps = sorted(df["reps"].astype(int).unique())
        if len(reps) > 1:
            print(f"NOTE: {label} has mixed repetition counts: {reps}.")
    combined = pd.concat(frames, ignore_index=True)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "pingpong_model_summary.csv", index=False)
    combined.to_csv(OUT / "pingpong_measurements_combined.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    for label, group in combined.groupby("dataset", sort=False):
        group = group.sort_values("bytes")
        ax.plot(group["message_size_MB"], group["t_oneway_s"], marker="o", label=label)
        row = summary.loc[summary["dataset"] == label].iloc[0]
        xs = np.linspace(group["message_size_MB"].min(), group["message_size_MB"].max(), 100)
        ys = row["startup_latency_s"] + xs * 1_000_000 * row["transfer_time_s_per_byte"]
        ax.plot(xs, ys, linestyle="--", linewidth=1)
    ax.set_xlabel("Message size (MB, decimal)")
    ax.set_ylabel("One-way communication time (s)")
    ax.set_title("MPI Ping-Pong: measured time and linear model")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "pingpong_time_model.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    for label, group in combined.groupby("dataset", sort=False):
        group = group.sort_values("bytes")
        ax.plot(group["message_size_MB"], group["bandwidth_MBps"], marker="o", label=label)
    ax.set_xlabel("Message size (MB, decimal)")
    ax.set_ylabel("Measured bandwidth (MB/s)")
    ax.set_title("MPI Ping-Pong: measured bandwidth comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "pingpong_bandwidth_comparison.png", dpi=180)
    plt.close(fig)

    print("\nMODEL SUMMARY")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nCreated in results/analysis/:")
    print("  pingpong_model_summary.csv")
    print("  pingpong_measurements_combined.csv")
    print("  pingpong_time_model.png")
    print("  pingpong_bandwidth_comparison.png")
    print("\nReview R-squared and plots; a single linear fit may not describe all sizes well.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
