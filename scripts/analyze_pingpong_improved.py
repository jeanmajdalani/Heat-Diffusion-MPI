#!/usr/bin/env python3
"""
Improved analysis for the five MPI ping-pong datasets.

Inputs (in ../results relative to this script):
  pingpong_same_vm_shm.csv
  pingpong_same_vm_tcp.csv
  pingpong_node1_node2.csv
  pingpong_node1_jean.csv
  pingpong_node2_jean.csv

Outputs:
  results/analysis/pingpong_improved_summary.csv
  results/analysis/pingpong_measurements_combined.csv
  results/analysis/pingpong_latency_log.png
  results/analysis/pingpong_bandwidth_log.png
  results/analysis/pingpong_large_message_model.png

The small-message latency is reported as the median measured one-way time
for 16 B through 4 KiB. Large-message bandwidth is the median reported
bandwidth for messages >= 1 MiB.

A separate linear fit is made only over 64 KiB through 4 MiB. Its intercept
is explicitly labelled an extrapolated large-regime intercept, not a direct
measurement of the small-message startup latency. Fit quality (R^2) is
reported so poor fits are visible.
"""
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

SMALL_MIN_BYTES = 16
SMALL_MAX_BYTES = 4096
LARGE_MIN_BYTES = 1_048_576
FIT_MIN_BYTES = 65_536
FIT_MAX_BYTES = 4_194_304
REQUIRED = {"bytes", "reps", "t_oneway_s", "bandwidth_MBps"}


def load_dataset(label, filename):
    path = RESULTS / filename
    if not path.is_file():
        raise FileNotFoundError(f"Missing input CSV: {path}")
    df = pd.read_csv(path)
    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"{filename} missing columns: {sorted(missing)}")
    for col in REQUIRED:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(REQUIRED))
    df = df[(df["bytes"] > 0) & (df["reps"] > 0) &
            (df["t_oneway_s"] > 0)]
    df = df.sort_values("bytes").reset_index(drop=True)
    if df.empty:
        raise ValueError(f"No valid measurements in {filename}")
    df["dataset"] = label
    df["message_size_bytes"] = df["bytes"]
    return df


def linear_fit(group):
    x = group["bytes"].to_numpy(dtype=float)
    y = group["t_oneway_s"].to_numpy(dtype=float)
    if len(group) < 3 or np.unique(x).size < 2:
        return (np.nan, np.nan, np.nan, len(group))
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    bw = (1.0 / slope) / 1_000_000 if slope > 0 else np.nan
    return intercept, slope, r2, len(group), bw


def main():
    frames = []
    summary_rows = []

    for label, filename in DATASETS.items():
        df = load_dataset(label, filename)
        frames.append(df)

        small = df[(df["bytes"] >= SMALL_MIN_BYTES) &
                   (df["bytes"] <= SMALL_MAX_BYTES)]
        large = df[df["bytes"] >= LARGE_MIN_BYTES]
        fit_group = df[(df["bytes"] >= FIT_MIN_BYTES) &
                       (df["bytes"] <= FIT_MAX_BYTES)]

        if small.empty:
            raise ValueError(f"{label}: no measurements in small-message range")
        if large.empty:
            raise ValueError(f"{label}: no measurements >= 1 MiB")

        fit_result = linear_fit(fit_group)
        intercept, slope, r2, fit_n, fit_bw = fit_result

        summary_rows.append({
            "dataset": label,
            "measurement_count": len(df),
            "repetition_counts": ", ".join(
                str(int(v)) for v in sorted(df["reps"].unique())
            ),
            "small_latency_range_bytes": f"{SMALL_MIN_BYTES}-{SMALL_MAX_BYTES}",
            "small_latency_median_ms": small["t_oneway_s"].median() * 1000,
            "small_latency_min_ms": small["t_oneway_s"].min() * 1000,
            "small_latency_max_ms": small["t_oneway_s"].max() * 1000,
            "large_bandwidth_range_bytes": f">={LARGE_MIN_BYTES}",
            "large_bandwidth_median_MBps": large["bandwidth_MBps"].median(),
            "large_bandwidth_min_MBps": large["bandwidth_MBps"].min(),
            "large_bandwidth_max_MBps": large["bandwidth_MBps"].max(),
            "large_fit_range_bytes": f"{FIT_MIN_BYTES}-{FIT_MAX_BYTES}",
            "large_fit_points": fit_n,
            "large_fit_intercept_ms_extrapolated": intercept * 1000,
            "large_fit_transfer_ns_per_byte": slope * 1e9,
            "large_fit_bandwidth_MBps": fit_bw,
            "large_fit_r_squared": r2,
        })

        if len(np.unique(df["reps"])) > 1:
            print(f"NOTE: {label} has mixed repetition counts: "
                  f"{sorted(int(v) for v in df['reps'].unique())}")

    combined = pd.concat(frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "pingpong_improved_summary.csv", index=False)
    combined.to_csv(OUT / "pingpong_measurements_combined.csv", index=False)

    # Plot 1: measured one-way time, log-log axes so local and network
    # communication are simultaneously visible.
    fig, ax = plt.subplots(figsize=(10, 6))
    for label, group in combined.groupby("dataset", sort=False):
        group = group.sort_values("bytes")
        ax.plot(group["bytes"], group["t_oneway_s"],
                marker="o", markersize=4, linewidth=1.4, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Message size (bytes, log scale)")
    ax.set_ylabel("Measured one-way time (seconds, log scale)")
    ax.set_title("MPI Ping-Pong: one-way latency across message sizes")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "pingpong_latency_log.png", dpi=180)
    plt.close(fig)

    # Plot 2: measured bandwidth, log-log axes.
    fig, ax = plt.subplots(figsize=(10, 6))
    for label, group in combined.groupby("dataset", sort=False):
        group = group[group["bandwidth_MBps"] > 0].sort_values("bytes")
        ax.plot(group["bytes"], group["bandwidth_MBps"],
                marker="o", markersize=4, linewidth=1.4, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Message size (bytes, log scale)")
    ax.set_ylabel("Measured bandwidth (MB/s, log scale)")
    ax.set_title("MPI Ping-Pong: bandwidth across message sizes")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "pingpong_bandwidth_log.png", dpi=180)
    plt.close(fig)

    # Plot 3: measured data and separate large-message fit, 64 KiB-4 MiB.
    fig, ax = plt.subplots(figsize=(10, 6))
    for label, group in combined.groupby("dataset", sort=False):
        group = group.sort_values("bytes")
        ax.scatter(group["bytes"], group["t_oneway_s"],
                   s=24, label=f"{label} measured")
        fit_group = group[(group["bytes"] >= FIT_MIN_BYTES) &
                          (group["bytes"] <= FIT_MAX_BYTES)]
        if len(fit_group) >= 3:
            slope, intercept = np.polyfit(
                fit_group["bytes"].to_numpy(float),
                fit_group["t_oneway_s"].to_numpy(float), 1
            )
            xs = np.linspace(FIT_MIN_BYTES, FIT_MAX_BYTES, 100)
            ys = intercept + slope * xs
            ax.plot(xs, ys, linestyle="--", linewidth=1.2,
                    label=f"{label} fit (64 KiB-4 MiB)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Message size (bytes, log scale)")
    ax.set_ylabel("Measured one-way time (seconds, log scale)")
    ax.set_title("MPI Ping-Pong: large-message linear fits")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "pingpong_large_message_model.png", dpi=180)
    plt.close(fig)

    cols = [
        "dataset", "measurement_count", "repetition_counts",
        "small_latency_median_ms", "large_bandwidth_median_MBps",
        "large_fit_intercept_ms_extrapolated",
        "large_fit_transfer_ns_per_byte",
        "large_fit_bandwidth_MBps", "large_fit_r_squared",
    ]
    print("\nIMPROVED SUMMARY")
    print(summary[cols].to_string(index=False,
          float_format=lambda v: f"{v:.6g}"))
    print("\nCreated in results/analysis/:")
    for name in [
        "pingpong_improved_summary.csv",
        "pingpong_measurements_combined.csv",
        "pingpong_latency_log.png",
        "pingpong_bandwidth_log.png",
        "pingpong_large_message_model.png",
    ]:
        print(f"  {name}")
    print("\nIMPORTANT: the large-fit intercept is an extrapolation within the")
    print("large-message regime, not the directly measured small-message latency.")
    print("Use the small-message median as the measured latency metric.")
    print("Check R-squared and the plots before interpreting the model.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
