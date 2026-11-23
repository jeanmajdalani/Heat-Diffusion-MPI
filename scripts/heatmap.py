#!/usr/bin/env python3
"""Heat map of a solver output file ((N+1)^2 little-endian doubles, row by row).
Usage: python3 scripts/heatmap.py file.bin --N 128 [--title "text"] [--out results/heatmap.png]
Several files: python3 scripts/heatmap.py a.bin b.bin --N 128 --title "A;B" --out results/heatmaps.png"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("files", nargs="+")
ap.add_argument("--N", type=int, required=True, help="grid size N (file holds (N+1)^2 values)")
ap.add_argument("--title", default=None, help="one title per file, separated by ';'")
ap.add_argument("--out", default="results/heatmap.png")
a = ap.parse_args()

titles = a.title.split(";") if a.title else [os.path.basename(f) for f in a.files]
fig, axes = plt.subplots(1, len(a.files), figsize=(5.2 * len(a.files), 4.6), squeeze=False)
for ax, f, t in zip(axes[0], a.files, titles):
    u = np.fromfile(f, dtype="<f8")
    if u.size != (a.N + 1) ** 2:
        sys.exit(f"{f}: {u.size} values, expected {(a.N + 1) ** 2} for N={a.N}")
    u = u.reshape(a.N + 1, a.N + 1)  # [j, i]
    im = ax.imshow(u, origin="lower", cmap="inferno", extent=(0, 1, 0, 1),
                   vmin=min(0, u.min()), vmax=max(1, u.max()))
    ax.set_title(t); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, fraction=0.046, label="temperature u")
fig.tight_layout()
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, dpi=150)
print("wrote", a.out)
