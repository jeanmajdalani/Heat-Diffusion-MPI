#!/usr/bin/env python3
"""Weak-scaling analysis for the Heat MPI project.

Weak scaling: the work PER PROCESS is kept (almost) constant, so N grows with P
(here (N-1)^2 / P ~ 1.05e6 interior nodes per process).  Ideally the run time
stays constant.

Reads the CSV files written by scripts/benchmark.sh
(columns: decomp,map_by,N,P,rep,steps,r,t_total,t_comp,t_comm),
one file per (decomposition, P), takes the MEDIAN of the repetitions and
computes

    weak efficiency  E_w(P) = T(1) / T(P)       (1.0 is ideal)
    comm share       = median over runs of t_comm / t_total

Usage (from the repository root):
    python3 scripts/weak.py 'results/weak_jean_*.csv' \
        --label "Jean VM (WSL2, 8 physical cores)" --out results/weak_jean

Writes <out>.png (figure) and <out>_summary.csv (the numbers behind it).
"""
import argparse
import csv
import glob
import os
import statistics as st
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scaling import COLORS, INK, INK_SOFT, GRID, STYLE, place_labels  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def read_rows(pattern):
    files = sorted(f for f in glob.glob(pattern)
                   if not f.endswith((".log", "_summary.csv")))
    if not files:
        raise SystemExit(f"No file matches {pattern!r}")
    rows = []
    for f in files:
        with open(f, newline="") as fh:
            lines = [ln for ln in fh if not ln.startswith("#")]
        for r in csv.DictReader(lines):
            rows.append(dict(decomp=r["decomp"], N=int(r["N"]), P=int(r["P"]),
                             t=float(r["t_total"]), comm=float(r["t_comm"])))
    return files, rows


def summarise(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["decomp"], r["P"])].append(r)
    out = {}
    for (d, P), rs in groups.items():
        Ns = {r["N"] for r in rs}
        if len(Ns) != 1:
            raise SystemExit(f"{d} P={P}: several N values {sorted(Ns)} - "
                             "weak scaling needs one N per P")
        N = Ns.pop()
        ts = [r["t"] for r in rs]
        out[(d, P)] = dict(
            N=N, n=len(rs), t_med=st.median(ts), t_min=min(ts), t_max=max(ts),
            nodes_pp=(N - 1) ** 2 / P,
            comm=st.median(r["comm"] / r["t"] for r in rs if r["t"] > 0))
    for (d, P), s in out.items():
        t1 = out.get((d, 1))
        s["eff"] = t1["t_med"] / s["t_med"] if t1 else float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern", help="glob of benchmark CSV files, quoted")
    ap.add_argument("--out", default="results/weak_scaling")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    files, rows = read_rows(a.pattern)
    S = summarise(rows)
    decomps = sorted({k[0] for k in S}, key=["strips", "blocks"].index)
    Ps = sorted({k[1] for k in S})

    with open(a.out + "_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["decomp", "P", "N", "nodes_per_process", "reps",
                    "t_median_s", "t_min_s", "t_max_s",
                    "weak_efficiency", "comm_share"])
        for (d, P) in sorted(S, key=lambda k: (k[0], k[1])):
            s = S[(d, P)]
            w.writerow([d, P, s["N"], f"{s['nodes_pp']:.0f}", s["n"],
                        f"{s['t_med']:.4f}", f"{s['t_min']:.4f}",
                        f"{s['t_max']:.4f}", f"{s['eff']:.3f}",
                        f"{s['comm']:.3f}"])
    print(f"read {len(rows)} runs from {len(files)} file(s)")
    print(f"{'decomp':7} {'P':>3} {'N':>5} {'nodes/proc':>11} "
          f"{'T med (s)':>10} {'min-max':>15} {'E_w':>5} {'comm':>6}")
    for (d, P) in sorted(S, key=lambda k: (k[0], k[1])):
        s = S[(d, P)]
        print(f"{d:7} {P:3d} {s['N']:5d} {s['nodes_pp']:11.0f} "
              f"{s['t_med']:10.3f} {s['t_min']:7.3f}-{s['t_max']:<7.3f} "
              f"{s['eff']:5.2f} {100*s['comm']:5.1f}%")

    color = {d: COLORS[i] for i, d in enumerate(decomps)}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    panels = [("t_med", "Run time (500 steps), work per process fixed", (0, None)),
              ("eff", "Weak efficiency  E_w(P) = T(1) / T(P)", (0, 1.15)),
              ("comm", "Communication share of run time", (0, None))]
    xlabels = []
    for P in Ps:
        Ns = sorted({S[k]["N"] for k in S if k[1] == P})
        xlabels.append(f"{P}\nN={'/'.join(str(n) for n in Ns)}")

    for ax, (key, title, ylim) in zip(axes, panels):
        ax.set_title(title, fontsize=11, color=INK, loc="left")
        ax.set_xscale("log", base=2)
        ax.set_xticks(Ps)
        ax.set_xticklabels(xlabels, fontsize=9)
        ax.minorticks_off()
        ax.set_xlabel("processes P", color=INK_SOFT)
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color("#b3b3b3")
        ax.tick_params(colors=INK_SOFT)

        if key == "eff":
            ax.axhline(1.0, color=INK_SOFT, lw=1.2, ls=":", zorder=1)
        ends = []
        for d in decomps:
            pts = [(P, S[(d, P)]) for P in Ps if (d, P) in S]
            if not pts:
                continue
            xs = [p for p, _ in pts]
            ys = [(100 * s[key] if key == "comm" else s[key]) for _, s in pts]
            if key == "t_med" and (d, 1) in S:      # ideal = flat at T(1)
                ax.axhline(S[(d, 1)]["t_med"], color=INK_SOFT, lw=1.0,
                           ls=":", zorder=1)
            st_ = STYLE[d]
            ax.plot(xs, ys, color=color[d], lw=2, ms=7, ls=st_["ls"],
                    marker=st_["marker"], markeredgecolor="white",
                    markeredgewidth=1.2, zorder=3)
            ends.append((ys[-1], xs[-1], d.capitalize()))
        if key == "comm":
            ax.set_ylabel("% of run time", color=INK_SOFT)
        if key == "t_med":
            ax.set_ylabel("seconds", color=INK_SOFT)
        if ylim:
            ax.set_ylim(*ylim)
        ax.margins(x=0.14)
        place_labels(fig, ax, ends)

    handles = [Line2D([0], [0], color=color[d], lw=2, ls=STYLE[d]["ls"],
                      marker=STYLE[d]["marker"], markeredgecolor="white", ms=7,
                      label=d.capitalize()) for d in decomps]
    handles.append(Line2D([0], [0], color=INK_SOFT, lw=1.2, ls=":",
                          label="ideal (flat / 1.0)"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=10, labelcolor=INK)
    sup = "Weak scaling" + (f" - {a.label}" if a.label else "")
    fig.suptitle(sup, x=0.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(a.out + ".png", dpi=170, facecolor="white")
    print(f"wrote {a.out}.png and {a.out}_summary.csv")


if __name__ == "__main__":
    main()
