#!/usr/bin/env python3
"""Strong-scaling analysis for the Heat MPI project.

Reads the CSV files written by scripts/benchmark.sh
(columns: decomp,map_by,N,P,rep,steps,r,t_total,t_comp,t_comm),
takes the MEDIAN of the repetitions of every (decomp, N, P), and computes

    speed-up    S(P) = T(1) / T(P)
    efficiency  E(P) = S(P) / P
    comm share  = median over runs of t_comm / t_total

T(1) is the P=1 median of the SAME decomposition and the same N.

Usage (from the repository root):
    python3 scripts/scaling.py 'results/strong_jean_*.csv' \
        --label "Jean VM (WSL2, 16 vCPU)" --out results/strong_jean

Writes <out>.png (figure) and <out>_summary.csv (the numbers behind it).
"""
import argparse
import csv
import glob
import statistics as st
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Categorical slots 1-3 of the validated palette (all-pairs CVD safe).
# Three of these sit below 3:1 contrast on white, so every line is also
# labelled directly at its right end, and a summary CSV is written.
COLORS = ["#257de0", "#f47230", "#00b894"]
INK = "#191919"        # text
INK_SOFT = "#666666"   # secondary text, ideal line
GRID = "#e6e6e6"
STYLE = {  # secondary encoding: line style + marker, so colour is never alone
    "strips": dict(ls="-", marker="o"),
    "blocks": dict(ls="--", marker="s"),
}


def place_labels(fig, ax, ends, gap=0.075):
    """Direct labels at the right of the axes, pushed apart so they never overlap."""
    if not ends:
        return
    fig.canvas.draw()
    to_axes = ax.transAxes.inverted()
    items = []
    for y, P, text in ends:
        px = ax.transData.transform((P, y))
        items.append([to_axes.transform(px)[1], to_axes.transform(px)[0], text])
    items.sort(key=lambda t: t[0])
    orig = [it[0] for it in items]
    for i in range(1, len(items)):          # push apart (upwards)
        items[i][0] = max(items[i][0], items[i - 1][0] + gap)
    # re-centre the block on the true positions, then keep it inside the axes
    shift = sum(it[0] for it in items) / len(items) - sum(orig) / len(orig)
    for it in items:
        it[0] -= shift
    low = min(it[0] for it in items)
    if low < 0.03:
        for it in items:
            it[0] += 0.03 - low
    for y_ax, x_ax, text in items:
        ax.annotate(text, xy=(x_ax, y_ax), xycoords="axes fraction",
                    xytext=(8, 0), textcoords="offset points",
                    fontsize=9, color=INK, va="center", annotation_clip=False)


def read_rows(pattern):
    rows = []
    files = sorted(f for f in glob.glob(pattern)
                   if not f.endswith((".log", "_summary.csv")))
    if not files:
        raise SystemExit(f"No file matches {pattern!r}")
    for f in files:
        with open(f, newline="") as fh:
            lines = [ln for ln in fh if not ln.startswith("#")]
        for r in csv.DictReader(lines):
            rows.append(dict(
                decomp=r["decomp"], N=int(r["N"]), P=int(r["P"]),
                t=float(r["t_total"]), comm=float(r["t_comm"]),
            ))
    return files, rows


def summarise(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["decomp"], r["N"], r["P"])].append(r)
    out = {}
    for key, rs in groups.items():
        ts = [r["t"] for r in rs]
        out[key] = dict(
            n=len(rs), t_med=st.median(ts), t_min=min(ts), t_max=max(ts),
            comm=st.median(r["comm"] / r["t"] for r in rs if r["t"] > 0),
        )
    for (d, N, P), s in out.items():
        t1 = out.get((d, N, 1))
        s["speedup"] = t1["t_med"] / s["t_med"] if t1 else float("nan")
        s["eff"] = s["speedup"] / P
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern", help="glob of benchmark CSV files, quoted")
    ap.add_argument("--out", default="results/strong_scaling")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    files, rows = read_rows(a.pattern)
    S = summarise(rows)
    decomps = sorted({k[0] for k in S}, key=["strips", "blocks"].index)
    Ns = sorted({k[1] for k in S})
    Ps = sorted({k[2] for k in S})

    # ---- numbers (the "table view") -------------------------------------
    with open(a.out + "_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["decomp", "N", "P", "reps", "t_median_s", "t_min_s",
                    "t_max_s", "speedup", "efficiency", "comm_share"])
        for (d, N, P) in sorted(S, key=lambda k: (k[0], k[1], k[2])):
            s = S[(d, N, P)]
            w.writerow([d, N, P, s["n"], f"{s['t_med']:.4f}",
                        f"{s['t_min']:.4f}", f"{s['t_max']:.4f}",
                        f"{s['speedup']:.3f}", f"{s['eff']:.3f}",
                        f"{s['comm']:.3f}"])
    print(f"read {len(rows)} runs from {len(files)} file(s)")
    print(f"{'decomp':7} {'N':>5} {'P':>3} {'T med (s)':>10} "
          f"{'S':>6} {'E':>5} {'comm':>6}")
    for (d, N, P) in sorted(S, key=lambda k: (k[0], k[1], k[2])):
        s = S[(d, N, P)]
        print(f"{d:7} {N:5d} {P:3d} {s['t_med']:10.3f} "
              f"{s['speedup']:6.2f} {s['eff']:5.2f} {100*s['comm']:5.1f}%")

    # ---- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    panels = [
        ("speedup", "Speed-up  S(P) = T(1) / T(P)", None),
        ("eff", "Efficiency  E(P) = S(P) / P", (0, 1.15)),
        ("comm", "Communication share of run time", (0, None)),
    ]
    for ax, (key, title, ylim) in zip(axes, panels):
        ends = []
        ax.set_title(title, fontsize=11, color=INK, loc="left")
        ax.set_xscale("log", base=2)
        ax.set_xticks(Ps)
        ax.set_xticklabels([str(p) for p in Ps])
        ax.minorticks_off()
        ax.set_xlabel("processes P", color=INK_SOFT)
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color("#b3b3b3")
        ax.tick_params(colors=INK_SOFT)

        if key == "speedup":
            ax.plot(Ps, Ps, color=INK_SOFT, lw=1.2, ls=":", zorder=1)
            ax.text(Ps[-1], Ps[-1] * 1.06, "ideal", color=INK_SOFT,
                    fontsize=9, ha="right", va="bottom")
        if key == "eff":
            ax.axhline(1.0, color=INK_SOFT, lw=1.2, ls=":", zorder=1)

        for i, N in enumerate(Ns):
            for d in decomps:
                pts = [(P, S[(d, N, P)]) for P in Ps if (d, N, P) in S]
                if not pts:
                    continue
                xs = [p for p, _ in pts]
                ys = [(100 * s[key] if key == "comm" else s[key])
                      for _, s in pts]
                st_ = STYLE[d]
                ax.plot(xs, ys, color=COLORS[i % 3], lw=2, ms=7,
                        ls=st_["ls"], marker=st_["marker"],
                        markeredgecolor="white", markeredgewidth=1.2,
                        zorder=3)
            # remember the right end of the first curve for a direct label
            d0 = decomps[0]
            last = [(P, S[(d0, N, P)]) for P in Ps if (d0, N, P) in S]
            if last:
                P, s = last[-1]
                y = 100 * s[key] if key == "comm" else s[key]
                ends.append((y, P, f"N={N}"))
        if key == "comm":
            ax.set_ylabel("% of run time", color=INK_SOFT)
        if ylim:
            ax.set_ylim(*ylim)
        ax.margins(x=0.12)
        place_labels(fig, ax, ends)

    handles = [Line2D([0], [0], color=COLORS[i % 3], lw=2, label=f"N = {N}")
               for i, N in enumerate(Ns)]
    handles += [Line2D([0], [0], color=INK_SOFT, lw=2, marker=STYLE[d]["marker"],
                       ls=STYLE[d]["ls"], markeredgecolor="white", ms=7,
                       label=d.capitalize()) for d in decomps]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=10, labelcolor=INK)
    sup = "Strong scaling" + (f" - {a.label}" if a.label else "")
    fig.suptitle(sup, x=0.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(a.out + ".png", dpi=170, facecolor="white")
    print(f"wrote {a.out}.png and {a.out}_summary.csv")


if __name__ == "__main__":
    main()
