import csv, math, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEAT = ROOT / "build" / "heat"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

rows = []
for N in [16, 32, 64, 128, 256]:
    out = subprocess.run(
        [str(HEAT), "--case", "V", "--N", str(N), "--r", "0.2", "--T", "0.05"],
        capture_output=True, text=True, check=True).stdout
    d = dict(t.split("=", 1) for t in out.split())
    rows.append((N, int(d["steps"]), float(d["err_inf"]), float(d["err_l2_h"])))

with open(OUT / "convergence.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["N", "h", "steps", "err_inf", "err_l2_h", "p_obs"])
    for k, (N, steps, einf, el2) in enumerate(rows):
        p = "" if k == 0 else f"{math.log2(rows[k-1][3] / el2):.4f}"
        w.writerow([N, 1.0 / N, steps, einf, el2, p])

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    hs = [1.0 / r[0] for r in rows]
    plt.loglog(hs, [r[3] for r in rows], "o-", label="||e||_2,h")
    plt.loglog(hs, [r[2] for r in rows], "s-", label="||e||_inf")
    plt.loglog(hs, [rows[0][3] * (h / hs[0]) ** 2 for h in hs], "k--", label="slope 2")
    plt.xlabel("h"); plt.ylabel("error"); plt.legend(); plt.grid(True, which="both")
    plt.savefig(OUT / "convergence.png", dpi=150)
except ImportError:
    print("matplotlib missing: CSV written, plot skipped", file=sys.stderr)