#!/usr/bin/env python3
"""Figure pour la Task 6.7 : instabilite a r = 1.05 * 0.25 (r > 1/4).
Lit results/unstable/ (cree par tools/unstable_demo). Sortie : results/unstable/instability.png"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 64
outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/unstable")
r_bad, r_ok = "0.2625", "0.24"
snaps = [100, 130, 160, 200]

def load(r, step):
    a = np.fromfile(outdir / f"snap_r{r}_step{step}.bin", dtype="<f8")
    return a.reshape(N + 1, N + 1)          # [j, i]

fig, ax = plt.subplots(2, 3, figsize=(13, 8))
for k, s in enumerate(snaps[:3]):
    u = load(r_bad, s)
    m = np.abs(u[1:-1, 1:-1]).max()
    im = ax.flat[k].imshow(u, origin="lower", cmap="RdBu_r", vmin=-m, vmax=m)
    ax.flat[k].set_title(f"r = {r_bad}, pas {s}  (max|u| = {m:.2g})")
    fig.colorbar(im, ax=ax.flat[k], fraction=0.046)

# zoom sur un coin au dernier pas : on voit l'alternance +/-
u = load(r_bad, snaps[-1]); z = u[1:17, 1:17]; m = np.abs(z).max()
im = ax[1, 0].imshow(z, origin="lower", cmap="RdBu_r", vmin=-m, vmax=m, extent=(0.5, 16.5, 0.5, 16.5))
ax[1, 0].set_title(f"Zoom 16x16 pres du coin, pas {snaps[-1]}\n(signe alterne d'un noeud a l'autre)")
ax[1, 0].set_xlabel("i"); ax[1, 0].set_ylabel("j")

# coupe 1D le long d'une ligne
row = u[N // 2, :]
ax[1, 1].plot(range(N + 1), row, marker="o", ms=3, lw=0.8)
ax[1, 1].axhline(0, color="k", lw=0.5)
ax[1, 1].set_title(f"Coupe j = {N//2}, pas {snaps[-1]}")
ax[1, 1].set_xlabel("i"); ax[1, 1].set_ylabel("u")

# croissance de max|u| pour r stable et instable
for r, c, lab in [(r_ok, "tab:green", "r = 0.24 (stable)"), (r_bad, "tab:red", "r = 0.2625 (1.05 x limite)")]:
    d = np.loadtxt(outdir / f"maxabs_r{r}.csv", delimiter=",", skiprows=1)
    ax[1, 2].semilogy(d[:, 0], np.maximum(d[:, 2], 1e-16), c, label=lab)
n = np.arange(0, 201)
G = 1 - 8 * float(r_bad)
n0 = 110  # repere : meme pente que |G|^n, recale pour passer par 1 au pas n0
ax[1, 2].semilogy(n, abs(G) ** (n - n0), "k--", lw=0.8, label=f"pente |G|^n, G = 1-8r = {G:.3f}")
ax[1, 2].set_ylim(1e-3, 1e6)
ax[1, 2].set_title("max|u| a l'interieur")
ax[1, 2].set_xlabel("pas de temps"); ax[1, 2].legend(fontsize=8)
fig.suptitle(f"Task 6.7 : schema explicite 2D, N = {N}, cas 1, r > 1/4 : instabilite en damier")
fig.tight_layout()
out = outdir / "instability.png"
fig.savefig(out, dpi=130)
print("ecrit", out)
