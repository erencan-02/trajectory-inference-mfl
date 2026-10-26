"""Figures for the handbook.

Run from the repository root:  python handbook/figures/make_figures.py
Creates the conceptual figures from scratch and extracts three result figures
from the executed notebook (trajectory_inference_mfl.ipynb).
"""
import base64
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "mathtext.fontset": "cm",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 9, "svg.fonttype": "path",
})


def lse(a, axis):
    m = a.max(axis=axis, keepdims=True)
    return np.squeeze(m + np.log(np.exp(a - m).sum(axis=axis, keepdims=True)), axis)


def entropic_plan(x, y, a, b, eps, n_iter=5000):
    # log-domain Sinkhorn for 1-D histograms, cost |x - y|^2 / 2
    C = 0.5 * (x[:, None] - y[None, :]) ** 2
    la, lb = np.log(a), np.log(b)
    f, g = np.zeros_like(a), np.zeros_like(b)
    for _ in range(n_iter):
        f = -eps * lse((g[None] - C) / eps + lb[None], 1)
        g = -eps * lse((f[:, None] - C) / eps + la[:, None], 0)
    return np.exp((f[:, None] + g[None] - C) / eps + la[:, None] + lb[None])


def mixture(x, comps):
    p = sum(w * np.exp(-(x - m) ** 2 / (2 * s ** 2)) / s for w, m, s in comps)
    return p / p.sum()


grid = np.linspace(-3, 3, 240)
mu = mixture(grid, [(0.5, -1.6, 0.25), (0.5, 0.4, 0.35)])
nu = mixture(grid, [(0.35, -0.6, 0.3), (0.65, 1.7, 0.25)])

# ---------------------------------------------------------------- plans
fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.2))
for ax, eps in zip(axes, [0.002, 0.05, 0.5]):
    P = entropic_plan(grid, grid, mu, nu, eps)
    ax.imshow(np.log(P.T + 1e-12), origin="lower", extent=[-3, 3, -3, 3], cmap="Blues",
              vmin=np.log(P.max()) - 8, vmax=np.log(P.max()))
    ax.set(title=rf"$\varepsilon = {eps}$", xlabel=r"$x \sim \mu$", xticks=[-2, 0, 2], yticks=[-2, 0, 2])
axes[0].set_ylabel(r"$y \sim \nu$")
fig.tight_layout()
fig.savefig(HERE / "eot_plans.svg")
plt.close(fig)

# ------------------------------------------- displacement vs entropic interpolation
ts = np.linspace(0, 1, 120)
cdf_mu, cdf_nu = np.cumsum(mu), np.cumsum(nu)
qs = np.linspace(0.001, 0.999, 4000)
Qmu, Qnu = np.interp(qs, cdf_mu, grid), np.interp(qs, cdf_nu, grid)


def displacement(t):
    # W2 geodesic in 1-D: interpolate quantile functions
    pts = (1 - t) * Qmu + t * Qnu
    h, _ = np.histogram(pts, bins=np.append(grid - (grid[1] - grid[0]) / 2, grid[-1] + (grid[1] - grid[0]) / 2))
    return h / h.sum()


def entropic(t, P, eps):
    # rho_t = int N((1-t)x + t y, eps t (1-t)) dpi(x, y): Brownian bridges glued to the static plan
    i, j = np.nonzero(P > 1e-10)
    w = P[i, j]
    m = (1 - t) * grid[i] + t * grid[j]
    s2 = max(eps * t * (1 - t), 1e-6)
    dens = (w[None] * np.exp(-(grid[:, None] - m[None]) ** 2 / (2 * s2))).sum(1)
    if t in (0.0, 1.0):
        return mu if t == 0 else nu
    return dens / dens.sum()


fig, axes = plt.subplots(1, 3, figsize=(6.6, 2.5), sharey=True)
panels = [("OT: displacement interpolation", None), (r"Schrödinger bridge, $\varepsilon = 0.05$", 0.05),
          (r"Schrödinger bridge, $\varepsilon = 0.5$", 0.5)]
rng = np.random.default_rng(0)
for ax, (title, eps) in zip(axes, panels):
    if eps is None:
        D = np.stack([displacement(t) for t in ts])
    else:
        P = entropic_plan(grid, grid, mu, nu, eps)
        D = np.stack([entropic(t, P, eps) for t in ts])
    ax.imshow(D.T, origin="lower", aspect="auto", extent=[0, 1, -3, 3], cmap="magma_r",
              vmax=np.percentile(D, 99.5))
    # a few sample trajectories
    for _ in range(7):
        k = rng.choice(len(qs))
        if eps is None:
            ax.plot(ts, (1 - ts) * Qmu[k] + ts * Qnu[k], color="C0", lw=0.8)
        else:
            flat = P.ravel() / P.sum()
            ij = rng.choice(flat.size, p=flat)
            x0, x1 = grid[ij // len(grid)], grid[ij % len(grid)]
            B = np.concatenate([[0], np.cumsum(rng.standard_normal(len(ts) - 1) * np.sqrt(eps * np.diff(ts)))])
            ax.plot(ts, (1 - ts) * x0 + ts * x1 + B - ts * B[-1], color="C0", lw=0.6)
    ax.set(title=title, xlabel="time $t$")
axes[0].set_ylabel("space $x$")
fig.tight_layout()
fig.savefig(HERE / "interpolations.svg")
plt.close(fig)

# ------------------------------------------- figures extracted from the notebook
nb = json.loads((ROOT / "trajectory_inference_mfl.ipynb").read_text(encoding="utf-8"))
wanted = {
    "mfl_dynamics.png": "snaps, trace = mfl.run(2500",   # cell after the run plots dynamics
    "paths.png": "tt, paths = mfl.sample_paths(150",
    "branching.png": 'for name, g in [("MFL + branching"',
    "few_samples.png": "obs1 = make_snapshots(",
}
cells = nb["cells"]
for name, marker in wanted.items():
    idx = next(i for i, c in enumerate(cells) if c["cell_type"] == "code" and marker in "".join(c["source"]))
    if name == "mfl_dynamics.png":        # the plot lives in the cell right after the run
        idx += 1
    pngs = [o["data"]["image/png"] for o in cells[idx].get("outputs", []) if "image/png" in o.get("data", {})]
    (HERE / name).write_bytes(base64.b64decode(pngs[0]))
print("figures written to", HERE)
