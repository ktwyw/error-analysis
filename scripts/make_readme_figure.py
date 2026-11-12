"""Draw the overview figure shown in README.md (docs/images/overview.png).

Four results from the book, each recomputed here in a few lines:
(a) Chapter 15, (b) Chapter 25, (c) Chapter 17, (d) Chapter 26.

Run from the repository root:  python scripts/make_readme_figure.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import lfilter

try:
    import errprop as ep
except ImportError:          # running from a clone without `pip install -e .`
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import errprop as ep

ep.use_style()
rng = np.random.default_rng(2026)
fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))

# (a) Chapter 15: a constant level fitted to a histogram of 20 bins
ax = axes[0, 0]
mus = np.array([2, 3, 5, 10, 20, 50, 100])
neyman, pearson, likelihood = [], [], []
for mu in mus:
    n = rng.poisson(mu, (20000, 20))
    neyman.append(np.mean(20 / np.sum(1 / np.maximum(n, 1), axis=1)) - mu)   # weights 1/n: harmonic mean
    pearson.append(np.mean(np.sqrt(np.mean(n**2, axis=1))) - mu)            # weights 1/mu: rms
    likelihood.append(np.mean(n.mean(axis=1)) - mu)                         # Poisson likelihood: mean
ax.semilogx(mus, neyman, "o-", color="C3", label=r"$\chi^2$, weights from data")
ax.semilogx(mus, pearson, "s-", color="C1", label=r"$\chi^2$, weights from model")
ax.semilogx(mus, likelihood, "^-", color="C0", label="Poisson likelihood")
ax.axhline(0, color="0.5", lw=0.6)
ax.set_xlabel("true counts per bin"); ax.set_ylabel("bias of fitted level (counts per bin)")
ax.set_ylim(-1.5, 1.6); ax.legend(loc="upper right", fontsize=8.5)
ax.set_title("Fitting a histogram", fontsize=10); ep.panel_label(ax, "(a)", x=0.05, y=0.1)

# (b) Chapter 25: probability of detection near the detection limit
ax = axes[0, 1]
a, b, sigma, m = 2.0, 1.5, 0.8, 3
xs = np.repeat([0.0, 1, 2, 4, 6, 8, 10], 3)
N, nu, xbar = xs.size, xs.size - 2, xs.mean()
Sxx = np.sum((xs - xbar)**2)
k0 = np.sqrt(1 / m + 1 / N + xbar**2 / Sxx)
t95 = stats.t.ppf(0.95, nu)
levels = np.linspace(0, 3, 13)
p_currie, p_rule = [], []
for x in levels:
    Y = a + b * xs + rng.normal(0, sigma, (20000, N))                       # a new calibration each time
    slope = ((xs - xbar) * (Y - Y.mean(axis=1, keepdims=True))).sum(axis=1) / Sxx
    icpt = Y.mean(axis=1) - slope * xbar
    s = np.sqrt(((Y - icpt[:, None] - slope[:, None] * xs)**2).sum(axis=1) / nu)
    x0 = (a + b * x + rng.normal(0, sigma / np.sqrt(m), slope.size) - icpt) / slope
    p_currie.append(np.mean(x0 > t95 * s / slope * k0))
    p_rule.append(np.mean(x0 > 3 * s / slope))
ax.plot(levels, p_currie, "o-", color="C0", label="Currie decision threshold")
ax.plot(levels, p_rule, "s-", color="C3", label='"3s" rule of thumb')
ax.axvline(3 * sigma / b, color="C3", ls=":", lw=1)
ax.axhline(0.95, color="0.5", lw=0.6)
ax.set_xlabel("true concentration (µg/L)"); ax.set_ylabel("probability of detection")
ax.set_ylim(0, 1.05); ax.legend(loc="lower right", fontsize=8.5)
ax.set_title("A detection limit", fontsize=10); ep.panel_label(ax, "(b)", x=0.05)

# (c) Chapter 17: the error bar of a correlated time series
ax = axes[1, 0]
phi, n_steps = 0.95, 2**16
series = lfilter([1.0], [1.0, -phi], rng.normal(0, 1, n_steps + 2000))[2000:]
sizes = 2**np.arange(0, 12)
sem = [ep.block_means(series, int(bs)).std(ddof=1) / np.sqrt(n_steps // bs) for bs in sizes]
sem_true = np.sqrt((1 + phi) / (1 - phi) / (1 - phi**2) / n_steps)
ax.semilogx(sizes, sem, "o-", color="C0", label="block-averaged estimate")
ax.axhline(sem[0], color="C3", ls="--", label=r"naive $s/\sqrt{N}$")
ax.axhline(sem_true, color="k", ls=":", label="true standard error")
ax.set_xlabel("block size (steps)"); ax.set_ylabel("standard error of the mean")
ax.set_ylim(0, 1.25 * max(max(sem), sem_true)); ax.legend(loc="center right", fontsize=8.5)
ax.set_title("Correlated simulation data", fontsize=10); ep.panel_label(ax, "(c)", x=0.05)

# (d) Chapter 26: tolerance stack-up of a six-part assembly
ax = axes[1, 1]
tol = np.array([0.10, 0.04, 0.05, 0.05, 0.04, 0.03])
sign = np.array([1, -1, -1, -1, -1, -1])
g0, wc, rss = 0.30, tol.sum(), np.sqrt(np.sum(tol**2))
bins = np.linspace(-0.05, 0.65, 141)
ax.hist(g0 + rng.normal(0, tol / 3, (400000, 6)) @ sign, bins, density=True, histtype="step", lw=1.4,
        color="C0", label="normal, centred")
ax.hist(g0 + rng.uniform(-tol, tol, (400000, 6)) @ sign, bins, density=True, histtype="step", lw=1.4,
        color="C1", label="uniform in tolerance")
ax.axvspan(0.10, 0.50, color="0.9", zorder=0)
ax.hlines(10.4, g0 - wc, g0 + wc, color="k", lw=4); ax.text(g0 + wc + 0.01, 10.4, "worst case", fontsize=8, va="center")
ax.hlines(9.7, g0 - rss, g0 + rss, color="C0", lw=4); ax.text(g0 + rss + 0.01, 9.7, "RSS", fontsize=8, va="center", color="C0")
ax.set_xlabel("assembly gap (mm)"); ax.set_ylabel("probability density (1/mm)")
ax.set_xlim(-0.05, 0.78); ax.set_ylim(0, 11.2); ax.legend(loc="center right", fontsize=8.5)
ax.set_title("A tolerance stack-up", fontsize=10); ep.panel_label(ax, "(d)", x=0.05, y=0.6)

fig.savefig("docs/images/overview.png", dpi=160)
print("wrote docs/images/overview.png")
