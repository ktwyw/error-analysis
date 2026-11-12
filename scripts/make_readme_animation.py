"""Animation at the top of README.md: do error bars cover the truth as often as they claim?

An experiment (the mean of 25 readings, reported as mean ± standard error) is repeated 2000
times.  Left: the latest results, with the error bars that miss the true value in red.
Right: the distribution of the results against the prediction of the standard error.
Output: docs/images/coverage.gif (about 1 MB, 80 frames).

Run from the repository root:  python scripts/make_readme_animation.py   (needs Pillow)
"""
import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

try:
    import errprop as ep
except ImportError:          # running from a clone without `pip install -e .`
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import errprop as ep

ep.use_style()
rng = np.random.default_rng(68)
truth, sigma, n_read, n_exp, n_show = 9.81, 0.25, 25, 2000, 30
data = rng.normal(truth, sigma, (n_exp, n_read))
mean = data.mean(axis=1)
sem = data.std(axis=1, ddof=1) / np.sqrt(n_read)
hit = np.abs(mean - truth) < sem
sem_true = sigma / np.sqrt(n_read)

counts = np.unique(np.round(np.r_[np.arange(1, n_show + 1), np.geomspace(n_show, n_exp, 46)]).astype(int))
bins = np.linspace(truth - 4 * sem_true, truth + 4 * sem_true, 41)
xg = np.linspace(bins[0], bins[-1], 300)
gauss = np.exp(-0.5 * ((xg - truth) / sem_true)**2) / (sem_true * np.sqrt(2 * np.pi))

frames = []
for n in counts:
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.9))
    idx = np.arange(max(0, n - n_show), n)
    for row, i in enumerate(idx):
        c = "C0" if hit[i] else "C3"
        a1.errorbar(mean[i], row, xerr=sem[i], fmt="o", ms=3.5, color=c, lw=1.3, capsize=2.5)
    a1.axvline(truth, color="k", ls=":", lw=1)
    a1.set_xlim(bins[0], bins[-1]); a1.set_ylim(-1.5, n_show + 0.5); a1.set_yticks([])
    a1.set_xlabel(r"result of one experiment, $\bar x \pm s/\sqrt{N}$")
    a1.set_ylabel(f"latest {n_show} experiments")
    a1.set_title(f"experiment {n} of {n_exp}", fontsize=10)

    a2.hist(mean[:n], bins, density=True, color="C0", alpha=0.55)
    a2.plot(xg, gauss, "k", lw=1.2, label="predicted by the standard error")
    a2.axvline(truth, color="k", ls=":", lw=1)
    a2.set_xlim(bins[0], bins[-1]); a2.set_ylim(0, 1.45 * gauss.max())
    a2.set_xlabel("result"); a2.set_ylabel("probability density")
    a2.legend(loc="upper left", fontsize=8.5)
    a2.set_title(f"error bars that contain the true value: {100 * hit[:n].mean():.1f} %   (claimed: 68 %)", fontsize=10)

    buf = io.BytesIO(); fig.savefig(buf, dpi=95); plt.close(fig)
    frames.append(Image.open(buf).convert("RGB"))

frames += [frames[-1]] * 25                                   # hold the final frame
pal = frames[-1].quantize(colors=64)
frames = [f.quantize(colors=64, palette=pal, dither=Image.Dither.NONE) for f in frames]
frames[0].save("docs/images/coverage.gif", save_all=True, append_images=frames[1:], duration=110, loop=0, optimize=True)
import os
print(f"wrote docs/images/coverage.gif, {len(frames)} frames, {os.path.getsize('docs/images/coverage.gif') / 1e6:.2f} MB;"
      f" final coverage {hit.mean():.3f}")
