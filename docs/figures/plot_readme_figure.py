#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the two metric functions of the noncommutative dirty black hole.

    f(x) = 1 - erf(mu x)/x + 2 mu exp(-mu^2 x^2)/sqrt(pi)
    g(x) = (mu/2) [1 - erf(mu x) + 2 mu x exp(-mu^2 x^2)/sqrt(pi)]

as written in matrixassembler_nonext.mpl and matrixassembler_ext.mpl. The extremal mass is found
as matrixassembler_ext.mpl finds it, by solving f = 0 and df/dx = 0 together. From this folder,
    python3 plot_readme_figure.py
writes dirtybh-metric.png (Python 3 with NumPy, SciPy and Matplotlib). The style is the stamped
copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import brentq, fsolve
from scipy.special import erf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402


def f(x, mu):
    return 1 - erf(mu * x) / x + 2 * mu * np.exp(-(mu * x) ** 2) / np.sqrt(np.pi)


def dfdx(x, mu):
    return erf(mu * x) / x**2 - 2 * mu * np.exp(-(mu * x) ** 2) / (np.sqrt(np.pi) * x) \
        - 4 * mu**3 * x * np.exp(-(mu * x) ** 2) / np.sqrt(np.pi)


def g(x, mu):
    return mu / 2 * (1 - erf(mu * x) + 2 * mu * x * np.exp(-(mu * x) ** 2) / np.sqrt(np.pi))


xe, mu_e = fsolve(lambda v: [f(v[0], v[1]), dfdx(v[0], v[1])], (0.79, 1.90), xtol=1e-14)
print(f"extremal black hole: mu_e = {mu_e:.10f}, x_e = {xe:.10f}")
masses = [(1.3, "no horizon"), (mu_e, "extremal"), (2.3, "two horizons")]

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.75))
x = np.linspace(1e-3, 4.0, 4000)
for k, (mu, words) in enumerate(masses):
    style = hs.series(k)
    style["marker"] = None
    label = rf"$\mu = {mu:.5g}$ ({words})"
    left.plot(x, f(x, mu), lw=1.1, label=label, **style)
    right.plot(x, np.exp(-2 * g(x, mu)), lw=1.1, label=label, **style)
xh = brentq(lambda t: f(t, 2.3), 0.9, 1.1, xtol=1e-14)
print(f"mu = 2.3: event horizon x_h = {xh:.10f}")
for xr, k in ((xe, 1), (xh, 2)):
    left.plot([xr], [0.0], linestyle="none", marker="o", markersize=4.2,
              markerfacecolor="white", markeredgecolor=hs.SERIES[k], markeredgewidth=0.9)
left.plot(x[x > 0.5], 1 - 1 / x[x > 0.5], color=hs.SLATE, lw=0.8, linestyle=(0, (5.0, 1.5, 1.5, 1.5)),
          label="Schwarzschild, $1 - 1/x$")
left.axhline(0.0, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
left.set_xlim(0, 4.0)
left.set_ylim(-0.3, 1.05)
left.set_xlabel(r"$x = r/(2M)$")
left.set_ylabel(r"$f(x)$")
hs.panel(left, r"(a) the function $f$, whose zeros are the horizons")
right.set_xlim(0, 4.0)
right.set_ylim(0, 1.05)
right.set_xlabel(r"$x = r/(2M)$")
right.set_ylabel(r"$e^{-2g(x)}$")
hs.panel(right, r"(b) the factor $e^{-2g}$ of the dirtiness")
hs.legend_outside(fig, where="lower", ncol=4)

out = HERE / "dirtybh-metric.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
