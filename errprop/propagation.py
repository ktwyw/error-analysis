"""Propagation of uncertainty: first-order (linear) formula and Monte Carlo.

Conventions
-----------
A model is an ordinary Python function ``f(x1, x2, ...)`` of its input
quantities.  Inputs are given as a sequence of best estimates ``x`` together
with either standard uncertainties ``u`` (independent inputs) or a full
covariance matrix ``cov``.

The functions are written to be read.  For production work on long
expressions, the ``uncertainties`` package does the same job with automatic
differentiation.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "jacobian",
    "propagate",
    "budget",
    "format_budget",
    "monte_carlo",
    "summarize",
]


def _covariance(u, cov, n):
    """Return the n-by-n input covariance matrix from ``u`` or ``cov``."""
    if cov is not None:
        cov = np.asarray(cov, dtype=float)
        if cov.shape != (n, n):
            raise ValueError(f"cov must have shape ({n}, {n}), got {cov.shape}")
        return cov
    if u is None:
        raise ValueError("give either the standard uncertainties u or cov")
    u = np.asarray(u, dtype=float)
    if u.shape != (n,):
        raise ValueError(f"u must have length {n}, got shape {u.shape}")
    return np.diag(u**2)


def jacobian(f, x, u=None, rel_step=1e-5):
    """Sensitivity coefficients df/dx_i by central differences.

    Parameters
    ----------
    f : callable
        Model ``f(x1, x2, ...)`` returning a scalar or a 1-D array.
    x : sequence of float
        Point at which the derivatives are evaluated.
    u : sequence of float, optional
        Standard uncertainties; used only to choose a step for inputs whose
        best estimate is exactly zero.
    rel_step : float
        Step relative to ``|x_i|``.

    Returns
    -------
    J : ndarray, shape (n_outputs, n_inputs)
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    y0 = np.atleast_1d(np.asarray(f(*x), dtype=float))
    J = np.empty((y0.size, n))
    for i in range(n):
        if x[i] != 0.0:
            h = rel_step * abs(x[i])
        elif u is not None and np.asarray(u, dtype=float)[i] > 0:
            h = rel_step * float(np.asarray(u, dtype=float)[i])
        else:
            h = rel_step
        xp, xm = x.copy(), x.copy()
        xp[i] += h
        xm[i] -= h
        yp = np.atleast_1d(np.asarray(f(*xp), dtype=float))
        ym = np.atleast_1d(np.asarray(f(*xm), dtype=float))
        J[:, i] = (yp - ym) / (2.0 * h)
    return J


def propagate(f, x, u=None, cov=None):
    """First-order propagation of uncertainty, ``C_y = J C_x J^T``.

    Returns
    -------
    y, u_y : float, float
        For a scalar model: the output estimate and its standard uncertainty.
    y, C_y : ndarray, ndarray
        For a vector model: the output estimates and their covariance matrix.
    """
    x = np.asarray(x, dtype=float)
    C = _covariance(u, cov, x.size)
    y = np.asarray(f(*x), dtype=float)
    J = jacobian(f, x, u=u)
    Cy = J @ C @ J.T
    if y.ndim == 0:
        return float(y), float(np.sqrt(Cy[0, 0]))
    return y, Cy


def budget(f, x, u, names=None):
    """Uncertainty budget for a scalar model with independent inputs.

    Each row lists the input, its standard uncertainty ``u_i``, the
    sensitivity coefficient ``c_i = df/dx_i``, the contribution ``|c_i| u_i``
    to the output uncertainty, and the fraction of the output variance that
    the input is responsible for.

    Returns
    -------
    rows : list of dict
    y : float
    u_y : float
    """
    x = np.asarray(x, dtype=float)
    u = np.asarray(u, dtype=float)
    if names is None:
        names = [f"x{i + 1}" for i in range(x.size)]
    y = float(f(*x))
    c = jacobian(f, x, u=u)[0]
    contrib = np.abs(c) * u
    var = float(np.sum(contrib**2))
    rows = []
    for i in range(x.size):
        rows.append(
            dict(
                name=names[i],
                value=float(x[i]),
                u=float(u[i]),
                sensitivity=float(c[i]),
                contribution=float(contrib[i]),
                variance_fraction=float(contrib[i] ** 2 / var) if var > 0 else 0.0,
            )
        )
    return rows, y, float(np.sqrt(var))


def format_budget(rows, y=None, u_y=None):
    """Plain-text table of an uncertainty budget returned by :func:`budget`."""
    head = f"{'input':<10}{'value':>13}{'u':>12}{'rel. u':>10}{'|c| u':>13}{'share':>9}"
    lines = [head, "-" * len(head)]
    for r in rows:
        rel = abs(r["u"] / r["value"]) if r["value"] != 0 else float("nan")
        lines.append(
            f"{r['name']:<10}{r['value']:>13.5g}{r['u']:>12.2g}{rel:>10.2%}"
            f"{r['contribution']:>13.3g}{r['variance_fraction']:>9.1%}"
        )
    if y is not None and u_y is not None:
        lines.append("-" * len(head))
        lines.append(
            f"{'output':<10}{y:>13.5g}{u_y:>12.2g}{abs(u_y / y):>10.2%}"
        )
    return "\n".join(lines)


def monte_carlo(f, x, u=None, cov=None, n=200_000, dist=None, rng=None,
                vectorized=True):
    """Propagate distributions through a model by random sampling.

    Parameters
    ----------
    f : callable
        Model ``f(x1, x2, ...)``.  With ``vectorized=True`` it must accept
        NumPy arrays of samples for every argument.
    x, u, cov
        Best estimates and either independent standard uncertainties or a
        covariance matrix (which implies a joint normal distribution).
    n : int
        Number of samples.
    dist : sequence of {"normal", "uniform"}, optional
        Distribution of each independent input.  A "uniform" input has
        half-width ``u * sqrt(3)`` so that its standard deviation is ``u``.
    rng : numpy.random.Generator or int, optional

    Returns
    -------
    samples : ndarray
        ``n`` samples of the output (shape ``(n,)`` for a scalar model).
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(x, dtype=float)
    m = x.size
    if cov is not None:
        if dist is not None:
            raise ValueError("dist cannot be combined with cov")
        X = rng.multivariate_normal(x, _covariance(None, cov, m), size=n)
    else:
        u = np.asarray(u, dtype=float)
        if u.shape != (m,):
            raise ValueError(f"u must have length {m}")
        if dist is None:
            dist = ["normal"] * m
        X = np.empty((n, m))
        for i, d in enumerate(dist):
            if d == "normal":
                X[:, i] = rng.normal(x[i], u[i], size=n)
            elif d == "uniform":
                a = u[i] * np.sqrt(3.0)
                X[:, i] = rng.uniform(x[i] - a, x[i] + a, size=n)
            else:
                raise ValueError(f"unknown distribution {d!r}")
    if vectorized:
        return np.asarray(f(*X.T), dtype=float)
    return np.array([f(*row) for row in X], dtype=float)


def summarize(samples, level=0.6827):
    """Mean, standard deviation, median and a central probability interval."""
    s = np.asarray(samples, dtype=float)
    lo, hi = np.quantile(s, [(1 - level) / 2, (1 + level) / 2])
    return dict(
        mean=float(s.mean()),
        std=float(s.std(ddof=1)),
        median=float(np.median(s)),
        low=float(lo),
        high=float(hi),
        level=level,
    )
