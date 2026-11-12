"""Sensitivity analysis: local (derivative-based) and global (variance-based).

Local measures answer "how strongly does the output respond to a small change
of this input?"; variance-based (Sobol) indices answer "what fraction of the
output variance is due to this input?" over the full range of the input
distributions, including nonlinearity and interactions.
"""

from __future__ import annotations

import numpy as np

from .propagation import jacobian

__all__ = ["elasticities", "sobol_indices", "ishigami"]


def elasticities(f, x, u=None):
    """Normalised (logarithmic) sensitivity coefficients of a scalar model.

    ``E_i = (df/dx_i) * x_i / f``: the percentage change of the output per
    percentage change of input ``i``.  For a product of powers
    ``f = C x1**p1 x2**p2 ...`` the elasticities are the exponents.
    """
    x = np.asarray(x, dtype=float)
    y = float(f(*x))
    c = jacobian(f, x, u=u)[0]
    return c * x / y


def sobol_indices(f, sample, n=2**14, rng=None, n_boot=200):
    """First-order and total Sobol indices by Monte Carlo.

    Uses the Saltelli (2010) estimator for the first-order index and the
    Jansen (1999) estimator for the total index, which together need
    ``n * (k + 2)`` model evaluations for ``k`` inputs.

    Parameters
    ----------
    f : callable
        Vectorised model ``f(x1, x2, ..., xk)`` accepting arrays of samples.
    sample : callable
        ``sample(n, rng)`` returning an array of shape ``(n, k)`` of
        independent draws of the inputs.
    n : int
        Base sample size.
    rng : numpy.random.Generator or int, optional
    n_boot : int
        Bootstrap replicas used for the standard errors (0 to skip).

    Returns
    -------
    dict with keys
        ``S`` (first-order), ``ST`` (total), ``S_se``, ``ST_se`` (bootstrap
        standard errors), ``var`` (output variance) and ``n_evals``.
    """
    rng = np.random.default_rng(rng)
    A = np.asarray(sample(n, rng), dtype=float)
    B = np.asarray(sample(n, rng), dtype=float)
    k = A.shape[1]
    fA = np.asarray(f(*A.T), dtype=float)
    fB = np.asarray(f(*B.T), dtype=float)
    fAB = np.empty((k, n))
    for i in range(k):
        ABi = A.copy()
        ABi[:, i] = B[:, i]
        fAB[i] = f(*ABi.T)

    def estimate(idx):
        a, b, ab = fA[idx], fB[idx], fAB[:, idx]
        var = np.var(np.concatenate([a, b]), ddof=1)
        S = np.mean(b * (ab - a), axis=1) / var
        ST = 0.5 * np.mean((a - ab) ** 2, axis=1) / var
        return S, ST, var

    S, ST, var = estimate(np.arange(n))
    out = dict(S=S, ST=ST, var=float(var), n_evals=n * (k + 2))
    if n_boot:
        reps = [estimate(rng.integers(0, n, n)) for _ in range(n_boot)]
        out["S_se"] = np.std([r[0] for r in reps], axis=0, ddof=1)
        out["ST_se"] = np.std([r[1] for r in reps], axis=0, ddof=1)
    return out


def ishigami(x1, x2, x3, a=7.0, b=0.1):
    """Ishigami test function, a standard benchmark for Sobol indices.

    With ``x_i`` uniform on ``[-pi, pi]`` the exact indices are known in
    closed form; see :func:`ishigami_exact`.
    """
    return np.sin(x1) + a * np.sin(x2) ** 2 + b * x3**4 * np.sin(x1)


def ishigami_exact(a=7.0, b=0.1):
    """Exact first-order and total Sobol indices of the Ishigami function."""
    pi = np.pi
    V1 = 0.5 * (1 + b * pi**4 / 5) ** 2
    V2 = a**2 / 8
    V13 = b**2 * pi**8 * (1 / 18 - 1 / 50)
    V = V1 + V2 + V13
    S = np.array([V1, V2, 0.0]) / V
    ST = np.array([V1 + V13, V2, V13]) / V
    return S, ST


__all__.append("ishigami_exact")
