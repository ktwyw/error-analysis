"""Counting intervals, robust estimators, nonlinear fitting and numerical tools."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

__all__ = [
    "poisson_interval",
    "mad",
    "huber_location",
    "NonlinearFit",
    "nonlinear_fit",
    "richardson",
    "observed_order",
]


# ------------------------------------------------------------------ counting
def poisson_interval(n, cl=0.6827):
    """Exact (Garwood) central confidence interval for a Poisson mean.

    Given an observed count ``n``, returns ``(lower, upper)`` such that the
    interval covers the true mean with probability at least ``cl``. Uses the
    relation between the Poisson and chi-square distributions.
    """
    n = np.asarray(n, dtype=float)
    a = 1.0 - cl
    lower = np.where(n > 0, 0.5 * stats.chi2.ppf(a / 2, 2 * n), 0.0)
    upper = 0.5 * stats.chi2.ppf(1 - a / 2, 2 * (n + 1))
    return lower, upper


# ------------------------------------------------------------------ robust statistics
def mad(x, scale="normal"):
    """Median absolute deviation, scaled to estimate sigma for normal data."""
    x = np.asarray(x, dtype=float)
    m = np.median(np.abs(x - np.median(x)))
    return float(m * 1.482602218505602) if scale == "normal" else float(m)


def huber_location(x, c=1.345, tol=1e-10, max_iter=200):
    """Huber M-estimate of location with a MAD scale (iteratively reweighted mean).

    Points within ``c`` robust standard deviations of the estimate get full
    weight; points further out are down-weighted in proportion to their
    distance. ``c = 1.345`` gives 95 % efficiency for normal data.
    """
    x = np.asarray(x, dtype=float)
    s = mad(x)
    if s == 0:
        return float(np.median(x))
    mu = np.median(x)
    for _ in range(max_iter):
        r = (x - mu) / s
        w = np.minimum(1.0, c / np.maximum(np.abs(r), 1e-300))
        mu_new = np.sum(w * x) / np.sum(w)
        if abs(mu_new - mu) < tol * s:
            mu = mu_new
            break
        mu = mu_new
    return float(mu)


# ------------------------------------------------------------------ nonlinear least squares
@dataclass
class NonlinearFit:
    """Result of :func:`nonlinear_fit`."""

    params: np.ndarray
    cov: np.ndarray
    chi2: float
    dof: int
    residuals: np.ndarray       # (y - model) / sigma
    success: bool
    nfev: int

    @property
    def errors(self):
        return np.sqrt(np.diag(self.cov))

    @property
    def chi2_red(self):
        return self.chi2 / self.dof

    @property
    def p_value(self):
        return float(stats.chi2.sf(self.chi2, self.dof))

    @property
    def correlation(self):
        e = self.errors
        return self.cov / np.outer(e, e)


def nonlinear_fit(model, x, y, sigma, p0, scale_cov=False, **kwargs):
    """Weighted nonlinear least squares by the Levenberg-Marquardt method.

    Parameters
    ----------
    model : callable
        ``model(x, *params)``.
    x, y, sigma : array_like
        Data and the standard uncertainties of ``y``.
    p0 : sequence
        Starting values.
    scale_cov : bool
        Multiply the covariance by ``chi2/dof`` (use when ``sigma`` gives only
        relative weights).

    The covariance is ``(J^T W J)^{-1}`` evaluated at the minimum, where
    ``J`` is the Jacobian of the model with respect to the parameters.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), y.shape)

    def res(p):
        return (y - model(x, *p)) / sigma

    kw = dict(method="lm", x_scale="jac")
    kw.update(kwargs)
    sol = optimize.least_squares(res, np.asarray(p0, dtype=float), **kw)
    J = sol.jac                              # d(res)/dp = -(1/sigma) d(model)/dp
    cov = np.linalg.inv(J.T @ J)
    r = res(sol.x)
    chi2 = float(r @ r)
    dof = y.size - sol.x.size
    if scale_cov:
        cov = cov * chi2 / dof
    return NonlinearFit(params=sol.x, cov=cov, chi2=chi2, dof=dof, residuals=r,
                        success=bool(sol.success), nfev=int(sol.nfev))


# ------------------------------------------------------------------ numerics
def observed_order(f1, f2, f3, ratio=2.0):
    """Observed order of convergence from results on three grids.

    ``f1, f2, f3`` are computed with step sizes ``h, h*ratio, h*ratio**2``
    (fine to coarse).
    """
    return float(np.log(abs((f3 - f2) / (f2 - f1))) / np.log(ratio))


def richardson(f_fine, f_coarse, p, ratio=2.0):
    """Richardson extrapolation and an error estimate for the fine result.

    Returns ``(f_extrapolated, error_estimate_of_f_fine)`` for a method of
    order ``p`` whose results ``f_fine`` (step ``h``) and ``f_coarse`` (step
    ``ratio*h``) are given.
    """
    delta = (f_fine - f_coarse) / (ratio**p - 1.0)
    return f_fine + delta, abs(delta)
