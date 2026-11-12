"""Least-squares fitting with parameter covariance.

Only linear least squares is implemented here, in closed form, so that every
step can be followed.  For nonlinear models use ``scipy.optimize.curve_fit``,
which returns the same kind of covariance matrix.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["general_linear_fit", "linear_fit", "LineFit"]


def general_linear_fit(A, y, sigma=None):
    """Weighted linear least squares for a model ``y = A @ beta``.

    Parameters
    ----------
    A : ndarray, shape (N, p)
        Design matrix; column ``j`` is the j-th basis function evaluated at
        the data points.
    y : ndarray, shape (N,)
    sigma : ndarray, shape (N,), optional
        Standard uncertainties of ``y``.  If omitted, all points get equal
        weight and the covariance is scaled by the residual variance.

    Returns
    -------
    beta : ndarray, shape (p,)
    cov : ndarray, shape (p, p)
    chi2 : float
        ``sum(((y - A beta) / sigma)**2)``, or the residual sum of squares
        when ``sigma`` is omitted.
    dof : int
    """
    A = np.asarray(A, dtype=float)
    y = np.asarray(y, dtype=float)
    N, p = A.shape
    dof = N - p
    if sigma is None:
        w = np.ones(N)
    else:
        w = 1.0 / np.asarray(sigma, dtype=float) ** 2
    AtWA = A.T @ (w[:, None] * A)
    cov = np.linalg.inv(AtWA)
    beta = cov @ (A.T @ (w * y))
    r = y - A @ beta
    chi2 = float(np.sum(w * r**2))
    if sigma is None:
        cov = cov * chi2 / dof
    return beta, cov, chi2, dof


@dataclass
class LineFit:
    """Result of a straight-line fit ``y = intercept + slope * x``."""

    intercept: float
    slope: float
    cov: np.ndarray  # order: [intercept, slope]
    chi2: float
    dof: int
    weighted: bool
    x: np.ndarray
    y: np.ndarray

    @property
    def u_intercept(self):
        return float(np.sqrt(self.cov[0, 0]))

    @property
    def u_slope(self):
        return float(np.sqrt(self.cov[1, 1]))

    @property
    def correlation(self):
        """Correlation coefficient between intercept and slope."""
        return float(self.cov[0, 1] / np.sqrt(self.cov[0, 0] * self.cov[1, 1]))

    @property
    def chi2_red(self):
        """Reduced chi-square (meaningful only for a weighted fit)."""
        return self.chi2 / self.dof if self.weighted else float("nan")

    @property
    def residual_std(self):
        """Scatter of the points about the line, ``sqrt(RSS / (N - 2))``.

        For an unweighted fit this is the estimate of the y-uncertainty of a
        single point."""
        r = self.residuals
        return float(np.sqrt(np.sum(r**2) / self.dof))

    @property
    def residuals(self):
        return self.y - (self.intercept + self.slope * self.x)

    def predict(self, x0):
        """Fitted value and its standard uncertainty at ``x0``.

        The uncertainty includes the covariance between intercept and slope,
        ``u^2 = C_aa + 2 x0 C_ab + x0^2 C_bb``.
        """
        x0 = np.asarray(x0, dtype=float)
        y0 = self.intercept + self.slope * x0
        var = self.cov[0, 0] + 2 * x0 * self.cov[0, 1] + x0**2 * self.cov[1, 1]
        return y0, np.sqrt(var)

    def __str__(self):
        s = (
            f"intercept = {self.intercept:.6g} +/- {self.u_intercept:.2g}\n"
            f"slope     = {self.slope:.6g} +/- {self.u_slope:.2g}\n"
            f"correlation(intercept, slope) = {self.correlation:+.3f}\n"
        )
        if self.weighted:
            s += f"chi2 / dof = {self.chi2:.2f} / {self.dof} = {self.chi2_red:.2f}"
        else:
            s += f"residual standard deviation = {self.residual_std:.3g}  (dof = {self.dof})"
        return s


def linear_fit(x, y, sigma=None, scale_cov=False):
    """Straight-line least-squares fit ``y = a + b x``.

    Parameters
    ----------
    x, y : array_like
    sigma : array_like or float, optional
        Standard uncertainties of ``y``.  Without them the fit is unweighted
        and the parameter uncertainties are based on the scatter of the
        points about the line.
    scale_cov : bool
        For a weighted fit, multiply the covariance by ``chi2 / dof``.  Use
        this when ``sigma`` gives only the relative weights of the points.

    Returns
    -------
    LineFit
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if sigma is not None:
        sigma = np.broadcast_to(np.asarray(sigma, dtype=float), y.shape)
    A = np.column_stack([np.ones_like(x), x])
    beta, cov, chi2, dof = general_linear_fit(A, y, sigma)
    if sigma is not None and scale_cov:
        cov = cov * chi2 / dof
    return LineFit(
        intercept=float(beta[0]),
        slope=float(beta[1]),
        cov=cov,
        chi2=chi2,
        dof=dof,
        weighted=sigma is not None,
        x=x,
        y=y,
    )
