"""Statistics of repeated and correlated data.

Standard error of the mean, weighted mean, autocorrelation analysis, block
averaging, bootstrap and jackknife.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "mean_sem",
    "weighted_mean",
    "autocorrelation",
    "integrated_autocorrelation_time",
    "sem_correlated",
    "block_means",
    "block_sem",
    "blocking_curve",
    "bootstrap",
    "jackknife",
]


def mean_sem(x):
    """Mean, sample standard deviation and standard error of the mean.

    The standard error ``s / sqrt(N)`` is valid for independent readings.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 2:
        raise ValueError("need at least two readings")
    s = x.std(ddof=1)
    return float(x.mean()), float(s), float(s / np.sqrt(n))


def weighted_mean(x, u):
    """Inverse-variance weighted mean of independent results.

    Returns
    -------
    mean : float
    u_mean : float
        ``1 / sqrt(sum(1 / u_i**2))``.
    chi2_red : float
        Reduced chi-square of the results about their weighted mean.  A value
        far above one says the results scatter more than their quoted
        uncertainties allow.
    """
    x = np.asarray(x, dtype=float)
    u = np.asarray(u, dtype=float)
    w = 1.0 / u**2
    mean = np.sum(w * x) / np.sum(w)
    u_mean = 1.0 / np.sqrt(np.sum(w))
    chi2_red = np.sum(w * (x - mean) ** 2) / (x.size - 1) if x.size > 1 else np.nan
    return float(mean), float(u_mean), float(chi2_red)


def autocorrelation(x, max_lag=None):
    """Normalised autocorrelation function ``rho(t)`` of a time series.

    Uses the FFT and the usual (biased) estimator that divides every lag by
    ``N``.  ``rho[0] == 1``.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    dx = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    fx = np.fft.rfft(dx, nfft)
    acov = np.fft.irfft(fx * np.conj(fx), nfft)[:n] / n
    rho = acov / acov[0]
    if max_lag is not None:
        rho = rho[: max_lag + 1]
    return rho


def integrated_autocorrelation_time(x, c=5.0):
    """Integrated autocorrelation time ``tau = 1/2 + sum_{t>=1} rho(t)``.

    The sum is cut at the smallest window ``M`` for which ``M >= c * tau(M)``
    (the automatic windowing rule of Sokal), which keeps the noise in the tail
    of ``rho`` out of the estimate.

    With this definition the variance of the mean of ``N`` correlated values
    is ``2 * tau * s**2 / N``; independent data have ``tau = 1/2``.

    Returns
    -------
    tau : float
        In units of the sampling interval.
    window : int
    """
    rho = autocorrelation(x)
    tau = 0.5 + np.cumsum(rho[1:])
    M = np.arange(1, tau.size + 1)
    ok = M >= c * tau
    i = int(np.argmax(ok)) if ok.any() else tau.size - 1
    return float(max(tau[i], 0.5)), int(M[i])


def sem_correlated(x, c=5.0):
    """Standard error of the mean of a correlated series, from ``tau``."""
    x = np.asarray(x, dtype=float)
    tau, _ = integrated_autocorrelation_time(x, c=c)
    return float(x.std(ddof=1) * np.sqrt(2.0 * tau / x.size))


def block_means(x, block_size):
    """Means of consecutive non-overlapping blocks (the remainder is dropped)."""
    x = np.asarray(x, dtype=float)
    nb = x.size // block_size
    if nb < 1:
        raise ValueError("block_size is longer than the series")
    return x[: nb * block_size].reshape(nb, block_size).mean(axis=1)


def block_sem(x, n_blocks=20):
    """Standard error of the mean from ``n_blocks`` block averages.

    Reliable only when each block is much longer than the correlation time.
    """
    x = np.asarray(x, dtype=float)
    m = block_means(x, x.size // n_blocks)
    return float(m.std(ddof=1) / np.sqrt(m.size))


def blocking_curve(x, min_blocks=16):
    """Blocking analysis in the manner of Flyvbjerg and Petersen.

    The series is cut into blocks of size 1, 2, 4, ... and the standard error
    of the mean is estimated from the block averages at every size.  For a
    long enough series the estimate rises and then levels off; the plateau is
    the standard error.

    Returns
    -------
    block_size : ndarray of int
    sem : ndarray
        Standard-error estimate at each block size.
    sem_err : ndarray
        Statistical uncertainty of that estimate, ``sem / sqrt(2 (n_b - 1))``.
    """
    x = np.asarray(x, dtype=float)
    sizes, sems, errs = [], [], []
    b = 1
    while x.size // b >= min_blocks:
        m = block_means(x, b)
        sem = m.std(ddof=1) / np.sqrt(m.size)
        sizes.append(b)
        sems.append(sem)
        errs.append(sem / np.sqrt(2.0 * (m.size - 1)))
        b *= 2
    return np.array(sizes), np.array(sems), np.array(errs)


def bootstrap(data, stat, n_boot=2000, rng=None):
    """Bootstrap standard error of ``stat(data)`` for independent rows.

    Returns
    -------
    estimate : float or ndarray
        ``stat`` evaluated on the original data.
    se : float or ndarray
        Standard deviation of the bootstrap replicas.
    replicas : ndarray
    """
    rng = np.random.default_rng(rng)
    data = np.asarray(data)
    n = data.shape[0]
    reps = np.array([stat(data[rng.integers(0, n, n)]) for _ in range(n_boot)])
    return stat(data), reps.std(axis=0, ddof=1), reps


def jackknife(data, stat):
    """Leave-one-out jackknife standard error of ``stat(data)``.

    ``data`` holds independent rows; pass block averages to treat a correlated
    series.

    Returns
    -------
    estimate : float or ndarray
    se : float or ndarray
    """
    data = np.asarray(data)
    n = data.shape[0]
    idx = np.arange(n)
    theta = np.array([stat(data[idx != i]) for i in range(n)])
    se = np.sqrt((n - 1) / n * np.sum((theta - theta.mean(axis=0)) ** 2, axis=0))
    return stat(data), se
