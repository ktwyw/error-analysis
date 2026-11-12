"""Rounding and formatting of a result with its uncertainty."""

from __future__ import annotations

import math

__all__ = ["round_to_uncertainty", "format_result", "format_concise"]


def round_to_uncertainty(value, u, sig=2):
    """Round ``u`` to ``sig`` significant figures and ``value`` to match.

    Returns
    -------
    value_rounded, u_rounded : float
    decimals : int
        Number of decimal places kept (negative means rounding to tens,
        hundreds, ...).
    """
    if not (u > 0) or not math.isfinite(u):
        raise ValueError("the uncertainty must be positive and finite")
    exponent = math.floor(math.log10(u))
    decimals = sig - 1 - exponent
    u_r = round(u, decimals)
    if math.floor(math.log10(u_r)) != exponent:  # e.g. 0.0996 -> 0.10
        exponent += 1
        decimals = sig - 1 - exponent
        u_r = round(u, decimals)
    return round(value, decimals), u_r, decimals


def format_result(value, u, unit="", sig=2, sci=None):
    """Format ``value ± u`` with consistent rounding.

    >>> format_result(9.81234, 0.0234, "m/s^2")
    '9.812 ± 0.023 m/s^2'
    >>> format_result(1.3456e-4, 2.1e-6, "1/s")
    '(1.346 ± 0.021) × 10^-4 1/s'
    """
    v, ur, dec = round_to_uncertainty(value, u, sig)
    lead = v if v != 0 else ur
    e = math.floor(math.log10(abs(lead)))
    if sci is None:
        sci = e < -3 or e >= 5
    if sci:
        d = max(dec + e, 0)
        scale = 10.0**e
        s = f"({v / scale:.{d}f} ± {ur / scale:.{d}f}) × 10^{e}"
    else:
        d = max(dec, 0)
        s = f"{v:.{d}f} ± {ur:.{d}f}"
    return f"{s} {unit}".rstrip()


def format_concise(value, u, unit="", sig=2):
    """Concise notation with the uncertainty in brackets, e.g. ``9.812(23)``.

    The bracketed digits are the uncertainty in units of the last quoted
    decimal place.
    """
    v, ur, dec = round_to_uncertainty(value, u, sig)
    if dec <= 0:
        return format_result(value, u, unit, sig, sci=False)
    digits = round(ur * 10**dec)
    return f"{v:.{dec}f}({digits:d}) {unit}".rstrip()
