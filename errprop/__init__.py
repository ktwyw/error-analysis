"""errprop: small, readable tools for error analysis and propagation.

The module accompanies the notebooks of the repository.  It favours clarity
over generality: every function is short enough to read in one sitting.
"""

from .fitting import LineFit, general_linear_fit, linear_fit
from .plotting import panel_label, use_style
from .propagation import (
    budget,
    format_budget,
    jacobian,
    monte_carlo,
    propagate,
    summarize,
)
from .reporting import format_concise, format_result, round_to_uncertainty
from .sensitivity import elasticities, ishigami, ishigami_exact, sobol_indices
from .extra import (
    NonlinearFit,
    huber_location,
    mad,
    nonlinear_fit,
    observed_order,
    poisson_interval,
    richardson,
)
from .stats import (
    autocorrelation,
    block_means,
    block_sem,
    blocking_curve,
    bootstrap,
    integrated_autocorrelation_time,
    jackknife,
    mean_sem,
    sem_correlated,
    weighted_mean,
)

__version__ = "0.3.0"

__all__ = [
    "NonlinearFit",
    "huber_location",
    "mad",
    "nonlinear_fit",
    "observed_order",
    "poisson_interval",
    "richardson",
    "LineFit",
    "autocorrelation",
    "block_means",
    "block_sem",
    "blocking_curve",
    "bootstrap",
    "budget",
    "elasticities",
    "format_budget",
    "format_concise",
    "format_result",
    "general_linear_fit",
    "integrated_autocorrelation_time",
    "ishigami",
    "ishigami_exact",
    "jackknife",
    "jacobian",
    "linear_fit",
    "mean_sem",
    "monte_carlo",
    "panel_label",
    "propagate",
    "round_to_uncertainty",
    "sem_correlated",
    "sobol_indices",
    "summarize",
    "use_style",
    "weighted_mean",
]
