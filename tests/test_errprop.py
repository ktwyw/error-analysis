"""Tests of the errprop module against results known in closed form."""

import numpy as np
import pytest
from scipy.signal import lfilter
from scipy import stats

import errprop as ep


# ----------------------------------------------------------------- propagation
def test_sum_and_difference_add_in_quadrature():
    for f in (lambda a, b: a + b, lambda a, b: a - b):
        _, u = ep.propagate(f, [10.0, 4.0], [0.3, 0.4])
        assert u == pytest.approx(0.5, rel=1e-8)


def test_product_and_power_rule():
    # y = a * b**2 / c : relative uncertainties 1 %, 2 x 2 %, 3 %
    f = lambda a, b, c: a * b**2 / c
    x = [2.0, 3.0, 5.0]
    u = [0.02, 0.06, 0.15]
    y, uy = ep.propagate(f, x, u)
    assert y == pytest.approx(3.6)
    assert uy / y == pytest.approx(np.sqrt(0.01**2 + 0.04**2 + 0.03**2), rel=1e-6)


def test_fully_correlated_inputs_cancel_in_a_difference():
    cov = 0.2**2 * np.ones((2, 2))
    _, u = ep.propagate(lambda a, b: a - b, [5.0, 3.0], cov=cov)
    assert u == pytest.approx(0.0, abs=1e-9)
    _, u = ep.propagate(lambda a, b: a + b, [5.0, 3.0], cov=cov)
    assert u == pytest.approx(0.4, rel=1e-8)


def test_vector_output_covariance():
    f = lambda a, b: np.array([a + b, a - b])
    y, C = ep.propagate(f, [1.0, 2.0], [0.3, 0.4])
    assert np.allclose(y, [3.0, -1.0])
    assert np.allclose(C, [[0.25, 0.09 - 0.16], [0.09 - 0.16, 0.25]])


def test_budget_shares_sum_to_one():
    f = lambda m, d, h: 4 * m / (np.pi * d**2 * h)
    rows, y, uy = ep.budget(f, [23.45e-3, 9.98e-3, 35.10e-3],
                            [0.01e-3, 0.02e-3, 0.05e-3], ["m", "d", "h"])
    assert sum(r["variance_fraction"] for r in rows) == pytest.approx(1.0)
    assert max(rows, key=lambda r: r["variance_fraction"])["name"] == "d"
    _, uy2 = ep.propagate(f, [23.45e-3, 9.98e-3, 35.10e-3],
                          [0.01e-3, 0.02e-3, 0.05e-3])
    assert uy == pytest.approx(uy2)


def test_monte_carlo_agrees_with_linear_formula_for_a_linear_model():
    f = lambda a, b: 2 * a - 3 * b
    s = ep.monte_carlo(f, [1.0, 1.0], [0.1, 0.2], n=400_000, rng=1)
    assert s.std(ddof=1) == pytest.approx(np.sqrt(0.2**2 + 0.6**2), rel=0.01)


def test_monte_carlo_uniform_input_has_requested_standard_deviation():
    s = ep.monte_carlo(lambda a: a, [0.0], [0.5], dist=["uniform"],
                       n=400_000, rng=2)
    assert s.std(ddof=1) == pytest.approx(0.5, rel=0.01)
    assert np.abs(s).max() <= 0.5 * np.sqrt(3) + 1e-12


def test_monte_carlo_lognormal_mean():
    s = ep.monte_carlo(np.exp, [0.0], [0.5], n=400_000, rng=3)
    assert s.mean() == pytest.approx(np.exp(0.125), rel=0.01)


# ------------------------------------------------------------------ statistics
def test_mean_sem():
    m, s, sem = ep.mean_sem([1.0, 2.0, 3.0, 4.0, 5.0])
    assert m == pytest.approx(3.0)
    assert s == pytest.approx(np.sqrt(2.5))
    assert sem == pytest.approx(np.sqrt(2.5 / 5))


def test_weighted_mean():
    m, u, _ = ep.weighted_mean([10.0, 12.0], [1.0, 1.0])
    assert m == pytest.approx(11.0)
    assert u == pytest.approx(1 / np.sqrt(2))
    m, u, _ = ep.weighted_mean([10.0, 12.0], [1.0, 2.0])
    assert m == pytest.approx(10.4)
    assert u == pytest.approx(np.sqrt(0.8))


def _ar1(phi, n, rng, burn=2000):
    xi = rng.normal(size=n + burn) * np.sqrt(1 - phi**2)
    return lfilter([1.0], [1.0, -phi], xi)[burn:]


def test_autocorrelation_time_of_ar1():
    phi = 0.9
    x = _ar1(phi, 400_000, np.random.default_rng(4))
    rho = ep.autocorrelation(x, max_lag=5)
    assert rho[0] == pytest.approx(1.0)
    assert np.allclose(rho[1:], phi ** np.arange(1, 6), atol=0.02)
    tau, _ = ep.integrated_autocorrelation_time(x)
    assert tau == pytest.approx(0.5 * (1 + phi) / (1 - phi), rel=0.1)


def test_white_noise_has_tau_one_half():
    x = np.random.default_rng(5).normal(size=200_000)
    tau, _ = ep.integrated_autocorrelation_time(x)
    assert tau == pytest.approx(0.5, abs=0.05)
    assert ep.sem_correlated(x) == pytest.approx(ep.mean_sem(x)[2], rel=0.1)


def test_blocking_plateau_matches_exact_standard_error():
    phi, n = 0.9, 2**18
    x = _ar1(phi, n, np.random.default_rng(6))
    exact = np.sqrt((1 + phi) / (1 - phi) / n)
    sizes, sem, err = ep.blocking_curve(x)
    assert sem[0] == pytest.approx(exact / np.sqrt(19), rel=0.05)  # naive
    plateau = sem[sizes >= 512]
    assert np.all(np.abs(plateau - exact) < 4 * err[sizes >= 512])
    assert ep.block_sem(x, 64) == pytest.approx(exact, rel=0.35)


def test_block_means():
    assert np.allclose(ep.block_means(np.arange(10.0), 4), [1.5, 5.5])


def test_bootstrap_and_jackknife_reproduce_sem_for_the_mean():
    x = np.random.default_rng(7).normal(size=400)
    sem = ep.mean_sem(x)[2]
    _, se_j = ep.jackknife(x, np.mean)
    assert se_j == pytest.approx(sem, rel=1e-10)
    _, se_b, reps = ep.bootstrap(x, np.mean, n_boot=4000, rng=8)
    assert reps.shape == (4000,)
    assert se_b == pytest.approx(sem, rel=0.06)


# --------------------------------------------------------------------- fitting
def test_unweighted_line_fit_matches_textbook_formulas():
    rng = np.random.default_rng(9)
    x = np.linspace(1, 10, 12)
    y = 2.0 + 0.5 * x + rng.normal(0, 0.1, x.size)
    fit = ep.linear_fit(x, y)
    Sxx = np.sum((x - x.mean()) ** 2)
    b = np.sum((x - x.mean()) * (y - y.mean())) / Sxx
    a = y.mean() - b * x.mean()
    s2 = np.sum((y - a - b * x) ** 2) / (x.size - 2)
    assert fit.slope == pytest.approx(b)
    assert fit.intercept == pytest.approx(a)
    assert fit.u_slope == pytest.approx(np.sqrt(s2 / Sxx))
    assert fit.u_intercept == pytest.approx(
        np.sqrt(s2 * (1 / x.size + x.mean() ** 2 / Sxx)))
    assert fit.cov[0, 1] == pytest.approx(-x.mean() * s2 / Sxx)
    assert np.allclose(np.polyfit(x, y, 1), [fit.slope, fit.intercept])


def test_weighted_fit_covariance_and_prediction():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([1.0, 3.1, 4.9, 7.0])
    fit = ep.linear_fit(x, y, sigma=0.1)
    Sxx = np.sum((x - x.mean()) ** 2)
    assert fit.u_slope == pytest.approx(0.1 / np.sqrt(Sxx))
    # at the centroid the uncertainty of the fitted value is sigma / sqrt(N)
    _, u0 = fit.predict(x.mean())
    assert u0 == pytest.approx(0.1 / 2)
    assert fit.dof == 2 and fit.weighted


def test_fit_uncertainty_has_correct_coverage():
    rng = np.random.default_rng(10)
    x = np.linspace(0, 5, 8)
    slopes, us = [], []
    for _ in range(2000):
        y = 1.0 + 2.0 * x + rng.normal(0, 0.3, x.size)
        fit = ep.linear_fit(x, y, sigma=0.3)
        slopes.append(fit.slope)
        us.append(fit.u_slope)
    assert np.std(slopes) == pytest.approx(us[0], rel=0.05)


def test_general_linear_fit_recovers_a_quadratic():
    x = np.linspace(-1, 1, 9)
    y = 1 - 2 * x + 3 * x**2
    A = np.column_stack([np.ones_like(x), x, x**2])
    beta, cov, chi2, dof = ep.general_linear_fit(A, y, sigma=np.full(9, 0.1))
    assert np.allclose(beta, [1, -2, 3])
    assert chi2 == pytest.approx(0.0, abs=1e-18)
    assert dof == 6 and cov.shape == (3, 3)


# ------------------------------------------------------------------- reporting
@pytest.mark.parametrize(
    "value, u, kwargs, expected",
    [
        (9.81234, 0.0234, dict(unit="m/s^2"), "9.812 ± 0.023 m/s^2"),
        (9.81234, 0.0234, dict(sig=1), "9.81 ± 0.02"),
        (8540.3, 36.4, dict(unit="kg/m^3"), "8540 ± 36 kg/m^3"),
        (1234.5, 123.0, dict(sig=1), "1200 ± 100"),
        (0.50321, 0.0996, dict(sig=1), "0.5 ± 0.1"),
        (1.3456e-4, 2.1e-6, dict(unit="1/s"), "(1.346 ± 0.021) × 10^-4 1/s"),
    ],
)
def test_format_result(value, u, kwargs, expected):
    assert ep.format_result(value, u, **kwargs) == expected


def test_format_concise():
    assert ep.format_concise(9.81234, 0.0234) == "9.812(23)"
    assert ep.format_concise(1.20, 0.05, sig=1) == "1.20(5)"


def test_zero_uncertainty_is_rejected():
    with pytest.raises(ValueError):
        ep.format_result(1.0, 0.0)


# ----------------------------------------------------------------- sensitivity
def test_elasticities_of_a_power_law_are_the_exponents():
    f = lambda a, b, c: 3.0 * a**2 * b**-0.5 / c
    assert np.allclose(ep.elasticities(f, [1.3, 2.7, 0.8]), [2, -0.5, -1], atol=1e-8)


def test_sobol_indices_of_the_ishigami_function():
    sample = lambda n, rng: rng.uniform(-np.pi, np.pi, (n, 3))
    res = ep.sobol_indices(ep.ishigami, sample, n=2**15, rng=11, n_boot=50)
    S, ST = ep.ishigami_exact()
    assert np.allclose(res["S"], S, atol=0.03)
    assert np.allclose(res["ST"], ST, atol=0.03)
    assert np.all(res["S_se"] < 0.03)


def test_sobol_indices_of_a_linear_model_equal_the_budget_shares():
    a = np.array([1.0, -2.0, 0.5])
    s = np.array([0.3, 0.2, 1.0])
    f = lambda x1, x2, x3: a[0] * x1 + a[1] * x2 + a[2] * x3
    sample = lambda n, rng: rng.normal(0, s, (n, 3))
    res = ep.sobol_indices(f, sample, n=2**15, rng=12, n_boot=0)
    shares = (a * s) ** 2 / np.sum((a * s) ** 2)
    assert np.allclose(res["S"], shares, atol=0.02)
    assert np.allclose(res["ST"], shares, atol=0.02)


# ---------------------------------------------------------------- extra tools
def test_poisson_interval_known_values_and_coverage():
    lo, hi = ep.poisson_interval(0, cl=0.90)
    assert lo == 0 and hi == pytest.approx(2.9957, abs=1e-3)      # 5 % one-sided upper limit for n = 0
    lo, hi = ep.poisson_interval(10)
    assert lo < 10 < hi
    mu = 3.7
    n = np.arange(0, 60)
    lo, hi = ep.poisson_interval(n)
    cover = stats.poisson.pmf(n, mu)[(lo <= mu) & (mu <= hi)].sum()
    assert cover >= 0.6827


def test_robust_estimators():
    x = np.random.default_rng(13).normal(5.0, 2.0, 20_000)
    assert ep.mad(x) == pytest.approx(2.0, rel=0.03)
    assert ep.huber_location(x) == pytest.approx(5.0, abs=0.05)
    y = np.concatenate([x[:1000], np.full(100, 100.0)])
    assert abs(ep.huber_location(y) - 5.0) < abs(y.mean() - 5.0) / 10


def test_nonlinear_fit_agrees_with_linear_fit_and_curve_fit():
    from scipy.optimize import curve_fit
    rng = np.random.default_rng(14)
    x = np.linspace(0, 4, 15)
    y = 2.0 + 1.5 * x + rng.normal(0, 0.2, x.size)
    nl = ep.nonlinear_fit(lambda x, a, b: a + b * x, x, y, 0.2, [0, 0])
    lin = ep.linear_fit(x, y, sigma=0.2)
    assert np.allclose(nl.params, [lin.intercept, lin.slope], rtol=1e-6)
    assert np.allclose(nl.cov, lin.cov, rtol=1e-5)
    f = lambda x, A, k: A * np.exp(-k * x)
    y2 = f(x, 5.0, 0.7) + rng.normal(0, 0.1, x.size)
    nl2 = ep.nonlinear_fit(f, x, y2, 0.1, [4, 1])
    popt, pcov = curve_fit(f, x, y2, p0=[4, 1], sigma=np.full(x.size, 0.1), absolute_sigma=True)
    assert np.allclose(nl2.params, popt, rtol=1e-5)
    assert np.allclose(nl2.cov, pcov, rtol=1e-3)


def test_richardson_and_observed_order():
    exact = np.e - 1
    def trap(n):
        x = np.linspace(0, 1, n + 1)
        y = np.exp(x)
        return (x[1] - x[0]) * (y.sum() - 0.5 * (y[0] + y[-1]))
    f1, f2, f3 = trap(64), trap(32), trap(16)
    assert ep.observed_order(f1, f2, f3) == pytest.approx(2.0, abs=0.01)
    fx, err = ep.richardson(f1, f2, p=2)
    assert abs(fx - exact) < 1e-3 * abs(f1 - exact)
    assert err == pytest.approx(abs(f1 - exact), rel=0.01)
