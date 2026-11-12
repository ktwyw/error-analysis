# Error analysis and propagation

*Evaluating, propagating and reporting uncertainty, from the first ± to the
covariance matrix of a fit, with examples across disciplines.*

Every number in science and engineering comes with an uncertainty, yet the
subject is usually taught once, as a set of recipes in a first-year
laboratory, and seldom carried into the places where it is needed later:
curve fitting, counting experiments, correlated simulation data, numerical
computation, and the systematic effects that no formula for random error
will reveal. This book treats error analysis as one connected subject.

Each chapter is an executed Jupyter notebook. Every idea is derived,
demonstrated on simulated data for which the true answer is known, and
checked by simulation, so that every error bar can be compared with the
truth. The notebooks share a small Python module, `errprop`, written to be
read as much as to be used.

* **Core** chapters assume partial derivatives and a first look at Python.
* **Advanced** chapters are for graduate students and researchers who fit
  data, run simulations or publish numbers with error bars.

The [reading guide](reading_guide.md) maps every chapter to the standard
textbooks and to the GUM, and suggests routes through the material.

### Part I. Measurement and uncertainty

| Chapter | Track | Topics |
|---|---|---|
| [1. Measurements, errors and how to report a result](../notebooks/part1_measurement/01_measurement_and_uncertainty.ipynb) | core | error vs uncertainty; random and systematic effects; Type A and Type B evaluation; quadrature; rounding and reporting; comparing results |
| [2. Measuring instruments](../notebooks/part1_measurement/02_measuring_instruments.ipynb) | core | static characteristics; resolution, noise and dithering; calibration and accuracy specifications; first-order dynamic response and lag correction |
| [3. Propagation of uncertainty](../notebooks/part1_measurement/03_propagation_of_uncertainty.ipynb) | core | the propagation law; rules for sums, products and powers; uncertainty budgets; Monte Carlo check; differences of nearly equal numbers |
| [4. Repeated measurements](../notebooks/part1_measurement/04_repeated_measurements.ipynb) | core | standard deviation and standard error; the uncertainty of the uncertainty; Student's *t*; what averaging does not cure; weighted mean |

### Part II. Probability and statistics

| Chapter | Track | Topics |
|---|---|---|
| [5. Probability distributions and the central limit theorem](../notebooks/part2_probability_statistics/05_probability_distributions.ipynb) | core | binomial, Poisson, normal, uniform, exponential and Cauchy distributions; the central limit theorem; χ², *t* and *F* |
| [6. Counting statistics](../notebooks/part2_probability_statistics/06_counting_statistics.ipynb) | core | the square-root rule and the dispersion test; background subtraction and counting-time split; Garwood intervals for small counts; dead time |
| [7. Estimators and maximum likelihood](../notebooks/part2_probability_statistics/07_estimators_and_maximum_likelihood.ipynb) | advanced | bias, variance and mean-square error; efficient and inefficient estimators; maximum likelihood; likelihood intervals; Cramér–Rao bound |
| [8. Hypothesis tests and goodness of fit](../notebooks/part2_probability_statistics/08_hypothesis_tests.ipynb) | core | p-values, power and errors of both kinds; *z*, *t* and Welch tests; χ² test of a distribution; *F* test; correlation coefficient; many tests at once |
| [9. Outliers and the rejection of data](../notebooks/part2_probability_statistics/09_outliers_and_data_rejection.ipynb) | core | Chauvenet's criterion tested; robust estimators (median, trimmed mean, Huber); the bootstrap; a policy for suspect data |
| [10. Monte Carlo methods](../notebooks/part2_probability_statistics/10_monte_carlo_methods.ipynb) | advanced | sampling by inverse transform and rejection; correlated inputs; the accuracy of a Monte Carlo result; Latin hypercube sampling |

### Part III. Data analysis

| Chapter | Track | Topics |
|---|---|---|
| [11. Fitting a straight line](../notebooks/part3_data_analysis/11_fitting_a_straight_line.ipynb) | core | least squares and its parameter uncertainties; residuals; weighted fits and χ²; coverage check; pitfalls |
| [12. Covariance and the limits of linearisation](../notebooks/part3_data_analysis/12_covariance_and_limits_of_linearisation.ipynb) | advanced | the general law **C**<sub>y</sub> = **J C**<sub>x</sub> **J**<sup>T</sup>; shared errors; staged calculations; second-order bias; Monte Carlo propagation |
| [13. Fit parameters and their covariance](../notebooks/part3_data_analysis/13_fit_parameter_covariance.ipynb) | advanced | matrix form of least squares; correlated slope and intercept; centring; transformed data and weights; absolute vs relative error bars; errors in both variables |
| [14. Nonlinear least squares](../notebooks/part3_data_analysis/14_nonlinear_least_squares.ipynb) | advanced | Gauss–Newton and Levenberg–Marquardt; local minima and starting values; covariance at the minimum; profile Δχ² intervals; derived quantities; coverage |
| [15. Fitting counts: maximum likelihood](../notebooks/part3_data_analysis/15_maximum_likelihood_fitting.ipynb) | advanced | Neyman vs Pearson χ² vs Poisson likelihood for histograms; goodness of fit with few counts; unbinned likelihood with a finite window |
| [16. Rules of thumb, tested](../notebooks/part3_data_analysis/16_rules_of_thumb_tested.ipynb) | core | nine rules from course guides tested by simulation: worst case vs quadrature, overlapping error bars, limits of linearisation, the dependency trap, measuring many at once, dropping small terms, significant figures, offsets and gains, cancelling errors |

### Part IV. Advanced topics

| Chapter | Track | Topics |
|---|---|---|
| [17. Correlated data and block averaging](../notebooks/part4_advanced/17_correlated_data_and_block_averaging.ipynb) | advanced | variance of the mean of a time series; integrated autocorrelation time; blocking analysis; coverage test; jackknife for derived quantities |
| [18. Sensitivity analysis and uncertainty budgets](../notebooks/part4_advanced/18_sensitivity_analysis_and_budgets.ipynb) | advanced | sensitivity, elasticity and importance; budgets across operating points; bottom-up vs top-down; Sobol indices; aleatory vs epistemic uncertainty; model-form error |
| [19. Numerical error](../notebooks/part4_advanced/19_numerical_error.ipynb) | advanced | floating point and cancellation; Welford's algorithm; conditioning; step size of numerical derivatives; Richardson extrapolation; long-time integration |

### Part V. Worked examples

| Chapter | Track | Topics |
|---|---|---|
| [20. Physics: the acceleration of gravity from a pendulum](../notebooks/part5_examples/20_pendulum_g.ipynb) | core | physics: a systematic error that shows up as a trend; a fit that absorbs it; a correction common to all points |
| [21. Chemical engineering: viscosity from a falling sphere](../notebooks/part5_examples/21_falling_sphere_viscosity.ipynb) | core | chemical engineering: Type B evaluation; a full budget; a model error fourteen times the budget; defining the measurand |
| [22. Chemistry: activation energy from an Arrhenius plot](../notebooks/part5_examples/22_arrhenius_activation_energy.ipynb) | advanced | chemistry: parameters with correlation 0.999; prediction with and without covariance; the statistical compensation effect |
| [23. Molecular simulation: pressure and diffusion in a Lennard-Jones liquid](../notebooks/part5_examples/23_md_lennard_jones.ipynb) | advanced | molecular simulation: correlated time series; block averaging on real MD data; transport coefficients; what an error bar leaves out |
| [24. Nuclear physics: a half-life from decay counts](../notebooks/part5_examples/24_radioactive_decay_counts.ipynb) | advanced | nuclear physics: why the log-linear fit is biased; background in the likelihood; checking four methods by simulation |
| [25. Analytical chemistry: a calibration curve and its detection limit](../notebooks/part5_examples/25_calibration_curve_detection_limit.ipynb) | core | analytical chemistry: inverse prediction and its uncertainty; decision threshold and detection limit (Currie); the 3s rule tested |
| [26. Mechanical engineering: tolerance stack-up and inspection](../notebooks/part5_examples/26_tolerance_stack_up.ipynb) | core | mechanical engineering: worst case vs RSS vs Monte Carlo tolerancing; yield and C<sub>pk</sub>; the Abbe error; guard bands at inspection |

## Licence and citation

Code: MIT. Text, figures and data: CC BY 4.0. Please cite the material
using the metadata in `CITATION.cff` in the repository.
