# Reading guide

Each chapter of this book is self-contained, but none replaces a good
textbook. The table maps every chapter to the corresponding chapters or
sections of four textbooks and to the clauses of the *Guide to the
Expression of Uncertainty in Measurement* (GUM). A blank cell means the
source does not treat the topic. The chapters themselves list further
books and articles; all of them are collected in [References](references.md).

| Chapter | Taylor | Bevington & Robinson | Fornasini | Rouaud | GUM / JCGM |
|---|---|---|---|---|---|
| [1. Measurements, errors and how to report a result](../notebooks/part1_measurement/01_measurement_and_uncertainty.ipynb) | 1–2 | 1 | 1, 4 | I.A–C, I.I | 2–4 |
| [2. Measuring instruments](../notebooks/part1_measurement/02_measuring_instruments.ipynb) | 1.5 |  | 3, 4.2, D.1, E.7 |  | F.2.2 |
| [3. Propagation of uncertainty](../notebooks/part1_measurement/03_propagation_of_uncertainty.ipynb) | 3 | 3 | 8.1–8.3 | II.B | 5.1 |
| [4. Repeated measurements](../notebooks/part1_measurement/04_repeated_measurements.ipynb) | 4, 5.7, 7 | 1.3–1.4, 4.1–4.2 | 7, 9.2 | I.E | 4.2, G |
| [5. Probability distributions and the central limit theorem](../notebooks/part2_probability_statistics/05_probability_distributions.ipynb) | 5, 10, 11 | 2 | 5, 6 | I.D–F, III | C |
| [6. Counting statistics](../notebooks/part2_probability_statistics/06_counting_statistics.ipynb) | 3.2, 11 | 2.2 | 1.5, 6.4 | III.A |  |
| [7. Estimators and maximum likelihood](../notebooks/part2_probability_statistics/07_estimators_and_maximum_likelihood.ipynb) | 5.5, 7 | 4.1, 10 | 7.3 | IV |  |
| [8. Hypothesis tests and goodness of fit](../notebooks/part2_probability_statistics/08_hypothesis_tests.ipynb) | 5.8, 9, 12 | 4.4, 11 | 9–11 | I.G–H |  |
| [9. Outliers and the rejection of data](../notebooks/part2_probability_statistics/09_outliers_and_data_rejection.ipynb) | 6 |  | 9.3 |  |  |
| [10. Monte Carlo methods](../notebooks/part2_probability_statistics/10_monte_carlo_methods.ipynb) |  | 5 |  | III.D | JCGM 101 |
| [11. Fitting a straight line](../notebooks/part3_data_analysis/11_fitting_a_straight_line.ipynb) | 8 | 6 | 10.3–10.4 | II.C |  |
| [12. Covariance and the limits of linearisation](../notebooks/part3_data_analysis/12_covariance_and_limits_of_linearisation.ipynb) | 3.11, 9 | 3.2 | 8.4 | II.B | 5.2; JCGM 101, 102 |
| [13. Fit parameters and their covariance](../notebooks/part3_data_analysis/13_fit_parameter_covariance.ipynb) |  | 6–7 | 10.4 | II.C–D | H.3 |
| [14. Nonlinear least squares](../notebooks/part3_data_analysis/14_nonlinear_least_squares.ipynb) |  | 8, 9, 11.5 |  | II.D |  |
| [15. Fitting counts: maximum likelihood](../notebooks/part3_data_analysis/15_maximum_likelihood_fitting.ipynb) |  | 10 |  | IV.B |  |
| [16. Rules of thumb, tested](../notebooks/part3_data_analysis/16_rules_of_thumb_tested.ipynb) | 2–3 | 3 |  |  |  |
| [17. Correlated data and block averaging](../notebooks/part4_advanced/17_correlated_data_and_block_averaging.ipynb) |  |  |  |  | 4.2.7 |
| [18. Sensitivity analysis and uncertainty budgets](../notebooks/part4_advanced/18_sensitivity_analysis_and_budgets.ipynb) |  | 3 |  |  | 5.1.3–5.1.4 |
| [19. Numerical error](../notebooks/part4_advanced/19_numerical_error.ipynb) |  | A.3 |  |  |  |
| [20. Physics: the acceleration of gravity from a pendulum](../notebooks/part5_examples/20_pendulum_g.ipynb) | 3, 8 |  | E.2, E.5–E.6 |  |  |
| [21. Chemical engineering: viscosity from a falling sphere](../notebooks/part5_examples/21_falling_sphere_viscosity.ipynb) | 3 |  |  |  | Annex H |
| [22. Chemistry: activation energy from an Arrhenius plot](../notebooks/part5_examples/22_arrhenius_activation_energy.ipynb) | 8.6 | 6 |  | II.C |  |
| [23. Molecular simulation: pressure and diffusion in a Lennard-Jones liquid](../notebooks/part5_examples/23_md_lennard_jones.ipynb) |  |  |  |  |  |
| [24. Nuclear physics: a half-life from decay counts](../notebooks/part5_examples/24_radioactive_decay_counts.ipynb) | 8.6, 11.4 | 8 (Example 8.1) | 6.4 |  |  |
| [25. Analytical chemistry: a calibration curve and its detection limit](../notebooks/part5_examples/25_calibration_curve_detection_limit.ipynb) | 8 | 6 |  | II.C | H.3 |
| [26. Mechanical engineering: tolerance stack-up and inspection](../notebooks/part5_examples/26_tolerance_stack_up.ipynb) | 3.3–3.6 |  | E.1 | V | 5.1, 6 |

Numbers are chapters or sections (Taylor, Bevington & Robinson, Fornasini),
chapter and section letters (Rouaud), or clauses and annexes (GUM).

## The four textbooks

* **J. R. Taylor**, *An Introduction to Error Analysis: The Study of
  Uncertainties in Physical Measurements*, 2nd ed. (University Science Books,
  1997). The standard first course: clear, gentle and full of worked
  examples. Best companion to Parts I and II.
* **P. R. Bevington and D. K. Robinson**, *Data Reduction and Error Analysis
  for the Physical Sciences*, 3rd ed. (McGraw-Hill, 2003). Distributions,
  propagation and, above all, fitting: linear, nonlinear and maximum
  likelihood. Best companion to Part III.
* **P. Fornasini**, *The Uncertainty in Physical Measurements: An
  Introduction to Data Analysis in the Physics Laboratory* (Springer, 2008).
  Strong on instruments, measurement units and the probabilistic foundations,
  with laboratory experiments in its appendices.
* **M. Rouaud**, *Probability, Statistics and Estimation: Propagation of
  Uncertainties in Experimental Measurement* (2017; freely available under a
  CC BY-NC 4.0 licence). A compact text with many exercises. It is cited here
  for further reading only; none of its material is reproduced.
* **JCGM 100:2008**, *Evaluation of measurement data — Guide to the
  expression of uncertainty in measurement* (GUM), with its supplements
  JCGM 101:2008 (Monte Carlo propagation) and JCGM 102:2011 (multivariate
  models). The international standard for terminology and procedure; free
  from the BIPM.

## Suggested routes

* **A first course** (undergraduate laboratory): Chapters 1–4, 11 and 16,
  with Examples 20 and 25. Read alongside Taylor.
* **Statistics for experimenters**: Part II (Chapters 5–10), then
  Chapters 14 and 15, with Example 24.
* **Fitting and data analysis** (graduate): Chapters 11–15 and 17, with
  Examples 22 and 24. Read alongside Bevington & Robinson.
* **Engineering uncertainty analysis**: Chapters 2, 3, 12, 18 and 19, with
  Examples 21 and 26.
* **Computational science**: Chapters 10, 17 and 19, with Example 23.
