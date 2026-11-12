# Changelog

## 0.3.0 (2026-11-12)

The repository is reorganised as a book in five parts with 26 chapters, and
grows from 13 to 26 notebooks.

### New chapters

- **2. Measuring instruments**: static characteristics; resolution and
  dithering; calibration and accuracy specifications; first-order dynamic
  response and lag correction.
- **Part II, Probability and statistics** (Chapters 5–10): distributions and
  the central limit theorem; counting statistics; estimators and maximum
  likelihood; hypothesis tests; outliers and robust estimation; Monte Carlo
  methods.
- **14. Nonlinear least squares** and **15. Fitting counts: maximum
  likelihood**.
- **19. Numerical error**.
- Worked examples **24. Half-life from decay counts**, **25. Calibration curve
  and detection limit** and **26. Tolerance stack-up and inspection**.

### Other changes

- Every chapter carries a *Track* (core or advanced) and a *Reading* line
  pointing to Taylor, Bevington & Robinson, Fornasini, Rouaud and the GUM.
  New pages `docs/reading_guide.md` and `docs/references.md` collect them.
- References now cite books and published articles first. Chapter 16 (the
  former Notebook 5) maps each rule to its textbook source and keeps course
  notes only for rules specific to them; the Sokal lectures are cited as the
  published book chapter.
- Navigation links (previous/next) at the end of every notebook.
- New README in the style of the author's other repositories, with an
  animation, an overview figure and a gallery (`scripts/make_readme_*.py`).
- The book can be built with Jupyter Book 2 (`myst.yml`) and is published to
  GitHub Pages by `.github/workflows/book.yml`.
- `errprop` 0.3.0 adds `nonlinear_fit`, `poisson_interval`, `mad`,
  `huber_location`, `observed_order` and `richardson`, with tests.

### Renamed notebooks

| 0.2 | 0.3 |
|---|---|
| `core/01_measurement_and_uncertainty` | `part1_measurement/01_measurement_and_uncertainty` |
| `core/02_propagation_of_uncertainty` | `part1_measurement/03_propagation_of_uncertainty` |
| `core/03_repeated_measurements` | `part1_measurement/04_repeated_measurements` |
| `core/04_fitting_a_straight_line` | `part3_data_analysis/11_fitting_a_straight_line` |
| `core/05_rules_of_thumb_tested` | `part3_data_analysis/16_rules_of_thumb_tested` |
| `advanced/A1_covariance_and_limits_of_linearisation` | `part3_data_analysis/12_covariance_and_limits_of_linearisation` |
| `advanced/A3_fit_parameter_covariance` | `part3_data_analysis/13_fit_parameter_covariance` |
| `advanced/A2_correlated_data_and_block_averaging` | `part4_advanced/17_correlated_data_and_block_averaging` |
| `advanced/A4_sensitivity_analysis_and_budgets` | `part4_advanced/18_sensitivity_analysis_and_budgets` |
| `examples/E1_pendulum_g` | `part5_examples/20_pendulum_g` |
| `examples/E2_falling_sphere_viscosity` | `part5_examples/21_falling_sphere_viscosity` |
| `examples/E3_arrhenius_activation_energy` | `part5_examples/22_arrhenius_activation_energy` |
| `examples/E4_md_lennard_jones` | `part5_examples/23_md_lennard_jones` |

Section numbers inside the reused chapters follow the new chapter numbers.
The numerical results of the reused chapters are unchanged.

## 0.2.0

- Added "Rules of thumb, tested" (Notebook 5) and "Sensitivity analysis and
  uncertainty budgets" (A4); `errprop.sensitivity`.

## 0.1.0

- First release: core notebooks 1–4, advanced A1–A3, examples E1–E4 and the
  `errprop` module.
