# Changelog

## v0.4.0 — 2026-09-21
### Changed
- **A1C mu bounds.** Before, the bounds came from the reference interval (4.0–5.6), which gave (1.6, 8.0). Patients whose A1C setpoint was above 8% could not be represented. The bounds are now the empirical 0.5th/99.5th percentiles of isolated A1C measurements:

  | sex | (min_mu, max_mu) | n measurements |
  |-----|------------------|----------------|
  | ALL | (4.5, 12.8) | 247,655 |
  | F   | (4.5, 12.8) | 116,364 |
  | M   | (4.5, 12.9) | 131,291 |

  With grid size 201, the mu grid spacing goes from about 0.032 to about 0.0415.
- **A1C `log_lambda_`**, re-optimized under the new bounds: ALL 1.4827 (was 0.9950), F 1.5791 (was 1.0735), M 1.2502 (was 0.9690). bsi fits A1C unstratified, so use `sex="ALL"` to match it.
- A1C `min_sigma`/`max_sigma` (0.001 / 0.408), the grid size and the flat prior are unchanged. All other markers are unchanged.

### Added
- `perri/data/mu_bound_overrides.csv` is a per-marker, per-sex table of mu bounds. Its bounds are used directly instead of the RI-derived ones in `bayesian_hyperparameters.csv`. Each row records its provenance (method, percentiles, n, data version, source). Rows are validated on load: missing, non-finite or non-increasing bounds, unknown markers and duplicate rows raise `ValueError`.
- A pytest suite (`tests/`) using only synthetic data. Install with `pip install -e ".[dev]"`.
