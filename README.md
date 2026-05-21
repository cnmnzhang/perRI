# PerRI — Personalized Reference Intervals

Code accompanying the paper: *Personalized clinical reference intervals for routine precision medical care*

Built on University of Washington Medicine (UWM) data. All default hyperparameters, units, and reference ranges are calibrated to UWM laboratory standards. 

> **Note on units and reference ranges:** Default hyperparameters assume UWM laboratory units and reference ranges ([UW Lab Test Guide](https://testguide.labmed.uw.edu/)). If your institution uses different units or reference ranges, update `bayesian_hyperparameters.csv` to reflect own lab's reference ranges. 

---

## What does it do?

Standard clinical reference intervals are derived from a healthy reference population and applied uniformly to all patients. PerRI takes a different approach: for any laboratory marker that follows homeostatic dynamics, it fits a Bayesian model to an individual's longitudinal lab measurements to estimate that patient's personal physiological setpoint as a normal distribution N(μ, σ).

The model starts from a population-level prior derived from the standard reference range for the marker, then refines its estimate with each new observation. At any time point, the personalized reference interval is the 95% credible interval: **μ ± 1.96σ**. The result is a reference interval specific to the individual — narrower for patients with stable physiology, and appropriately wider for those with greater intra-individual variability.

---

### Inputs

The model takes an ordered sequence of measurement values and timestamps for a single patient and marker. Five hyperparameters control the search space and temporal weighting:

| Parameter | Description |
|-----------|-------------|
| `log_lambda_` | Log of the temporal decay rate. Higher values give more weight to recent measurements over older ones. Optimized per marker, procedure mentioned in the paper |
| `min_mu`, `max_mu` | Search grid bounds for the setpoint μ. Typically set as: **population RI center ± 2.5 × population RI range**. |
| `min_sigma`, `max_sigma` | Search grid bounds for within-person SD σ. |

Pre-optimized values for all 43 supported markers are in `perri/data/bayesian_hyperparameters.csv`, stratified by sex (`ALL`, `M`, `F`).

### Supported Markers

43 common lab markers:

`['A1C', 'ALB', 'ALK', 'ALT', 'AST', 'BIL', 'BILD', 'BUN', 'CA', 'CHOL', 'CL', 'CO2', 'CRE', 'FER', 'GLU', 'HB', 'HCT', 'HDL', 'HSCRP', 'IGAP', 'K', 'LD', 'LDL', 'LYMPH', 'MCH', 'MCHC', 'MCV', 'MG', 'MONOC', 'NA', 'NONHDL', 'P', 'PLT', 'PROINR', 'PROPAT', 'RBC', 'RDWCV', 'TNEUT', 'TP', 'TRIG', 'TSH', 'VITDT', 'WBC']`

### Measurement Selection

By default, measurements within 90 days of any neighbour are filtered out before fitting to avoid clusters of repeat draws (e.g. serial phlebotomy, acute illness) biasing the setpoint estimate. This can be disabled by passing `filter_isolated_measurements=False` to `fit_patient()` or `fit_batch()`, and the gap threshold can be adjusted via `min_gap_days`.

### Scoring
We report RMSE, KS statistic, and KS quintile range. To enable comparison across markers with different scales, RMSE is also reported as RMSE_Scaled, which is RMSE normalized by the median intra-patient standard deviation. 
The theoretical lower bound for RMSE Scaled is 1. These are pre-calculated from the UW cohort and stored in `perri/data/marker_intra_patient_std.csv`.

---

## Requirements

Python 3.9+.

To install the package in editable mode with all dependencies from `pyproject.toml`:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

---

## Running the Examples

Sample data (synthetic, 3 patients × 43 markers × 12 visits) is bundled at `perri/data/sample_data.csv`.

**Single-patient fit:**

```bash
python perri/example_single.py --test_code HB --patient_id 3
```

Outputs the estimated setpoint, within-person SD, personalized 95% reference interval, and saves a plot to the `perri/` directory.

**Batch fit + cohort evaluation:**

```bash
python perri/example_batch.py --test_code HB
```

Outputs RMSE, KS statistic, and KS quintile range across all patients for the specified marker.

Optional flags for both scripts:

| Flag | Default | Description |
|------|---------|-------------|
| `--test_code` | `HB` | Lab marker (see supported list above) |
| `--patient_id` | `3` | Patient ID — single-patient script only |
| `--sex` | `ALL` | `ALL`, `M`, or `F` |

---

## Custom Hyperparameters

To override the optimized defaults, pass a `params` dict to `fit_patient()` or `fit_batch()`:

```python
import perri as pri

# Start from the defaults and tweak one value
params = pri.get_default_params("HB")
params["log_lambda_"] = 0.5  # stronger recency weighting

result = pri.fit_patient(
    values=patient_df["result_value"].values,
    timestamps=patient_df["result_date"].values,
    params=params,
)
```

Or specify all five parameters from scratch:

```python
result = pri.fit_patient(
    values=patient_df["result_value"].values,
    timestamps=patient_df["result_date"].values,
    params={
        "log_lambda_": 0.5,
        "min_mu": 5.0,
        "max_mu": 20.0,
        "min_sigma": 0.01,
        "max_sigma": 2.0,
    },
)
```

When `params` is provided, `test_code` is not required for parameter loading (but can still be passed for display purposes).

---

## Interactive Streamlit App

Requires `pip install streamlit plotly`.

```bash
pip install -e ".[app]"
streamlit run perri/app.py
```
