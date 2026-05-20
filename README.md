# PerRI — Personalized Reference Intervals

Code accompanying the paper. Fits a Bayesian model to an individual's longitudinal lab measurements to produce a personalized reference interval (setpoint) for 43 common lab markers:
`['A1C', 'ALB', 'ALK', 'ALT', 'AST', 'BIL', 'BILD', 'BUN', 'CA', 'CHOL', 'CL', 'CO2', 'CRE', 'FER', 'GLU', 'HB', 'HCT', 'HDL', 'HSCRP', 'IGAP', 'K', 'LD', 'LDL', 'LYMPH', 'MCH', 'MCHC', 'MCV', 'MG', 'MONOC', 'NA', 'NONHDL', 'P', 'PLT', 'PROINR', 'PROPAT', 'RBC', 'RDWCV', 'TNEUT', 'TP', 'TRIG', 'TSH', 'VITDT', 'WBC']`

By default, measurements within 90 days of any neighbour are filtered out before fitting to avoid clusters of repeat draws (e.g. serial phlebotomy, acute illness) biasing the setpoint estimate. This can be disabled by passing `filter_isolated_measurements=False` to `fit_patient()` or `fit_batch()`, and the gap threshold adjusted via `min_gap_days`.

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
| `--test_code` | `HB` | Lab marker (see supported list below) |
| `--patient_id` | `3` | Patient ID — single-patient script only |
| `--sex` | `ALL` | `ALL`, `M`, or `F` |

---

## Custom Hyperparameters

To override the optimised defaults, pass a `params` dict to `fit_patient()` or `fit_batch()`:

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

**Interactive Streamlit app** (requires `pip install streamlit plotly`):

If you also want the Streamlit app:
```bash     
pip install -e ".[app]"
streamlit run perri/app.py
```