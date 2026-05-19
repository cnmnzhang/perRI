# PerRI — Personalized Reference Intervals

Code accompanying the paper. Fits a Bayesian model to an individual's longitudinal lab measurements to produce a personalized reference interval (setpoint) for 43 common lab markers.

---

## Requirements

Python 3.9+. Install dependencies:

Method 1. Installs the package in editable mode with all dependencies from `pyproject.toml`. 
```bash
  python3 -m venv .venv                                                         
  source .venv/bin/activate
  pip install -e .                                                              
```          
If you also want the Streamlit app:
```bash                                                                
  pip install -e ".[app]"
```

Or without installing:

```bash
pip install numpy pandas scipy numba matplotlib
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

**List all supported markers:**

```bash
python -c "import perri as pri; print(pri.list_supported_markers())"
```

**Interactive Streamlit app** (requires `pip install streamlit`):

```bash
streamlit run perri/app.py
```

---

## Data Format

To run on your own data, provide a CSV with columns:

```
patient_id, test_code, result_value, result_date
```

where `result_date` is ISO 8601 (`YYYY-MM-DD`) and `test_code` is one of the 43 supported marker codes.
