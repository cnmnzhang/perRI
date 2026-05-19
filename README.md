# PerRI — Personalized Reference Intervals

Code accompanying the paper. Fits a Bayesian model to an individual's longitudinal lab measurements to produce a personalized reference interval (setpoint) for 43 common lab markers.

---

## Requirements

Python 3.9+. 

To install the package in editable mode with all dependencies from `pyproject.toml`:
```bash
  python3 -m venv .venv                                                         
  source .venv/bin/activate
  pip install -e .                                                              
```          
If you also want the Streamlit app:
```bash                                                                
  pip install -e ".[app]"
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

**Interactive Streamlit app** (requires `pip install streamlit plotly`):

```bash
streamlit run perri/app.py
```
'
Supports
`['A1C', 'ALB', 'ALK', 'ALT', 'AST', 'BIL', 'BILD', 'BUN', 'CA', 'CHOL', 'CL', 'CO2', 'CRE', 'FER', 'GLU', 'HB', 'HCT', 'HDL', 'HSCRP', 'IGAP', 'K', 'LD', 'LDL', 'LYMPH', 'MCH', 'MCHC', 'MCV', 'MG', 'MONOC', 'NA', 'NONHDL', 'P', 'PLT', 'PROINR', 'PROPAT', 'RBC', 'RDWCV', 'TNEUT', 'TP', 'TRIG', 'TSH', 'VITDT', 'WBC']`