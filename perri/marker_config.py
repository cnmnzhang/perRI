"""
Standalone marker metadata for personalized_reference_intervals.

Covers all 43 supported markers: battery grouping, display units, human-readable
name, and sex-specific population reference interval (pop_ri).

pop_ri values are (lower, upper) tuples. HDL has a one-sided lower bound
(upper = float("inf")). ALK's age-stratified RI is collapsed to the adult
range spanning ages 18–75.

No imports from the parent research repo.
"""

MARKER_CONFIG: dict = {
    # ── CBC ─────────────────────────────────────────────────────────────────
    "HB":    {"battery": "CBC",     "units": "g/dL",      "full_name": "Hemoglobin",
              "pop_ri": {"F": (11.5, 15.5), "M": (13.0, 18.0)}},
    "HCT":   {"battery": "CBC",     "units": "%",         "full_name": "Hematocrit",
              "pop_ri": {"F": (36.0, 45.0), "M": (38.0, 50.0)}},
    "RBC":   {"battery": "CBC",     "units": "x10⁶/μL", "full_name": "Red Blood Cell Count",
              "pop_ri": {"F": (3.8, 6.0),  "M": (4.4, 5.6)}},
    "MCH":   {"battery": "CBC",     "units": "pg",        "full_name": "Mean Corpuscular Hemoglobin",
              "pop_ri": {"F": (27.3, 33.6), "M": (27.3, 33.6)}},
    "MCHC":  {"battery": "CBC",     "units": "g/dL",      "full_name": "Mean Corpuscular Hemoglobin Concentration",
              "pop_ri": {"F": (32.2, 36.5), "M": (32.2, 36.5)}},
    "MCV":   {"battery": "CBC",     "units": "fL",        "full_name": "Mean Corpuscular Volume",
              "pop_ri": {"F": (81.0, 98.0), "M": (81.0, 98.0)}},
    "RDWCV": {"battery": "CBC",     "units": "%",         "full_name": "Red Cell Distribution Width (CV)",
              "pop_ri": {"F": (10.0, 14.5), "M": (10.0, 14.5)}},
    "PLT":   {"battery": "CBC",     "units": "x10³/μL", "full_name": "Platelets",
              "pop_ri": {"F": (150.0, 400.0), "M": (150.0, 400.0)}},
    "WBC":   {"battery": "CBC",     "units": "x10³/μL", "full_name": "White Blood Cell Count",
              "pop_ri": {"F": (4.3, 10.0), "M": (4.3, 10.0)}},

    # ── WBC Differential ────────────────────────────────────────────────────
    "TNEUT": {"battery": "WBC diff", "units": "x10³/μL", "full_name": "Total Neutrophils",
              "pop_ri": {"F": (1.80, 7.0), "M": (1.80, 7.0)}},
    "LYMPH": {"battery": "WBC diff", "units": "x10³/μL", "full_name": "Lymphocytes",
              "pop_ri": {"F": (1.00, 4.80), "M": (1.00, 4.80)}},
    "MONOC": {"battery": "WBC diff", "units": "x10³/μL", "full_name": "Monocytes",
              "pop_ri": {"F": (0.00, 0.80), "M": (0.00, 0.80)}},

    # ── Basic Metabolic Panel ────────────────────────────────────────────────
    "NA":    {"battery": "BMP",     "units": "mEq/L",     "full_name": "Sodium",
              "pop_ri": {"F": (135.0, 145.0), "M": (135.0, 145.0)}},
    "K":     {"battery": "BMP",     "units": "mEq/L",     "full_name": "Potassium",
              "pop_ri": {"F": (3.6, 5.2), "M": (3.6, 5.2)}},
    "CL":    {"battery": "BMP",     "units": "mEq/L",     "full_name": "Chloride",
              "pop_ri": {"F": (98.0, 108.0), "M": (98.0, 108.0)}},
    "CO2":   {"battery": "BMP",     "units": "mEq/L",     "full_name": "Carbon Dioxide (Bicarbonate)",
              "pop_ri": {"F": (22.0, 32.0), "M": (22.0, 32.0)}},
    "IGAP":  {"battery": "BMP",     "units": "mEq/L",     "full_name": "Anion Gap",
              "pop_ri": {"F": (4.0, 12.0), "M": (4.0, 12.0)}},
    "GLU":   {"battery": "BMP",     "units": "mg/dL",     "full_name": "Glucose",
              "pop_ri": {"F": (62.0, 125.0), "M": (62.0, 125.0)}},
    "BUN":   {"battery": "BMP",     "units": "mg/dL",     "full_name": "Blood Urea Nitrogen",
              "pop_ri": {"F": (8.0, 21.0), "M": (8.0, 21.0)}},
    "CRE":   {"battery": "BMP",     "units": "mg/dL",     "full_name": "Creatinine",
              "pop_ri": {"F": (0.38, 1.02), "M": (0.51, 1.18)}},
    "CA":    {"battery": "BMP",     "units": "mg/dL",     "full_name": "Calcium",
              "pop_ri": {"F": (8.9, 10.2), "M": (8.9, 10.2)}},

    # ── Hepatic Panel ────────────────────────────────────────────────────────
    "ALB":   {"battery": "Hepatic", "units": "g/dL",      "full_name": "Albumin",
              "pop_ri": {"F": (3.5, 5.2), "M": (3.5, 5.2)}},
    "TP":    {"battery": "Hepatic", "units": "g/dL",      "full_name": "Total Protein",
              "pop_ri": {"F": (6.0, 8.2), "M": (6.0, 8.2)}},
    "AST":   {"battery": "Hepatic", "units": "U/L",       "full_name": "Aspartate Aminotransferase",
              "pop_ri": {"F": (9.0, 33.0), "M": (9.0, 33.0)}},
    "ALT":   {"battery": "Hepatic", "units": "U/L",       "full_name": "Alanine Aminotransferase",
              "pop_ri": {"F": (7.0, 33.0), "M": (10.0, 64.0)}},
    # ALK pop_ri collapsed from age-stratified adult (18–75) brackets:
    #   F: min lower = 25, max upper = 172  |  M: min lower = 35, max upper = 161
    "ALK":   {"battery": "Hepatic", "units": "U/L",       "full_name": "Alkaline Phosphatase",
              "pop_ri": {"F": (25.0, 172.0), "M": (35.0, 161.0)}},
    "BIL":   {"battery": "Hepatic", "units": "mg/dL",     "full_name": "Total Bilirubin",
              "pop_ri": {"F": (0.2, 1.3), "M": (0.2, 1.3)}},
    "BILD":  {"battery": "Hepatic", "units": "mg/dL",     "full_name": "Direct Bilirubin",
              "pop_ri": {"F": (0.0, 0.3), "M": (0.0, 0.3)}},

    # ── Lipid Panel ──────────────────────────────────────────────────────────
    "CHOL":   {"battery": "Lipid",  "units": "mg/dL",     "full_name": "Total Cholesterol",
               "pop_ri": {"F": (0.0, 210.0), "M": (0.0, 210.0)}},
    "TRIG":   {"battery": "Lipid",  "units": "mg/dL",     "full_name": "Triglycerides",
               "pop_ri": {"F": (0.0, 175.0), "M": (0.0, 175.0)}},
    # HDL: one-sided lower bound (higher is better); upper = inf
    "HDL":    {"battery": "Lipid",  "units": "mg/dL",     "full_name": "HDL Cholesterol",
               "pop_ri": {"F": (49.0, float("inf")), "M": (39.0, float("inf"))}},
    "NONHDL": {"battery": "Lipid",  "units": "mg/dL",     "full_name": "Non-HDL Cholesterol",
               "pop_ri": {"F": (0.0, 190.0), "M": (0.0, 190.0)}},
    "LDL":    {"battery": "Lipid",  "units": "mg/dL",     "full_name": "LDL Cholesterol",
               "pop_ri": {"F": (0.0, 130.0), "M": (0.0, 130.0)}},

    # ── Coagulation ──────────────────────────────────────────────────────────
    "PROPAT": {"battery": "Coag",   "units": "sec",       "full_name": "Prothrombin Time",
               "pop_ri": {"F": (10.7, 15.6), "M": (10.7, 15.6)}},
    "PROINR": {"battery": "Coag",   "units": "ratio",     "full_name": "INR",
               "pop_ri": {"F": (0.8, 1.3), "M": (0.8, 1.3)}},

    # ── Miscellaneous ────────────────────────────────────────────────────────
    "HSCRP":  {"battery": "Misc",   "units": "mg/L",      "full_name": "High-Sensitivity C-Reactive Protein",
               "pop_ri": {"F": (0.0, 10.0), "M": (0.0, 10.0)}},
    "LD":     {"battery": "Misc",   "units": "U/L",       "full_name": "Lactate Dehydrogenase",
               "pop_ri": {"F": (0.0, 210.0), "M": (0.0, 210.0)}},
    "MG":     {"battery": "Misc",   "units": "mg/dL",     "full_name": "Magnesium",
               "pop_ri": {"F": (1.8, 2.4), "M": (1.8, 2.4)}},
    "P":      {"battery": "Misc",   "units": "mg/dL",     "full_name": "Phosphate",
               "pop_ri": {"F": (2.5, 4.5), "M": (2.5, 4.5)}},
    "TSH":    {"battery": "Misc",   "units": "mIU/L",     "full_name": "Thyroid-Stimulating Hormone",
               "pop_ri": {"F": (0.4, 5.0), "M": (0.4, 5.0)}},
    "FER":    {"battery": "Misc",   "units": "ng/mL",     "full_name": "Ferritin",
               "pop_ri": {"F": (10.0, 180.0), "M": (20.0, 230.0)}},
    "VITDT":  {"battery": "Misc",   "units": "ng/mL",     "full_name": "Vitamin D (25-OH)",
               "pop_ri": {"F": (20.1, 50.0), "M": (20.1, 50.0)}},
    "A1C":    {"battery": "Misc",   "units": "%",         "full_name": "Hemoglobin A1c",
               "pop_ri": {"F": (4.0, 5.6), "M": (4.0, 5.6)}},
}

# ── Derived lookup tables ────────────────────────────────────────────────────

BATTERY2TESTCODE: dict = {}
for _tc, _cfg in MARKER_CONFIG.items():
    BATTERY2TESTCODE.setdefault(_cfg["battery"], []).append(_tc)

MARKER_UNITS: dict = {_tc: _cfg["units"] for _tc, _cfg in MARKER_CONFIG.items()}
MARKER_FULL_NAMES: dict = {_tc: _cfg["full_name"] for _tc, _cfg in MARKER_CONFIG.items()}


def get_population_ri(test_code: str, sex: str = "ALL") -> tuple:
    """
    Return the population reference interval (lower, upper) for a marker.

    Parameters
    ----------
    test_code : str
        Lab marker code, e.g. "HB", "PLT".
    sex : str
        "ALL", "F", or "M".
        - "ALL": returns the union of F and M bounds (min lower, max upper).
        - "F" / "M": returns the sex-specific tuple.
        Falls back to the union if sex-specific bounds are not available.

    Returns
    -------
    tuple (lower, upper)
        upper may be float("inf") for one-sided markers like HDL.

    Raises
    ------
    KeyError
        If test_code is not in MARKER_CONFIG.
    """
    cfg = MARKER_CONFIG[test_code]
    ri = cfg["pop_ri"]
    if sex in ri:
        return ri[sex]
    # Fall back to union across all available sexes
    lo = min(v[0] for v in ri.values())
    hi = max(v[1] for v in ri.values())
    return (lo, hi)
