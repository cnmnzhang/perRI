"""
Batch fit and evaluation example.

Fits the model for all patients of a given marker and reports cohort-level
error metrics (RMSE normalized by intra-patient SD, KS, KS quintile range).

Run from the repo root:
    python personalized_reference_intervals/batch_example.py
    python personalized_reference_intervals/batch_example.py --test_code TSH
    python personalized_reference_intervals/batch_example.py --test_code HB --sex M
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd

import perri as pri

DATA_DIR = Path(__file__).parent / "data"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--test_code", default="HB",  help="Lab marker code (default: HB)")
    parser.add_argument("--sex",       default="ALL", choices=["ALL", "M", "F"])
    args = parser.parse_args()

    df = pd.read_csv(DATA_DIR / "sample_data.csv", keep_default_na=False)
    full_name = pri.MARKER_FULL_NAMES.get(args.test_code, args.test_code)
    units     = pri.MARKER_UNITS.get(args.test_code, "")

    marker_df = df[df["test_code"] == args.test_code]
    if marker_df.empty:
        print(f"No data for test_code={args.test_code!r}. Supported: {pri.list_supported_markers()}")
        sys.exit(1)

    print(f"\n{full_name} ({args.test_code})  [sex={args.sex}]")
    print(f"  Input rows: {len(marker_df)}  across {marker_df['patient_id'].nunique()} patients")

    batch_df = pri.fit_batch(
        marker_df,
        value_col="result_value",
        timestamp_col="result_date",
        patient_id_col="patient_id",
        test_code=args.test_code,
        sex=args.sex,
    )

    if batch_df.empty:
        print("  No patients passed the isolation filter threshold.")
        sys.exit(1)

    scores = pri.evaluate_batch(batch_df)

    try:
        intra_std = pri.get_intra_patient_std(args.test_code, args.sex)
    except ValueError:
        intra_std = None

    print(f"\n  Patients fit:        {scores['n_patients']}")
    print(f"  Predictions:         {scores['n_predictions']}")
    print(f"  RMSE:                {scores['rmse']:.4f} {units}")
    if intra_std is not None:
        print(f"  RMSE Scaled*:        {scores['rmse'] / intra_std:.3f}")
    print(f"  KS statistic:        {scores['ks']:.4f}")
    print(f"  KS quintile range:   {scores['ks_quintile_range']:.4f}")

    # Per-patient breakdown
    print(f"\n  Per-patient summary:")
    for pid, grp in batch_df.groupby("patient_id"):
        mus   = grp.sort_values("measurement_index")["mu"].values
        sigs  = grp.sort_values("measurement_index")["sigma"].values
        vals  = grp.sort_values("measurement_index")["value"].values
        mu_f, sig_f = mus[-1], sigs[-1]
        from scipy import stats as _stats
        z95 = _stats.norm.ppf(0.975)
        lo, hi = mu_f - z95 * sig_f, mu_f + z95 * sig_f
        cv     = sig_f / mu_f * 100
        p_rmse = pri.compute_rmse([mus], [vals])
        p_ks   = pri.compute_ks([mus], [sigs], [vals])
        print(f"    Patient {pid}:  μ={mu_f:.2f}  SD={sig_f:.2f}  CV={cv:.1f}%  "
              f"95% CI=[{lo:.2f}, {hi:.2f}]  "
              f"RMSE={p_rmse:.4f}  KS={p_ks:.4f}")

    if intra_std is not None:
        print(f"\n  * RMSE normalized by median intra-patient SD; enables comparison across markers. 1 is the theoretical lower bound.")


if __name__ == "__main__":
    main()
