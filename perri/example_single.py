"""
Single-patient fit example.

Run from the repo root:
    python personalized_reference_intervals/single_example.py
    python personalized_reference_intervals/single_example.py --test_code PLT --patient_id 2
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

import perri as pri

DATA_DIR = Path(__file__).parent / "data"
Z95 = stats.norm.ppf(0.975)  # 1.96


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--test_code",  default="HB",  help="Lab marker code (default: HB)")
    parser.add_argument("--patient_id", default=3, type=int, help="Patient ID in sample data (1, 2, or 3; default: 3)")
    parser.add_argument("--sex",        default="ALL", choices=["ALL", "M", "F"])
    args = parser.parse_args()

    df = pd.read_csv(DATA_DIR / "sample_data.csv", keep_default_na=False)
    full_name = pri.MARKER_FULL_NAMES.get(args.test_code, args.test_code)
    units     = pri.MARKER_UNITS.get(args.test_code, "")

    patient_df = df[(df["test_code"] == args.test_code) & (df["patient_id"] == args.patient_id)]
    if patient_df.empty:
        print(f"No data for test_code={args.test_code!r}, patient_id={args.patient_id}")
        sys.exit(1)

    result = pri.fit_patient(
        values=patient_df["result_value"].values,
        timestamps=patient_df["result_date"].values,
        test_code=args.test_code,
        sex=args.sex,
    )
    if result is None:
        print("Not enough isolated measurements to fit the model.")
        sys.exit(1)

    mu, sigma = result.mu_history[-1], result.sigma_history[-1]
    lo, hi    = mu - Z95 * sigma, mu + Z95 * sigma
    cv        = sigma / mu * 100

    print(f"\n{full_name} ({args.test_code})  — Patient {args.patient_id}  [sex={args.sex}]")
    print(f"  Measurements used:   {len(result.values)}")
    print(f"  Setpoint (μ):        {mu:.2f} {units}")
    print(f"  Setpoint SD:         {sigma:.2f} {units}")
    print(f"  Setpoint CV:         {cv:.1f}%")
    print(f"  Personalized 95% CI: [{lo:.2f}, {hi:.2f}] {units}")

    # Population reference interval comparison
    try:
        pop_lo, pop_hi = pri.get_population_ri(args.test_code, args.sex)
        pop_str = f"≥ {pop_lo:.2f}" if pop_hi == float("inf") else f"[{pop_lo:.2f}, {pop_hi:.2f}]"
        print(f"  Population RI:       {pop_str} {units}")
    except Exception:
        pop_lo, pop_hi = None, None

    # Error metrics (require ≥ 2 measurements)
    if len(result.values) >= 2:
        rmse  = pri.compute_rmse([result.mu_history], [result.values])
        ks    = pri.compute_ks([result.mu_history], [result.sigma_history], [result.values])
        ks_qr = pri.compute_ks_quintile_range([result.mu_history], [result.sigma_history], [result.values])
        print(f"\n  RMSE:                {rmse:.4f} {units}")
        try:
            intra_std = pri.get_intra_patient_std(args.test_code, args.sex)
            print(f"  RMSE Scaled*:        {rmse / intra_std:.3f}")
        except ValueError:
            pass
        print(f"  KS statistic:        {ks:.4f}")
        print(f"  KS quintile range:   {ks_qr:.4f}")
        print(f"\n  * RMSE normalized by median intra-patient SD, enables comparison across markers. 1 is the theoretical lower bound.")

    # Plot
    fig, ax = pri.plot_fit(
        result.values,
        result.timestamps,
        result.mu_history,
        result.sigma_history,
        pop_lo=pop_lo,
        pop_hi=pop_hi,
        title=f"{full_name} — Patient {args.patient_id}",
        ylabel=f"{full_name} ({units})",
    )
    plt.tight_layout()
    figures_dir = Path(__file__).parent.parent / "figures"
    figures_dir.mkdir(exist_ok=True)
    out_path = figures_dir / f"patient{args.patient_id}_{args.test_code.lower()}_fit.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\n  Plot saved → {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
