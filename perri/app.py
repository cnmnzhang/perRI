"""
Personalized Reference Intervals — Streamlit App

Run with:
    streamlit run personalized_reference_intervals/app.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

import perri as pri

# ---------------------------------------------------------------------------
# Battery groupings
# ---------------------------------------------------------------------------
_SUPPORTED = set(pri.list_supported_markers())
_BATTERIES_RAW = {
    "CBC":      ["HB", "HCT", "RBC", "MCH", "MCHC", "MCV", "RDWCV", "PLT", "WBC"],
    "WBC diff": ["TNEUT", "LYMPH", "MONOC"],
    "BMP":      ["NA", "K", "CL", "CO2", "IGAP", "GLU", "BUN", "CRE", "CA"],
    "Hepatic":  ["ALB", "ALT", "AST", "ALK", "BIL", "BILD", "TP"],
    "Lipid":    ["CHOL", "LDL", "HDL", "NONHDL", "TRIG"],
    "Coag":     ["PROINR", "PROPAT"],
    "Misc":     ["A1C", "FER", "HSCRP", "LD", "MG", "P", "TSH", "VITDT"],
}
BATTERIES = {k: [tc for tc in v if tc in _SUPPORTED] for k, v in _BATTERIES_RAW.items()}
BATTERIES = {k: v for k, v in BATTERIES.items() if v}

# Color palette for up to 8 patients
_PALETTE = [
    ("steelblue",    "rgba(70,130,180,0.15)"),
    ("coral",        "rgba(255,127,80,0.15)"),
    ("seagreen",     "rgba(46,139,87,0.15)"),
    ("mediumpurple", "rgba(147,112,219,0.15)"),
    ("darkorange",   "rgba(255,140,0,0.15)"),
    ("crimson",      "rgba(220,20,60,0.15)"),
    ("teal",         "rgba(0,128,128,0.15)"),
    ("chocolate",    "rgba(210,105,30,0.15)"),
]

def _color(i):
    return _PALETTE[i % len(_PALETTE)]

# ---------------------------------------------------------------------------
# App config
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Personalized Reference Intervals", layout="wide")
st.title("Personalized Reference Intervals")
st.caption("Infers each patient's personal lab reference range from a time series of repeated measurements.")

# ---------------------------------------------------------------------------
# Sample data (loaded once at startup)
# ---------------------------------------------------------------------------
@st.cache_data
def _load_sample_data() -> pd.DataFrame:
    path = Path(__file__).parent / "data" / "sample_data.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, keep_default_na=False)

_SAMPLE_DATA = _load_sample_data()

# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------
def _init():
    defaults = {
        "results": {},
        "input_data": None,
        "excluded_by_id": {},
        "paste_seed": 0,
        "loaded_paste": "",
        "_autorun": False,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v
_init()

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
col_cfg, col_main, col_scores = st.columns([1, 2.2, 1])

# ============================================================
# LEFT: Configuration
# ============================================================
with col_cfg:
    st.header("Configuration")

    battery   = st.selectbox("Battery", list(BATTERIES.keys()), key="battery")
    test_code = st.selectbox(
        "Test Code",
        BATTERIES[battery],
        format_func=lambda tc: f"{pri.MARKER_FULL_NAMES.get(tc, tc)} ({tc})",
        key="test_code",
    )

    full_name = pri.MARKER_FULL_NAMES.get(test_code, test_code)
    units_str = pri.MARKER_UNITS.get(test_code, "")

    sex = st.radio("Sex", ["ALL", "M", "F"], horizontal=True, key="sex")

    st.markdown("---")
    with st.expander("Model Parameters", expanded=False):
        use_defaults = st.toggle("Use optimized defaults", value=True, key="use_defaults")
        defaults_p = pri.get_default_params(test_code, sex)
        if use_defaults:
            log_lam = defaults_p["log_lambda_"]
            st.caption(f"log(λ) = {log_lam:.3f}  →  λ = {np.exp(log_lam):.3f}")
            st.caption(f"μ grid: [{defaults_p['min_mu']:.3f}, {defaults_p['max_mu']:.3f}]")
            st.caption(f"σ grid: [{defaults_p['min_sigma']:.4f}, {defaults_p['max_sigma']:.3f}]")
            params = None
        else:
            log_lambda = st.slider(
                "log(λ)  —  temporal decay",
                min_value=-3.0, max_value=4.0,
                value=float(defaults_p["log_lambda_"]),
                step=0.05,
                help=(
                    "λ controls how quickly older measurements are down-weighted. "
                    "Optimized per marker over log(λ) ∈ [−3, 4], corresponding to λ ∈ [0.05, 54.6]. "
                    "Higher → recent measurements dominate; lower → all measurements weighted equally."
                ),
                key="log_lambda_slider",
            )
            st.caption(f"λ = {np.exp(log_lambda):.3f}")
            params = {**defaults_p, "log_lambda_": log_lambda}

    st.markdown("---")
    with st.expander("Isolation Filter", expanded=False):
        st.caption(
            "Removes measurements taken within a minimum gap of any neighbor. "
            "This prevents clusters of serial draws — e.g. daily labs during a hospitalization — "
            "from pulling the estimated setpoint toward an acutely abnormal value."
        )
        filter_iso = st.toggle("Enable isolation filter", value=True, key="filter_iso")
        min_gap = st.slider(
            "Min gap (days)", 0, 365, 90, 10,
            disabled=not filter_iso,
            help="Keep only measurements ≥ this many days from all neighbors.",
            key="min_gap",
        )
        min_meas = st.number_input("Min measurements required", min_value=1, max_value=20, value=5, key="min_meas")

    st.markdown("---")
    with st.expander("Plot Settings", expanded=False):
        ci = st.slider("Confidence interval", 0.5, 0.99, 0.95, 0.01, key="ci")

# ============================================================
# MIDDLE: Data entry + plot
# ============================================================
with col_main:
    st.header("Measurements")

    st.caption(
        "One measurement per line: `id, value, date` (YYYY-MM-DD). "
        "Dates are optional — omit to auto-space 6 months apart. "
        "Separate patients with a blank line."
    )

    load_col, clear_col = st.columns([1, 1])
    with load_col:
        if st.button("Load sample data", use_container_width=True):
            marker_rows = _SAMPLE_DATA[_SAMPLE_DATA["test_code"] == test_code]
            if not marker_rows.empty:
                lines = [
                    f"{row.patient_id}, {row.result_value}, {row.result_date}"
                    for _, row in marker_rows.iterrows()
                ]
                st.session_state["loaded_paste"] = "\n".join(lines)
                st.session_state["paste_seed"] += 1
                st.session_state["_autorun"] = True
            else:
                st.warning(f"No sample data available for {test_code}.")
    with clear_col:
        if st.button("Clear", use_container_width=True):
            st.session_state["loaded_paste"] = ""
            st.session_state["paste_seed"] += 1
            st.session_state["results"] = {}
            st.session_state["input_data"] = None

    paste_text = st.text_area(
        "Paste measurements",
        value=st.session_state.get("loaded_paste", ""),
        height=220,
        placeholder=(
            "1, 14.2, 2015-03-01\n1, 13.9, 2015-10-14\n1, 14.5, 2016-05-08\n"
            "\n"
            "2, 15.1, 2016-01-10\n2, 15.3, 2016-08-22\n2, 14.8, 2017-03-15"
        ),
        key=f"paste_text_{st.session_state['paste_seed']}",
        label_visibility="collapsed",
    )

    # ── Parser ────────────────────────────────────────────────
    def _fill_dates(dates, base="2010-01-01"):
        base_ts = pd.Timestamp(base)
        filled, counter = [], 0
        for d in dates:
            if not d or str(d).strip() == "":
                filled.append((base_ts + pd.DateOffset(months=6 * counter)).strftime("%Y-%m-%d"))
                counter += 1
            else:
                filled.append(str(d).strip())
                counter += 1
        return filled

    def parse_paste(text: str) -> dict:
        """Returns {patient_id: (values_list, dates_list)}."""
        if not text or not text.strip():
            return {}
        grouped: dict = {}
        auto_id = 1
        current_vals, current_dates = [], []

        def _flush():
            nonlocal auto_id
            if current_vals:
                grouped[auto_id] = (list(current_vals), _fill_dates(list(current_dates)))
                auto_id += 1

        for line in text.strip().splitlines():
            line = line.strip()
            if not line:
                _flush()
                current_vals.clear()
                current_dates.clear()
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 3:
                try:
                    pid = int(parts[0])
                    val = float(parts[1])
                    date = parts[2] if len(parts) > 2 else ""
                    if pid not in grouped:
                        grouped[pid] = ([], [])
                    grouped[pid][0].append(val)
                    grouped[pid][1].append(date)
                    continue
                except ValueError:
                    pass
            try:
                val = float(parts[0])
                date = parts[1] if len(parts) >= 2 else ""
                current_vals.append(val)
                current_dates.append(date)
            except (ValueError, IndexError):
                continue

        _flush()
        return {pid: (vals, _fill_dates(dates)) for pid, (vals, dates) in grouped.items()}

    input_data = parse_paste(paste_text or "")
    n_patients = len(input_data)
    total_rows = sum(len(v) for v, _ in input_data.values())

    if input_data:
        st.caption(f"{n_patients} patient{'s' if n_patients != 1 else ''}, {total_rows} total measurements.")

    # ── Run button ────────────────────────────────────────────
    run = st.button("▶ Run", type="primary")

    if run or st.session_state.get("_autorun"):
        st.session_state["_autorun"] = False
        if not input_data:
            st.warning("Enter at least one measurement.")
        else:
            results = {}
            excluded_by_id = {}
            with st.spinner(f"Fitting model for {n_patients} patient(s)..."):
                for pid, (vals, dates) in input_data.items():
                    try:
                        res = pri.fit_patient(
                            values=vals,
                            timestamps=dates,
                            test_code=test_code,
                            sex=sex,
                            params=None if use_defaults else params,
                            filter_isolated_measurements=filter_iso,
                            min_gap_days=min_gap,
                            min_measurements=int(min_meas),
                        )
                        if res is not None:
                            results[pid] = res
                            ts_input = pd.to_datetime(dates)
                            ts_kept  = set(t.date() for t in res.timestamps)
                            excluded_by_id[pid] = [t.date() not in ts_kept for t in ts_input]
                        else:
                            st.warning(f"Patient {pid}: not enough isolated measurements (need ≥ {int(min_meas)}).")
                    except Exception as e:
                        st.error(f"Patient {pid}: {e}")
            st.session_state.results = results
            st.session_state.excluded_by_id = excluded_by_id
            st.session_state.input_data = input_data

    # ── Plot ──────────────────────────────────────────────────
    results = st.session_state.results
    if results:
        z = stats.norm.ppf(1 - (1 - ci) / 2)
        fig = go.Figure()
        input_data_saved = st.session_state.get("input_data") or {}
        excluded_by_id   = st.session_state.get("excluded_by_id", {})

        # Population reference interval — gray background shading
        try:
            pop_lo, pop_hi = pri.get_population_ri(test_code, sex)
            if pop_hi == float("inf"):
                fig.add_hline(
                    y=pop_lo,
                    line_dash="dash",
                    line_color="lightgray",
                    annotation_text=f"Pop RI lower: {pop_lo} {units_str}",
                    annotation_position="bottom right",
                )
            else:
                fig.add_hrect(
                    y0=pop_lo, y1=pop_hi,
                    fillcolor="lightgray",
                    opacity=0.18,
                    layer="below",
                    line_width=0,
                    annotation_text=f"Pop RI",
                    annotation_position="top left",
                )
        except Exception:
            pass

        for i, (pid, res) in enumerate(sorted(results.items())):
            color_line, color_fill = _color(i)
            mus  = res.mu_history
            sigs = res.sigma_history
            ts   = res.timestamps
            ts_fwd = list(ts[1:]) + [ts[-1] + pd.Timedelta(days=90)]
            grp  = f"Patient {pid}"

            fig.add_trace(go.Scatter(
                x=list(ts_fwd) + list(reversed(ts_fwd)),
                y=list(mus + z * sigs) + list(reversed(mus - z * sigs)),
                fill="toself",
                fillcolor=color_fill,
                line=dict(color="rgba(0,0,0,0)"),
                hoverinfo="skip",
                name=f"{int(ci*100)}% PerRI",
                legendgroup=grp,
                showlegend=True,
            ))

            fig.add_trace(go.Scatter(
                x=ts_fwd, y=mus,
                mode="lines",
                name=f"Setpoint",
                line=dict(color=color_line, width=2),
                legendgroup=grp,
                legendgrouptitle_text=grp,
                showlegend=True,
                hoverinfo="skip",
            ))

            if pid in excluded_by_id and pid in input_data_saved:
                raw_vals, raw_dates = input_data_saved[pid]
                excl = excluded_by_id[pid]
                raw_ts = pd.to_datetime(raw_dates)
                ex_x = [d for d, e in zip(raw_ts, excl) if e]
                ex_y = [v for v, e in zip(raw_vals, excl) if e]
                if ex_x:
                    fig.add_trace(go.Scatter(
                        x=ex_x, y=ex_y,
                        mode="markers",
                        name=f"non-isolated",
                        marker=dict(size=7, color="lightgray", symbol="x"),
                        legendgroup=grp,
                        showlegend=True,
                    ))

            fig.add_trace(go.Scatter(
                x=ts, y=res.values,
                mode="markers+lines",
                name=f"measurements",
                marker=dict(size=7, color=color_line),
                line=dict(width=1, color=color_line, dash="dot"),
                legendgroup=grp,
                showlegend=True,
                customdata=res.mu_history,
                hovertemplate=(
                    f"<b>Patient {pid}</b><br>"
                    "Date: %{x|%Y-%m-%d}<br>"
                    f"Value: %{{y:.2f}} {units_str}<br>"
                    f"Setpoint (μ): %{{customdata:.2f}} {units_str}"
                    "<extra></extra>"
                ),
            ))

        fig.update_layout(
            height=440,
            # title=dict(text=f"{full_name}", x=0.5),
            xaxis_title="Date",
            yaxis_title=f"{test_code} {units_str}",
            legend=dict(
                groupclick="toggleitem",
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
            ),
            hovermode="closest",
        )
        st.plotly_chart(fig, use_container_width=True)

    elif not run:
        st.info("Enter measurements above and click **▶ Run** to fit the model.")

# ============================================================
# RIGHT: Results
# ============================================================
with col_scores:
    st.header("Results")

    results = st.session_state.results
    if not results:
        st.caption("Run the model to see results.")
    else:
        z = stats.norm.ppf(1 - (1 - ci) / 2)

        # Per-patient summaries
        for i, (pid, res) in enumerate(sorted(results.items())):
            color_line, _ = _color(i)
            mu_f  = res.mu_history[-1]
            sig_f = res.sigma_history[-1]
            lo    = mu_f - z * sig_f
            hi    = mu_f + z * sig_f

            header = f"Patient {pid}" if len(results) > 1 else "Summary"
            st.markdown(
                f"<span style='color:{color_line}; font-weight:bold; font-size:1.05em'>"
                f"{header}</span>",
                unsafe_allow_html=True,
            )
            st.markdown(
                f"**Personalized range:** {lo:.2f} – {hi:.2f} {units_str}"
            )
            st.markdown(f"**Setpoint:** {mu_f:.2f} {units_str}")
            if i < len(results) - 1:
                st.markdown("---")

        st.markdown("---")

        # All metrics in one collapsed expander
        with st.expander("Model diagnostics", expanded=False):
            mus_l  = [r.mu_history for r in results.values()]
            sigs_l = [r.sigma_history for r in results.values()]
            meas_l = [r.values for r in results.values()]

            # Normalize RMSE by median intra-patient SD for this marker
            try:
                intra_std = pri.get_intra_patient_std(test_code, sex)
                rmse_label = "RMSE / σ_intra"
                rmse_help = (
                    f"RMSE normalized by median intra-patient σ ({intra_std:.3f} {units_str}). "
                    "Values < 1 mean the setpoint predicts the next measurement within "
                    "typical within-person noise."
                )
                def _norm_rmse(raw: float) -> str:
                    return f"{raw / intra_std:.3f}"
            except Exception:
                intra_std = None
                rmse_label = "RMSE"
                rmse_help = "Root mean squared prediction error."
                def _norm_rmse(raw: float) -> str:
                    return f"{raw:.4f}"

            if len(results) > 1:
                agg_rmse  = pri.compute_rmse(mus_l, meas_l)
                agg_ks    = pri.compute_ks(mus_l, sigs_l, meas_l)
                agg_ks_qr = pri.compute_ks_quintile_range(mus_l, sigs_l, meas_l)
                st.markdown("**Aggregate**")
                st.metric(rmse_label, _norm_rmse(agg_rmse), help=rmse_help)
                st.metric("KS", f"{agg_ks:.4f}", help="KS distance from Uniform(0,1) on PIT values.")
                st.metric("KS quintile range", f"{agg_ks_qr:.4f}", help="Max − min KS across quintile bins; measures calibration uniformity.")
                st.markdown("---")

            for i, (pid, res) in enumerate(sorted(results.items())):
                color_line, _ = _color(i)
                label = f"Patient {pid}" if len(results) > 1 else "Per-measurement"
                st.markdown(
                    f"<span style='color:{color_line}; font-weight:bold'>{label}</span>",
                    unsafe_allow_html=True,
                )

                if len(res.values) >= 2:
                    p_rmse  = pri.compute_rmse([res.mu_history], [res.values])
                    p_ks    = pri.compute_ks([res.mu_history], [res.sigma_history], [res.values])
                    p_ks_qr = pri.compute_ks_quintile_range([res.mu_history], [res.sigma_history], [res.values])
                    col_a, col_b = st.columns(2)
                    col_a.metric(rmse_label, _norm_rmse(p_rmse), help=rmse_help)
                    col_b.metric("KS",   f"{p_ks:.4f}")
                    st.metric("KS quintile range", f"{p_ks_qr:.4f}")

                pit = [None] + list(stats.norm.cdf(
                    (res.values[1:] - res.mu_history[:-1])
                    / np.clip(res.sigma_history[:-1], 1e-6, None)
                ))
                tbl = pd.DataFrame({
                    "Date":  [str(t.date()) if hasattr(t, "date") else str(t) for t in res.timestamps],
                    "Value": np.round(res.values, 2),
                    "μ":     np.round(res.mu_history, 2),
                    "σ":     np.round(res.sigma_history, 2),
                    "Pct":   [f"{p:.2f}" if p is not None else "—" for p in pit],
                })
                st.dataframe(tbl, width='stretch', hide_index=True)
                if i < len(results) - 1:
                    st.markdown("---")
