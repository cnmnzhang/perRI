"""
Matplotlib visualization of the Bayesian setpoint model fit.

Follows the style of plot_single_patient_panel() in the research repo:
  - x-axis in years from first measurement
  - Black o- markers for observations
  - #3B5CCB setpoint line + shaded PerRI band
  - Optional gray PopRI background band
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Match research repo constants (config/fig_config.py)
_COLOR      = "#3B5CCB"
_LINEWIDTH  = 1.25
_PER_ALPHA  = 0.6
_POP_ALPHA  = 0.3
_TITLE_FS   = 8


def _year_since_baseline(timestamps, extend_days: int = 90):
    """
    Convert timestamps to years from the first measurement.
    Returns (years_obs, years_est) where years_est shifts the setpoint
    trajectory forward by one step (mu[t] predicts the next observation)
    and appends a 90-day projection beyond the last measurement.
    """
    ts = pd.to_datetime(timestamps)
    ts = pd.Series(ts) if not isinstance(ts, pd.Series) else ts
    start = ts.min()

    new_point = ts.iloc[-1] + pd.Timedelta(days=extend_days)
    ts_est = pd.to_datetime(np.concatenate([ts.iloc[1:].to_numpy(), [new_point]]))

    years     = (ts     - start) / pd.Timedelta(days=365.25)
    years_est = (ts_est - start) / pd.Timedelta(days=365.25)
    return years.to_numpy(), years_est.to_numpy()


def plot_fit(
    values,
    timestamps,
    mu_history,
    sigma_history,
    confidence_interval: float = 0.95,
    pop_lo: float = None,
    pop_hi: float = None,
    title: str = None,
    ylabel: str = "Value",
    ax=None,
):
    """
    Plot observed measurements alongside the fitted setpoint trajectory.

    Matches the visual style of plot_single_patient_panel() in the research repo.
    The x-axis shows years from the patient's first measurement, making plots
    comparable across patients with different enrollment dates.

    Parameters
    ----------
    values : array-like of float
        Observed measurements (same length as mu_history).
    timestamps : array-like
        Corresponding dates (datetime-like or parseable strings).
    mu_history : array-like of float
        Posterior mean at each time step.
    sigma_history : array-like of float
        Posterior SD at each time step.
    confidence_interval : float
        Width of the PerRI band. Default 0.95 (95% CI).
    pop_lo, pop_hi : float, optional
        Population reference interval bounds. If both are provided, a gray
        background band is drawn. If pop_hi is inf (e.g. HDL), only pop_lo
        is drawn as a dashed line.
    title : str, optional
        Axes title (e.g. marker name).
    ylabel : str
        y-axis label. Defaults to "Value"; pass "HB (g/dL)" for full context.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. Creates a new figure if None.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    """
    values = np.asarray(values, dtype=float)
    mus    = np.asarray(mu_history, dtype=float)
    sigs   = np.asarray(sigma_history, dtype=float)

    years, years_est = _year_since_baseline(timestamps)
    z = stats.norm.ppf(0.5 + confidence_interval / 2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3))
    else:
        fig = ax.figure

    x_max = float(years_est[-1])
    ax.set_xlim(0, x_max)

    # Population RI background (drawn first, behind everything)
    if pop_lo is not None and pop_hi is not None:
        if np.isfinite(pop_hi):
            ax.fill_between(
                [0, x_max], pop_lo, pop_hi,
                alpha=_POP_ALPHA, color="gray", linewidth=0, label="Pop RI",
            )
        else:
            ax.axhline(pop_lo, color="gray", linewidth=_LINEWIDTH, linestyle="--",
                       alpha=0.6, label="Pop RI lower bound")

    # PerRI band
    ax.fill_between(
        years_est,
        mus - z * sigs,
        mus + z * sigs,
        alpha=_PER_ALPHA,
        color=_COLOR,
        linewidth=0,
        label=f"PerRI {int(confidence_interval * 100)}% CI",
    )

    # Setpoint trajectory
    ax.plot(
        years_est, mus,
        color=_COLOR, linewidth=_LINEWIDTH * 1.5,
        zorder=2, label="Setpoint",
    )

    # Observed measurements
    ax.plot(
        years, values,
        "o-",
        color="black",
        markersize=_LINEWIDTH + 1,
        linewidth=_LINEWIDTH,
        zorder=3,
        label="Observation",
    )

    ax.set_xlabel("Years from baseline")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=_TITLE_FS)
    ax.legend(fontsize=7)

    return fig, ax
