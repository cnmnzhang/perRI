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


def _time_coordinates(timestamps, extend_days: int = 90, projection_timestamp=None):
    """
    Convert timestamps to years from the first measurement.

    Returns:
    - years_obs: observed measurement times
    - years_fit: setpoint trajectory times, shifted forward by one step
      (mu[t] predicts the next observation) through the last fitted
      90-day estimate
    - fit_year: the final fitted estimate time (last measurement + extend_days)
    - projection_year: the requested projection timestamp, clipped to be at
      or after fit_year. Equals fit_year when projection_timestamp is None.
    """
    ts = pd.to_datetime(timestamps)
    ts = pd.Series(ts) if not isinstance(ts, pd.Series) else ts
    start = ts.min()

    fit_point = ts.iloc[-1] + pd.Timedelta(days=extend_days)
    projection_point = pd.to_datetime(projection_timestamp, errors="coerce")
    if pd.isna(projection_point) or projection_point < fit_point:
        projection_point = fit_point

    ts_fit = pd.to_datetime(np.concatenate([ts.iloc[1:].to_numpy(), [fit_point]]))

    years_obs = (ts - start) / pd.Timedelta(days=365.25)
    years_fit = (ts_fit - start) / pd.Timedelta(days=365.25)
    fit_year = float((fit_point - start) / pd.Timedelta(days=365.25))
    projection_year = float((projection_point - start) / pd.Timedelta(days=365.25))
    return years_obs.to_numpy(), years_fit.to_numpy(), fit_year, projection_year


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
    projection_timestamp=None,
    add_legend: bool = True,
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
    projection_timestamp : optional
        If given, extends the setpoint trajectory as a flat projection from
        the last fitted estimate out to this date (clipped to be at or after
        the normal 90-day fit horizon). Useful for showing "where the
        setpoint band would sit today" beyond the last measurement.
    add_legend : bool
        Whether to draw the legend. Default True.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    """
    values = np.asarray(values, dtype=float)
    mus    = np.asarray(mu_history, dtype=float)
    sigs   = np.asarray(sigma_history, dtype=float)

    years, years_fit, fit_year, projection_year = _time_coordinates(
        timestamps,
        projection_timestamp=projection_timestamp,
    )
    z = stats.norm.ppf(0.5 + confidence_interval / 2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3))
    else:
        fig = ax.figure

    x_max = projection_year
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
        years_fit,
        mus - z * sigs,
        mus + z * sigs,
        alpha=_PER_ALPHA,
        color=_COLOR,
        linewidth=0,
        label=f"PerRI {int(confidence_interval * 100)}% CI",
    )

    if projection_year > fit_year:
        proj_x = np.array([fit_year, projection_year], dtype=float)
        proj_lo = np.array([mus[-1] - z * sigs[-1], mus[-1] - z * sigs[-1]], dtype=float)
        proj_hi = np.array([mus[-1] + z * sigs[-1], mus[-1] + z * sigs[-1]], dtype=float)
        ax.fill_between(
            proj_x,
            proj_lo,
            proj_hi,
            alpha=_PER_ALPHA * 0.75,
            color=_COLOR,
            linewidth=0,
        )

    # Setpoint trajectory
    ax.plot(
        years_fit, mus,
        color=_COLOR, linewidth=_LINEWIDTH * 1.5,
        zorder=2, label="Setpoint",
    )

    if projection_year > fit_year:
        ax.plot(
            [fit_year, projection_year],
            [mus[-1], mus[-1]],
            color=_COLOR,
            linewidth=_LINEWIDTH * 1.5,
            linestyle="-",
            zorder=2,
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
    if add_legend:
        ax.legend(fontsize=7)

    return fig, ax
