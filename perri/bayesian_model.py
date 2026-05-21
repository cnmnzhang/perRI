"""
Core Bayesian setpoint model.

Algorithm overview
------------------
At each time step t, the model maintains a posterior distribution over
(mu, sigma) on a 2-D grid. Older measurements are down-weighted via an
exponential decay controlled by log_lambda_:

    weight(t) = exp(-lambda * (t_max - t) / t_max)

The posterior is proportional to exp(weighted_log_likelihood) * prior.
Marginal posteriors for mu and sigma are obtained by summing the joint over
the other axis. The reported estimate at each step is the posterior mean.

"""

import numpy as np
from numba import njit

DEFAULT_GRID_SIZE = 201


@njit
def fast_weight(log_lambda_, n):
    """Exponential decay weights for measurements 1..n (JIT-compiled)."""
    i_max = np.max(n)
    lambda_ = np.exp(log_lambda_)
    return np.exp(-lambda_ * (i_max - n) / i_max)


@njit
def fast_custom_log_likelihood(x, mu_grid, sig_grid, n, log_lambda_):
    """Time-weighted log-likelihood over the (mu, sigma) grid (JIT-compiled)."""
    T = len(x)
    M = len(mu_grid)
    S = len(sig_grid)
    ll = np.zeros((M, S))
    w = fast_weight(log_lambda_, n)

    for m in range(M):
        mu = mu_grid[m]
        for s in range(S):
            sigma = sig_grid[s]
            logp = 0.0
            for t in range(T):
                diff = (x[t] - mu) / sigma
                logp += -0.5 * w[t] * diff**2 - np.log(sigma) - 0.5 * np.log(2 * np.pi)
            ll[m, s] = logp
    return ll


def build_prior_flat(mu_grid, sig_grid, eps: float = 1e-12):
    """Uniform joint prior over (mu, sigma). Used as default."""
    mu_pdf = np.ones(mu_grid.size) / (mu_grid[-1] - mu_grid[0])
    sig_pdf = np.ones(sig_grid.size) / (sig_grid[-1] - sig_grid[0])
    prior = np.outer(mu_pdf, sig_pdf) + eps
    prior /= np.sum(prior)
    return prior


def bayesian(
    x,
    log_lambda_,
    min_mu,
    max_mu,
    min_sigma,
    max_sigma,
    grid_size=None,
):
    """
    Fit the Bayesian setpoint model to a sequence of measurements.

    Parameters
    ----------
    x : array-like of float
        Ordered measurement values (already filtered to isolated tests).
    log_lambda_ : float
        Log of the temporal decay rate. Higher values → recent measurements
        dominate. Optimized per-marker in the research repo.
    min_mu, max_mu : float
        Grid bounds for mu. Should bracket the physiological range for the marker.
    min_sigma, max_sigma : float
        Grid bounds for sigma (within-person SD).
    grid_size : int, optional
        Number of grid points per axis. Default 201.

    Returns
    -------
    mu_history : np.ndarray, shape (len(x),)
        Posterior mean of mu after observing each measurement.
    sigma_history : np.ndarray, shape (len(x),)
        Posterior mean of sigma after observing each measurement.
    """
    if grid_size is None:
        grid_size = DEFAULT_GRID_SIZE

    if int(grid_size) <= 0:
        return None, None, None

    x = np.array(x)
    grid_size = int(grid_size)

    mu_grid = np.linspace(min_mu, max_mu, grid_size)
    sig_grid = np.linspace(min_sigma, max_sigma, grid_size)

    prior0 = build_prior_flat(mu_grid, sig_grid)

    len_x = len(x)
    mu_history = np.zeros(len_x)
    sigma_history = np.zeros(len_x)

    for t in range(len_x):
        n = np.arange(1, t + 2)
        ll_values = fast_custom_log_likelihood(x[: t + 1], mu_grid, sig_grid, n, log_lambda_)

        joint = np.exp(ll_values - ll_values.max()) * prior0

        mu_post = np.sum(joint, axis=1)
        mu_post_sum = np.sum(mu_post)
        if mu_post_sum > 0:
            mu_post /= mu_post_sum
        else:
            mu_post = np.ones(grid_size) / grid_size

        sig_post = np.sum(joint, axis=0)
        sig_post_sum = np.sum(sig_post)
        if sig_post_sum > 0:
            sig_post /= sig_post_sum
        else:
            sig_post = np.ones(grid_size) / grid_size

        mu_history[t] = np.dot(mu_post, mu_grid)
        sigma_history[t] = np.dot(sig_post, sig_grid)

    return mu_history, sigma_history
