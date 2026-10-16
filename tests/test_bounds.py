"""Marker bound policies resolve against the site's isolated values."""

import pytest

from perri.bounds import grid_from_interval, resolve_grid_bounds
from perri.marker_config import MARKER_CONFIG


def test_a1c_empirical_mu_keeps_clinical_sigma():
    policy = MARKER_CONFIG["A1C"]["mu_bound_policy"]
    bounds = resolve_grid_bounds((4.0, 5.6), policy, [4.5, 5.0, 6.0, 10.0])
    assert bounds["min_mu"] == pytest.approx(4.5075)
    assert bounds["max_mu"] == pytest.approx(9.94)
    assert bounds["max_sigma"] == grid_from_interval(4.0, 5.6)["max_sigma"]


def test_hdl_keeps_finite_clinical_lower_bound():
    policy = MARKER_CONFIG["HDL"]["mu_bound_policy"]
    bounds = resolve_grid_bounds((49, float("inf")), policy, [30, 40, 70, 90])
    expected = grid_from_interval(49, 88.5)
    assert bounds["min_mu"] == expected["min_mu"]
    assert bounds["max_mu"] == expected["max_mu"]
    assert bounds["max_sigma"] == pytest.approx((expected["max_mu"] - expected["min_mu"]) / 4 / 3.92)


@pytest.mark.parametrize("values", [None, [], [5, 5]])
def test_empirical_policy_rejects_missing_or_degenerate_values(values):
    with pytest.raises(ValueError):
        resolve_grid_bounds((4.0, 5.6), {"method": "empirical_mu", "percentiles": (0.5, 99.5)}, values)
