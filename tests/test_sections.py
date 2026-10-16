import math
import numpy as np
import pytest
from gaussian_gabor import GaussianPlan, Resolution, ResourceLimitError, choose_resolution


def walnut(density, y, phase, length, radius=10):
    """Independent time-domain Gaussian product sum for rectangular lattices."""
    a, b = 1 / math.sqrt(density * y), math.sqrt(y / density)
    indices = np.arange(-length, length + 1)
    out = np.zeros((len(indices), len(indices)))
    for j in indices:
        tj = a * phase + j / b
        for m in range(-radius, radius + 1):
            if not -length <= j - m <= length:
                continue
            center = round((tj - m / (2 * b)) / a)
            k = np.arange(center - 10, center + 11)
            product = math.sqrt(2) * np.exp(
                -math.pi * ((tj - k * a) ** 2 + (tj - k * a - m / b) ** 2)
            )
            out[j + length, j - m + length] = a * np.sum(product)
    return out


@pytest.mark.parametrize("density", [math.sqrt(2), 1.5, 2.0])
def test_independent_walnut_and_periodic_phase(density):
    plan = GaussianPlan(density, y=1.2, resolution=Resolution(10, 8, 8))
    np.testing.assert_allclose(
        plan.matrix(0.173), walnut(density, 1.2, 0.173, 8), atol=5e-13, rtol=5e-13
    )
    np.testing.assert_allclose(plan.matrix(0.173), plan.matrix(1.173), atol=5e-13)


def test_shear_covariance_and_banded_vs_dense():
    plan = GaussianPlan(math.sqrt(2), x=0.41, y=1.1, resolution=Resolution(6, 12, 12))
    a, b = plan.matrix(0.173 + math.sqrt(2)), plan.matrix(0.173)
    np.testing.assert_allclose(a[:-1, :-1], b[1:, 1:], atol=3e-13, rtol=3e-13)
    dense, band = plan.bounds(backend="dense"), plan.bounds()
    for key in ("a_interval_ignoring_roundoff", "b_interval_ignoring_roundoff"):
        np.testing.assert_allclose(dense[key], band[key], atol=5e-13, rtol=5e-13)
    assert band["configuration"]["roundoff_certified"] is False


def test_density_two_theta_endpoints_and_refinement():
    plan = GaussianPlan(2.0, resolution=Resolution(7, 16, 32))
    result = plan.bounds()
    n = np.arange(-12, 13)
    weights = np.exp(-math.pi * n * n)
    exact = [np.sum((-1.0) ** n * weights) ** 2, np.sum(weights) ** 2]
    for key, value in zip(("a_interval_ignoring_roundoff", "b_interval_ignoring_roundoff"), exact):
        lo, hi = result[key]
        assert lo <= value <= hi
    assert abs(exact[1] / exact[0] - math.sqrt(2)) < 1e-12
    fine = GaussianPlan(2.0, resolution=Resolution(7, 32, 64)).bounds()
    assert sum(fine["errors"].values()) < sum(result["errors"].values())


def test_resolution_selection_and_resource_failures():
    resolution = choose_resolution(1.5, x=0.2, y=1.1, analytic_error=0.05)
    plan = GaussianPlan(1.5, x=0.2, y=1.1, resolution=resolution)
    assert sum(plan.errors.values()) <= 0.05
    with pytest.raises(ResourceLimitError):
        choose_resolution(1.5, analytic_error=1e-12, max_length=8)
    with pytest.raises(ResourceLimitError):
        GaussianPlan(2.0, max_workspace_bytes=1)
    with pytest.raises(ResourceLimitError):
        plan.matrix(0, max_dense_bytes=1)
    with pytest.raises(ResourceLimitError):
        GaussianPlan(2.0, max_phase_work=1)


@pytest.mark.parametrize(
    "args",
    [(1.0,), (float("inf"),), (2.0, float("nan")), (2.0, 0.0, 0.0), (2.0, 0.0, -1.0), (True,)],
)
def test_invalid_parameters(args):
    with pytest.raises(ValueError):
        GaussianPlan(*args)


@pytest.mark.parametrize("args", [(0, 4, 8), (2, 1, 8), (2, 4, 0), (2.0, 4, 8), (True, 4, 8)])
def test_invalid_resolutions(args):
    with pytest.raises(ValueError):
        Resolution(*args)


def test_returned_arrays_and_configuration_do_not_mutate_plan():
    plan = GaussianPlan(2.0, resolution=Resolution(4, 8, 8))
    before = plan.matrix(0.2)
    band = plan.band(0.2)
    band[:] = 0
    config = plan.configuration
    config["density"] = 10
    np.testing.assert_array_equal(before, plan.matrix(0.2))
    with pytest.raises(ValueError):
        plan.band(float("nan"))
    with pytest.raises(ValueError):
        plan.bounds(backend="guess")


def test_excessive_spectral_work_is_rejected_before_large_allocation():
    with pytest.raises(ResourceLimitError):
        GaussianPlan(2.0, resolution=Resolution(1, 100000, 1))
