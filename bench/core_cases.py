"""Matched normalization, truncation, phase grid and requested extreme intervals."""

import numpy as np
from gaussian_gabor import GaussianPlan, Resolution


def make(case, variant):
    length, grid = {"small": (12, 12), "medium": (32, 24), "large": (64, 48)}[case]
    resolution = Resolution(6, length, grid)
    reference = GaussianPlan(1.5, x=0.25, y=1.1, resolution=resolution).bounds(backend="dense")

    def call():
        # Include Fourier setup in both variants, not only warm kernel reuse.
        plan = GaussianPlan(1.5, x=0.25, y=1.1, resolution=resolution)
        return plan.bounds(backend="dense" if variant == "reference" else "banded")

    def validate(result):
        for key in ("a_interval_ignoring_roundoff", "b_interval_ignoring_roundoff"):
            if not np.allclose(result[key], reference[key], atol=2e-12, rtol=2e-12):
                raise ArithmeticError("Dense and banded spectral endpoints disagree")
        return {
            "passed": True,
            "endpoint_atol": 2e-12,
            "endpoint_rtol": 2e-12,
            "analytic_errors": result["errors"],
            "roundoff_certified": False,
        }

    return (
        call,
        validate,
        dict(
            density=1.5,
            x=0.25,
            y=1.1,
            radius=6,
            length=length,
            grid=grid,
            variant=variant,
            cost_scope="plan construction plus bounds",
        ),
    )
