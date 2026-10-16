"""Banded finite sections of the normalized Janssen operator.

All computed endpoints ignore binary64 roundoff and input representation error.
The analytic tail, localization and phase-grid terms are reported individually.
"""

from dataclasses import dataclass, asdict
import math
from numbers import Integral, Real

import numpy as np
from scipy.linalg import eig_banded


class ResourceLimitError(RuntimeError):
    """Requested finite-section work or workspace exceeds an explicit limit."""


def positive_integer(x, name, low=1):
    if isinstance(x, bool) or not isinstance(x, Integral) or x < low:
        raise ValueError(f"{name} must be an integer >= {low}")
    return int(x)


def finite_real(x, name):
    if isinstance(x, bool) or not isinstance(x, Real) or not math.isfinite(x):
        raise ValueError(f"{name} must be a finite real scalar")
    return float(x)


@dataclass(frozen=True)
class Resolution:
    radius: int = 6
    length: int = 64
    grid: int = 128

    def __post_init__(self):
        for name, low in [("radius", 1), ("length", 2), ("grid", 1)]:
            positive_integer(getattr(self, name), name, low)


class GaussianPlan:
    """Reusable immutable-configuration finite sections; no global thread state.

    Density >1; shape (x,y), y>0; Gaussian g(t)=2^(1/4) exp(-pi t^2).
    Returned spectral endpoints are A/density and B/density. The plan stores
    O(radius^2 * length) complex Fourier data, with explicit workspace admission.
    """

    def __init__(
        self,
        density,
        x=0.0,
        y=1.0,
        *,
        resolution=Resolution(),
        max_workspace_bytes=128 * 1024 * 1024,
        max_phase_work=200000000,
    ):
        density, x, y = (finite_real(density, "density"), finite_real(x, "x"), finite_real(y, "y"))
        if density <= 1 or y <= 0:
            raise ValueError("Require density >1 and y>0")
        if not isinstance(resolution, Resolution):
            raise TypeError("resolution must be a Resolution")
        budget = positive_integer(max_workspace_bytes, "max_workspace_bytes")
        work = positive_integer(max_phase_work, "max_phase_work")
        N, L, G = resolution.radius, resolution.length, resolution.grid
        size, width = 2 * L + 1, min(N, 2 * L) + 1
        # Conservative admission includes cached Fourier data and solver work.
        estimate = (
            16 * (width * (2 * N + 1) * size + 12 * width * size + 8 * size) + 64 * (2 * N + 1) ** 2
        )
        if estimate > budget or max(width * (2 * N + 1) * size * G, G * size * size * width) > work:
            raise ResourceLimitError("Finite-section workspace/work exceeds configured admission")
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            try:
                h = math.sqrt(density / y) * np.array([[1.0, x], [0.0, y]])
                gram = h.T @ h
                # det(gram)/lambda_max avoids subtracting nearly equal eigenvalues.
                largest = float(np.linalg.eigvalsh(gram)[-1])
                mu = density * density / largest
                c = math.pi * mu / 2
                if not math.isfinite(c) or c <= 0:
                    raise ValueError("Shape is outside the stable binary64 parameter range")
                n = np.arange(-N, N + 1)
                mm, nn = np.meshgrid(np.arange(-N, N + 1), n, indexing="ij")
                norms = (density / y) * ((mm + x * nn) ** 2 + (y * nn) ** 2)
                weights = np.exp(-math.pi * norms / 2)
                tail = (
                    2
                    * (1 + math.sqrt(math.pi / c))
                    * math.sqrt(math.pi / c)
                    * math.erfc(N * math.sqrt(c))
                )
                moment = float(np.sum(mm * mm * weights))
                lipschitz = 2 * math.pi * float(np.sum(np.abs(nn) * weights))
                fourier = np.zeros((width, 2 * N + 1, size), dtype=complex)
                for m in range(width):
                    cols = np.arange(size - m)
                    angles = n[:, None] * (density * (cols - L + m / 2))[None, :]
                    fourier[m, :, : size - m] = weights[N + m, :, None] * np.exp(
                        2j * math.pi * np.remainder(angles, 1.0)
                    )
            except (FloatingPointError, OverflowError) as exc:
                raise ValueError("Parameters exceed stable binary64 construction range") from exc
        if not math.isfinite(tail + moment + lipschitz):
            raise ValueError("Error indicators are nonfinite for this shape")
        n.flags.writeable = fourier.flags.writeable = False
        self._density, self._x, self._y = density, x, y
        self._resolution, self._n, self._fourier = resolution, n, fourier
        self._tail, self._moment, self._lipschitz = tail, moment, lipschitz
        self._workspace, self._max_phase_work = estimate, work

    @property
    def configuration(self):
        return dict(
            density=self._density,
            x=self._x,
            y=self._y,
            **asdict(self._resolution),
            estimated_workspace_bytes=self._workspace,
            normalization="A/density, B/density",
            arithmetic="binary64",
            roundoff_certified=False,
        )

    @property
    def errors(self):
        r = self._resolution
        return {
            "gaussian_tail": self._tail,
            "localization": math.pi**2 * self._moment / (8 * r.length**2),
            "phase_grid": self._lipschitz / (2 * r.grid),
        }

    def band(self, phase):
        """Lower Hermitian band in SciPy convention; caller owns the returned array."""
        phase = finite_real(phase, "phase") % 1.0
        factor = np.exp(2j * math.pi * np.remainder(self._n * phase, 1.0))
        result = np.einsum("n,mnj->mj", factor, self._fourier, optimize=False)
        # The diagonal is mathematically real; its rounding asymmetry is removed.
        result[0] = result[0].real
        return result

    def matrix(self, phase, *, max_dense_bytes=64 * 1024 * 1024):
        """Explicit small-case reference view; refuse unadmitted dense allocation."""
        n = 2 * self._resolution.length + 1
        if 16 * n * n > positive_integer(max_dense_bytes, "max_dense_bytes"):
            raise ResourceLimitError("Dense reference matrix exceeds max_dense_bytes")
        band = self.band(phase)
        result = np.zeros((n, n), dtype=complex)
        for m in range(len(band)):
            j = np.arange(n - m)
            result[j + m, j] = band[m, : n - m]
            if m:
                result[j, j + m] = band[m, : n - m].conj()
        return result

    def bounds(self, *, backend="banded"):
        """Sample a uniform phase grid and return analytic-error-inflated intervals.

        `backend='dense'` is an independent eigensolver control at the same
        truncation/resolution, not an algorithm for larger problem sizes.
        """
        if backend not in ("banded", "dense"):
            raise ValueError("backend must be banded or dense")
        r = self._resolution
        lows, highs = [], []
        for phase in np.arange(r.grid) / r.grid:
            if backend == "banded":
                band = self.band(float(phase))
                # A single band reduction beat two separate endpoint searches
                # on the profiled cases. Preserve both extremes at the same precision.
                eig = eig_banded(band, lower=True, eigvals_only=True, check_finite=True)
                low, high = eig[0], eig[-1]
            else:
                eig = np.linalg.eigvalsh(self.matrix(float(phase)))
                low, high = eig[0], eig[-1]
            lows.append(float(low))
            highs.append(float(high))
        a, b = min(lows), max(highs)
        err = self.errors
        pad = sum(err.values())
        tail = err["gaussian_tail"]
        al, au, bl, bu = a - pad, a + tail, b - tail, b + pad
        return {
            "configuration": self.configuration,
            "errors": err,
            "a_interval_ignoring_roundoff": [al, au],
            "b_interval_ignoring_roundoff": [bl, bu],
            "condition_interval_ignoring_roundoff": [bl / au, bu / al] if al > 0 else None,
            "a_sample_phase": int(np.argmin(lows)) / r.grid,
            "b_sample_phase": int(np.argmax(highs)) / r.grid,
            "positive_lower_bound_resolved_ignoring_roundoff": al > 0,
            "status": "NUMERICAL EVIDENCE; roundoff and input error uncertified",
        }


def choose_resolution(
    density, x=0.0, y=1.0, *, analytic_error=0.05, max_length=512, max_grid=4096, max_radius=24
):
    """Choose a finite resolution for the *analytic* indicators, without solving.

    This target does not control floating-point error. Fail rather than silently
    relax the target. Binary64 parameter and plan work limits remain in force.
    """
    target = finite_real(analytic_error, "analytic_error")
    if target <= 0:
        raise ValueError("analytic_error must be positive")
    for name, value in [
        ("max_length", max_length),
        ("max_grid", max_grid),
        ("max_radius", max_radius),
    ]:
        positive_integer(value, name)
    for N in range(1, max_radius + 1):
        probe = GaussianPlan(density, x, y, resolution=Resolution(N, 2, 1))
        if probe.errors["gaussian_tail"] <= target / 3:
            length = max(2, math.ceil(math.sqrt(3 * math.pi**2 * probe._moment / (8 * target))))
            grid = max(1, math.ceil(3 * probe._lipschitz / (2 * target)))
            if length > max_length or grid > max_grid:
                raise ResourceLimitError("Requested analytic error exceeds length/grid limits")
            result = Resolution(N, length, grid)
            GaussianPlan(density, x, y, resolution=result)
            return result
    raise ResourceLimitError("Requested Gaussian tail exceeds max_radius")
