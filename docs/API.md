# API and numerical contract

`Resolution(radius=6, length=64, grid=128)` fixes integer Gaussian cutoff N>=1,
section half-length L>=2 (matrix dimension 2L+1), and uniform phase count G>=1.
`GaussianPlan(density, x=0, y=1, resolution=..., max_workspace_bytes=128*1024**2,
max_phase_work=200000000)` validates finite density>1, finite x and finite y>0.
Inputs describe the unit-energy Gaussian g(t)=2^(1/4) exp(−pi t²) and lattice
shape H=sqrt(density/y)[[1,x],[0,y]] in the adjoint-coordinate representation.
The spectral endpoints are **A/density and B/density**, not A and B themselves.

- `plan.configuration`: detached effective settings, arithmetic, normalization,
  workspace estimate and `roundoff_certified=False`.
- `plan.errors`: detached Gaussian-tail, localization and phase-grid terms.
- `plan.band(phase)`: a newly owned lower Hermitian band in SciPy storage format;
  phase is taken modulo one. Cached Fourier arrays remain unchanged.
- `plan.matrix(phase, max_dense_bytes=64*1024**2)`: bounded dense control view.
- `plan.bounds(backend='banded')`: two intervals, sampled extremum phases, explicit
  errors, optional condition-number interval and evidence status. The banded
  backend computes one full band spectrum internally and retains its extremes.
  `backend='dense'` uses NumPy's dense eigensolver on the same finite operator.
- `choose_resolution(..., analytic_error=.05, max_length=512, max_grid=4096,
  max_radius=24)`: choose N,L,G so each analytic term is <=target/3, then apply
  normal plan admission. It performs no spectral solves, does not control
  roundoff, and never relaxes the target silently.

All arithmetic is binary64/complex128. The effective input is its represented
binary value; no irrational input is represented exactly. Endpoint phases are
samples, not exact optimizing phases. Underflow, cancellation, extreme shape
conditioning and large phases can invalidate a rigorous enclosure; no interval
or correctly rounded result is claimed. Obvious nonfinite/overflow constructions
are rejected. The caller should perform resolution and precision comparisons for
scientific use; a future validated-arithmetic backend would be a separate guarantee.

Storage is O(N² L) for cached Fourier data plus band solver work, rather than an
unbounded phase-by-matrix tensor. Workspace and phase-work caps are admission
estimates, not an OS sandbox; dense controls have their own cap. Repeated calls
reuse only immutable-configuration data. No RNG or persistent disk cache is used.
Invalid parameters raise ValueError; unadmitted resolutions raise ResourceLimitError.
A nonpositive lower endpoint is returned as unresolved, not called a failure or
proof of non-frame behavior. The condition interval is None in that case.

After installation run `python examples/frame_bounds.py`. Development uses
`.[dev]`, `python -m pytest` and `python -m build`.
