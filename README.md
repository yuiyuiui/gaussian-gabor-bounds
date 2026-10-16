# Gaussian Gabor Bounds

Explore Gaussian Gabor frame-bound endpoints using **denominator-independent
finite sections**, Hermitian banded eigensolvers, and explicit Gaussian-tail,
localization and phase-grid error indicators. Matrix size depends on the
requested section length, not a rational approximation's denominator.

```sh
python -m pip install 'git+https://github.com/yuiyuiui/gaussian-gabor-bounds.git'
```

```python
from gaussian_gabor import GaussianPlan, choose_resolution

resolution = choose_resolution(1.5, x=.2, y=1.1, analytic_error=.05)
plan = GaussianPlan(1.5, x=.2, y=1.1, resolution=resolution)
result = plan.bounds()
print(result['a_interval_ignoring_roundoff'])  # Normalized lower endpoint A/density.
print(result['b_interval_ignoring_roundoff'])  # Normalized upper endpoint B/density.
print(result['errors'])
```

**These are numerical evidence, not certified enclosures.** All endpoints and
analytic indicators are computed in binary64; roundoff and parameter representation
error are unbounded. The resolution selector controls only the reported analytic
indicators. A lower interval endpoint <=0 leaves positivity unresolved.

## Intended use

Researchers studying Gaussian frames can compare shapes and oversampling densities,
refine finite-section/phase resolutions and inspect the distinct approximation
errors, including represented irrational densities. This is a small spectral
component, not a signal-analysis toolkit, an interval-arithmetic proof engine,
or a proof of hexagonal optimality. For finite periodic Gabor transforms and
frame bounds see [LTFAT](https://ltfat.org/doc/gabor/gabframebounds.html); its
finite-periodic problem differs from the infinite-lattice indicators here.

Python >=3.11, NumPy and SciPy are required. Plan setup is reusable and exposes
its effective configuration/workspace estimate. Dense reference allocation is
explicitly bounded. No input arrays, global precision or thread settings are mutated.

See [API and normalization](docs/API.md), [method and provenance](docs/METHOD.md),
[time/accuracy/memory benchmarks](bench/README.md), [readiness](docs/RELEASE_READINESS.md),
[contributing](CONTRIBUTING.md) and [citation](CITATION.cff).

Code, tests and documentation were developed with AI assistance under Kaiwen
Jin's direction. The finite controls do not constitute independent human review,
formal verification or an optimal-lattice theorem.
