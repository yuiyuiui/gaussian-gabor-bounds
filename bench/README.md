# Reproducible core benchmarks

[Raw results](results/core_v0.1.0.json) include five warm samples per variant,
first-call and setup times, source hashes, profiles, validation, settings and
dependency versions. macOS arm64, CPython 3.13.3, one numerical worker and
BLAS/OpenMP thread. Variants alternate execution order between cases and run in
separate fresh processes. No random input seed is hidden.

| Case | Reference median s | Production median s | Ratio | Python peak B (ref/prod) | Process RSS B (ref/prod) |
|---|---:|---:|---:|---:|---:|
| small | 0.00058 | 0.000489375 | 1.19x | 73110 / 63802 | 54689792 / 54624256 |
| medium | 0.00645208 | 0.00343825 | 1.88x | 269174 / 151482 | 54558720 / 54771712 |
| large | 0.042504 | 0.0227791 | 1.87x | 796150 / 291770 | 55410688 / 55001088 |

Both paths construct the same Fourier plan and return the same two spectral intervals at the same density, shape, truncation, section and phase resolution. The reference uses NumPy's dense Hermitian eigensolver; production uses SciPy's Hermitian banded eigensolver. Setup is included in each measured call. An initial profile found two endpoint searches spending most time repeating the band reduction; one full band eigensolve per phase removed that duplication. The final comparison retains both numerical controls and the unchanged analytic error budget. Banded storage reduces the matrix representation, but cached Fourier data and overall process RSS must be considered separately. The error indicators are not a roundoff-certified accuracy bound.

Validation is excluded from timed calls, repeated after every sample, and must
agree with the independent reference. Separate profiling and tracemalloc calls
do not contribute to timing samples. Setup includes imports, case construction
and the independent reference; it is not a pure production startup measurement.
Process RSS includes these reference/validation dependencies. Python-traced peak
measures one separate call and excludes untraced native memory. Small RSS changes
are noise; no universal process-memory reduction is claimed. These finite cases
do not establish asymptotic speedups or performance on arbitrary hardware.

```sh
python -m pip install '.[dev]'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python bench/run.py \
  --repeats 5 --output /tmp/core-comparison.json
```

The output path must be new. The driver runs bounded small/medium/large cases;
performance thresholds are intentionally absent from shared CI. No network,
paid backend, GPU or hidden input data is needed. Each JSON `settings` object
specifies its requested output and cost boundary. Keep accuracy settings fixed
when comparing variants.

## Work and analytic error

| Case | Section dimension | Phases | Tail | Localization | Phase grid |
|---|---:|---:|---:|---:|---:|
| small | 25 | 12 | 1.27e-29 | 0.00246 | 0.0475 |
| medium | 65 | 24 | 1.27e-29 | 0.000346 | 0.0238 |
| large | 129 | 48 | 1.27e-29 | 8.66e-05 | 0.0119 |

These analytic terms shrink under refinement; they exclude arithmetic and input error.
