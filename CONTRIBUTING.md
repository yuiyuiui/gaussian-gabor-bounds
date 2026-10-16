# Contributing and support

Use GitHub issues with a minimal reproducer, effective parameters, dependency
versions and the full exception. Contributions should include independent
numerical/exact checks and document error, precision, cost and memory effects.
Do not loosen tolerances to improve speed. Install `.[dev]`, run `python -m pytest`,
`python -m ruff check src tests bench examples`, and `python -m build`.
Test a wheel in a separate environment outside this checkout. Timings are not
correctness gates and no noisy speed threshold is enforced in shared CI.

Patch versions preserve the public contract. During 0.x, incompatible API changes
use a minor version and changelog entry. No third-party private inputs are required.
Security reports should avoid private data; use GitHub's private reporting if
available, otherwise report a minimal non-sensitive reproducer.

The implementation, tests and documentation were developed with AI assistance
under Kaiwen Jin's direction. This does not imply independent human review,
formal verification or community adoption. Credit original algorithms and
contributors; citation requests do not add license restrictions.
