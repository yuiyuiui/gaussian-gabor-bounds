# 0.1.0 readiness evidence

| Gate | Status | Evidence |
|---|---|---|
| Value/scope | READY | README alternatives and bounded user contract; focused reusable component |
| Reliability | READY | 19 tests with independent controls; limitations in docs/API.md |
| Performance | READY | bench/ raw samples, profiles, error checks and time/memory scope |
| Usability/distribution | READY | Fresh isolated wheel outside checkout; full suite and quickstart passed; sdist-to-wheel rebuild and Twine passed |
| Community/credit | READY | MIT, actual contributor, algorithm provenance, CITATION.cff and support guide |
| Hosted identity/CI | NEEDS_WORK | Local candidate; actual publication/CI results are recorded separately |
| DOI / paper / GPU | NOT_APPLICABLE | Not requested and not needed for this CPU component |

Local checks use macOS arm64 and CPython 3.13.3. Configured CI covers Linux and
macOS only after actual jobs pass. No formal verification or external human
review is claimed. Numerical evidence and exact finite results retain their
separate meanings. No private research/control files are included.
