# Finite sections and error accounting

For integer hopping/phase frequencies m,n, define

    a[m,n] = exp(-pi * |H (m,n)|² / 2).
    T_theta[j,j-m] = sum_n a[m,n] exp(2 pi i n (theta+density*j-density*m/2)).

Truncate |m|,|n|<=N and restrict j to −L..L. The resulting matrix is Hermitian.
The implementation stores its lower band and precomputes weighted Fourier
factors for repeated phases. It retains one O(N² L) plan rather than storing
G dense matrices. A single band eigensolve per phase avoids repeating the same
band-to-tridiagonal reduction for the two endpoint searches.

Let c=pi*lambda_min(H^T H)/2. The analytic indicators used in the AIM-449
computational handoff are

    tau = 2(1+sqrt(pi/c)) sqrt(pi/c) erfc(N sqrt(c)),
    eta = pi² sum_{|m|,|n|<=N} m² a[m,n] / (8 L²),
    D = 2 pi sum_{|m|,|n|<=N} |n| a[m,n], grid_error = D/(2G).

The Gaussian radial tail bounds omitted coefficients; the localization term
controls the open finite section; the Lipschitz term controls gaps between
uniform phase samples. If a and b are the smallest/largest sampled eigenvalues,
the reported endpoint intervals are [a−tau−eta−grid_error,a+tau] and
[b−tau,b+tau+eta+grid_error]. Their applicability as rigorous intervals would
also require bounding arithmetic/input error, which this implementation **does
not do**. The bound formulas are reused from the prior handoff, not a new proof
certification or a solution of the originating optimality problem.

Independent controls compare the rectangular matrix with a time-domain Walnut
Gaussian product sum, use the density-two square theta endpoints, check shear
phase covariance, compare banded/dense eigensolvers and refine error indicators.
All are finite software controls at stated tolerances. The normalized density-two
square endpoint ratio is sqrt(2); reproducing it does not validate other lattices.

The code was adapted from the user's AIM-449 finite-section computational seed;
no third-party implementation was copied. SciPy supplies the established
[Hermitian banded eigensolver](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig_banded.html).
[LTFAT's Gabor frame bounds](https://ltfat.org/doc/gabor/gabframebounds.html)
address finite periodic systems and are a complementary existing tool, not a
matched infinite-lattice benchmark target. No new algorithmic priority is claimed.

## Direct-use result

For represented density delta>1, shape (x,y) with finite y>0, Gaussian cutoff
N>=1, section half-length L>=2 and uniform phase count G>=1, the API reports
[a−tau−eta−D/(2G),a+tau] for A/delta and
[b−tau,b+tau+eta+D/(2G)] for B/delta. Here a,b are the extreme sampled finite-
section eigenvalues; tau is the computed Gaussian tail indicator, eta the
computed localization indicator and D the computed phase Lipschitz indicator,
as defined above. A,B denote the Gaussian Gabor lower and upper frame endpoints.
These are **non-certified numerical intervals**: binary64 roundoff and input
representation error are excluded. A nonpositive lower endpoint leaves a
positive bound unresolved; no optimal-lattice conclusion follows.
