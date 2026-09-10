# Computational certificate for Lemma 4.1

This directory contains the reproducible computation used in Lemma 4.1. The
program performs two separate tasks:

1. it proves, using outward-rounded interval arithmetic, that
   \[
   \max_{(r,m,\lambda)\in\mathcal D}
   \min\{\widehat C_\eta(r,m,\lambda),
          \widehat C_{1/3}(r,m,\lambda)\}<\frac{79}{25};
   \]
2. it searches for a feasible point with a large objective value, certifies
   the value of that witness from below, and maps the point back to the
   variables of the CVRP reduction.

Here
\[
\mathcal D=\left\{(r,m,\lambda):
\frac{33}{50}\le r\le1,\quad m,\lambda\ge0,\quad
m+\lambda\le\frac3{25}\right\}.
\]
The statement of Lemma 4.1 still uses the old names \((u,g,\delta)\). The
identification is exactly
\[
u=r,\qquad g=m,\qquad \delta=\lambda.
\]

## Functions checked

After the normalization \(\mathrm{OPT}=1\), the program evaluates the two
expressions displayed after Lemma 4.1:
\[
\begin{aligned}
\widehat C_\eta(r,m,\lambda)
={}&\frac52-r+2r\int_0^1\min\left\{1,
\frac{1-(\frac{16}{25}+2m+3\lambda)t}
{2r-(2r+\frac{16}{25}+2m+4\lambda)t}\right\}\,dt,
\\
\widehat C_{1/3}(r,m,\lambda)
={}&3-\frac32r+\left(\frac32r+\frac4{25}+m\right)
\int_0^1\min\left\{1,
\frac{1-(\frac{16}{25}+2m+3\lambda)t}
{\frac32r-\frac1{50}+2m+\frac32\lambda
 -(\frac32r+\frac{16}{25}+\frac52m+3\lambda)t},
\right.\\[-2mm]
&\hspace{57mm}\left.
\frac{1-(\frac{16}{25}+2m+3\lambda)t}
{\frac32r+\frac4{25}+m
 -(\frac32r+\frac45+3m+4\lambda)t}
\right\}\,dt.
\end{aligned}
\]
As stipulated in the paper, a fraction is treated as infinity when it becomes
negative. Since the pointwise minimum already contains the constant branch
\(1\), the implementation equivalently represents an inadmissible quotient
branch by \(1\). A quotient is evaluated only where its numerator and
denominator are positive.

## Mapping of a witness

For every reported witness, the program also prints the following normalized
quantities:
\[
\begin{array}{rclcrcl}
R_1(V_0^\eta)&=&1-r,
&&R_1(V_\eta^{1/3})&=&r-\frac{14}{25}+\lambda,\\
R_1(V_{\mathrm{single}})&=&\frac8{25}+m+\lambda,
&&R_1(V_{\mathrm{double}})&=&\frac6{25}-m-2\lambda,\\
R_0(V_{\mathrm{single}})&=&\frac{16}{25}+2m+3\lambda,
&&R_0(V_{\mathrm{double}})&=&\frac{18}{25}-4m-6\lambda.
\end{array}
\]
It also reports
\[
\sigma_\eta(V_\eta^1)=2r,
\qquad
\sigma_{1/3}(V_\eta^1)=\frac32r+\frac4{25}+m.
\]

## Rigorous upper certificate

On a parameter box, each quotient has the form
\[
q(t)=\frac{1-pt}{a-bt}.
\]
If \(p_-\), \(a_-\), and \(b_+\) are outward-rounded coefficient bounds,
then, whenever all relevant numerators and denominators are positive,
\[
q(t)\le \frac{1-p_-t}{a_--b_+t}.
\]
The program integrates this upper function analytically, using
\[
\int_l^r\frac{1-pt}{a-bt}\,dt
=\frac pb(r-l)+\frac{b-ap}{b^2}
  \log\!\left(\frac{a-bl}{a-br}\right).
\]
Every arithmetic operation and the logarithm are evaluated with directed
rounding by `Boost.Numeric.Interval`. Adaptive subdivision in \(t\) resolves
validity thresholds and changes of the active branch. A multithreaded,
best-first subdivision of \(\mathcal D\) then certifies every parameter box at
or below the downward-rounded representation of \(79/25\). Consequently, the
accepted bound is strictly below the exact rational target.

## Second implementation in Python

`verify.py` proves a sharper variant of the upper bound, on top of
[python-flint](https://github.com/flintlib/python-flint). It uses the
coefficient \(\frac18\) instead of \(\frac14\) in front of
\(R_0(V_{\mathrm{double}})\) in \(\widehat C_{1/3}\): the constant term of the
first denominator of \(\widehat C_{1/3}\) is
\[
\sigma_{1/3}(V_\eta^1)-\frac18R_0(V_{\mathrm{double}})
=\frac32r+\frac7{100}+\frac32m+\frac34\lambda
\]
instead of
\(\sigma_{1/3}(V_\eta^1)-\frac14R_0(V_{\mathrm{double}})
=\frac32r-\frac1{50}+2m+\frac32\lambda\) as displayed above. The slope of
that denominator, the second quotient and \(\widehat C_\eta\) are unchanged.
For these functions it proves
\[
\max_{(r,m,\lambda)\in\mathcal D}
\min\{\widehat C_\eta(r,m,\lambda),
       \widehat C_{1/3}(r,m,\lambda)\}<\frac{3159}{1000}=3.159.
\]
This does not contradict the witness value \(3.1597523\) reported below: the
witness is for the functions with the coefficient \(\frac14\), which the C++
and Rust programs check against \(79/25\).

`verify.py` is standalone: no Rust or C++ toolchain, and no file other than
`verify.py` has to be read in order to check the proof.

It differs from the C++ and Rust programs in where the rigour comes from. Every
coefficient of the problem is an affine function, with rational coefficients, of
the endpoints of the parameter box, and bisecting a rational interval yields
rational endpoints. So `verify.py` keeps the box endpoints, the branch data
\(p_-,a_-,b_+\), the cut points and the constraint \(m+\lambda\le3/25\) as
**exact rationals** (`flint.fmpq`); every decision the search makes -- branch
admissibility, acceptance of a box, the final comparison against
\(3159/1000\) -- is an exact rational comparison, with no rounding to analyse.
The only quantity that is not rational is the logarithm in the closed form
above; it alone is evaluated in Arb ball arithmetic (`flint.arb`, rigorous by
construction) and its ball upper endpoint is converted back to an exact
rational. The whole numerical-soundness argument is therefore confined to the
functions `upper_fmpq` and `segment_upper_bound`. Because the target
\(3159/1000\) is an exact rational and each box bound is compared against it
with a strict `<`, no "downward-rounded target" argument is needed.

Parallelism is a static decomposition: the initial box is bisected into
\(2^{12}\) task boxes whose union is exactly the initial box, and each task is
proved independently in a worker process by the pure function `prove`. Running
with `--jobs 1` performs the identical computation in one process.

```bash
source /path/to/venv/bin/activate     # provides python-flint
python3 verify.py                     # ~35 s on 14 threads
python3 verify.py --jobs 1            # same result, one process, ~100 s
```

Exit status is `0` when the bound is proved and `2` when a box could not be
certified. A reference run reports

```
PROVED: max min{C_eta, C_1/3} <= 3.158999999759877 < 3159/1000
boxes: 2783916, max depth: 24, jobs: 14, elapsed: 33.9s
```

The companion `selftest.py` cross-checks `verify.py` against independently
written, deliberately naive implementations -- the closed-form antiderivative
against a rigorous interval Riemann sum, the per-box bound against a plain-float
re-transcription of the formulas above (with the \(\frac18\) coefficient), the
bound itself against a naive
Riemann sum of the envelope at interior points, and the task decomposition
against the domain \(\mathcal D\). It is a test of the transcription, not part
of the proof, and `verify.py` does not refer to it:

```bash
python3 selftest.py
```

## Feasible-witness search

The lower search begins with a grid in \(r\) and the triangular
\((m,\lambda)\)-domain. It then combines:

- best-first branch-and-bound, using the same rigorous box upper bounds as the
  certificate;
- feasible samples from the interiors and slanted boundaries of promising
  boxes;
- local mesh refinement, including refinement along
  \(m+\lambda=3/25\) and balancing the two cost bounds.

At a fixed point, all branch changes are zeros of linear functions, so the
search value is obtained from analytic antiderivatives. Finally, the best
witness is evaluated again with a lower Riemann sum whose operations are all
outward-rounded intervals. The line `Certified feasible lower bound` is
therefore a rigorous lower bound on the maximum in Lemma 4.1. The more precise
line `Numerical value at witness` and the reported coordinates locate the
apparent worst case; they are not used in the proof of the upper bound.
The line `Rigorous search upper bound` is the largest bound among the remaining
branch-and-bound boxes (and the safely pruned region), so the search mode also
returns a rigorous, generally coarser upper bound on the same maximum.

With the default witness resolution, a reference run gives the certified
two-sided conclusion
\[
3.1597523
<\max_{(r,m,\lambda)\in\mathcal D}
  \min\{\widehat C_\eta,\widehat C_{1/3}\}
<3.16.
\]
The numerical maximizer is approximately
\[
(r,m,\lambda)=(0.829526,\ 0.110834,\ 0.009166),
\]
where \(m+\lambda=3/25\) and the two cost bounds agree to the displayed
precision. The corresponding original normalized variables are approximately
\[
\begin{array}{rclcrcl}
R_1(V_0^\eta)&=&0.170474,
&&R_1(V_\eta^{1/3})&=&0.278692,\\
R_1(V_{\mathrm{single}})&=&0.440000,
&&R_1(V_{\mathrm{double}})&=&0.110834,\\
R_0(V_{\mathrm{single}})&=&0.889166,
&&R_0(V_{\mathrm{double}})&=&0.221668.
\end{array}
\]

## Build and run

The only nonstandard dependency is Boost (specifically
`Boost.Numeric.Interval`). With a C++20 compiler:

```bash
cmake -S interval_checker -B interval_checker/build -DCMAKE_BUILD_TYPE=Release
cmake --build interval_checker/build -j
./interval_checker/build/cvrp_interval_check
```

The two computations can be run separately:

```bash
./interval_checker/build/cvrp_interval_check --certificate-only
./interval_checker/build/cvrp_interval_check --search-only
```

## Visualizing fixed-r slices

The companion program `visualize_search_space.py` evaluates a triangular grid
in \((m,\lambda)\) for several fixed values of \(r\). It generates one colored
heatmap per slice and one combined multi-panel picture:

```bash
python3 interval_checker/visualize_search_space.py
```

By default, the pictures are written to `interval_checker/figures/`. The color
shows \(\min\{\widehat C_\eta,\widehat C_{1/3}\}\), the dashed white curve marks
where the two bounds agree, and the blue star marks the largest grid value in
that slice. The black diagonal is the boundary
\(m+\lambda=3/25\).

The fixed values of \(r\), grid resolution, and output directory are
configurable. For example:

```bash
python3 interval_checker/visualize_search_space.py \
  --r-values 0.75 0.82 0.829526 0.86 0.95 \
  --grid-size 1001 \
  --output-dir interval_checker/figures
```

This visualization uses the same analytic antiderivatives as the numerical
witness search, but ordinary floating-point arithmetic. It illustrates the
search landscape; only `cvrp_interval_check` provides rigorous bounds. The
plotting program additionally requires NumPy and Matplotlib.

Run `./interval_checker/build/cvrp_interval_check --help` for the subdivision,
search, threading, and reporting options. Increasing `--search-boxes` improves
the branch-and-bound search; increasing `--witness-bins` tightens the certified
lower value at the final witness.

The certificate assumes an IEEE-754 environment in which the C floating-point
rounding modes are honored by arithmetic operations and `log`, as required by
`Boost.Numeric.Interval`. CMake enables `-frounding-math`, disables
floating-point contraction, and forbids fast-math on GNU and Clang compilers.

## Exit status

- `0`: every requested computation completed, and, if requested, the upper
  bound was proved;
- `1`: invalid arguments or a runtime error;
- `2`: the configured certificate depth or evaluation limit was reached, so
  the rigorous upper check is inconclusive.
