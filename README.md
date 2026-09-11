# Computational certificate for the reduction to a three-dimensional prism

This directory contains a rigorous, reproducible proof of the hypothesis of the
lemma "Reduction to a three-dimensional prism" (`New_Reduction.txt`): with
\(\rho=\frac{3159}{1000}=3.159\),
\[
\max_{(r,m,\lambda)\in\mathcal D}
\min\{\widehat C_\eta(r,m,\lambda),\widehat C_{1/3}(r,m,\lambda)\}<\rho
\]
on the triangular prism
\[
\mathcal D=\left\{(r,m,\lambda):
\frac{659}{1000}\le r\le1,\quad m,\lambda\ge0,\quad
m+\lambda\le\frac{91}{750}\right\}.
\]

- `verify.py` is the proof. It is standalone: no file other than `verify.py`
  has to be read in order to check it.
- `selftest.py` cross-checks the transcription of the formulas in `verify.py`.
- `visualize_search_space.py` plots the search space; it has not been updated
  to this reduction yet (see the last section).

## The reduction

After the normalization \(\mathrm{OPT}=1\), the cost bounds are evaluated with
\[
\begin{array}{rclcrcl}
R_1(V_0^\eta)&=&1-r,
&&R_1(V_\eta^{1/3})&=&r-\frac{841}{1500}+\lambda,\\
R_1(V_{\mathrm{single}})&=&\frac{159}{500}+m+\lambda,
&&R_1(V_{\mathrm{double}})&=&\frac{91}{375}-m-2\lambda,\\
R_0(V_{\mathrm{single}})&=&\frac{159}{250}+2m+3\lambda,
&&R_0(V_{\mathrm{double}})&=&\frac{91}{125}-4m-6\lambda.
\end{array}
\]
Hence \(R_1(V_\eta^1)=1-R_1(V_0^\eta)=r\) and
\[
\begin{aligned}
\sigma_{1/3}(V_\eta^1)
&=\frac32R_1(V_\eta^{1/3})
 +3\bigl(R_1(V_{\mathrm{single}})+R_1(V_{\mathrm{double}})\bigr)
 -\frac12\bigl(R_0(V_{\mathrm{single}})+R_0(V_{\mathrm{double}})\bigr)\\
&=\frac32r+\frac{159}{1000}+m.
\end{aligned}
\]

## Functions checked

With \(\alpha=\frac32\) and the substitution \(x=2R_1(V_\eta^1)\,t=2rt\), the
definition of \(\widehat C_\eta\) reads
\[
\widehat C_\eta
=\alpha+R_1(V_0^\eta)+2r\int_0^1\min\left\{1,
\frac{1-R_0(V_{\mathrm{single}})\,t}
{2r-\bigl(2r+2R_0(V_{\mathrm{single}})-2R_1(V_{\mathrm{single}})\bigr)t}
\right\}dt .
\]
Likewise, with \(\sigma=\sigma_{1/3}(V_\eta^1)\), the substitution
\(x=\sigma t\), \(\sigma_{1/3}(V_0^\eta)=\frac32R_1(V_0^\eta)\) and
\(\Delta=R_0(V_{\mathrm{single}})-R_1(V_{\mathrm{single}})\),
\[
\widehat C_{1/3}
=\alpha+\frac32R_1(V_0^\eta)+\sigma\int_0^1\min\left\{1,
\frac{1-R_0(V_{\mathrm{single}})\,t}
{\sigma-\frac18R_0(V_{\mathrm{double}})-(\sigma+\frac32\Delta)t},
\frac{1-R_0(V_{\mathrm{single}})\,t}
{\sigma-(\sigma+2\Delta)t}
\right\}dt .
\]
Substituting the values above gives the two expressions that `verify.py`
checks:
\[
\begin{aligned}
\widehat C_\eta(r,m,\lambda)
={}&\frac52-r+2r\int_0^1\min\left\{1,
\frac{1-(\frac{159}{250}+2m+3\lambda)t}
{2r-(2r+\frac{159}{250}+2m+4\lambda)t}\right\}\,dt,
\\
\widehat C_{1/3}(r,m,\lambda)
={}&3-\frac32r+\left(\frac32r+\frac{159}{1000}+m\right)
\int_0^1\min\left\{1,
\frac{1-(\frac{159}{250}+2m+3\lambda)t}
{\frac32r+\frac{17}{250}+\frac32m+\frac34\lambda
 -(\frac32r+\frac{159}{250}+\frac52m+3\lambda)t},
\right.\\[-2mm]
&\hspace{57mm}\left.
\frac{1-(\frac{159}{250}+2m+3\lambda)t}
{\frac32r+\frac{159}{1000}+m
 -(\frac32r+\frac{159}{200}+3m+4\lambda)t}
\right\}\,dt.
\end{aligned}
\]
As stipulated in the paper, a fraction is treated as infinity when it becomes
negative. Since the pointwise minimum already contains the constant branch
\(1\), the implementation equivalently represents an inadmissible quotient
branch by \(1\). A quotient is evaluated only where its numerator and
denominator are positive.

## The rigorous certificate: `verify.py`

On a parameter box, each quotient has the form
\[
q(t)=\frac{1-pt}{a-bt}.
\]
If \(p\ge p_-\), \(a\ge a_-\) and \(b\le b_+\) throughout the box, then,
wherever \(1-p_+t>0\) and \(a_--b_+t>0\),
\[
q(t)\le \frac{1-p_-t}{a_--b_+t},
\]
and this upper function is integrated in closed form,
\[
\int_l^r\frac{1-pt}{a-bt}\,dt
=\frac pb(r-l)+\frac{b-ap}{b^2}
  \log\!\left(\frac{a-bl}{a-br}\right).
\]
The interval \([0,1]\) is cut where consecutive quotients swap and where the
last one meets the constant \(1\); on each segment the bound uses the quotient
that attains the minimum there, or \(1\). The domain \(\mathcal D\) is covered
by adaptive bisection, and a box is accepted as soon as its bound is below
\(\rho\).

Every coefficient of the problem is an affine function, with rational
coefficients, of the endpoints of the parameter box, and bisecting a rational
interval yields rational endpoints. So `verify.py` keeps the box endpoints, the
branch data \(p_-,a_-,b_+\), the cut points and the constraint
\(m+\lambda\le91/750\) as **exact rationals**
([python-flint](https://github.com/flintlib/python-flint)'s `fmpq`); every
decision the search makes -- branch admissibility, acceptance of a box, the
final comparison against \(\rho=3159/1000\) -- is an exact rational comparison,
with no rounding to analyse. The only quantity that is not rational is the
logarithm in the closed form above; it alone is evaluated in Arb ball
arithmetic (`flint.arb`, rigorous by construction) and its ball upper endpoint
is converted back to an exact rational. The whole numerical-soundness argument
is therefore confined to the functions `upper_fmpq` and `segment_upper_bound`.
Because the target is an exact rational and each box bound is compared against
it with a strict `<`, no rounding of the target has to be considered.

Parallelism is a static decomposition: the initial box is bisected into
\(2^{12}\) task boxes whose union is exactly the initial box, and each task is
proved independently in a worker process by the pure function `prove`. Running
with `--jobs 1` performs the identical computation in one process.

```bash
source /path/to/venv/bin/activate     # provides python-flint
python3 verify.py                     # ~80 s on 14 threads
python3 verify.py --jobs 1            # same result, one process, ~3 min
```

Exit status is `0` when the bound is proved and `2` when a box could not be
certified; in that case the box and its bound are printed. `verify.py` refuses
to run under `python3 -O`, because it relies on its assertions. A reference run
reports

```
PROVED: max min{C_eta, C_1/3} <= 3.158999999711591 < 3159/1000
boxes: 5582748, max depth: 25, jobs: 14, elapsed: 80.4s
```

## Cross-checks: `selftest.py`

`selftest.py` cross-checks `verify.py` against independently written,
deliberately naive implementations -- the closed-form antiderivative against a
rigorous interval Riemann sum, the per-box bound against a plain-float
re-transcription of the formulas above, the bound itself against a naive
Riemann sum of the envelope at interior points, and the task decomposition
against the domain \(\mathcal D\). It is a test of the transcription, not part
of the proof, and `verify.py` does not refer to it:

```bash
python3 selftest.py
```

## Visualizing fixed-r slices (not yet updated)

`visualize_search_space.py` evaluates \(\min\{\widehat C_\eta,\widehat C_{1/3}\}\)
on a triangular grid in \((m,\lambda)\) for several fixed values of \(r\). It
writes one heatmap per slice and one combined multi-panel picture to `figures/`
next to the script:

```bash
python3 visualize_search_space.py
python3 visualize_search_space.py --r-values 0.75 0.82 0.86 0.95 --grid-size 1001
```

The color shows \(\min\{\widehat C_\eta,\widehat C_{1/3}\}\), the dashed white
curve marks where the two bounds agree, the blue star marks the largest grid
value in the slice, and the black diagonal is the boundary of the triangle.

**The script still evaluates the previous reduction**: the domain
\(r\ge\frac{33}{50}\), \(m+\lambda\le\frac3{25}\), the previous constants and the
coefficient \(\frac14\) in front of \(R_0(V_{\mathrm{double}})\). It uses
ordinary floating-point arithmetic, illustrates the landscape only, and is not
part of the proof. It requires NumPy and Matplotlib.
