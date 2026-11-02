# Interval arithmetic for the α + 1.659 approximation for CVRP

`verify.py` machine-checks the hypothesis of the lemma *Reduction to a
three-dimensional prism*: that

```
    max_Δ min{ G_η , G_1/3 }  <  β = 1659/1000
```

over the triangular prism Δ, with `G_η` and `G_1/3` evaluated under that lemma's
variable substitution and `Opt` normalised to 1. In terms of the numerator
`N` and the denominators `D_η`, `D_1`, `D_2` defined in the paper,

```
    G_η    =        R₁(V_0^η)  +  2R₁(V_η¹) · ∫₀¹ min{ 1, N/D_η } dt
    G_1/3  =  (3/2)·R₁(V_0^η)  +  σ_1/3(V_η¹) · ∫₀¹ min{ 1, N/D_1, N/D_2 } dt
```

The implementation iteratively splits the domain into small boxes and computes
an upper bound on each box.
Every computation is exact rational arithmetic, except the per-segment
antiderivative (that contains a logarithm) which is evaluated in Arb
ball arithmetic and then rounded outward to a rational upper bound.


## Running it

The only dependency is [`python-flint`](https://pypi.org/project/python-flint/)
(Python ≥ 3.10).

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python3 verify.py          # ~80 s on 14 threads
```

Exit status is `0` when the bound is proved and `2` when a box could not be
certified, in which case that box and its bound are printed. `verify.py` refuses
to run under `-O`, because it relies on its assertions. A run reports

```
PROVED: max min{G_eta, G_1/3} <= 1.658999999711591 < 1659/1000
boxes: 5582748, max depth per task: 25, jobs: 14, elapsed: 82.3s
```

## Implementation Details

### Coefficient form

With `Opt = 1` each cost density is a ratio of two affine functions of `t`, so
the code stores it as `(1 − p·t) / (a − b·t)`. Writing `σ := σ_1/3(V_η¹)` and,
as a shorthand of this repository only, `Θ := R₀(V_single) − R₁(V_single)`:

| density | `p` | `a` | `b` |
| --- | --- | --- | --- |
| `N/D_η` | `R₀(V_single)` | `2R₁(V_η¹)` | `2R₁(V_η¹) + 2Θ` |
| `N/D_1` | `R₀(V_single)` | `σ` | `σ + 2Θ` |
| `N/D_2` | `R₀(V_single)` | `σ − R₀(V_double)/8` | `σ + (3/2)·Θ` |

`D_η` uses `σ_η(V_η¹)` raised to `2R₁(V_η¹)`, as the paper does for `C_η`.


### The box relaxation

On a box of `(r, w, λ)` the code over-estimates each density by

```
    (1 − p·t) / (a − b·t)  ≤  (1 − p₋·t) / (a₋ − b₊·t) ,
```

where `p₋`, `a₋` bound `p`, `a` from below and `b₊` bounds `b` from above.
Besides `a₋ − b₊·t > 0`, this needs `1 − p₋·t ≥ 0` — for a negative numerator
the worse case would be the *largest* denominator, not the smallest — and that
holds throughout Δ, since `w + λ ≤ 91/750` gives

```
    p₋ = 159/250 + 2w + 3λ  ≤  159/250 + 3·(91/750)  =  1 .
```

Each of `R₁(V_0^η)`, `R₁(V_η¹)`, `σ`, `Θ` and `R₀(V_double)` is affine in
`(r, w, λ)` with coefficients of one sign, so each of those bounds is attained
at a corner of the box; `verify.py` names the corner at every use.

A density that is inadmissible on a segment is represented there by the constant
`1` rather than by `+∞`: the pointwise minimum already contains `1`, so the
envelope is unchanged. The integral of the envelope is bounded segment by
segment in closed form, and a box whose bound does not fall below `β` is
bisected.


