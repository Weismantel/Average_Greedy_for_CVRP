# Interval Arithmetic for the 3.159-approximation for CVRP

This repository contains the machine check of the conditions of the lemma
"Reduction to a three-dimensional prism" in the paper, implemented in
`verify.py`. Every computation is done in exact rational arithmetic, except the
per-segment antiderivative — the one place a logarithm appears — which is
evaluated in Arb ball arithmetic and then rounded outward to a rational upper
bound, so rounding errors cannot invalidate the result.

A second script, `worst_case_search.py`, searches the opposite direction for a
worst-case instance. It is a convenience, not part of the certificate; see
[Looking for the worst case](#looking-for-the-worst-case) below.


## Problem statement

We prove

```
For every (r, w, λ) satisfying

    659/1000 ≤ r ≤ 1        w ≥ 0        λ ≥ 0        w + λ ≤ 91/750

it holds that

    min{ Ĉ_η(r, w, λ) − α·Opt ,  Ĉ_1/3(r, w, λ) − α·Opt }  <  1659/1000 .
```

Here `Ĉ_η` and `Ĉ_1/3` are the two cost bounds of the paper (see definitions below),
evaluated on an instance normalised to `Opt = 1` under the variable substitution
of the reduction lemma (see below). Because the target of the lemma is `ρ = α + 1.659`, this *is* the lemma's
hypothesis

```
    min{ Ĉ_η , Ĉ_1/3 }  <  ρ        on the prism.
```


### Notation

We need the following notation from the paper

| symbol | meaning |
| --- | --- |
| `Opt` | cost of an optimal solution. The instance is scaled so that `Opt = 1`. |
| `α` | the approximation ratio for TSP |
| `ρ = α + 1.659` | the target of the lemma. The paper writes the second summand as `β`. |
| `R₁(·)`, `R₀(·)` | the two quantities of the paper's variable list, on a set of customers. |
| `V_0^η`, `V_η^{1/3}`, `V_single`, `V_double` | the customer sets the substitution is stated for. |
| `V_η¹` | the customer set the integrals run over. In general `R₁(V_η¹) = R₁(V) − R₁(V_0^η)`; this is `Opt − R₁(V_0^η)` only on the face `R₁(V) = Opt`, which is where the reduction lemma puts every point and where the whole prism lies. |
| `σ_1/3(V_η¹)` | the `1/3`-potential of `V_η¹`. Abbreviated `σ` below. |
| `Δ` | shorthand used here for `R₀(V_single) − R₁(V_single)`. **Not** the paper's `Δ`; see the next table. |
| `r`, `w`, `λ` | the three coordinates of the prism. |
| `Ĉ_η`, `Ĉ_1/3` | the two cost bounds, with their error terms omitted. |

### Symbols that differ between the paper and this repository

| paper | here | meaning |
| --- | --- | --- |
| `Δ` | `D` | the triangular prism itself |
| — | `Δ`, `R0_minus_R1_single()` | `R₀(V_single) − R₁(V_single)` |
| `w` | `w` | the second prism coordinate. The paper's equation labels (`eq:m-lambda-…`) still spell it `m`, as did earlier versions of this code. |
| `β = 1.659` | `RHO_MINUS_ALPHA` | the target, `ρ − α` |
| `M` | — | the Matching Savings bound. Not needed for the certificate; `worst_case_search.py` does evaluate it. |

### The cost bound `Ĉ_η`

```
Ĉ_η  =  α·Opt  +  R₁(V_0^η)  +  ∫₀^{2·R₁(V_η¹)} min{ 1 , N(s) / D(s) } ds

    N(s) =  Opt  −  ( R₀(V_single) / (2·R₁(V_η¹)) ) · s

    D(s) =  2·R₁(V_η¹)  −  ( 1 + ( 2·R₀(V_single) − 2·R₁(V_single) ) / (2·R₁(V_η¹)) ) · s
```

Substituting `s = 2·R₁(V_η¹)·t` maps the integral onto `[0, 1]` and clears the
nested fractions. With `Δ = R₀(V_single) − R₁(V_single)`:

```
Ĉ_η  =  α·Opt  +  R₁(V_0^η)  +  2·R₁(V_η¹) · ∫₀¹ min{ 1 , q_η(t) } dt

                        Opt − R₀(V_single)·t
        q_η(t)  =  ──────────────────────────────────
                    2R₁(V_η¹) − (2R₁(V_η¹) + 2Δ)·t
```

### The cost bound `Ĉ_1/3`

The paper's integral runs over `[0, σ]` with `σ = σ_1/3(V_η¹)`. The substitution
`s = σ·t` maps it onto `[0, 1]`, which is the form `verify.py` uses:

```
Ĉ_1/3  =  α·Opt  +  (3/2)·R₁(V_0^η)  +  σ · ∫₀¹ min{ 1 , Φ⁽¹⁾(t) , Φ⁽²⁾(t) } dt
```

where the two cost density functions, written directly in `t`, are

```
              Opt − R₀(V_single)·t                                   Opt − R₀(V_single)·t
Φ⁽¹⁾(t) = ─────────────────────────────       Φ⁽²⁾(t) = ───────────────────────────────────────────────
               (1 − t)·σ  −  2·t·Δ                      (1 − t)·σ  −  (3/2)·t·Δ  −  (1/8)·R₀(V_double)
```

### The potential `σ_1/3(V_η^1)`

```
σ_1/3(V_η¹) = (3/2)·R₁(V_η^{1/3}) + 3·( R₁(V_single) + R₁(V_double) )
                                  − (1/2)·( R₀(V_single) + R₀(V_double) )
```

### The prism and the variable substitution

The lemma's region is the triangular prism

```
    659/1000 ≤ r ≤ 1        w ≥ 0        λ ≥ 0        w + λ ≤ 91/750
```

on which the cost bounds are evaluated with (after the normalisation `Opt = 1`)

| quantity | value on the prism |
| --- | --- |
| `R₁(V_0^η)` | `1 − r` |
| `R₁(V_η^{1/3})` | `r − 841/1500 + λ` |
| `R₁(V_single)` | `159/500 + w + λ` |
| `R₁(V_double)` | `91/375 − w − 2λ` |
| `R₀(V_single)` | `159/250 + 2w + 3λ` |
| `R₀(V_double)` | `91/125 − 4w − 6λ` |

The three quantities that the cost bounds need on top of these follow:

| derived quantity | definition | value on the prism |
| --- | --- | --- |
| `R₁(V_η¹)` | `R₁(V) − R₁(V_0^η)`, and `R₁(V) = Opt` here | `r` |
| `Δ` | `R₀(V_single) − R₁(V_single)` | `159/500 + w + 2λ` |
| `σ = σ_1/3(V_η¹)` | see above | `(3/2)·r + 159/1000 + w` |



### Quotient conventions

Every one of the three quotients above has the shape `(1 − p·t) / (a − b·t)`
once `Opt = 1`. Collecting the `t` terms of the denominators:

```
(1 − t)·σ − 2·t·Δ                              =  σ                     −  ( σ + 2Δ )·t
(1 − t)·σ − (3/2)·t·Δ − R₀(V_double)/8         =  σ − R₀(V_double)/8    −  ( σ + (3/2)Δ )·t
```

so the three branches are

| branch | `p` (numerator slope) | `a` (denominator constant) | `b` (denominator slope) |
| --- | --- | --- | --- |
| `Ĉ_η` | `R₀(V_single)` | `2·R₁(V_η¹)` | `2·R₁(V_η¹) + 2Δ` |
| `Φ⁽¹⁾` | `R₀(V_single)` | `σ` | `σ + 2Δ` |
| `Φ⁽²⁾` | `R₀(V_single)` | `σ − R₀(V_double)/8` | `σ + (3/2)·Δ` |

As stipulated in the paper, a fraction is `+∞` exactly when its denominator is
negative — the sign of the numerator does not matter, so a quotient with a
negative numerator is simply a negative value, and the pointwise minimum takes
it. Since that minimum already contains the constant branch `1`, the
implementation equivalently represents an inadmissible quotient by `1`. A
denominator of exactly `0` occurs at a single `t`, so it changes no integral;
`verify.py` declines the quotient there.

On a parameter box `verify.py` bounds `p` from below by `p₋` and uses
`(1 − p₋·t) / (a₋ − b₊·t)`. Besides a positive denominator, that relaxation
needs `1 − p₋·t ≥ 0` — for a negative numerator the worse case would be the
*largest* denominator, not the smallest. It holds automatically, because
`w + λ ≤ 91/750` gives

```
p₋ = 159/250 + 2w + 3λ  ≤  159/250 + 3·(91/750)  =  159/250 + 91/250  =  1
```

and hence `1 − p₋·t ≥ 1 − t ≥ 0` on `[0, 1]`.

## Running it

The only dependency is [`python-flint`](https://pypi.org/project/python-flint/)
(Python ≥ 3.10). Create the virtual environment once:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Then run the verification:

```bash
source .venv/bin/activate
python3 verify.py # ~80 s on 14 threads
```

Exit status is `0` when the bound is proved and `2` when a box could not be
certified; in that case the box and its bound are printed. `verify.py` refuses
to run under `python3 -O`, because it relies on its assertions. A run reports

```
PROVED: max min{C_eta - alpha, C_1/3 - alpha} <= 1.658999999711591 < 1659/1000
boxes: 5582748, max depth per task: 25, jobs: 14, elapsed: 82.3s
```

The bound, the box count and the depth are deterministic: the subdivision is
driven only by exact rational comparisons and the tasks share no state, so
`--jobs` does not affect them and a re-run reproduces all three exactly. Only
`elapsed` depends on the machine. Note that the depth is counted *below a task
box*, so the deepest box of a run sits `TASK_SPLITS + 25` bisections below the
initial box.

## Looking for the worst case

`worst_case_search.py` goes the other way: it uses ordinary floating-point
arithmetic to hunt for a feasible instance that pushes
`min{C_M, C_η, C_1/3} − α·Opt` as close to `1.659` as it can, over the six
original variables rather than the prism. It proves nothing — it only shows how
much slack the analysis has left, and it is the one place in the repository
where the Matching Savings bound `C_M` is evaluated.

```bash
python3 worst_case_search.py            # no dependencies beyond the standard library
```
