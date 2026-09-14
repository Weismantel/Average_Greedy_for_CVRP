# Interval Arithmetic for the 3.159-approximation for CVRP

This repository contains a python script to check the conditions of the lemma
"Reduction to a three-dimensional prism" in the paper. This is implemented in
`verify.py`. Every computation is done using rational arithmetic, except for
the logarithm, which is done using Arb balls to rule out rounding errors.
The file `selftest.py` only does redundant checks to catch
implementation bugs, and is not part of the certificate.


## Problem statement

We prove

```
For every (r, m, λ) satisfying

    659/1000 ≤ r ≤ 1        m ≥ 0        λ ≥ 0        m + λ ≤ 91/750

it holds that

    min{ Ĉ_η(r, m, λ) − α·Opt ,  Ĉ_1/3(r, m, λ) − α·Opt }  <  1659/1000 .
```

Here `Ĉ_η` and `Ĉ_1/3` are the two cost bounds of the paper (see definitions below),
evaluated on an instance normalised to `Opt = 1` under the variable substitution
of the reduction lemma (Section 4). Because the target of the lemma is `ρ = α + 1.659`, this *is* the lemma's
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
| `ρ = α + 1.659` | the target of the lemma. |
| `R₁(·)`, `R₀(·)` | the two quantities of the paper's variable list, on a set of customers. |
| `V_0^η`, `V_η^{1/3}`, `V_single`, `V_double` | the customer sets the substitution is stated for. |
| `V_η¹` | the customer set the integrals run over; `R₁(V_η¹) = Opt − R₁(V_0^η)`. |
| `σ_1/3(V_η¹)` | the `1/3`-potential of `V_η¹`. Abbreviated `σ` below. |
| `Δ` | shorthand used here for `R₀(V_single) − R₁(V_single)`. |
| `r`, `m`, `λ` | the three coordinates of the prism. |
| `Ĉ_η`, `Ĉ_1/3` | the two cost bounds, with their error terms omitted. |

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

```
Ĉ_1/3  =  α·Opt  +  (3/2)·R₁(V_0^η)  +  ∫₀^σ min{ Φ⁽¹⁾ , Φ⁽²⁾ , 1 } ds ,    σ = σ_1/3(V_η¹)
```

where, with `t = s / σ`, the two cost density functions are

```
             Opt − R₀(V_single)·t                          Opt − R₀(V_single)·t
Φ⁽¹⁾ = ─────────────────────────────       Φ⁽²⁾ = ───────────────────────────────────────────────
         (1 − t)·σ  −  2·t·Δ                        (1 − t)·σ  −  (3/2)·t·Δ  −  (1/8)·R₀(V_double)
```

The same substitution `t = s/σ` maps this integral onto `[0, 1]` as well:

```
Ĉ_1/3  =  α·Opt  +  (3/2)·R₁(V_0^η)  +  σ · ∫₀¹ min{ Φ⁽¹⁾ , Φ⁽²⁾ , 1 } dt
```

### The potential `σ_1/3(V_η^1)`

```
σ_1/3(V_η¹) = (3/2)·R₁(V_η^{1/3}) + 3·( R₁(V_single) + R₁(V_double) )
                                  − (1/2)·( R₀(V_single) + R₀(V_double) )
```

### The prism and the variable substitution

The lemma's region is the triangular prism

```
    659/1000 ≤ r ≤ 1        m ≥ 0        λ ≥ 0        m + λ ≤ 91/750
```

on which the cost bounds are evaluated with (after the normalisation `Opt = 1`)

| quantity | value on the prism |
| --- | --- |
| `R₁(V_0^η)` | `1 − r` |
| `R₁(V_η^{1/3})` | `r − 841/1500 + λ` |
| `R₁(V_single)` | `159/500 + m + λ` |
| `R₁(V_double)` | `91/375 − m − 2λ` |
| `R₀(V_single)` | `159/250 + 2m + 3λ` |
| `R₀(V_double)` | `91/125 − 4m − 6λ` |

The three quantities that the cost bounds need on top of these follow:

| derived quantity | definition | value on the prism |
| --- | --- | --- |
| `R₁(V_η¹)` | `Opt − R₁(V_0^η)` | `r` |
| `Δ` | `R₀(V_single) − R₁(V_single)` | `159/500 + m + 2λ` |
| `σ = σ_1/3(V_η¹)` | see 3.3 | `(3/2)·r + 159/1000 + m` |



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

As stipulated in the paper, a fraction is treated as `+∞` whenever its
denominator is non-positive. Since the pointwise minimum already contains the
constant branch `1`, the implementation equivalently represents an inadmissible
quotient by `1`.

## 8. Running it

```bash
source /path/to/venv/bin/activate     # provides python-flint; $VENV in the dev container
python3 verify.py                     # ~80 s on 14 threads
python3 verify.py --jobs 1            # same result, one process, ~3.5 min
```

Exit status is `0` when the bound is proved and `2` when a box could not be
certified; in that case the box and its bound are printed. `verify.py` refuses
to run under `python3 -O`, because it relies on its assertions. A reference run
reports

```
PROVED: max min{C_eta - alpha, C_1/3 - alpha} <= 1.658999999711591 < 1659/1000
boxes: 5582748, max depth: 25, jobs: 14, elapsed: 82.3s
```

