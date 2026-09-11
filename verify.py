#!/usr/bin/env python3
"""Rigorous verification of the certificate of the reduction to a
three-dimensional prism (see README.md):

    max_{(r,m,lambda) in D} min{C_eta(r,m,lambda), C_1/3(r,m,lambda)} < 3159/1000,
    D = {659/1000 <= r <= 1, m >= 0, lambda >= 0, m + lambda <= 91/750},

with the functions as displayed in README.md.

How soundness is organized
--------------------------
Every coefficient of the problem is an affine function, with rational
coefficients, of the endpoints of the parameter box, and bisecting a rational
interval produces rational endpoints.  Therefore *all* of the following are
computed exactly, as `fmpq`, with no rounding whatsoever:

  * the box endpoints and the constraint m + lambda <= 91/750,
  * the branch data (p_lo, p_hi, a_lo, b_hi, s) and the cut points,
  * every comparison that steers the algorithm (branch admissibility,
    acceptance of a box, the final comparison against 3159/1000).

The one quantity that is not rational is the logarithm in the closed form of
the branch integral.  It -- and only it -- is evaluated in Arb ball arithmetic,
which is rigorous by construction, and the ball's upper endpoint is converted
back to an exact rational.  So the whole numerical-soundness argument is
confined to `upper_fmpq` and `segment_upper_bound` below; everything else is
exact arithmetic and needs no rounding analysis.

Usage:  python3 verify.py [--jobs N]      (exit 0 proved, 2 inconclusive)
"""

import argparse
import multiprocessing
import sys
import time

from flint import arb, ctx, fmpq

# Working precision of the ball arithmetic.  This influences only how tight the
# computed bounds are, never whether they are valid: an Arb ball encloses the
# exact value at any precision.
ctx.prec = 128

# The exact rational target rho = 3159/1000 of the reduction lemma.  Bounds
# are compared against it with a strict "<" in exact rational arithmetic, so a
# successful run establishes the strict inequality literally, with no rounding
# argument.
# TARGET = fmpq(79, 25)
TARGET = fmpq(3159, 1000)

# Fails only if a box cannot be subdivided far enough; 2^90 is far beyond what
# the certificate needs.
MAX_DEPTH = 90

# The domain is first split into 2^TASK_SPLITS boxes, which are the units of
# parallel work.  Purely a scheduling parameter: their union is the initial box.
TASK_SPLITS = 12

# A box is the 6-tuple of exact rational endpoints (r_lo, r_hi, m_lo, m_hi,
# lambda_lo, lambda_hi).  A superset of D: the constraint m + lambda <= 91/750
# is imposed in restrict(), the other bounds are exact.
CAP = fmpq(91, 750)
INITIAL_BOX = (fmpq(659, 1000), fmpq(1), fmpq(0), CAP, fmpq(0), CAP)


class Inconclusive(Exception):
    """A box that could not be certified within MAX_DEPTH subdivisions."""

    def __init__(self, box):
        # The box is the only argument, so the exception survives the pickling
        # that carries it from a worker process back to the parent.
        super().__init__(box)
        self.box = box


def upper_fmpq(x):
    """The exact rational value of the upper endpoint of the ball `x`.

    This is the only place where a ball is turned back into a number.
    `arb.upper()` is Arb's rigorous upper bound for the ball (rounded upwards),
    returned as an exact ball, and `man_exp()` reads off its exact value
    mantissa * 2^exponent.  Hence the result is >= every point of `x`, and in
    particular >= the exact quantity that `x` was built to enclose.
    """
    u = x.upper()
    assert u.is_exact(), "arb.upper() returned an inexact ball"
    man, exp = u.man_exp()
    exp = int(exp)
    two_pow = fmpq(2) ** abs(exp)
    return fmpq(man) * two_pow if exp >= 0 else fmpq(man) / two_pow


# ---------------------------------------------------------------------------
# Upper bound for the integral of the envelope on one parameter box
# ---------------------------------------------------------------------------
#
# A branch is the pair (a_lo, b_hi) of exact rational coefficient bounds of one
# quotient, already reduced to the form that bounds it on the whole parameter
# box.  With p in [p_lo, p_hi], a >= a_lo and b <= b_hi, the true quotient
# satisfies, for t >= 0,
#
#     (1 - p t) / (a - b t) <= (1 - p_lo t) / (a_lo - b_hi t)
#
# as soon as the right-hand denominator is positive: the numerator only grows
# and the denominator only shrinks, while staying positive.


def segment_upper_bound(p_lo, p_hi, a_lo, b_hi, l, r):
    """Certified upper bound, as an exact rational, for the integral over
    [l, r] of the true integrand

        F(t) = min{1, quotients of the branches that are admissible at t},

    using only this branch.  Returns the segment width -- always a valid bound,
    since F <= 1 -- unless the branch provably bounds F on all of [l, r].
    """
    width = r - l

    # The quotient may replace F only where numerator and denominator are
    # positive for *every* parameter of the box, i.e. where 1 - p_hi t > 0 and
    # a_lo - b_hi t > 0.  Note the asymmetry: admissibility is decided with
    # p_hi, the bound itself uses p_lo.  Both expressions are affine in t, so
    # positivity at the two endpoints certifies the whole segment.  These are
    # exact rational comparisons.
    if b_hi <= 0:
        return width
    for t in (l, r):
        if not (1 - p_hi * t > 0 and a_lo - b_hi * t > 0):
            return width

    # int_l^r (1 - p t) / (a - b t) dt
    #     = (p/b) (r - l) + ((b - a p) / b^2) log((a - b l) / (a - b r)).
    # Evaluated in ball arithmetic from exact rational inputs, so `value`
    # encloses the exact integral of the bounding quotient.
    p, a, b = arb(p_lo), arb(a_lo), arb(b_hi)
    lo, hi = arb(l), arb(r)
    value = (p / b) * (hi - lo) + ((b - a * p) / (b * b)) * (
        (a - b * lo) / (a - b * hi)
    ).log()

    # Non-finite only if the enclosure degenerated; fall back to the constant
    # branch then.  (The positivity checks above make this unreachable in
    # practice; it is kept so that no path can return an invalid number.)
    if not value.is_finite():
        return width

    bound = upper_fmpq(value)
    return bound if bound < width else width


def cut_points(p_lo, branches):
    """Cut points for the partition of [0, 1], as exact rationals.

    The branches are listed in the order in which they attain the minimum, and
    after the last of them the constant 1 does, so the envelope needs one cut
    per branch: consecutive branches swap where their denominators agree (they
    share their numerator, so the smaller quotient is the one with the larger
    denominator), and the last branch meets 1 where 1 - p_lo t = a - b t.  No
    cut is needed where a quotient stops being admissible: that happens after it
    has passed 1, where the constant branch is the minimum anyway.

    Correctness does not depend on any of this.  Whatever the cuts are, the
    bound below stays valid, because on every segment it takes the minimum of
    the constant 1 and of branches that are bounds for the whole segment; a
    misplaced cut costs tightness only.
    """
    cuts = [fmpq(0), fmpq(1)]

    def add(numerator, denominator):
        if denominator != 0:
            t = numerator / denominator
            if 0 < t < 1:
                cuts.append(t)

    for (a0, b0), (a1, b1) in zip(branches, branches[1:]):
        add(a0 - a1, b0 - b1)
    a_last, b_last = branches[-1]
    add(1 - a_last, p_lo - b_last)

    cuts.sort()
    return cuts


def integral_upper(p_lo, p_hi, branches):
    """Certified upper bound, as an exact rational, for int_0^1 F(t) dt on the
    parameter box described by p_lo, p_hi and the branches."""
    cuts = cut_points(p_lo, branches)
    total = fmpq(0)
    for j in range(len(cuts) - 1):
        l, r = cuts[j], cuts[j + 1]
        # The constant branch 1 bounds F on every segment.
        best = r - l
        # The cuts are placed so that segment j is where branch j attains the
        # minimum; the remaining segment belongs to the constant branch.  Only
        # this one branch is evaluated -- the others cannot improve the bound
        # here, and evaluating them would only cost logarithms.
        if j < len(branches):
            a_lo, b_hi = branches[j]
            candidate = segment_upper_bound(p_lo, p_hi, a_lo, b_hi, l, r)
            if candidate < best:
                best = candidate
        total += best
    assert total >= 0, "integral bound is negative"
    # F <= 1 on [0, 1], so the integral never exceeds 1.
    return total if total < 1 else fmpq(1)


# ---------------------------------------------------------------------------
# Subdivision of the parameter domain
# ---------------------------------------------------------------------------


def upper_bound(box, target=TARGET):
    """Certified upper bound, as an exact rational, for min{C_eta, C_1/3} on the
    box.  The eta bound alone is returned as soon as it settles the box, since
    the minimum of the two is bounded by either one."""
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = box

    # p = R_0(V_single) = 159/250 + 2m + 3lambda.
    p_lo = fmpq(159, 250) + 2 * m_lo + 3 * lam_lo
    p_hi = fmpq(159, 250) + 2 * m_hi + 3 * lam_hi

    # C_eta = 5/2 - r + 2r int_0^1 min{1, (1-pt) / (2r - bt)} dt with the slope
    # b = 2r + 2R_0(V_single) - 2R_1(V_single) = 2r + 159/250 + 2m + 4lambda.
    eta = (2 * r_lo, 2 * r_hi + fmpq(159, 250) + 2 * m_hi + 4 * lam_hi)
    eta_integral = integral_upper(p_lo, p_hi, [eta])
    eta_cost = fmpq(5, 2) - r_lo + 2 * r_hi * eta_integral
    if eta_cost < target:
        return eta_cost

    # C_1/3 = 3 - 3r/2 + s int_0^1 min{1, two quotients} dt with
    # s = sigma_1/3(V_eta^1) = 3r/2 + 159/1000 + m.  The constant term of the
    # second denominator is s - R_0(V_double)/8, with R_0(V_double) =
    # 91/125 - 4m - 6lambda; the branches are listed in the order in which
    # they attain the minimum: at t = 0 the second denominator is
    # smaller by 91/1000 - m/2 - 3lambda/4 >= 0 on D (since m/2 + 3lambda/4 <=
    # 3(m + lambda)/4 <= 91/1000), and it decreases more slowly, by
    # 159/1000 + m/2 + lambda per unit of t, so the two swap exactly once.
    s_lo = fmpq(3, 2) * r_lo + fmpq(159, 1000) + m_lo
    s_hi = fmpq(3, 2) * r_hi + fmpq(159, 1000) + m_hi
    branches = [
        (s_lo, fmpq(3, 2) * r_hi + fmpq(159, 200) + 3 * m_hi + 4 * lam_hi),
        (
            fmpq(3, 2) * r_lo + fmpq(17, 250) + fmpq(3, 2) * m_lo + fmpq(3, 4) * lam_lo,
            fmpq(3, 2) * r_hi + fmpq(159, 250) + fmpq(5, 2) * m_hi + 3 * lam_hi,
        ),
    ]
    one_third_integral = integral_upper(p_lo, p_hi, branches)
    one_third_cost = 3 - fmpq(3, 2) * r_lo + s_hi * one_third_integral

    return one_third_cost if one_third_cost < eta_cost else eta_cost


def restrict(box):
    """Drops boxes that are disjoint from D and shrinks the remaining ones to
    the constraint m + lambda <= 91/750.  Returns None if the box is disjoint
    from D."""
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = box

    # Every point of the box violates the constraint if already the smallest
    # sum in it exceeds 91/750.
    if m_lo + lam_lo > CAP:
        return None

    # Every feasible point of the box has m <= 91/750 - lambda_lo and
    # lambda <= 91/750 - m_lo, so these two contractions keep the whole
    # intersection of the box with D.  Neither can fall below the corresponding
    # lower endpoint, because m_lo + lambda_lo <= 91/750 was just checked.
    m_cap = CAP - lam_lo
    if m_cap < m_hi:
        m_hi = m_cap
    lam_cap = CAP - m_lo
    if lam_cap < lam_hi:
        lam_hi = lam_cap
    return (r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi)


# Reciprocals of the extents of D: 341/1000 in r, 91/750 in m and lambda.
_R_SCALE = fmpq(1000, 341)
_ML_SCALE = fmpq(750, 91)


def split(box):
    """Bisects the coordinate that is widest relative to its extent in D.  The
    two halves share the midpoint, so their union is exactly the parent box;
    which coordinate is chosen is a pure heuristic.  Midpoints are exact
    rationals, so a non-degenerate box always splits into two non-degenerate
    halves."""
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = box
    w_r = (r_hi - r_lo) * _R_SCALE
    w_m = (m_hi - m_lo) * _ML_SCALE
    w_lam = (lam_hi - lam_lo) * _ML_SCALE

    if w_r >= w_m and w_r >= w_lam:
        mid = (r_lo + r_hi) / 2
        return (
            (r_lo, mid, m_lo, m_hi, lam_lo, lam_hi),
            (mid, r_hi, m_lo, m_hi, lam_lo, lam_hi),
        )
    if w_m >= w_lam:
        mid = (m_lo + m_hi) / 2
        return (
            (r_lo, r_hi, m_lo, mid, lam_lo, lam_hi),
            (r_lo, r_hi, mid, m_hi, lam_lo, lam_hi),
        )
    mid = (lam_lo + lam_hi) / 2
    return (
        (r_lo, r_hi, m_lo, m_hi, lam_lo, mid),
        (r_lo, r_hi, m_lo, m_hi, mid, lam_hi),
    )


def prove(box):
    """Proves that min{C_eta, C_1/3} < TARGET on the whole intersection of `box`
    with D.  Returns (largest certified bound or None, boxes visited, maximum
    depth reached); raises Inconclusive if some sub-box could not be settled.

    The traversal keeps an explicit stack of (box, depth).  Its invariant is
    that the union of the stacked boxes together with the already accepted ones
    covers box ∩ D: a popped box is either dropped by restrict() -- which only
    ever discards points outside D -- or accepted, or replaced by the two halves
    returned by split(), whose union is the popped box.
    """
    stack = [(box, 0)]
    best = None
    visited = 0
    max_depth = 0
    while stack:
        current, depth = stack.pop()
        current = restrict(current)
        if current is None:
            continue
        visited += 1
        if depth > max_depth:
            max_depth = depth
        bound = upper_bound(current)
        if bound < TARGET:
            if best is None or bound > best:
                best = bound
            continue
        if depth >= MAX_DEPTH:
            raise Inconclusive(current)
        lower, upper = split(current)
        stack.append((lower, depth + 1))
        stack.append((upper, depth + 1))
    return best, visited, max_depth


def subdivide(box, splits):
    """The 2^splits boxes obtained by bisecting `box` repeatedly.  Since split()
    returns two halves whose union is its argument, the union of the returned
    boxes is exactly `box`."""
    boxes = [box]
    for _ in range(splits):
        halves = []
        for current in boxes:
            halves.extend(split(current))
        boxes = halves
    return boxes


def decimal_upper(value, digits=15):
    """A decimal string with `digits` fractional digits whose value is an exact
    rational >= `value`.  Used only for reporting, so that the printed number is
    itself a valid upper bound rather than a rounded float."""
    scaled = (value * 10**digits).ceil()
    sign = "-" if scaled < 0 else ""
    scaled = abs(int(scaled))
    text = str(scaled).rjust(digits + 1, "0")
    return "%s%s.%s" % (sign, text[:-digits], text[-digits:])


def main():
    # Several safety nets in this file are `assert`s, which python -O removes.
    if not __debug__:
        raise SystemExit("run without -O: this program relies on its assertions")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=multiprocessing.cpu_count(),
        help="worker processes (1 runs everything in this process)",
    )
    args = parser.parse_args()

    # The tasks are a partition of INITIAL_BOX into 2^TASK_SPLITS boxes; proving
    # the bound on each of them proves it on INITIAL_BOX, hence on D.  prove()
    # is a pure function of one box and shares no state, so the only thing the
    # parallel run must guarantee is that every task is actually proved: Pool.map
    # returns one result per task and re-raises any exception from a worker.
    tasks = subdivide(INITIAL_BOX, TASK_SPLITS)

    start = time.time()
    try:
        if args.jobs == 1:
            results = [prove(task) for task in tasks]
        else:
            with multiprocessing.Pool(args.jobs) as pool:
                results = pool.map(prove, tasks, chunksize=1)
    except Inconclusive as failure:
        elapsed = time.time() - start
        print("INCONCLUSIVE: a box could not be certified.", file=sys.stderr)
        for i, name in enumerate(("r", "m", "lambda")):
            lo, hi = failure.box[2 * i], failure.box[2 * i + 1]
            print("  %-6s = [%s, %s]" % (name, lo, hi), file=sys.stderr)
        print("  bound  = %s" % decimal_upper(upper_bound(failure.box)), file=sys.stderr)
        print("elapsed: %.1fs" % elapsed, file=sys.stderr)
        return 2
    elapsed = time.time() - start

    assert len(results) == len(tasks), "a worker did not return a result"
    bounds = [bound for bound, _, _ in results if bound is not None]
    assert bounds, "no box intersected D"
    bound = max(bounds)
    assert bound < TARGET, "accepted a bound that is not below the target"

    boxes = sum(visited for _, visited, _ in results)
    depth = max(d for _, _, d in results)
    print(
        "PROVED: max min{C_eta, C_1/3} <= %s < %s" % (decimal_upper(bound), TARGET)
    )
    print(
        "boxes: %d, max depth: %d, jobs: %d, elapsed: %.1fs"
        % (boxes, depth, args.jobs, elapsed)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
