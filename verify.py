#!/usr/bin/env python3
"""
See README.md for the problem statement.
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

# The target of the reduction lemma. Bounds are compared against it with
# a strict "<" in exact rational arithmetic, so a successful run establishes the
# strict inequality
RHO_MINUS_ALPHA = fmpq(1659, 1000)

# Fails only if a box cannot be subdivided far enough; 2^90 is far beyond what
# the certificate needs.
MAX_DEPTH = 90

# The domain is first split into 2^TASK_SPLITS boxes, which are the units of
# parallel work.  Purely a scheduling parameter: their union is the initial box.
TASK_SPLITS = 12

# The triangular prism D of the reduction lemma:
#     659/1000 <= r <= 1,  m >= 0,  lambda >= 0,  m + lambda <= 91/750.
R_MIN = fmpq(659, 1000)
R_MAX = fmpq(1)
M_PLUS_LAMBDA_CAP = fmpq(91, 750)

# A box is the 6-tuple of exact rational endpoints (r_lo, r_hi, m_lo, m_hi,
# lambda_lo, lambda_hi).  INITIAL_BOX is a superset of D: the constraint
# m + lambda <= 91/750 is imposed in restrict(), the other bounds are exact.
INITIAL_BOX = (R_MIN, R_MAX, fmpq(0), M_PLUS_LAMBDA_CAP, fmpq(0), M_PLUS_LAMBDA_CAP)


class Inconclusive(Exception):
    """A box that could not be certified within MAX_DEPTH subdivisions."""

    def __init__(self, box):
        # The box is the only argument, so the exception survives the pickling
        # that carries it from a worker process back to the parent.
        super().__init__(box)
        self.box = box

#
# The variable transformation
#

def R1_V0_eta(r, m, lam):
    return 1 - r


def R1_V_eta_third(r, m, lam):
    return r - fmpq(841, 1500) + lam


def R1_V_single(r, m, lam):
    return fmpq(159, 500) + m + lam


def R1_V_double(r, m, lam):
    return fmpq(91, 375) - m - 2 * lam


def R0_V_single(r, m, lam):
    return fmpq(159, 250) + 2 * m + 3 * lam


def R0_V_double(r, m, lam):
    return fmpq(91, 125) - 4 * m - 6 * lam


def R1_V_eta_one(r, m, lam):
    return 1 - R1_V0_eta(r, m, lam)


def sigma_third(r, m, lam):
    return (
        fmpq(3, 2) * R1_V_eta_third(r, m, lam)
        + 3 * (R1_V_single(r, m, lam) + R1_V_double(r, m, lam))
        - fmpq(1, 2) * (R0_V_single(r, m, lam) + R0_V_double(r, m, lam))
    )


def Delta(r, m, lam):
    return R0_V_single(r, m, lam) - R1_V_single(r, m, lam)


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
# box.  With p >= p_lo, a >= a_lo and b <= b_hi, the true quotient satisfies,
# for t >= 0,
#
#     (1 - p t) / (a - b t) <= (1 - p_lo t) / (a_lo - b_hi t)
#
# wherever a_lo - b_hi t > 0 and 1 - p_lo t >= 0: the numerator only grows and
# the denominator only shrinks, while both stay non-negative.


# Compute an upper bound for the integral of
#     F(t) = min{1, (1 - p t) / (a - b t)}
# on the interval [l, r] for p >= p_lo, a >= a_lo, b <= b_hi.
def segment_upper_bound(p_lo, a_lo, b_hi, l, r):
    width = r - l

    # As in the paper, a fraction is infinity where its denominator is
    # negative. So the quotient may be used wherever a_lo - b_hi t > 0,
    # plus 1 - p_lo t >= 0 (needed for the relaxation above). Both are
    # affine in t, so it suffices to check the endpoints
    if b_hi <= 0:
        return width
    for t in (l, r):
        if not (a_lo - b_hi * t > 0 and 1 - p_lo * t >= 0):
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


def integral_upper(p_lo, branches):
    """Certified upper bound, as an exact rational, for int_0^1 F(t) dt on the
    parameter box described by p_lo and the branches."""
    cuts = cut_points(p_lo, branches)
    total = fmpq(0)
    for j in range(len(cuts) - 1):
        l, r = cuts[j], cuts[j + 1]
        # The constant branch 1 bounds F on every segment.
        best = r - l
        # The cuts are placed so that segment j is where branch j attains the
        # minimum; the remaining segment belongs to the constant branch. Note
        # that choosing the wrong branch can only be pessimistic
        if j < len(branches):
            a_lo, b_hi = branches[j]
            candidate = segment_upper_bound(p_lo, a_lo, b_hi, l, r)
            if candidate < best:
                best = candidate
        total += best
    assert total >= 0, "integral bound is negative"
    # F <= 1 on [0, 1], so the integral never exceeds 1.
    return total if total < 1 else fmpq(1)


# ---------------------------------------------------------------------------
# Bounding and subdivision of the parameter domain
# ---------------------------------------------------------------------------


def upper_bound(box, target=RHO_MINUS_ALPHA):
    """Certified upper bound, as an exact rational, for
    min{C_eta - alpha * Opt, C_1/3 - alpha * Opt} on the box.
    """
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = box

    # Every quantity of the paper is affine in (r, m, lambda) with coefficients
    # of constant sign, so each of its extremes over the box is attained at one
    # of the two corners below.  Which corner is stated for every use.
    lo = (r_lo, m_lo, lam_lo)
    hi = (r_hi, m_hi, lam_hi)

    # p = R_0(V_single) increases in m and in lambda; only its minimum is
    # needed, since it appears in a numerator that has to be over-estimated.
    p_lo = R0_V_single(*lo)

    # Delta = R_0(V_single) - R_1(V_single) increases in m and in lambda.  It
    # enters only the denominator slopes, which have to be over-estimated, so
    # only its maximum is needed.
    delta_max = Delta(*hi)

    # R_1(V_0^eta) = Opt - r decreases in r, so the lower corner maximizes it.
    R1_V0_eta_max = R1_V0_eta(*lo)

    # R_1(V_eta^1) = r increases in r.  Its minimum gives the denominator
    # constant of the eta quotient, which has to be under-estimated; its maximum
    # gives both the denominator slope and the factor in front of the integral,
    # which have to be over-estimated.
    R1_V_eta_one_min = R1_V_eta_one(*lo)
    R1_V_eta_one_max = R1_V_eta_one(*hi)

    # C_eta - alpha * Opt
    #   = R_1(V_0^eta)
    #     + 2 R_1(V_eta^1) int_0^1 min{1, (1 - R_0(V_single) t) / (a - b t)} dt,
    #   a = 2 R_1(V_eta^1),   b = 2 R_1(V_eta^1) + 2 Delta.
    branch_eta = (2 * R1_V_eta_one_min, 2 * R1_V_eta_one_max + 2 * delta_max)
    integral_eta = integral_upper(p_lo, [branch_eta])
    C_hat_eta_minus_alpha = R1_V0_eta_max + 2 * R1_V_eta_one_max * integral_eta
    if C_hat_eta_minus_alpha < target:
        return C_hat_eta_minus_alpha

    # sigma_{1/3}(V_eta^1) = 3/2 r + 159/1000 + m increases in r and in m.
    sigma_min = sigma_third(*lo)
    sigma_max = sigma_third(*hi)

    # R_0(V_double) = 91/125 - 4m - 6 lambda *decreases* in m and in lambda, so
    # its maximum sits at the *lower* corner.  It is subtracted in the
    # denominator constant of Phi^(2), which therefore needs that maximum.
    R0_double_max = R0_V_double(*lo)

    # C_1/3 - alpha * Opt
    #   = 3/2 R_1(V_0^eta)
    #     + sigma int_0^1 min{1, Phi^(1), Phi^(2)} dt,
    #   Phi^(1): a = sigma,                      b = sigma + 2 Delta,
    #   Phi^(2): a = sigma - R_0(V_double) / 8,  b = sigma + 3/2 Delta.
    # The branches are listed in the order in which they attain the minimum: at
    # t = 0 the denominator of Phi^(2) is smaller by R_0(V_double)/8 =
    # 91/1000 - m/2 - 3lambda/4 >= 0 on D (since m/2 + 3lambda/4 <=
    # 3(m + lambda)/4 <= 91/1000), and it decreases more slowly, by
    # Delta/2 = 159/1000 + m/2 + lambda per unit of t, so the two swap exactly
    # once.
    branches_third = [
        # Phi^(1)
        (sigma_min, sigma_max + 2 * delta_max),
        # Phi^(2)
        (sigma_min - R0_double_max / 8, sigma_max + fmpq(3, 2) * delta_max),
    ]
    integral_third = integral_upper(p_lo, branches_third)
    C_hat_third_minus_alpha = (
        fmpq(3, 2) * R1_V0_eta_max + sigma_max * integral_third
    )

    if C_hat_third_minus_alpha < C_hat_eta_minus_alpha:
        return C_hat_third_minus_alpha
    return C_hat_eta_minus_alpha


def restrict(box):
    """Drops boxes that are disjoint from D and shrinks the remaining ones to
    the constraint m + lambda <= 91/750.  Returns None if the box is disjoint
    from D."""
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = box

    # Every point of the box violates the constraint if already the smallest
    # sum in it exceeds 91/750.
    if m_lo + lam_lo > M_PLUS_LAMBDA_CAP:
        return None

    # Every feasible point of the box has m <= 91/750 - lambda_lo and
    # lambda <= 91/750 - m_lo, so these two contractions keep the whole
    # intersection of the box with D.  Neither can fall below the corresponding
    # lower endpoint, because m_lo + lambda_lo <= 91/750 was just checked.
    m_cap = M_PLUS_LAMBDA_CAP - lam_lo
    if m_cap < m_hi:
        m_hi = m_cap
    lam_cap = M_PLUS_LAMBDA_CAP - m_lo
    if lam_cap < lam_hi:
        lam_hi = lam_cap
    return (r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi)


# Reciprocals of the extents of D: 341/1000 in r, 91/750 in m and lambda.
_R_SCALE = 1 / (R_MAX - R_MIN)
_ML_SCALE = 1 / M_PLUS_LAMBDA_CAP


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
    """Proves that min{C_eta - alpha * Opt, C_1/3 - alpha * Opt} <
    RHO_MINUS_ALPHA on the whole intersection of `box` with D.  Returns
    (largest certified bound or None, boxes visited, maximum depth reached);
    raises Inconclusive if some sub-box could not be settled.

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
        if bound < RHO_MINUS_ALPHA:
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
    assert bound < RHO_MINUS_ALPHA, "accepted a bound that is not below the target"

    boxes = sum(visited for _, visited, _ in results)
    depth = max(d for _, _, d in results)
    print(
        "PROVED: max min{C_eta - alpha, C_1/3 - alpha} <= %s < %s"
        % (decimal_upper(bound), RHO_MINUS_ALPHA)
    )
    print(
        "boxes: %d, max depth: %d, jobs: %d, elapsed: %.1fs"
        % (boxes, depth, args.jobs, elapsed)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
