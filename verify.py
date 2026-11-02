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

# The paper's beta.  Bounds are compared against it with a strict "<" in exact
# rational arithmetic, so a successful run establishes the strict inequality
BETA = fmpq(1659, 1000)

# Fails only if a box cannot be subdivided far enough.  Depth is counted from
# the task box, not from INITIAL_BOX, so the effective limit is
# TASK_SPLITS + MAX_DEPTH bisections; 2^90 below a task is far beyond what the
# certificate needs.
MAX_DEPTH = 90

# The domain is first split into 2^TASK_SPLITS boxes, which are the units of
# parallel work.  Purely a scheduling parameter: their union is the initial box.
TASK_SPLITS = 12

# The triangular prism Delta of the reduction lemma:
#     659/1000 <= r <= 1,  w >= 0,  lambda >= 0,  w + lambda <= 91/750.
R_MIN = fmpq(659, 1000)
R_MAX = fmpq(1)
W_PLUS_LAMBDA_CAP = fmpq(91, 750)

# A box is the 6-tuple of exact rational endpoints (r_lo, r_hi, w_lo, w_hi,
# lambda_lo, lambda_hi).  INITIAL_BOX is a superset of Delta: the constraint
# w + lambda <= 91/750 is imposed in restrict(), the other bounds are exact.
INITIAL_BOX = (R_MIN, R_MAX, fmpq(0), W_PLUS_LAMBDA_CAP, fmpq(0), W_PLUS_LAMBDA_CAP)


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

def R1_V0_eta(r, w, lam):
    return 1 - r


def R1_V_eta_third(r, w, lam):
    return r - fmpq(841, 1500) + lam


def R1_V_single(r, w, lam):
    return fmpq(159, 500) + w + lam


def R1_V_double(r, w, lam):
    return fmpq(91, 375) - w - 2 * lam


def R0_V_single(r, w, lam):
    return fmpq(159, 250) + 2 * w + 3 * lam


def R0_V_double(r, w, lam):
    return fmpq(91, 125) - 4 * w - 6 * lam


def R1_V_eta_one(r, w, lam):
    """R_1(V_eta^1) = R_1(V) - R_1(V_0^eta).

    Equal to 1 - R_1(V_0^eta) because the reduction lemma moves every point
    onto the face R_1(V) = Opt = 1, on which all of Delta lies."""
    return 1 - R1_V0_eta(r, w, lam)


def sigma_third(r, w, lam):
    return (
        fmpq(3, 2) * R1_V_eta_third(r, w, lam)
        + 3 * (R1_V_single(r, w, lam) + R1_V_double(r, w, lam))
        - fmpq(1, 2) * (R0_V_single(r, w, lam) + R0_V_double(r, w, lam))
    )


def Theta(r, w, lam):
    """Theta = R_0(V_single) - R_1(V_single).

    A shorthand of this repository only: the paper always writes it out, inside
    D_eta, D_1 and D_2.  Deliberately not called Delta, which is the prism."""
    return R0_V_single(r, w, lam) - R1_V_single(r, w, lam)


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

    # Not an admissibility test: a negative b_hi is a denominator that *grows*
    # in t, which the relaxation above handles fine.  The closed form below
    # divides by b_hi and b_hi^2, so this guard only keeps a division by zero
    # off every path.  It is unreachable for this certificate, where b_hi is
    # 2 R_1(V_eta^1) + 2 Theta or sigma + c Theta and hence positive.
    if b_hi <= 0:
        return width

    # As in the paper, a fraction is infinity where its denominator is
    # negative. So the quotient may be used wherever a_lo - b_hi t > 0,
    # plus 1 - p_lo t >= 0 (needed for the relaxation above). Both are
    # affine in t, so it suffices to check the endpoints
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
    # Each segment contributes at most its own width and the widths sum to 1,
    # so total <= 1 already; the cap cannot bite and only records that F <= 1.
    return total if total < 1 else fmpq(1)


# ---------------------------------------------------------------------------
# Bounding and subdivision of the parameter domain
# ---------------------------------------------------------------------------


def upper_bound(box, target=BETA):
    """Certified upper bound, as an exact rational, for the supremum over the
    box of min{G_eta, G_1/3}.

    G_eta alone is returned as soon as it settles the box: the minimum of the
    two is bounded by either one, so a value already below `target` certifies
    the box even where G_1/3 would have been smaller.
    """
    r_lo, r_hi, w_lo, w_hi, lam_lo, lam_hi = box

    # Every quantity of the paper is affine in (r, w, lambda) with coefficients
    # of constant sign, so each of its extremes over the box is attained at one
    # of the two corners below.  Which corner is stated for every use.
    lo = (r_lo, w_lo, lam_lo)
    hi = (r_hi, w_hi, lam_hi)

    # p = R_0(V_single) increases in w and in lambda; only its minimum is
    # needed, since it appears in a numerator that has to be over-estimated.
    p_lo = R0_V_single(*lo)

    # Theta = R_0(V_single) - R_1(V_single) increases in w and in lambda.  It
    # enters only the denominator slopes, which have to be over-estimated, so
    # only its maximum is needed.
    theta_max = Theta(*hi)

    # R_1(V_0^eta) = Opt - r decreases in r, so the lower corner maximizes it.
    R1_V0_eta_max = R1_V0_eta(*lo)

    # R_1(V_eta^1) = r increases in r.  Its minimum gives the denominator
    # constant of the eta quotient, which has to be under-estimated; its maximum
    # gives both the denominator slope and the factor in front of the integral,
    # which have to be over-estimated.
    R1_V_eta_one_min = R1_V_eta_one(*lo)
    R1_V_eta_one_max = R1_V_eta_one(*hi)

    # G_eta = R_1(V_0^eta) + 2 R_1(V_eta^1) int_0^1 min{1, N/D_eta} dt, where
    # N/D_eta = (1 - R_0(V_single) t) / (a - b t) with
    #   a = 2 R_1(V_eta^1),   b = 2 R_1(V_eta^1) + 2 Theta.
    branch_D_eta = (2 * R1_V_eta_one_min, 2 * R1_V_eta_one_max + 2 * theta_max)
    integral_eta = integral_upper(p_lo, [branch_D_eta])
    G_eta = R1_V0_eta_max + 2 * R1_V_eta_one_max * integral_eta
    if G_eta < target:
        return G_eta

    # sigma_{1/3}(V_eta^1) = 3/2 r + 159/1000 + w increases in r and in w.
    sigma_min = sigma_third(*lo)
    sigma_max = sigma_third(*hi)

    # R_0(V_double) = 91/125 - 4w - 6 lambda *decreases* in w and in lambda, so
    # its maximum sits at the *lower* corner.  It is subtracted in the
    # denominator constant of D_2, which therefore needs that maximum.
    R0_double_max = R0_V_double(*lo)

    # G_1/3 = 3/2 R_1(V_0^eta) + sigma int_0^1 min{1, N/D_1, N/D_2} dt, where
    #   D_1: a = sigma,                      b = sigma + 2 Theta,
    #   D_2: a = sigma - R_0(V_double) / 8,  b = sigma + 3/2 Theta.
    # The branches are listed in the order in which they attain the minimum: at
    # t = 0 the denominator D_2 is smaller than D_1 by R_0(V_double)/8 =
    # 91/1000 - w/2 - 3lambda/4 >= 0 on Delta (since w/2 + 3lambda/4 <=
    # 3(w + lambda)/4 <= 91/1000), and it decreases more slowly, by
    # Theta/2 = 159/1000 + w/2 + lambda per unit of t, so the two swap exactly
    # once.
    branches_third = [
        # N/D_1
        (sigma_min, sigma_max + 2 * theta_max),
        # N/D_2
        (sigma_min - R0_double_max / 8, sigma_max + fmpq(3, 2) * theta_max),
    ]
    integral_third = integral_upper(p_lo, branches_third)
    G_third = fmpq(3, 2) * R1_V0_eta_max + sigma_max * integral_third

    if G_third < G_eta:
        return G_third
    return G_eta


def restrict(box):
    """Shrinks a box to the constraint w + lambda <= 91/750, and returns None if
    that constraint already excludes every point of it.

    This one constraint is all that has to be tested: every box descends from
    INITIAL_BOX by bisection, so its r range is contained in Delta's and its w
    and lambda are non-negative by construction."""
    r_lo, r_hi, w_lo, w_hi, lam_lo, lam_hi = box

    # Every point of the box violates the constraint if already the smallest
    # sum in it exceeds 91/750.
    if w_lo + lam_lo > W_PLUS_LAMBDA_CAP:
        return None

    # Every feasible point of the box has w <= 91/750 - lambda_lo and
    # lambda <= 91/750 - w_lo, so these two contractions keep the whole
    # intersection of the box with Delta.  Neither can fall below the
    # corresponding lower endpoint, because w_lo + lam_lo <= 91/750 was checked.
    w_cap = W_PLUS_LAMBDA_CAP - lam_lo
    if w_cap < w_hi:
        w_hi = w_cap
    lam_cap = W_PLUS_LAMBDA_CAP - w_lo
    if lam_cap < lam_hi:
        lam_hi = lam_cap
    return (r_lo, r_hi, w_lo, w_hi, lam_lo, lam_hi)


# The extents of Delta are 341/1000 in r and 91/750 in w and lambda; these are
# their reciprocals, so split() can compare the three widths on a common scale.
_R_SCALE = 1 / (R_MAX - R_MIN)
_WL_SCALE = 1 / W_PLUS_LAMBDA_CAP


def split(box):
    """Bisects the coordinate that is widest relative to its extent in Delta.
    The two halves share the midpoint, so their union is exactly the parent box;
    which coordinate is chosen is a pure heuristic.  Midpoints are exact
    rationals, so a non-degenerate box always splits into two non-degenerate
    halves."""
    r_lo, r_hi, w_lo, w_hi, lam_lo, lam_hi = box
    width_r = (r_hi - r_lo) * _R_SCALE
    width_w = (w_hi - w_lo) * _WL_SCALE
    width_lam = (lam_hi - lam_lo) * _WL_SCALE

    if width_r >= width_w and width_r >= width_lam:
        mid = (r_lo + r_hi) / 2
        return (
            (r_lo, mid, w_lo, w_hi, lam_lo, lam_hi),
            (mid, r_hi, w_lo, w_hi, lam_lo, lam_hi),
        )
    if width_w >= width_lam:
        mid = (w_lo + w_hi) / 2
        return (
            (r_lo, r_hi, w_lo, mid, lam_lo, lam_hi),
            (r_lo, r_hi, mid, w_hi, lam_lo, lam_hi),
        )
    mid = (lam_lo + lam_hi) / 2
    return (
        (r_lo, r_hi, w_lo, w_hi, lam_lo, mid),
        (r_lo, r_hi, w_lo, w_hi, mid, lam_hi),
    )


def prove(box):
    """Proves that min{G_eta, G_1/3} < BETA on the whole intersection of `box`
    with Delta.  Returns (largest certified bound or None, boxes visited,
    maximum depth reached below `box`); raises Inconclusive if some sub-box
    could not be settled.

    The traversal keeps an explicit stack of (box, depth).  Its invariant is
    that the union of the stacked boxes together with the already accepted ones
    covers box ∩ Delta: a popped box is either dropped by restrict() -- which
    only ever discards points outside Delta -- or accepted, or replaced by the
    two halves returned by split(), whose union is the popped box.
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
        if bound < BETA:
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
    # the bound on each proves it on INITIAL_BOX, hence on Delta.  prove() is a
    # pure function of one box and shares no state, so the only thing the
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
        for i, name in enumerate(("r", "w", "lambda")):
            lo, hi = failure.box[2 * i], failure.box[2 * i + 1]
            print("  %-6s = [%s, %s]" % (name, lo, hi), file=sys.stderr)
        print("  bound  = %s" % decimal_upper(upper_bound(failure.box)), file=sys.stderr)
        print("elapsed: %.1fs" % elapsed, file=sys.stderr)
        return 2
    elapsed = time.time() - start

    assert len(results) == len(tasks), "a worker did not return a result"
    bounds = [bound for bound, _, _ in results if bound is not None]
    assert bounds, "no box intersected Delta"
    bound = max(bounds)
    assert bound < BETA, "accepted a bound that is not below the target"

    boxes = sum(visited for _, visited, _ in results)
    # Depth is counted below a task box, so the deepest box of the run sits
    # TASK_SPLITS + depth bisections below INITIAL_BOX.
    depth = max(d for _, _, d in results)
    print(
        "PROVED: max min{G_eta, G_1/3} <= %s < %s"
        % (decimal_upper(bound), BETA)
    )
    print(
        "boxes: %d, max depth per task: %d, jobs: %d, elapsed: %.1fs"
        % (boxes, depth, args.jobs, elapsed)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
