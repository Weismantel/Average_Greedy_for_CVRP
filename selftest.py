#!/usr/bin/env python3
"""Independent cross-checks for verify.py.

verify.py is exact wherever a decision is made and rigorous (Arb balls) where it
is not, so it cannot suffer from rounding errors.  What exactness cannot rule
out is a *transcription* error: a mistyped coefficient, a wrong endpoint, a
botched antiderivative.  This file checks exactly that, by comparing verify.py
against independently written, deliberately naive implementations:

  check_domain          that the searched region covers D and the target is
                        79/25, both spelled out again from the lemma;
  check_upper_fmpq      the one function that turns an Arb ball back into a
                        number, against the ball's own endpoints;
  check_antiderivative  the closed form of the branch integral, against a
                        rigorous interval Riemann sum that uses no
                        antiderivative at all;
  check_box_bound       the whole per-box bound, against a plain-float
                        re-transcription of the coefficients displayed in
                        README.md -- this is what catches a mistyped constant
                        or a swapped box endpoint, in either direction;
  check_upper_bound     the defining property, that the per-box bound really is
                        above min{C_eta, C_1/3} at points of the box, against a
                        naive Riemann sum of the envelope that shares no
                        structure with verify.py;
  check_partition       that the task boxes handed to the worker processes tile
                        the initial box, which is what makes the parallel run
                        equivalent to a serial one.

Nothing in verify.py depends on this file: reviewing the upper-bound proof means
reading verify.py alone.

Usage:  python3 selftest.py        (exit 0 if all checks passed)
"""

import math
import random
import sys

from flint import arb, ctx, fmpq

import verify

ctx.prec = 128


# The domain D of Lemma 4.1, written out here independently of verify.py.
R_MIN, R_MAX, SUM_MAX = fmpq(33, 50), fmpq(1), fmpq(3, 25)


def sample_point(rng):
    """A random point (r, m, lambda) of D, as exact rationals."""
    r = fmpq(rng.randrange(3300, 5000), 5000)
    m = fmpq(rng.randrange(0, 1201), 10000)
    lam = fmpq(rng.randrange(0, 1201 - int(m * 10000)), 10000)
    assert R_MIN <= r <= R_MAX and m >= 0 and lam >= 0 and m + lam <= SUM_MAX
    return r, m, lam


# ---------------------------------------------------------------------------
# Check 0: the domain and the target are the ones the lemma talks about
# ---------------------------------------------------------------------------


def check_domain(samples=2000, seed=20260907):
    """The searched region must cover D = {33/50 <= r <= 1, m, lambda >= 0,
    m + lambda <= 3/25} and the target must be 79/25."""
    assert verify.TARGET == fmpq(79, 25), "wrong target"
    assert verify.INITIAL_BOX[0] <= R_MIN and verify.INITIAL_BOX[1] >= R_MAX, "r range"
    assert verify.INITIAL_BOX[2] <= 0 and verify.INITIAL_BOX[4] <= 0, "m, lambda start above 0"
    assert (verify.INITIAL_BOX[3] >= SUM_MAX
            and verify.INITIAL_BOX[5] >= SUM_MAX), "m, lambda stop below 3/25"

    # restrict() may only remove points that are outside D: every point of D
    # that lies in a box must still lie in the restricted box.
    rng = random.Random(seed)
    for _ in range(samples):
        r, m, lam = sample_point(rng)
        pad = [fmpq(rng.randrange(0, 2**20), 2**20) / 8 for _ in range(6)]
        box = verify.restrict((r - pad[0], r + pad[1], m - pad[2], m + pad[3],
                               lam - pad[4], lam + pad[5]))
        assert box is not None, "restrict() dropped a box containing a point of D"
        assert (box[0] <= r <= box[1] and box[2] <= m <= box[3]
                and box[4] <= lam <= box[5]), (
            "restrict() cut away a point of D: r=%s m=%s lambda=%s -> %s"
            % (r, m, lam, box)
        )
    print("ok  domain: target 79/25, %d boxes keep their points of D" % samples)


# ---------------------------------------------------------------------------
# Check 1: upper_fmpq, the single bridge from balls back to exact numbers
# ---------------------------------------------------------------------------


def ball_endpoints(x):
    """The exact rational midpoint and radius of the ball `x`, read off its
    floating-point representation without going through arb.upper()."""
    man, exp = x.mid().man_exp()
    mid = fmpq(man) * fmpq(2) ** int(exp)
    man, exp = x.rad().man_exp()
    rad = fmpq(man) * fmpq(2) ** int(exp)
    return mid, rad


def check_upper_fmpq(samples=3000, seed=20260907):
    """upper_fmpq(x) must be >= every point of the ball x, i.e. >= mid + rad,
    at any working precision."""
    rng = random.Random(seed)
    for precision in (8, 20, 53, 128, 333):
        ctx.prec = precision
        for _ in range(samples):
            x = arb(fmpq(rng.randrange(-10**9, 10**9), rng.randrange(1, 10**9)))
            x = x * arb(fmpq(rng.randrange(1, 10**6), rng.randrange(1, 10**6)))
            mid, rad = ball_endpoints(x)
            assert verify.upper_fmpq(x) >= mid + rad, (
                "upper_fmpq is below the ball at precision %d: %s" % (precision, x)
            )
    ctx.prec = 128
    print("ok  upper_fmpq: %d balls at 5 precisions" % (5 * samples))


# ---------------------------------------------------------------------------
# Check 2: the closed-form antiderivative used by segment_upper_bound
# ---------------------------------------------------------------------------


def quotient_integral_enclosure(p, a, b, l, r, pieces=4000):
    """A rigorous enclosure, as an Arb ball, of

        int_l^r (1 - p t) / (a - b t) dt,

    computed without any antiderivative: the integrand is enclosed by ball
    arithmetic over each of `pieces` subintervals (an interval Riemann sum), so
    the result is guaranteed to contain the exact integral.
    """
    p, a, b = arb(p), arb(a), arb(b)
    total = arb(0)
    step = (r - l) / pieces
    for i in range(pieces):
        lo = arb(l + step * i)
        hi = arb(l + step * (i + 1))
        t = lo.union(hi)  # a ball containing the whole subinterval
        total += ((1 - p * t) / (a - b * t)) * (hi - lo)
    return total


def check_antiderivative(samples=60, seed=20260907):
    """segment_upper_bound returns min(integral of the quotient, segment width),
    so it must lie between the two endpoints of an independent enclosure of that
    integral."""
    rng = random.Random(seed)
    analytic_used = 0
    for _ in range(samples):
        # Parameters of the size the certificate actually meets, with a margin
        # that keeps numerator and denominator positive.
        p = fmpq(rng.randrange(30, 200), 100)
        a = fmpq(rng.randrange(120, 300), 100)
        b = fmpq(rng.randrange(10, 100), 100)
        l = fmpq(rng.randrange(0, 40), 100)
        r = l + fmpq(rng.randrange(5, 50), 100)
        if not (1 - p * r > 0 and a - b * r > 0):
            continue

        got = verify.segment_upper_bound(p, p, a, b, l, r)

        enclosure = quotient_integral_enclosure(p, a, b, l, r)
        low = -verify.upper_fmpq(-enclosure)  # exact lower endpoint of the ball
        high = verify.upper_fmpq(enclosure)
        width = r - l
        expected_low, expected_high = min(low, width), min(high, width)

        assert expected_low <= got <= expected_high, (
            "closed form disagrees with quadrature: p=%s a=%s b=%s on [%s, %s]: "
            "got %s, expected in [%s, %s]"
            % (p, a, b, l, r, got, expected_low, expected_high)
        )
        if got < width:
            analytic_used += 1

    assert analytic_used >= samples // 2, "too few samples exercised the closed form"
    print("ok  antiderivative: %d samples, %d exercised the closed form"
          % (samples, analytic_used))


# ---------------------------------------------------------------------------
# Check 3: the coefficients, re-transcribed from README.md in plain floats
# ---------------------------------------------------------------------------
#
# This mirrors the *algorithm* of verify.py deliberately -- it is a differential
# test of the numbers, not of the method.  Every constant below was read off the
# display in README.md, and every choice of box endpoint follows the rule stated
# there: the bounding quotient uses the smallest numerator coefficient p_lo, the
# smallest denominator constant a_lo and the largest denominator slope b_hi.


def float_integral_upper(p_lo, p_hi, branches):
    cuts = [0.0, 1.0]
    for (a0, b0), (a1, b1) in zip(branches, branches[1:]):
        if b0 != b1:
            t = (a0 - a1) / (b0 - b1)
            if 0.0 < t < 1.0:
                cuts.append(t)
    a_last, b_last = branches[-1]
    if p_lo != b_last:
        t = (1.0 - a_last) / (p_lo - b_last)
        if 0.0 < t < 1.0:
            cuts.append(t)
    cuts.sort()

    total = 0.0
    for j in range(len(cuts) - 1):
        l, r = cuts[j], cuts[j + 1]
        best = r - l
        if j < len(branches):
            a, b = branches[j]
            usable = b > 0.0 and all(
                1.0 - p_hi * t > 0.0 and a - b * t > 0.0 for t in (l, r)
            )
            if usable:
                value = (p_lo / b) * (r - l) + ((b - a * p_lo) / (b * b)) * math.log(
                    (a - b * l) / (a - b * r)
                )
                best = min(best, value)
        total += best
    return min(total, 1.0)


def float_upper_bound(box):
    """min{C_eta, C_1/3} bound for a box, in plain floats, from README.md."""
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = (float(x) for x in box)

    # p = R_0(V_single) = 16/25 + 2m + 3 lambda
    p_lo = 16 / 25 + 2 * m_lo + 3 * lam_lo
    p_hi = 16 / 25 + 2 * m_hi + 3 * lam_hi

    eta = (2 * r_lo, 2 * r_hi + 16 / 25 + 2 * m_hi + 4 * lam_hi)
    eta_cost = 5 / 2 - r_lo + 2 * r_hi * float_integral_upper(p_lo, p_hi, [eta])

    s_lo = 3 / 2 * r_lo + 4 / 25 + m_lo
    s_hi = 3 / 2 * r_hi + 4 / 25 + m_hi
    branches = [
        (s_lo, 3 / 2 * r_hi + 4 / 5 + 3 * m_hi + 4 * lam_hi),
        (3 / 2 * r_lo - 1 / 50 + 2 * m_lo + 3 / 2 * lam_lo,
         3 / 2 * r_hi + 16 / 25 + 5 / 2 * m_hi + 3 * lam_hi),
    ]
    one_third_cost = 3 - 3 / 2 * r_lo + s_hi * float_integral_upper(p_lo, p_hi, branches)

    return min(eta_cost, one_third_cost)


def check_box_bound(samples=400, seed=20260907):
    """verify.upper_bound must agree with the float re-transcription, on
    degenerate boxes as well as on boxes of every size that occurs during the
    subdivision.  A single mistyped constant or endpoint shows up here."""
    rng = random.Random(seed)
    worst = 0.0
    for i in range(samples):
        r, m, lam = sample_point(rng)
        # Widths from 0 (a point) up to the whole domain; every fourth sample is
        # anchored at the corner m = lambda = 0, where wide boxes make p_hi
        # exceed 1 and the admissibility test in t actually bites.
        width = fmpq(0) if i % 4 == 1 else fmpq(1, 2 ** rng.randrange(1, 30))
        if i % 4 == 0:
            m = lam = fmpq(0)
            r = R_MIN
        box = verify.restrict((r, r + width, m, m + width, lam, lam + width))
        if box is None:
            continue

        # A target no cost can reach disables the early return of the eta bound,
        # so that both cost functions are evaluated.
        got = float(verify.upper_bound(box, fmpq(-1)))
        expected = float_upper_bound(box)
        assert abs(got - expected) <= 1e-9 * abs(expected), (
            "per-box bound disagrees with the float transcription: box=%s "
            "got %.15f, expected %.15f" % (box, got, expected)
        )
        worst = max(worst, abs(got - expected))

    print("ok  coefficients: %d boxes, largest disagreement %.2e" % (samples, worst))


# ---------------------------------------------------------------------------
# Check 4: the bound really is an upper bound for the envelope
# ---------------------------------------------------------------------------
#
# Here nothing is shared with verify.py: the two cost functions are evaluated by
# a naive Riemann sum of the envelope min{1, quotients}, with no cut points, no
# branch ordering and no antiderivative.


def envelope_integral(quotients, pieces):
    """Midpoint Riemann sum of  t -> min{1, admissible quotients}  over [0, 1].

    Following README.md, a quotient counts only where its numerator and its
    denominator are positive; elsewhere the fraction is treated as infinity,
    which the constant branch 1 absorbs.
    """
    total = 0.0
    for i in range(pieces):
        t = (i + 0.5) / pieces
        value = 1.0
        for p, a, b in quotients:
            numerator = 1.0 - p * t
            denominator = a - b * t
            if numerator > 0.0 and denominator > 0.0:
                value = min(value, numerator / denominator)
        total += value
    return total / pieces


def cost_eta(r, m, lam, pieces=100000):
    """C_eta as displayed in README.md."""
    p = 16 / 25 + 2 * m + 3 * lam
    quotients = [(p, 2 * r, 2 * r + 16 / 25 + 2 * m + 4 * lam)]
    return 5 / 2 - r + 2 * r * envelope_integral(quotients, pieces)


def cost_one_third(r, m, lam, pieces=100000):
    """C_1/3 as displayed in README.md."""
    p = 16 / 25 + 2 * m + 3 * lam
    quotients = [
        (p, 3 / 2 * r - 1 / 50 + 2 * m + 3 / 2 * lam,
         3 / 2 * r + 16 / 25 + 5 / 2 * m + 3 * lam),
        (p, 3 / 2 * r + 4 / 25 + m,
         3 / 2 * r + 4 / 5 + 3 * m + 4 * lam),
    ]
    return 3 - 3 / 2 * r + (3 / 2 * r + 4 / 25 + m) * envelope_integral(quotients, pieces)


def check_upper_bound(samples=60, seed=20260907):
    """For a box and a point inside it, verify.upper_bound(box) must be at least
    min{C_eta, C_1/3} at that point -- with the real target as well as without
    the early return.  The Riemann sum is accurate to about 1e-5."""
    rng = random.Random(seed)
    tolerance = 1e-5
    worst = float("inf")
    for i in range(samples):
        r, m, lam = sample_point(rng)
        width = fmpq(1, 2 ** rng.randrange(4, 24))
        box = verify.restrict((r - width, r + width, m, m + width, lam, lam + width))
        if box is None or box[0] < verify.INITIAL_BOX[0]:
            continue

        point = min(cost_eta(float(r), float(m), float(lam)),
                    cost_one_third(float(r), float(m), float(lam)))
        for target in (verify.TARGET, fmpq(-1)):
            bound = float(verify.upper_bound(box, target))
            assert bound >= point - tolerance, (
                "bound is below the envelope at an interior point: r=%s m=%s "
                "lambda=%s target=%s bound=%.9f point=%.9f"
                % (r, m, lam, target, bound, point)
            )
        worst = min(worst, float(verify.upper_bound(box, fmpq(-1))) - point)

    assert math.isfinite(worst), "no sample produced a usable box"
    print("ok  upper bound: %d boxes, smallest margin over an interior point %.2e"
          % (samples, worst))


# ---------------------------------------------------------------------------
# Check 5: the decomposition the parallel run relies on
# ---------------------------------------------------------------------------


def volume(box):
    r_lo, r_hi, m_lo, m_hi, lam_lo, lam_hi = box
    return (r_hi - r_lo) * (m_hi - m_lo) * (lam_hi - lam_lo)


def check_partition(splits=7, points=2000, seed=20260907):
    """The task boxes must tile INITIAL_BOX: each is contained in it, their
    volumes add up exactly, and every sampled point of INITIAL_BOX lies in at
    least one of them."""
    tasks = verify.subdivide(verify.INITIAL_BOX, splits)
    assert len(tasks) == 2**splits

    total = fmpq(0)
    for box in tasks:
        for i in range(3):
            assert (verify.INITIAL_BOX[2 * i] <= box[2 * i] <= box[2 * i + 1]
                    <= verify.INITIAL_BOX[2 * i + 1]), "task box escapes the initial box"
        total += volume(box)
    assert total == volume(verify.INITIAL_BOX), "task volumes do not add up"

    rng = random.Random(seed)
    for _ in range(points):
        point = [
            verify.INITIAL_BOX[2 * i]
            + (verify.INITIAL_BOX[2 * i + 1] - verify.INITIAL_BOX[2 * i])
            * fmpq(rng.randrange(0, 10**6), 10**6)
            for i in range(3)
        ]
        covered = any(
            all(box[2 * i] <= point[i] <= box[2 * i + 1] for i in range(3))
            for box in tasks
        )
        assert covered, "point of the initial box lies in no task: %s" % (point,)

    print("ok  partition: %d tasks tile the initial box, %d points covered"
          % (len(tasks), points))


def main():
    check_domain()
    check_upper_fmpq()
    check_antiderivative()
    check_box_bound()
    check_upper_bound()
    check_partition()
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
