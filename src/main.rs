// Rigorous verification of the upper bound of Lemma 4.1:
//
//     max_{(r,m,lambda) in D} min{C_eta(r,m,lambda), C_1/3(r,m,lambda)} < 79/25,
//     D = {33/50 <= r <= 1, m >= 0, lambda >= 0, m + lambda <= 3/25},
//
// with the functions as displayed in README.md.  Every floating-point operation
// below is an interval operation of `inari` (IEEE 1788), so each printed or
// compared bound encloses the exact value.

use std::process::exit;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use inari::{const_interval, interval, Interval};

const ONE: Interval = const_interval!(1.0, 1.0);
const TWO: Interval = const_interval!(2.0, 2.0);
const THREE: Interval = const_interval!(3.0, 3.0);
const FOUR: Interval = const_interval!(4.0, 4.0);

/// Fails only if the domain box cannot be subdivided far enough; 2^90 is far
/// beyond what the certificate needs.
const MAX_DEPTH: u32 = 90;

static BOXES: AtomicU64 = AtomicU64::new(0);

/// The degenerate interval [x, x].
fn pt(x: f64) -> Interval {
    interval!(x, x).unwrap()
}

/// An enclosure of the exact rational number n/d.
fn rat(n: f64, d: f64) -> Interval {
    pt(n) / pt(d)
}

// ---------------------------------------------------------------------------
// Upper bound for the integral of the envelope on one parameter box
// ---------------------------------------------------------------------------

/// One quotient branch, already reduced to the form that bounds it on the whole
/// parameter box.  With p in [p_lo, p_hi], a >= a_lo and b <= b_hi, the true
/// quotient satisfies, for t >= 0,
///
///     (1 - p t) / (a - b t) <= (1 - p_lo t) / (a_lo - b_hi t)
///
/// as soon as the right-hand denominator is positive: the numerator only grows
/// and the denominator only shrinks, while staying positive.
#[derive(Clone, Copy)]
struct Branch {
    a_lo: f64,
    b_hi: f64,
}

/// Certified upper bound for the integral over [l, r] of the true integrand
///
///     F(t) = min{1, quotients of the branches that are admissible at t},
///
/// using only this branch.  Returns the segment width -- always a valid bound,
/// since F <= 1 -- unless the branch provably bounds F on all of [l, r].
fn branch_segment_upper(p_lo: f64, p_hi: f64, br: Branch, l: f64, r: f64) -> f64 {
    let width = (pt(r) - pt(l)).sup();
    let p = pt(p_lo);
    let a = pt(br.a_lo);
    let b = pt(br.b_hi);

    // The quotient may replace F only where numerator and denominator are
    // positive for *every* parameter of the box, i.e. where 1 - p_hi t > 0 and
    // a_lo - b_hi t > 0.  Note the asymmetry: admissibility is decided with
    // p_hi, the bound itself uses p_lo.  Both expressions are affine in t, so
    // positivity at the two endpoints certifies the whole segment.
    let admissible = |t: f64| (ONE - pt(p_hi) * pt(t)).inf() > 0.0 && (a - b * pt(t)).inf() > 0.0;
    if br.b_hi <= 0.0 || !admissible(l) || !admissible(r) {
        return width;
    }

    // int_l^r (1 - p t) / (a - b t) dt
    //     = (p/b) (r - l) + ((b - a p) / b^2) log((a - b l) / (a - b r)).
    let value = (p / b) * (pt(r) - pt(l))
        + ((b - a * p) / (b * b)) * (((a - b * pt(l)) / (a - b * pt(r))).ln());

    let bound = value.sup();
    // `bound` is NaN only if the enclosure came out empty; fall back to the
    // constant branch then.
    if bound.is_nan() {
        width
    } else {
        bound.min(width)
    }
}

/// Cut points for the partition of [0, 1], as exact floating-point numbers.
///
/// The branches are listed in the order in which they attain the minimum, and
/// after the last of them the constant 1 does, so the envelope needs one cut
/// per branch: consecutive branches swap where their denominators agree (they
/// share their numerator, so the smaller quotient is the one with the larger
/// denominator), and the last branch meets 1 where 1 - p_lo t = a - b t.  No
/// cut is needed where a quotient stops being admissible: that happens after it
/// has passed 1, where the constant branch is the minimum anyway.
///
/// Correctness does not depend on any of this.  Whatever the cuts are, the
/// bound below stays valid, because on every segment it takes the minimum of
/// the constant 1 and of branches that are bounds for the whole segment; a
/// misplaced cut costs tightness only.
fn cut_points(p_lo: f64, branches: &[Branch]) -> ([f64; 4], usize) {
    let mut cuts = [0.0f64; 4];
    let mut n = 1;
    let mut add = |t: f64| {
        if t > 0.0 && t < 1.0 {
            cuts[n] = t;
            n += 1;
        }
    };

    for pair in branches.windows(2) {
        add((pair[0].a_lo - pair[1].a_lo) / (pair[0].b_hi - pair[1].b_hi));
    }
    let last = branches[branches.len() - 1];
    add((1.0 - last.a_lo) / (p_lo - last.b_hi));

    cuts[n] = 1.0;
    n += 1;
    cuts[..n].sort_unstable_by(f64::total_cmp);
    (cuts, n)
}

/// Certified upper bound for int_0^1 F(t) dt on the parameter box described by
/// p_lo, p_hi and the branches.
fn integral_upper(p_lo: f64, p_hi: f64, branches: &[Branch]) -> f64 {
    let (cuts, n) = cut_points(p_lo, branches);
    let mut total = const_interval!(0.0, 0.0);
    for (j, w) in cuts[..n].windows(2).enumerate() {
        let (l, r) = (w[0], w[1]);
        // The constant branch 1 bounds F on every segment.
        let mut best = (pt(r) - pt(l)).sup();
        // The cuts are placed so that segment j is where branch j attains the
        // minimum; the remaining segment belongs to the constant branch.  Only
        // this one branch is evaluated -- the others cannot improve the bound
        // here, and evaluating them would only cost logarithms.
        if let Some(br) = branches.get(j) {
            best = best.min(branch_segment_upper(p_lo, p_hi, *br, l, r));
        }
        total += pt(best);
    }
    let bound = total.sup();
    assert!(
        bound.is_finite() && bound >= 0.0,
        "integral bound is not a number"
    );
    // F <= 1 on [0, 1], so the integral never exceeds 1.
    bound.min(1.0)
}

// ---------------------------------------------------------------------------
// Subdivision of the parameter domain
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Box3 {
    r: Interval,
    m: Interval,
    lambda: Interval,
}

/// Certified upper bound for min{C_eta, C_1/3} on the box.  The eta bound alone
/// is returned as soon as it settles the box, since the minimum of the two is
/// bounded by either one.
fn upper_bound(bx: &Box3, target: f64) -> f64 {
    let (r, m, lambda) = (bx.r, bx.m, bx.lambda);

    // p = R_0(V_single) = 16/25 + 2m + 3lambda.
    let p = rat(16.0, 25.0) + TWO * m + THREE * lambda;
    let (p_lo, p_hi) = (p.inf(), p.sup());

    // C_eta = 5/2 - r + 2r int_0^1 min{1, (1-pt) / (2r - (2r+16/25+2m+4lambda)t)} dt.
    let eta = Branch {
        a_lo: (TWO * r).inf(),
        b_hi: (TWO * r + rat(16.0, 25.0) + TWO * m + FOUR * lambda).sup(),
    };
    let eta_integral = integral_upper(p_lo, p_hi, &[eta]);
    let eta_cost = (rat(5.0, 2.0) - r + TWO * r * interval!(0.0, eta_integral).unwrap()).sup();
    if eta_cost <= target {
        return eta_cost;
    }

    // C_1/3 = 3 - 3r/2 + s int_0^1 min{1, two quotients} dt, with the two
    // denominators displayed after Lemma 4.1, listed in the order in which they
    // attain the minimum: at t = 0 the second denominator is smaller by
    // 9/50 - m - 3lambda/2 >= 0 on D, and it decreases more slowly, by
    // 4/25 + m/2 + lambda per unit of t, so the two swap exactly once.
    let s = rat(3.0, 2.0) * r + rat(4.0, 25.0) + m;
    let branches = [
        Branch {
            a_lo: s.inf(),
            b_hi: (rat(3.0, 2.0) * r + rat(4.0, 5.0) + THREE * m + FOUR * lambda).sup(),
        },
        Branch {
            a_lo: (rat(3.0, 2.0) * r - rat(1.0, 50.0) + TWO * m + rat(3.0, 2.0) * lambda).inf(),
            b_hi: (rat(3.0, 2.0) * r + rat(16.0, 25.0) + rat(5.0, 2.0) * m + THREE * lambda).sup(),
        },
    ];
    let one_third_integral = integral_upper(p_lo, p_hi, &branches);
    let one_third_cost =
        (THREE - rat(3.0, 2.0) * r + s * interval!(0.0, one_third_integral).unwrap()).sup();

    eta_cost.min(one_third_cost)
}

/// Drops boxes that are disjoint from D and shrinks the remaining ones to the
/// constraint m + lambda <= 3/25.
fn restrict(bx: Box3) -> Option<Box3> {
    let cap = rat(3.0, 25.0);

    // Every point of the box violates the constraint if already the smallest
    // sum in it exceeds 3/25.
    if (pt(bx.m.inf()) + pt(bx.lambda.inf())).inf() > cap.sup() {
        return None;
    }

    // Every feasible point of the box has m <= 3/25 - lambda_lo and
    // lambda <= 3/25 - m_lo, so these two contractions keep the whole
    // intersection of the box with D.  A contraction below the lower endpoint
    // means the box is disjoint from D.
    let m = interval!(
        bx.m.inf(),
        bx.m.sup().min((cap - pt(bx.lambda.inf())).sup())
    )
    .ok()?;
    let lambda = interval!(
        bx.lambda.inf(),
        bx.lambda.sup().min((cap - pt(bx.m.inf())).sup())
    )
    .ok()?;
    Some(Box3 { r: bx.r, m, lambda })
}

/// Bisects the coordinate that is widest relative to its extent in D.  The two
/// halves share the midpoint, so their union is the parent box; which
/// coordinate is chosen is a pure heuristic.
fn split(bx: &Box3) -> Option<(Box3, Box3)> {
    // Widths relative to the extents 17/50 in r and 3/25 in m and lambda.
    let widths = [bx.r.wid() / 0.34, bx.m.wid() / 0.12, bx.lambda.wid() / 0.12];
    let axis = if widths[0] >= widths[1] && widths[0] >= widths[2] {
        0
    } else if widths[1] >= widths[2] {
        1
    } else {
        2
    };

    let coordinate = [bx.r, bx.m, bx.lambda][axis];
    let middle = coordinate.mid();
    if !(middle > coordinate.inf() && middle < coordinate.sup()) {
        return None;
    }
    let lower = interval!(coordinate.inf(), middle).ok()?;
    let upper = interval!(middle, coordinate.sup()).ok()?;

    let mut left = *bx;
    let mut right = *bx;
    match axis {
        0 => (left.r, right.r) = (lower, upper),
        1 => (left.m, right.m) = (lower, upper),
        _ => (left.lambda, right.lambda) = (lower, upper),
    }
    Some((left, right))
}

/// Proves that min{C_eta, C_1/3} <= target on the whole intersection of the box
/// with D, and returns the largest certified bound; on failure it returns a box
/// that could not be settled.
fn prove(bx: Box3, depth: u32, target: f64) -> Result<f64, Box3> {
    let Some(bx) = restrict(bx) else {
        return Ok(f64::NEG_INFINITY);
    };

    BOXES.fetch_add(1, Ordering::Relaxed);
    let bound = upper_bound(&bx, target);
    if bound <= target {
        return Ok(bound);
    }
    if depth == 0 {
        return Err(bx);
    }

    let (left, right) = split(&bx).ok_or(bx)?;
    let (left, right) = rayon::join(
        || prove(left, depth - 1, target),
        || prove(right, depth - 1, target),
    );
    Ok(left?.max(right?))
}

fn main() {
    let target = rat(79.0, 25.0);
    // 79/25 is not a binary floating-point number, hence this enclosure is not
    // degenerate and its lower endpoint is strictly below 79/25.  A box
    // accepted at `target` is therefore strictly below the exact target.
    assert!(target.inf() < target.sup());
    let target = target.inf();

    // A superset of D: the constraint m + lambda <= 3/25 is imposed in
    // restrict(), the other bounds are widened outwards.
    let initial = Box3 {
        r: interval!(rat(33.0, 50.0).inf(), 1.0).unwrap(),
        m: interval!(0.0, rat(3.0, 25.0).sup()).unwrap(),
        lambda: interval!(0.0, rat(3.0, 25.0).sup()).unwrap(),
    };

    let start = Instant::now();
    let result = prove(initial, MAX_DEPTH, target);
    let elapsed = start.elapsed().as_secs_f64();
    let boxes = BOXES.load(Ordering::Relaxed);

    match result {
        Ok(bound) => {
            println!("PROVED: max min{{C_eta, C_1/3}} <= {bound:.15} < 79/25");
            println!("boxes: {boxes}, elapsed: {elapsed:.1}s");
        }
        Err(bx) => {
            eprintln!("INCONCLUSIVE: a box could not be certified.");
            eprintln!("  r      = [{}, {}]", bx.r.inf(), bx.r.sup());
            eprintln!("  m      = [{}, {}]", bx.m.inf(), bx.m.sup());
            eprintln!("  lambda = [{}, {}]", bx.lambda.inf(), bx.lambda.sup());
            eprintln!("  bound  = {}", upper_bound(&bx, target));
            eprintln!("boxes: {boxes}, elapsed: {elapsed:.1}s");
            exit(2);
        }
    }
}
