// Reproducible interval certificate and feasible-witness search for Lemma 4.1.
// The formulas, variable map, and proof obligations are documented in
// README.md.

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/rounded_transc.hpp>
#include <boost/numeric/interval/transc.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fenv.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace bn = boost::numeric;
namespace bil = boost::numeric::interval_lib;

// rounded_transc_std provides directed rounding for +,-,*,/, and log.  The
// save_state wrapper restores each worker thread's floating-point environment
// after every interval operation.
using Rounding = bil::save_state<bil::rounded_transc_std<double>>;
using Policies = bil::policies<Rounding, bil::checking_strict<double>>;
using Interval = bn::interval<double, Policies>;

namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

// ---------------------------------------------------------------------------
// Command-line handling and exact rational constants
// ---------------------------------------------------------------------------

Interval point(double x) { return Interval(x); }

Interval rational(int numerator, int denominator) {
  return point(static_cast<double>(numerator)) /
         point(static_cast<double>(denominator));
}

double add_up(double a, double b) { return (point(a) + point(b)).upper(); }

std::string format_duration(std::chrono::steady_clock::duration duration) {
  const auto seconds =
      std::chrono::duration_cast<std::chrono::seconds>(duration).count();
  const auto hours = seconds / 3600;
  const auto minutes = (seconds % 3600) / 60;
  const auto remainder = seconds % 60;
  std::ostringstream out;
  out << std::setfill('0');
  if (hours > 0) {
    out << std::setw(2) << hours << ':';
  }
  out << std::setw(2) << minutes << ':' << std::setw(2) << remainder;
  return out.str();
}

struct Options {
  unsigned threads = std::max(1u, std::thread::hardware_concurrency());
  std::size_t batch_per_thread = 16;
  double report_seconds = 2.0;
  std::uint32_t max_box_depth = 90;
  std::uint64_t max_evaluations = 100'000'000;
  int integral_max_depth = 36;
  unsigned grid_r_steps = 24;
  unsigned grid_triangle_steps = 48;
  std::uint64_t search_boxes = 250'000;
  double search_gap = 1e-10;
  unsigned witness_bins = 1'048'576;
  bool run_certificate = true;
  bool run_search = true;
  bool quiet = false;
};

void print_usage(const char *program) {
  std::cout
      << "Usage: " << program << " [options]\n\n"
      << "Certificate and feasible-witness search for Lemma 4.1:\n"
      << "  max min{C_eta(r,m,lambda), C_1/3(r,m,lambda)} < 79/25\n"
      << "on 33/50 <= r <= 1, m >= 0, lambda >= 0,"
         " m+lambda <= 3/25.\n"
      << "The current statement uses the aliases u=r, g=m, delta=lambda.\n\n"
      << "Options:\n"
      << "  --threads N              worker threads (default: hardware count)\n"
      << "  --batch-per-thread N     parents processed per worker and batch\n"
      << "  --report-seconds X       progress-report period (default: 2)\n"
      << "  --max-box-depth N        fail rather than split beyond N\n"
      << "  --max-evaluations N      safety limit on evaluated boxes\n"
      << "  --integral-max-depth N   maximum splits at branch crossings\n"
      << "  --grid-r-steps N         initial grid subdivisions in r\n"
      << "  --grid-triangle-steps N  initial simplex grid subdivisions\n"
      << "  --search-boxes N         branch-and-bound search budget\n"
      << "  --search-gap X           stop search when the rigorous gap <= X\n"
      << "  --witness-bins N         bins for the certified witness value\n"
      << "  --certificate-only       skip the feasible-witness search\n"
      << "  --search-only            skip the rigorous upper certificate\n"
      << "  --quiet                  print only the final result\n"
      << "  --help                   show this message\n";
}

template <class Integer>
Integer parse_integer(const std::string &text, const char *option) {
  std::size_t consumed = 0;
  unsigned long long value = 0;
  try {
    value = std::stoull(text, &consumed);
  } catch (const std::exception &) {
    throw std::runtime_error(std::string("invalid value for ") + option);
  }
  if (consumed != text.size() ||
      value > static_cast<unsigned long long>(
                  std::numeric_limits<Integer>::max())) {
    throw std::runtime_error(std::string("invalid value for ") + option);
  }
  return static_cast<Integer>(value);
}

double parse_real(const std::string &text, const char *option) {
  std::size_t consumed = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &consumed);
  } catch (const std::exception &) {
    throw std::runtime_error(std::string("invalid value for ") + option);
  }
  if (consumed != text.size() || !std::isfinite(value)) {
    throw std::runtime_error(std::string("invalid value for ") + option);
  }
  return value;
}

Options parse_options(int argc, char **argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    auto require_value = [&](const char *name) -> std::string {
      if (++i >= argc) {
        throw std::runtime_error(std::string("missing value for ") + name);
      }
      return argv[i];
    };

    if (argument == "--threads") {
      options.threads =
          parse_integer<unsigned>(require_value("--threads"), "--threads");
    } else if (argument == "--batch-per-thread") {
      options.batch_per_thread = parse_integer<std::size_t>(
          require_value("--batch-per-thread"), "--batch-per-thread");
    } else if (argument == "--report-seconds") {
      options.report_seconds =
          parse_real(require_value("--report-seconds"), "--report-seconds");
    } else if (argument == "--max-box-depth") {
      options.max_box_depth = parse_integer<std::uint32_t>(
          require_value("--max-box-depth"), "--max-box-depth");
    } else if (argument == "--max-evaluations") {
      options.max_evaluations = parse_integer<std::uint64_t>(
          require_value("--max-evaluations"), "--max-evaluations");
    } else if (argument == "--integral-max-depth") {
      options.integral_max_depth = parse_integer<int>(
          require_value("--integral-max-depth"), "--integral-max-depth");
    } else if (argument == "--grid-r-steps") {
      options.grid_r_steps = parse_integer<unsigned>(
          require_value("--grid-r-steps"), "--grid-r-steps");
    } else if (argument == "--grid-triangle-steps") {
      options.grid_triangle_steps = parse_integer<unsigned>(
          require_value("--grid-triangle-steps"), "--grid-triangle-steps");
    } else if (argument == "--search-boxes") {
      options.search_boxes = parse_integer<std::uint64_t>(
          require_value("--search-boxes"), "--search-boxes");
    } else if (argument == "--search-gap") {
      options.search_gap =
          parse_real(require_value("--search-gap"), "--search-gap");
    } else if (argument == "--witness-bins") {
      options.witness_bins = parse_integer<unsigned>(
          require_value("--witness-bins"), "--witness-bins");
    } else if (argument == "--certificate-only") {
      options.run_search = false;
    } else if (argument == "--search-only") {
      options.run_certificate = false;
    } else if (argument == "--quiet") {
      options.quiet = true;
    } else if (argument == "--help") {
      print_usage(argv[0]);
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + argument);
    }
  }

  if (options.threads == 0 || options.batch_per_thread == 0 ||
      options.report_seconds <= 0.0 || options.integral_max_depth < 0 ||
      options.grid_r_steps == 0 || options.grid_triangle_steps == 0 ||
      options.search_boxes == 0 || options.search_gap <= 0.0 ||
      options.witness_bins == 0) {
    throw std::runtime_error("option values must be positive");
  }
  if (!options.run_certificate && !options.run_search) {
    throw std::runtime_error("at least one run mode must be enabled");
  }
  return options;
}

struct Constants {
  Interval three_over_two = rational(3, 2);
  Interval five_over_two = rational(5, 2);
  Interval four_over_twenty_five = rational(4, 25);
  Interval one_over_fifty = rational(1, 50);
  Interval sixteen_over_twenty_five = rational(16, 25);
  Interval four_over_five = rational(4, 5);
  Interval thirty_three_over_fifty = rational(33, 50);
  Interval seventeen_over_fifty = rational(17, 50);
  Interval three_over_twenty_five = rational(3, 25);
  Interval target = rational(79, 25);
};

struct RationalBranch {
  // This branch is the scalar upper function
  //       (1 - p_lo*t) / (a_lo - b_hi*t).
  double a_lo = 0.0;
  double b_hi = 0.0;
};

// ---------------------------------------------------------------------------
// Rigorous upper bounds on parameter boxes
// ---------------------------------------------------------------------------

// Certified upper bound for the analytic integral of
//        (1-p*t)/(a-b*t)
// over [left,right].  All four arguments p,a,b,left,right are treated as exact
// binary floating-point numbers.  Boost interval arithmetic encloses every
// arithmetic operation and the logarithm.
double rational_integral_upper(double p, double a, double b, double left,
                               double right) {
  const Interval P = point(p);
  const Interval A = point(a);
  const Interval B = point(b);
  const Interval L = point(left);
  const Interval R = point(right);
  const Interval width = R - L;

  const Interval denominator_left = A - B * L;
  const Interval denominator_right = A - B * R;
  if (denominator_right.lower() <= 0.0 || b <= 0.0) {
    return width.upper();
  }

  // (1-p*t)/(a-b*t)
  //   = p/b + (b-a*p)/(b*(a-b*t)).
  const Interval value =
      (P / B) * width +
      ((B - A * P) / (B * B)) * bn::log(denominator_left / denominator_right);

  if (!std::isfinite(value.upper())) {
    return width.upper();
  }
  // The constant branch 1 is always available, so an integral bound larger
  // than the segment width is never useful.
  return std::clamp(value.upper(), 0.0, width.upper());
}

class AnalyticEnvelopeIntegrator {
public:
  AnalyticEnvelopeIntegrator(double p_lo, double p_hi,
                             std::vector<RationalBranch> branches,
                             int max_depth)
      : p_lo_(p_lo), p_hi_(p_hi), branches_(std::move(branches)),
        max_depth_(max_depth) {}

  double upper_bound() const { return refine(0.0, 1.0, 0); }

private:
  bool valid_at(const RationalBranch &branch, double t) const {
    // Once either expression becomes nonpositive, it stays nonpositive as
    // t grows.  Strict positivity is required before a quotient is used.
    const Interval T = point(t);
    const Interval numerator_lower = point(1.0) - point(p_hi_) * T;
    const Interval denominator_lower =
        point(branch.a_lo) - point(branch.b_hi) * T;
    return numerator_lower.lower() > 0.0 && denominator_lower.lower() > 0.0;
  }

  double coarse(double left, double right) const {
    double best = (point(right) - point(left)).upper();
    for (const RationalBranch &branch : branches_) {
      if (!valid_at(branch, right)) {
        continue;
      }
      best = std::min(best, rational_integral_upper(p_lo_, branch.a_lo,
                                                    branch.b_hi, left, right));
    }
    return best;
  }

  bool linear_nonnegative(const Interval &constant, const Interval &slope,
                          double left, double right) const {
    const Interval values = constant + slope * Interval(left, right);
    return values.lower() >= 0.0;
  }

  bool linear_nonpositive(const Interval &constant, const Interval &slope,
                          double left, double right) const {
    const Interval values = constant + slope * Interval(left, right);
    return values.upper() <= 0.0;
  }

  double refine(double left, double right, int depth) const {
    // A branch which is valid at left but not at right changes status in
    // this segment.  Split until that single threshold is isolated.
    std::vector<std::size_t> full_segment_branches;
    bool has_validity_change = false;
    for (std::size_t i = 0; i < branches_.size(); ++i) {
      const bool valid_left = valid_at(branches_[i], left);
      const bool valid_right = valid_at(branches_[i], right);
      if (valid_right) {
        full_segment_branches.push_back(i);
      } else if (valid_left) {
        has_validity_change = true;
      }
    }

    if (!has_validity_change) {
      if (full_segment_branches.empty()) {
        // No quotient is a uniform upper bound here.
        return (point(right) - point(left)).upper();
      }

      // The quotients share their positive numerator.  Thus the
      // smallest quotient is the one with the largest denominator
      // d_i(t)=a_i-b_i t.  Look for a denominator which dominates on
      // the complete segment.
      std::optional<std::size_t> active;
      for (const std::size_t i : full_segment_branches) {
        bool dominates = true;
        for (const std::size_t j : full_segment_branches) {
          const Interval constant =
              point(branches_[i].a_lo) - point(branches_[j].a_lo);
          const Interval slope =
              point(branches_[j].b_hi) - point(branches_[i].b_hi);
          if (!linear_nonnegative(constant, slope, left, right)) {
            dominates = false;
            break;
          }
        }
        if (dominates) {
          active = i;
          break;
        }
      }

      if (active.has_value()) {
        const RationalBranch &branch = branches_[*active];

        // d(t)-(1-p_lo t) decides whether the active quotient or the
        // constant branch 1 is the lower upper bound.
        const Interval constant = point(branch.a_lo) - point(1.0);
        const Interval slope = point(p_lo_) - point(branch.b_hi);
        if (linear_nonnegative(constant, slope, left, right)) {
          return rational_integral_upper(p_lo_, branch.a_lo, branch.b_hi, left,
                                         right);
        }
        if (linear_nonpositive(constant, slope, left, right)) {
          return (point(right) - point(left)).upper();
        }
      }
    }

    // The only unresolved possibilities are a validity threshold, a
    // crossing of two denominator lines, or a crossing with 1.  Each is
    // the zero of a linear function, so only the segment containing that
    // zero keeps splitting.
    const double parent = coarse(left, right);
    if (depth >= max_depth_) {
      return parent;
    }
    const double middle = std::midpoint(left, right);
    if (middle == left || middle == right) {
      return parent;
    }
    return std::min(parent, add_up(refine(left, middle, depth + 1),
                                   refine(middle, right, depth + 1)));
  }

  double p_lo_;
  double p_hi_;
  std::vector<RationalBranch> branches_;
  int max_depth_;
};

struct Box {
  double r_lo = 0.0;
  double r_hi = 0.0;
  double m_lo = 0.0;
  double m_hi = 0.0;
  double lambda_lo = 0.0;
  double lambda_hi = 0.0;
  double eta_cost_upper = kInfinity;
  double one_third_cost_upper = kInfinity;
  double upper = kInfinity;
  std::uint32_t depth = 0;
};

struct WorseUpperFirst {
  bool operator()(const Box &left, const Box &right) const {
    return left.upper < right.upper;
  }
};

class BoundEvaluator {
public:
  BoundEvaluator(const Constants &constants, const Options &options)
      : c_(constants), options_(options) {}

  bool tighten_to_domain(Box &box) const {
    // Discard a box only when its lower corner is provably above the exact
    // line m+lambda=3/25 from Lemma 4.1.
    const Interval lower_sum = point(box.m_lo) + point(box.lambda_lo);
    if (lower_sum.lower() > c_.three_over_twenty_five.upper()) {
      return false;
    }

    // These two safe upper contractions retain the complete intersection
    // of the box with the triangular domain.
    const double lambda_cap =
        (c_.three_over_twenty_five - point(box.m_lo)).upper();
    const double m_cap =
        (c_.three_over_twenty_five - point(box.lambda_lo)).upper();
    box.lambda_hi = std::min(box.lambda_hi, lambda_cap);
    box.m_hi = std::min(box.m_hi, m_cap);
    return box.r_lo <= box.r_hi && box.m_lo <= box.m_hi &&
           box.lambda_lo <= box.lambda_hi;
  }

  void evaluate(Box &box) const {
    const Interval R(box.r_lo, box.r_hi);
    const Interval M(box.m_lo, box.m_hi);
    const Interval Lambda(box.lambda_lo, box.lambda_hi);

    // The current statement of Lemma 4.1 calls these variables
    // (u,g,delta); the reduction proof calls them (r,m,lambda).
    const Interval P =
        c_.sixteen_over_twenty_five + point(2.0) * M + point(3.0) * Lambda;

    // sigma_eta branch:
    // (1-P*t) / (2r-(2r+16/25+2m+4lambda)t).
    const Interval eta_a = point(2.0) * R;
    const Interval eta_b = point(2.0) * R + c_.sixteen_over_twenty_five +
                           point(2.0) * M + point(4.0) * Lambda;
    const AnalyticEnvelopeIntegrator eta_integral(
        P.lower(), P.upper(), {{eta_a.lower(), eta_b.upper()}},
        options_.integral_max_depth);
    const double eta_integral_upper = eta_integral.upper_bound();
    const Interval eta_cost =
        c_.five_over_two - R +
        point(2.0) * R * Interval(0.0, eta_integral_upper);
    box.eta_cost_upper = eta_cost.upper();

    // sigma_1/3 branches.  For the first branch, the requested
    // -R_0(V_double)/4 term gives
    //   a_1 = S-(18/25-4m-6lambda)/4
    //       = 3r/2-1/50+2m+3lambda/2.
    // These are exactly the two denominators displayed after Lemma 4.1.
    const Interval S = c_.three_over_two * R + c_.four_over_twenty_five + M;
    const Interval first_a = c_.three_over_two * R - c_.one_over_fifty +
                             point(2.0) * M + c_.three_over_two * Lambda;
    const Interval first_b = c_.three_over_two * R +
                             c_.sixteen_over_twenty_five +
                             c_.five_over_two * M + point(3.0) * Lambda;
    const Interval second_a = S;
    const Interval second_b = c_.three_over_two * R + c_.four_over_five +
                              point(3.0) * M + point(4.0) * Lambda;

    const AnalyticEnvelopeIntegrator one_third_integral(
        P.lower(), P.upper(),
        {{first_a.lower(), first_b.upper()},
         {second_a.lower(), second_b.upper()}},
        options_.integral_max_depth);
    const double one_third_integral_upper = one_third_integral.upper_bound();
    const Interval one_third_cost = point(3.0) - c_.three_over_two * R +
                                    S * Interval(0.0, one_third_integral_upper);
    box.one_third_cost_upper = one_third_cost.upper();
    box.upper = std::min(box.eta_cost_upper, box.one_third_cost_upper);
  }

private:
  const Constants &c_;
  const Options &options_;
};

struct ReducedPoint {
  double r = 0.0;
  double m = 0.0;
  double lambda = 0.0;
};

// ---------------------------------------------------------------------------
// Point evaluation and certified feasible lower bounds
// ---------------------------------------------------------------------------

struct PointCosts {
  long double eta = 0.0L;
  long double one_third = 0.0L;
  long double objective = 0.0L;
};

struct ScalarBranch {
  long double a = 0.0L;
  long double b = 0.0L;
};

void add_breakpoint(std::vector<long double> &points, long double value) {
  if (std::isfinite(value) && value > 0.0L && value < 1.0L) {
    points.push_back(value);
  }
}

long double scalar_rational_integral(long double p, const ScalarBranch &branch,
                                     long double left, long double right) {
  const long double width = right - left;
  const long double coefficient = branch.b - branch.a * p;
  const long double scale =
      std::max({1.0L, std::abs(branch.b), std::abs(branch.a * p)});
  if (std::abs(coefficient) <=
      64.0L * std::numeric_limits<long double>::epsilon() * scale) {
    return std::clamp((p / branch.b) * width, 0.0L, width);
  }

  const long double denominator_left = branch.a - branch.b * left;
  const long double denominator_right = branch.a - branch.b * right;
  if (denominator_left <= 0.0L || denominator_right <= 0.0L) {
    throw std::runtime_error(
        "internal error: active quotient has a nonpositive denominator");
  }
  const long double value = (p / branch.b) * width +
                            (coefficient / (branch.b * branch.b)) *
                                std::log(denominator_left / denominator_right);
  return std::clamp(value, 0.0L, width);
}

// Numerical evaluation at one feasible point.  All changes of the active
// branch occur at zeros of linear functions; after inserting those zeros, the
// envelope can be integrated analytically on every remaining interval.
long double
scalar_envelope_integral(long double p,
                         const std::vector<ScalarBranch> &branches) {
  std::vector<long double> points{0.0L, 1.0L};
  add_breakpoint(points, 1.0L / p);
  for (const ScalarBranch &branch : branches) {
    add_breakpoint(points, branch.a / branch.b);
    if (p != branch.b) {
      // The quotient equals 1 when
      // (a-bt)-(1-pt)=(a-1)+(p-b)t vanishes.
      add_breakpoint(points, (1.0L - branch.a) / (p - branch.b));
    }
  }
  for (std::size_t i = 0; i < branches.size(); ++i) {
    for (std::size_t j = i + 1; j < branches.size(); ++j) {
      const long double slope = branches[i].b - branches[j].b;
      if (slope != 0.0L) {
        add_breakpoint(points, (branches[i].a - branches[j].a) / slope);
      }
    }
  }
  std::sort(points.begin(), points.end());
  points.erase(
      std::unique(points.begin(), points.end(),
                  [](long double left, long double right) {
                    const long double scale =
                        std::max({1.0L, std::abs(left), std::abs(right)});
                    return std::abs(left - right) <=
                           32.0L * std::numeric_limits<long double>::epsilon() *
                               scale;
                  }),
      points.end());

  long double integral = 0.0L;
  for (std::size_t i = 0; i + 1 < points.size(); ++i) {
    const long double left = points[i];
    const long double right = points[i + 1];
    const long double middle = std::midpoint(left, right);
    const long double numerator = 1.0L - p * middle;
    long double value = 1.0L;
    std::optional<std::size_t> active;
    if (numerator > 0.0L) {
      for (std::size_t j = 0; j < branches.size(); ++j) {
        const long double denominator = branches[j].a - branches[j].b * middle;
        if (denominator <= 0.0L) {
          continue;
        }
        const long double quotient = numerator / denominator;
        if (quotient < value) {
          value = quotient;
          active = j;
        }
      }
    }
    integral +=
        active.has_value()
            ? scalar_rational_integral(p, branches[*active], left, right)
            : right - left;
  }
  return integral;
}

PointCosts evaluate_point(const ReducedPoint &point_value) {
  const long double r = point_value.r;
  const long double m = point_value.m;
  const long double lambda = point_value.lambda;
  const long double p = 16.0L / 25.0L + 2.0L * m + 3.0L * lambda;

  const long double eta_integral = scalar_envelope_integral(
      p, {{2.0L * r, 2.0L * r + 16.0L / 25.0L + 2.0L * m + 4.0L * lambda}});
  const long double eta = 5.0L / 2.0L - r + 2.0L * r * eta_integral;

  const long double sigma_one_third = 3.0L * r / 2.0L + 4.0L / 25.0L + m;
  const long double one_third_integral = scalar_envelope_integral(
      p, {{3.0L * r / 2.0L - 1.0L / 50.0L + 2.0L * m + 3.0L * lambda / 2.0L,
           3.0L * r / 2.0L + 16.0L / 25.0L + 5.0L * m / 2.0L + 3.0L * lambda},
          {sigma_one_third,
           3.0L * r / 2.0L + 4.0L / 5.0L + 3.0L * m + 4.0L * lambda}});
  const long double one_third =
      3.0L - 3.0L * r / 2.0L + sigma_one_third * one_third_integral;
  return {eta, one_third, std::min(eta, one_third)};
}

double branch_extension_lower(const Interval &numerator,
                              const Interval &denominator) {
  if (numerator.lower() > 0.0 && denominator.lower() > 0.0) {
    const Interval quotient = numerator / denominator;
    return std::clamp(quotient.lower(), 0.0, 1.0);
  }
  if (numerator.upper() <= 0.0 || denominator.upper() <= 0.0) {
    return 1.0;
  }
  // This integration bin meets a validity threshold. Nonnegativity is the
  // safe lower bound on the quotient capped by the constant branch 1.
  return 0.0;
}

struct IntervalBranch {
  Interval a;
  Interval b;
};

// A certified lower Riemann sum for the envelope at a fixed, feasible binary
// floating-point point.  It is used only once, for the final witness; the grid
// and branch-and-bound search use the much faster analytic evaluator above.
double envelope_integral_lower(const Interval &p,
                               const std::vector<IntervalBranch> &branches,
                               unsigned bins) {
  Interval total = point(0.0);
  for (unsigned i = 0; i < bins; ++i) {
    const double left = static_cast<double>(i) / bins;
    const double right = static_cast<double>(i + 1) / bins;
    const Interval T(left, right);
    const Interval numerator = point(1.0) - p * T;
    double lower = 1.0;
    for (const IntervalBranch &branch : branches) {
      const Interval denominator = branch.a - branch.b * T;
      lower = std::min(lower, branch_extension_lower(numerator, denominator));
    }
    total += (point(right) - point(left)) * point(lower);
  }
  return total.lower();
}

PointCosts certified_point_lower(const ReducedPoint &point_value,
                                 unsigned bins) {
  const Interval R = point(point_value.r);
  const Interval M = point(point_value.m);
  const Interval Lambda = point(point_value.lambda);
  const Interval P = rational(16, 25) + point(2.0) * M + point(3.0) * Lambda;

  const double eta_integral = envelope_integral_lower(
      P,
      {{point(2.0) * R, point(2.0) * R + rational(16, 25) + point(2.0) * M +
                            point(4.0) * Lambda}},
      bins);
  const Interval eta =
      rational(5, 2) - R + point(2.0) * R * point(eta_integral);

  const Interval sigma_one_third = rational(3, 2) * R + rational(4, 25) + M;
  const double one_third_integral = envelope_integral_lower(
      P,
      {{rational(3, 2) * R - rational(1, 50) + point(2.0) * M +
            rational(3, 2) * Lambda,
        rational(3, 2) * R + rational(16, 25) + rational(5, 2) * M +
            point(3.0) * Lambda},
       {sigma_one_third, rational(3, 2) * R + rational(4, 5) + point(3.0) * M +
                             point(4.0) * Lambda}},
      bins);
  const Interval one_third = point(3.0) - rational(3, 2) * R +
                             sigma_one_third * point(one_third_integral);
  return {eta.lower(), one_third.lower(),
          std::min(eta.lower(), one_third.lower())};
}

class ParallelEvaluator {
public:
  ParallelEvaluator(unsigned thread_count, const BoundEvaluator &evaluator)
      : evaluator_(evaluator) {
    workers_.reserve(thread_count);
    for (unsigned i = 0; i < thread_count; ++i) {
      workers_.emplace_back([this] { worker_loop(); });
    }
  }

  ParallelEvaluator(const ParallelEvaluator &) = delete;
  ParallelEvaluator &operator=(const ParallelEvaluator &) = delete;

  ~ParallelEvaluator() {
    {
      std::lock_guard lock(mutex_);
      stop_ = true;
      ++generation_;
    }
    work_ready_.notify_all();
    for (std::thread &worker : workers_) {
      worker.join();
    }
  }

  void evaluate(std::vector<Box> &boxes) {
    if (boxes.empty()) {
      return;
    }
    {
      std::lock_guard lock(mutex_);
      boxes_ = &boxes;
      next_.store(0, std::memory_order_relaxed);
      workers_remaining_ = workers_.size();
      ++generation_;
    }
    work_ready_.notify_all();
    std::unique_lock lock(mutex_);
    work_done_.wait(lock, [this] { return workers_remaining_ == 0; });
    boxes_ = nullptr;
  }

private:
  void worker_loop() {
    std::size_t observed_generation = 0;
    while (true) {
      std::vector<Box> *local_boxes = nullptr;
      {
        std::unique_lock lock(mutex_);
        work_ready_.wait(lock, [this, &observed_generation] {
          return stop_ || generation_ != observed_generation;
        });
        if (stop_) {
          return;
        }
        observed_generation = generation_;
        local_boxes = boxes_;
      }

      while (true) {
        const std::size_t index = next_.fetch_add(1, std::memory_order_relaxed);
        if (index >= local_boxes->size()) {
          break;
        }
        evaluator_.evaluate((*local_boxes)[index]);
      }

      {
        std::lock_guard lock(mutex_);
        if (--workers_remaining_ == 0) {
          work_done_.notify_one();
        }
      }
    }
  }

  const BoundEvaluator &evaluator_;
  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable work_ready_;
  std::condition_variable work_done_;
  std::vector<Box> *boxes_ = nullptr;
  std::atomic<std::size_t> next_{0};
  std::size_t workers_remaining_ = 0;
  std::size_t generation_ = 0;
  bool stop_ = false;
};

// ---------------------------------------------------------------------------
// Parameter subdivision shared by the search and the certificate
// ---------------------------------------------------------------------------

std::optional<std::pair<Box, Box>> split_box(const Box &parent,
                                             const Constants &c) {
  const double r_scale = c.seventeen_over_fifty.upper();
  const double triangle_scale = c.three_over_twenty_five.upper();
  const double r_width = (parent.r_hi - parent.r_lo) / r_scale;
  const double m_width = (parent.m_hi - parent.m_lo) / triangle_scale;
  const double lambda_width =
      (parent.lambda_hi - parent.lambda_lo) / triangle_scale;

  Box left = parent;
  Box right = parent;
  left.depth = right.depth = parent.depth + 1;
  left.upper = right.upper = kInfinity;

  if (r_width >= m_width && r_width >= lambda_width) {
    const double middle = std::midpoint(parent.r_lo, parent.r_hi);
    if (middle == parent.r_lo || middle == parent.r_hi) {
      return std::nullopt;
    }
    left.r_hi = middle;
    right.r_lo = middle;
  } else if (m_width >= lambda_width) {
    const double middle = std::midpoint(parent.m_lo, parent.m_hi);
    if (middle == parent.m_lo || middle == parent.m_hi) {
      return std::nullopt;
    }
    left.m_hi = middle;
    right.m_lo = middle;
  } else {
    const double middle = std::midpoint(parent.lambda_lo, parent.lambda_hi);
    if (middle == parent.lambda_lo || middle == parent.lambda_hi) {
      return std::nullopt;
    }
    left.lambda_hi = middle;
    right.lambda_lo = middle;
  }
  return std::pair<Box, Box>{left, right};
}

Box initial_box(const Constants &constants) {
  Box initial;
  initial.r_lo = constants.thirty_three_over_fifty.lower();
  initial.r_hi = 1.0;
  initial.m_lo = 0.0;
  initial.m_hi = constants.three_over_twenty_five.upper();
  initial.lambda_lo = 0.0;
  initial.lambda_hi = constants.three_over_twenty_five.upper();
  return initial;
}

bool is_feasible_point(const ReducedPoint &candidate,
                       const Constants &constants) {
  return candidate.r >= constants.thirty_three_over_fifty.upper() &&
         candidate.r <= 1.0 && candidate.m >= 0.0 && candidate.lambda >= 0.0 &&
         (point(candidate.m) + point(candidate.lambda)).upper() <=
             constants.three_over_twenty_five.lower();
}

bool make_feasible(ReducedPoint &candidate, const Constants &constants) {
  candidate.r =
      std::clamp(candidate.r, constants.thirty_three_over_fifty.upper(), 1.0);
  candidate.m =
      std::clamp(candidate.m, 0.0, constants.three_over_twenty_five.lower());
  candidate.lambda = std::max(0.0, candidate.lambda);
  const double lambda_cap =
      (constants.three_over_twenty_five - point(candidate.m)).lower();
  candidate.lambda = std::min(candidate.lambda, lambda_cap);
  while ((point(candidate.m) + point(candidate.lambda)).upper() >
         constants.three_over_twenty_five.lower()) {
    if (candidate.lambda > 0.0) {
      candidate.lambda = std::nextafter(candidate.lambda, 0.0);
    } else if (candidate.m > 0.0) {
      candidate.m = std::nextafter(candidate.m, 0.0);
    } else {
      return false;
    }
  }
  return is_feasible_point(candidate, constants);
}

void consider_point(ReducedPoint candidate, const Constants &constants,
                    ReducedPoint &best_point, PointCosts &best_costs,
                    std::uint64_t &evaluated_points) {
  if (!make_feasible(candidate, constants)) {
    return;
  }
  const PointCosts costs = evaluate_point(candidate);
  ++evaluated_points;
  if (costs.objective > best_costs.objective) {
    best_point = candidate;
    best_costs = costs;
  }
}

void sample_box(const Box &box, const Constants &constants,
                ReducedPoint &best_point, PointCosts &best_costs,
                std::uint64_t &evaluated_points) {
  const double cap = constants.three_over_twenty_five.lower();
  const double r = std::midpoint(
      std::max(box.r_lo, constants.thirty_three_over_fifty.upper()), box.r_hi);

  const double m_hi = std::min(box.m_hi, cap - box.lambda_lo);
  if (m_hi >= box.m_lo) {
    const double m = std::midpoint(box.m_lo, m_hi);
    const double lambda_hi = std::min(box.lambda_hi, cap - m);
    if (lambda_hi >= box.lambda_lo) {
      consider_point({r, m, std::midpoint(box.lambda_lo, lambda_hi)}, constants,
                     best_point, best_costs, evaluated_points);
    }
  }

  // Also sample the intersection with m+lambda=3/25, since a maximizer may
  // lie on the slanted boundary of the feasible triangle.
  const double diagonal_m_lo = std::max(box.m_lo, cap - box.lambda_hi);
  const double diagonal_m_hi = std::min(box.m_hi, cap - box.lambda_lo);
  if (diagonal_m_lo <= diagonal_m_hi) {
    const double m = std::midpoint(diagonal_m_lo, diagonal_m_hi);
    consider_point({r, m, cap - m}, constants, best_point, best_costs,
                   evaluated_points);
  }
}

void refine_witness(const Constants &constants, const Options &options,
                    ReducedPoint &best_point, PointCosts &best_costs,
                    std::uint64_t &evaluated_points) {
  double r_step = constants.seventeen_over_fifty.upper() / options.grid_r_steps;
  double triangle_step =
      constants.three_over_twenty_five.upper() / options.grid_triangle_steps;
  for (unsigned iteration = 0; iteration < 160; ++iteration) {
    const ReducedPoint center = best_point;
    const long double old_value = best_costs.objective;
    for (int dr = -1; dr <= 1; ++dr) {
      for (int dm = -1; dm <= 1; ++dm) {
        for (int dlambda = -1; dlambda <= 1; ++dlambda) {
          if (dr == 0 && dm == 0 && dlambda == 0) {
            continue;
          }
          consider_point({center.r + dr * r_step, center.m + dm * triangle_step,
                          center.lambda + dlambda * triangle_step},
                         constants, best_point, best_costs, evaluated_points);
        }
      }
    }
    if (best_costs.objective <= old_value) {
      r_step /= 2.0;
      triangle_step /= 2.0;
      if (std::max(r_step, triangle_step) < 1e-13) {
        return;
      }
    }
  }
}

void refine_diagonal_witness(const Constants &constants, const Options &options,
                             ReducedPoint &best_point, PointCosts &best_costs,
                             std::uint64_t &evaluated_points) {
  const double cap = constants.three_over_twenty_five.lower();
  consider_point({best_point.r, best_point.m, cap - best_point.m}, constants,
                 best_point, best_costs, evaluated_points);
  consider_point({best_point.r, cap - best_point.lambda, best_point.lambda},
                 constants, best_point, best_costs, evaluated_points);

  double r_step = constants.seventeen_over_fifty.upper() / options.grid_r_steps;
  double m_step =
      constants.three_over_twenty_five.upper() / options.grid_triangle_steps;
  for (unsigned iteration = 0; iteration < 160; ++iteration) {
    const ReducedPoint center = best_point;
    const long double old_value = best_costs.objective;
    for (int dr = -1; dr <= 1; ++dr) {
      for (int dm = -1; dm <= 1; ++dm) {
        if (dr == 0 && dm == 0) {
          continue;
        }
        const double m = center.m + dm * m_step;
        consider_point({center.r + dr * r_step, m, cap - m}, constants,
                       best_point, best_costs, evaluated_points);
      }
    }
    if (best_costs.objective <= old_value) {
      r_step /= 2.0;
      m_step /= 2.0;
      if (std::max(r_step, m_step) < 1e-13) {
        return;
      }
    }
  }
}

void balance_diagonal_at_r(double r, const Constants &constants,
                           ReducedPoint &best_point, PointCosts &best_costs,
                           std::uint64_t &evaluated_points) {
  const double cap = constants.three_over_twenty_five.lower();
  constexpr unsigned scan_steps = 32;
  double previous_m = 0.0;
  long double previous_difference = 0.0L;
  bool have_previous = false;

  for (unsigned i = 0; i <= scan_steps; ++i) {
    const double m = cap * i / scan_steps;
    ReducedPoint candidate{r, m, cap - m};
    if (!make_feasible(candidate, constants)) {
      continue;
    }
    const PointCosts costs = evaluate_point(candidate);
    ++evaluated_points;
    if (costs.objective > best_costs.objective) {
      best_point = candidate;
      best_costs = costs;
    }
    const long double difference = costs.eta - costs.one_third;
    if (have_previous &&
        ((previous_difference <= 0.0L && difference >= 0.0L) ||
         (previous_difference >= 0.0L && difference <= 0.0L))) {
      double left = previous_m;
      double right = m;
      long double left_difference = previous_difference;
      for (unsigned iteration = 0; iteration < 64; ++iteration) {
        const double middle = std::midpoint(left, right);
        if (middle == left || middle == right) {
          break;
        }
        ReducedPoint middle_point{r, middle, cap - middle};
        if (!make_feasible(middle_point, constants)) {
          break;
        }
        const PointCosts middle_costs = evaluate_point(middle_point);
        ++evaluated_points;
        if (middle_costs.objective > best_costs.objective) {
          best_point = middle_point;
          best_costs = middle_costs;
        }
        const long double middle_difference =
            middle_costs.eta - middle_costs.one_third;
        if ((left_difference <= 0.0L && middle_difference >= 0.0L) ||
            (left_difference >= 0.0L && middle_difference <= 0.0L)) {
          right = middle;
        } else {
          left = middle;
          left_difference = middle_difference;
        }
      }
    }
    previous_m = m;
    previous_difference = difference;
    have_previous = true;
  }
}

void refine_balanced_diagonal(const Constants &constants,
                              const Options &options, ReducedPoint &best_point,
                              PointCosts &best_costs,
                              std::uint64_t &evaluated_points) {
  double r_step = constants.seventeen_over_fifty.upper() / options.grid_r_steps;
  for (unsigned iteration = 0; iteration < 100; ++iteration) {
    const double center = best_point.r;
    const long double old_value = best_costs.objective;
    balance_diagonal_at_r(center - r_step, constants, best_point, best_costs,
                          evaluated_points);
    balance_diagonal_at_r(center + r_step, constants, best_point, best_costs,
                          evaluated_points);
    if (best_costs.objective <= old_value) {
      r_step /= 2.0;
      if (r_step < 1e-13) {
        return;
      }
    }
  }
}

void print_witness(const ReducedPoint &witness, const PointCosts &numerical,
                   const PointCosts &certified) {
  const long double r = witness.r;
  const long double m = witness.m;
  const long double lambda = witness.lambda;
  std::cout << std::fixed << std::setprecision(15)
            << "Certified feasible lower bound : "
            << static_cast<double>(certified.objective) << '\n'
            << "Numerical value at witness      : "
            << static_cast<double>(numerical.objective)
            << "  (C_eta=" << static_cast<double>(numerical.eta)
            << ", C_1/3=" << static_cast<double>(numerical.one_third) << ")\n"
            << "Reduced variables               : r=" << static_cast<double>(r)
            << ", m=" << static_cast<double>(m)
            << ", lambda=" << static_cast<double>(lambda) << '\n'
            << "Statement aliases               : u=" << static_cast<double>(r)
            << ", g=" << static_cast<double>(m)
            << ", delta=" << static_cast<double>(lambda) << '\n'
            << "Original normalized variables (OPT=1):\n"
            << "  R_1(V_0^eta)       = " << static_cast<double>(1.0L - r)
            << '\n'
            << "  R_1(V_eta^(1/3))   = "
            << static_cast<double>(r - 14.0L / 25.0L + lambda) << '\n'
            << "  R_1(V_single)      = "
            << static_cast<double>(8.0L / 25.0L + m + lambda) << '\n'
            << "  R_1(V_double)      = "
            << static_cast<double>(6.0L / 25.0L - m - 2.0L * lambda) << '\n'
            << "  R_0(V_single)      = "
            << static_cast<double>(16.0L / 25.0L + 2.0L * m + 3.0L * lambda)
            << '\n'
            << "  R_0(V_double)      = "
            << static_cast<double>(18.0L / 25.0L - 4.0L * m - 6.0L * lambda)
            << '\n'
            << "  sigma_eta(V_eta^1) = " << static_cast<double>(2.0L * r)
            << '\n'
            << "  sigma_1/3(V_eta^1) = "
            << static_cast<double>(3.0L * r / 2.0L + 4.0L / 25.0L + m) << '\n';
}

// ---------------------------------------------------------------------------
// Grid and branch-and-bound search for a feasible lower-bound witness
// ---------------------------------------------------------------------------

int run_search(const Options &options) {
  const Constants constants;
  const BoundEvaluator bound_evaluator(constants, options);
  ReducedPoint best_point;
  PointCosts best_costs;
  std::uint64_t evaluated_points = 0;

  const double r_min = constants.thirty_three_over_fifty.upper();
  const double r_width = 1.0 - r_min;
  const double cap = constants.three_over_twenty_five.lower();
  for (unsigned i = 0; i <= options.grid_r_steps; ++i) {
    const double r = i == options.grid_r_steps
                         ? 1.0
                         : r_min + r_width * i / options.grid_r_steps;
    for (unsigned j = 0; j <= options.grid_triangle_steps; ++j) {
      const double m = cap * j / options.grid_triangle_steps;
      for (unsigned k = 0; k + j <= options.grid_triangle_steps; ++k) {
        const double lambda = cap * k / options.grid_triangle_steps;
        consider_point({r, m, lambda}, constants, best_point, best_costs,
                       evaluated_points);
      }
    }
  }

  if (evaluated_points == 0) {
    throw std::runtime_error("internal error: the witness grid is empty");
  }
  const ReducedPoint grid_point = best_point;
  const PointCosts grid_costs = best_costs;
  refine_witness(constants, options, best_point, best_costs, evaluated_points);
  refine_diagonal_witness(constants, options, best_point, best_costs,
                          evaluated_points);
  refine_balanced_diagonal(constants, options, best_point, best_costs,
                           evaluated_points);
  const unsigned pruning_bins = std::min(options.witness_bins, 16'384u);
  const PointCosts pruning_certified =
      certified_point_lower(best_point, pruning_bins);

  Box initial = initial_box(constants);
  if (!bound_evaluator.tighten_to_domain(initial)) {
    throw std::runtime_error("internal error: empty initial domain");
  }
  bound_evaluator.evaluate(initial);
  std::priority_queue<Box, std::vector<Box>, WorseUpperFirst> pending;
  pending.push(initial);
  std::uint64_t evaluated_boxes = 1;
  std::uint64_t pruned_boxes = 0;
  std::uint32_t deepest = 0;

  const auto start = std::chrono::steady_clock::now();
  auto last_report = start;
  while (!pending.empty() && evaluated_boxes + 2 <= options.search_boxes) {
    const Box parent = pending.top();
    if (parent.upper - static_cast<double>(pruning_certified.objective) <=
        options.search_gap) {
      break;
    }
    pending.pop();
    if (parent.depth >= options.max_box_depth) {
      pending.push(parent);
      break;
    }
    const auto split = split_box(parent, constants);
    if (!split) {
      pending.push(parent);
      break;
    }
    for (Box child : {split->first, split->second}) {
      deepest = std::max(deepest, child.depth);
      if (!bound_evaluator.tighten_to_domain(child)) {
        continue;
      }
      bound_evaluator.evaluate(child);
      ++evaluated_boxes;
      sample_box(child, constants, best_point, best_costs, evaluated_points);
      if (child.upper > static_cast<double>(pruning_certified.objective)) {
        pending.push(child);
      } else {
        ++pruned_boxes;
      }
    }

    const auto now = std::chrono::steady_clock::now();
    if (!options.quiet &&
        std::chrono::duration<double>(now - last_report).count() >=
            options.report_seconds) {
      std::cout << '[' << format_duration(now - start)
                << "] search witness=" << std::fixed << std::setprecision(12)
                << static_cast<double>(best_costs.objective) << "  queue-upper="
                << (pending.empty()
                        ? static_cast<double>(pruning_certified.objective)
                        : pending.top().upper)
                << "  pending=" << pending.size()
                << "  boxes=" << evaluated_boxes << '\n';
      last_report = now;
    }
  }

  refine_witness(constants, options, best_point, best_costs, evaluated_points);
  refine_diagonal_witness(constants, options, best_point, best_costs,
                          evaluated_points);
  refine_balanced_diagonal(constants, options, best_point, best_costs,
                           evaluated_points);
  const PointCosts certified =
      certified_point_lower(best_point, options.witness_bins);
  const double remaining_upper =
      pending.empty()
          ? static_cast<double>(pruning_certified.objective)
          : std::max(static_cast<double>(pruning_certified.objective),
                     pending.top().upper);
  std::cout << "\nFEASIBLE-WITNESS SEARCH\n";
  if (!options.quiet) {
    std::cout << std::fixed << std::setprecision(15)
              << "Coarse-grid maximum            : "
              << static_cast<double>(grid_costs.objective)
              << " at (r,m,lambda)=(" << grid_point.r << ',' << grid_point.m
              << ',' << grid_point.lambda << ")\n";
  }
  print_witness(best_point, best_costs, certified);
  if (!options.quiet) {
    std::cout << "Rigorous search upper bound     : " << std::fixed
              << std::setprecision(15) << remaining_upper << '\n'
              << "Search effort                   : " << evaluated_points
              << " points, " << evaluated_boxes << " boxes, " << pruned_boxes
              << " pruned boxes, depth " << deepest << ", "
              << format_duration(std::chrono::steady_clock::now() - start)
              << '\n';
  }
  return 0;
}

// ---------------------------------------------------------------------------
// Rigorous proof of the strict 79/25 upper bound
// ---------------------------------------------------------------------------

void print_box(const Box &box) {
  std::cerr << std::setprecision(17) << "  r      = [" << box.r_lo << ", "
            << box.r_hi << "]\n"
            << "  m      = [" << box.m_lo << ", " << box.m_hi << "]\n"
            << "  lambda = [" << box.lambda_lo << ", " << box.lambda_hi << "]\n"
            << "  C_eta upper = " << box.eta_cost_upper << "\n"
            << "  C_1/3 upper = " << box.one_third_cost_upper << "\n"
            << "  min upper   = " << box.upper << '\n';
}

int run_certificate(const Options &options) {
  const Constants constants;
  const BoundEvaluator bound_evaluator(constants, options);
  ParallelEvaluator parallel(options.threads, bound_evaluator);

  // Compare against a downward-rounded endpoint.  Therefore an accepted box
  // is certainly below the exact rational number 79/25.
  const double target_down = constants.target.lower();

  Box initial = initial_box(constants);
  if (!bound_evaluator.tighten_to_domain(initial)) {
    throw std::runtime_error("internal error: empty initial domain");
  }
  bound_evaluator.evaluate(initial);

  std::priority_queue<Box, std::vector<Box>, WorseUpperFirst> pending;
  double frozen_upper = -kInfinity;
  if (initial.upper <= target_down) {
    frozen_upper = initial.upper;
  } else {
    pending.push(initial);
  }

  std::uint64_t evaluated = 1;
  std::uint64_t certified = initial.upper <= target_down ? 1 : 0;
  std::uint64_t outside = 0;
  std::uint32_t deepest = 0;
  const auto start = std::chrono::steady_clock::now();
  auto last_report = start;

  if (!options.quiet) {
    std::cout << std::fixed << std::setprecision(12)
              << "Target                         : 79/25 = "
              << constants.target.upper() << '\n'
              << "Threads                        : " << options.threads << '\n'
              << "Initial certified upper bound : " << initial.upper
              << "  (eta=" << initial.eta_cost_upper
              << ", one-third=" << initial.one_third_cost_upper << ")\n";
  }

  while (!pending.empty()) {
    const std::size_t parents_to_process = std::min<std::size_t>(
        pending.size(), options.threads * options.batch_per_thread);
    std::vector<Box> children;
    children.reserve(2 * parents_to_process);

    for (std::size_t i = 0; i < parents_to_process; ++i) {
      const Box parent = pending.top();
      pending.pop();
      if (parent.depth >= options.max_box_depth) {
        std::cerr << "INCONCLUSIVE: maximum box depth reached.\n";
        print_box(parent);
        return 2;
      }
      const auto split = split_box(parent, constants);
      if (!split) {
        std::cerr << "INCONCLUSIVE: a box can no longer be split.\n";
        print_box(parent);
        return 2;
      }
      Box left = split->first;
      Box right = split->second;
      deepest = std::max(deepest, left.depth);
      if (bound_evaluator.tighten_to_domain(left)) {
        children.push_back(left);
      } else {
        ++outside;
      }
      if (bound_evaluator.tighten_to_domain(right)) {
        children.push_back(right);
      } else {
        ++outside;
      }
    }

    if (evaluated + children.size() > options.max_evaluations) {
      std::cerr << "INCONCLUSIVE: maximum number of box evaluations "
                   "reached.\n";
      if (!pending.empty()) {
        std::cerr << "Worst unresolved box:\n";
        print_box(pending.top());
      }
      return 2;
    }
    parallel.evaluate(children);
    evaluated += children.size();

    for (const Box &child : children) {
      if (!std::isfinite(child.upper)) {
        std::cerr << "INCONCLUSIVE: non-finite interval bound.\n";
        print_box(child);
        return 2;
      }
      if (child.upper <= target_down) {
        frozen_upper = std::max(frozen_upper, child.upper);
        ++certified;
      } else {
        pending.push(child);
      }
    }

    const auto now = std::chrono::steady_clock::now();
    if (!options.quiet &&
        std::chrono::duration<double>(now - last_report).count() >=
            options.report_seconds) {
      const double unresolved_upper =
          pending.empty() ? -kInfinity : pending.top().upper;
      const double global_upper = std::max(frozen_upper, unresolved_upper);
      const double elapsed_seconds =
          std::chrono::duration<double>(now - start).count();
      const double rate = evaluated / std::max(1e-9, elapsed_seconds);
      std::cout << '[' << format_duration(now - start)
                << "] upper=" << std::setprecision(12) << global_upper
                << "  gap=" << std::showpos
                << (global_upper - constants.target.upper()) << std::noshowpos
                << "  pending=" << pending.size() << "  evaluated=" << evaluated
                << "  certified=" << certified << "  outside=" << outside
                << "  rate=" << std::setprecision(0) << rate
                << " boxes/s  depth=" << deepest << '\n'
                << std::setprecision(12);
      last_report = now;
    }
  }

  const auto finish = std::chrono::steady_clock::now();
  std::cout << std::fixed << std::setprecision(15)
            << "PROVED: max min{C_eta,C_1/3} <= " << frozen_upper
            << " < 79/25.\n";
  if (!options.quiet) {
    std::cout << "Elapsed: " << format_duration(finish - start)
              << ", evaluated boxes: " << evaluated
              << ", certified leaves: " << certified
              << ", outside boxes: " << outside
              << ", maximum depth: " << deepest << '\n';
  }
  return 0;
}

void check_floating_point_environment() {
  // Boost's interval policy relies on the platform exposing IEEE directed
  // rounding through <fenv.h>.
  if (std::fesetround(FE_TONEAREST) != 0 || std::fegetround() != FE_TONEAREST) {
    throw std::runtime_error("the platform does not expose IEEE directed "
                             "rounding through <fenv.h>");
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Options options = parse_options(argc, argv);
    check_floating_point_environment();
    if (options.run_search) {
      const int status = run_search(options);
      if (status != 0) {
        return status;
      }
    }
    if (options.run_certificate) {
      if (options.run_search) {
        std::cout << "\nRIGOROUS UPPER CERTIFICATE\n";
      }
      return run_certificate(options);
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
