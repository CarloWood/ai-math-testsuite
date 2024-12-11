#include "sys.h"
#include "cairowindow/QuickGraph.h"
#include "math/CubicPolynomial.h"
#include "mpreal/mpreal.h"
#include <random>
#include <chrono>
#include <sstream>

using mpreal = mpfr::mpreal;

constexpr int number_of_cubics = 1000000;
constexpr mp_prec_t precision_in_bits = 256;
constexpr double log10_of_two = 0.30102999566;
constexpr int mpprecision = precision_in_bits * log10_of_two;
double const rel_err_threshold = std::pow(10.0, -0.5 * mpprecision);

#ifdef CWDEBUG
// This can be used from within gdb to print an mpreal value.
void print(mpreal const& val)
{
  std::cout << std::setprecision(mpprecision) << val << std::endl;
}
#endif

mpreal exact_root(double guess, math::CubicPolynomial<mpreal> cubic)
{
  DoutEntering(dc::cubic, "exact_root(" << guess << ", " << cubic << ")");
  mpreal u = guess;

  mpreal prev_u = 0;
  mpreal step = 0;
  do
  {
    prev_u = u;
    mpreal Q_u = cubic.evaluate(u);
    Dout(dc::cubic, "Q_u = " << Q_u);
    mpreal dQ_u = cubic.derivative(u);
    mpreal half_ddQ_u = cubic.half_second_derivative(u);
    step = Q_u * dQ_u / (utils::square(dQ_u) - Q_u * half_ddQ_u);
    Dout(dc::cubic, "step = " << step);
    u -= step;
    Dout(dc::cubic, "u = " << u);
  }
  while (abs((prev_u - u) / u) > rel_err_threshold);

  return u;
}

double halley_iteration(double u, math::CubicPolynomial<double> const& cubic)
{
  double Q_u = cubic.evaluate(u);
  double dQ_u = cubic.derivative(u);
  double half_ddQ_u = cubic.half_second_derivative(u);
  double step = Q_u * dQ_u / (utils::square(dQ_u) - Q_u * half_ddQ_u);
  u -= step;
  return u;
}

class Statistics
{
 private:
  double minimum_;
  double maximum_;
  double sum_{0};
  size_t count_{0};

 public:
  void add_value(double const value)
  {
    if (count_ == 0)
      minimum_ = maximum_ = value;
    else
    {
      minimum_ = std::min(minimum_, value);
      maximum_ = std::max(maximum_, value);
    }
    sum_ += value;
    ++count_;
  }

  // Accessors.
  double get_minimum() const { return minimum_; }
  double get_maximum() const { return maximum_; }
  double get_average() const { return sum_ / count_; }

  friend std::ostream& operator<<(std::ostream& os, Statistics const& stats)
  {
    return os << std::setprecision(6) << stats.minimum_ << '/' << stats.get_average() << '/' << stats.maximum_;
  }
};

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Handle random seed.
  unsigned int seed;
  if (argc > 1)
  {
    std::istringstream ss(argv[1]);
    ss >> seed;
  }
  else
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  Dout(dc::notice, "Seed: " << seed);

  // Use the mersenne twister engine.
  std::mt19937 engine(seed);

  // Generate random cubics.
  Dout(dc::notice|flush_cf, "Generating random cubic polynomials...");
  std::vector<math::CubicPolynomial<double>> cubics;
  std::uniform_real_distribution<double> coefficient_dist(-10.0, 10.0);
  for (int i = 0; i < number_of_cubics; ++i)
    cubics.emplace_back(coefficient_dist(engine), coefficient_dist(engine), coefficient_dist(engine), coefficient_dist(engine));
  Dout(dc::notice, "Done (" << cubics.size() << " cubics).");

  Statistics stats;
  std::array<int, 1001> buckets = {};

  int ci = 0;
  for (auto&& cubic : cubics)
  {
    Dout(dc::notice|flush_cf, std::setprecision(std::numeric_limits<double>::max_digits10) << ci << " : P(x) = " << cubic);

    // Calculate the, or one of, roots to very high precision.
    std::array<double, 3> roots;
    int n = cubic.get_roots(roots);

    // If there is more than one root pick the one with the smallest absolute derivative.
    if (n > 1)
    {
      int min_derivative_i = 0;
      double min_derivative = cubic.derivative(roots[0]);
      for (int i = 1; i < n; ++i)
      {
        double derivative = cubic.derivative(roots[i]);
        if (std::abs(derivative) < std::abs(min_derivative))
        {
          min_derivative_i = i;
          min_derivative = derivative;
        }
      }
      n = min_derivative_i;
    }
    else
      n = 0;

    math::CubicPolynomial<mpreal> mpcubic(cubic[0], cubic[1], cubic[2], cubic[3]);
    mpreal root = exact_root(roots[n], mpcubic);

    Dout(dc::notice, "P(" << root << ") = " << mpcubic(root));

    // Use a starting point with an absolute relative error between 0.001 and 0.002,
    // and calculate the result of using Halley's method by iterating one iterating
    // further than an absolute relative error of 1e-8.

    std::uniform_real_distribution<double> relative_error_dist(0.0, 0.002);
    double relative_error = relative_error_dist(engine);
    if (relative_error < 0.001)
      relative_error -= 0.002;
    // relative_error is now in [-0.002, -0.001) or [0.001, 0.002].
    // relative_error = (starting_point - root) / |root| -->
    // starting_point = root + relative_error * |root|.
    double x = root.toDouble() + relative_error * std::abs(root.toDouble());
    constexpr double relative_error_threshold = 1e-8;
    do
    {
      x = halley_iteration(x, cubic);
      relative_error = (x - root.toDouble()) / std::abs(root.toDouble());
    }
    while (relative_error >= relative_error_threshold);
    // Do one more after reaching a relative error of 1e-8.
    x = halley_iteration(x, cubic);

    // Find the probability distribution of the final relative error.
    relative_error = (x - root.toDouble()) / std::abs(root.toDouble());
    Dout(dc::notice, "relative_error = " << relative_error);
    stats.add_value(relative_error);

    // Plot the probability distribution.
    // Ignore outliers.
    if (-1e-14 < relative_error && relative_error < 1e-14)
    {
      // Use 1000 buckets.
      int bucket = (relative_error + 1e-14) * (1000.0 / 2e-14);
      ++buckets[bucket];
    }

    ++ci;
  }

  cairowindow::QuickGraph graph("Probability distribution of rel.err.", "bucket", "count", {0, buckets.size() - 1}, {0, 35000});

  for (int b = 0; b <= 1000; ++b)
  {
    std::cout << b << " : " << buckets[b] << '\n';
    graph.add_line(cairowindow::LinePiece{{static_cast<double>(b), 0}, {static_cast<double>(b), static_cast<double>(buckets[b])}});
  }

  graph.wait_for_keypress();

  Dout(dc::notice, "relative_error (min/avg/max) = " << stats);
}
