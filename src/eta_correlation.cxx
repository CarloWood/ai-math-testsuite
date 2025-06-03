#include "sys.h"
#include "cairowindow/QuickGraph.h"
#include "math/CubicPolynomial.h"
#include "mpreal/mpreal.h"
#include "utils/square.h"
#include <string>
#include <sstream>
#include <random>
#include <limits>
#include "utils/debug_ostream_operators.h"
#include "debug.h"

using namespace mpfr;

constexpr mp_prec_t precision_in_bits = 128;
constexpr double log10_of_two = 0.30102999566;
constexpr int mpprecision = precision_in_bits * log10_of_two;

// For gdb.
void print(mpreal const& val)
{
  std::cout << val << std::endl;
}

mpreal random_mpreal(mpreal const& min, mpreal const& max, gmp_randstate_t& state)
{
  mpreal result = urandom(state);       // Generates uniform random in [0,1].
  return min + result * (max - min);
}

double eta_to_machine_epsilons(mpreal eta)
{
  return eta.toDouble() / std::numeric_limits<double>::epsilon();
}

class Buckets
{
 private:
  double min_value_;
  double max_value_;
  std::vector<int> buckets_;
  int max_count_;

 public:
  Buckets(int number_of_buckets, double min_value, double max_value) :
    min_value_(min_value), max_value_(max_value), buckets_(number_of_buckets), max_count_(0) { }

  void insert(double value)
  {
    int count = buckets_[value_to_bucket(value)] += 1;
    max_count_ = std::max(max_count_, count);
  }

  // Accessors.
  double min_value() const { return min_value_; }
  double max_value() const { return max_value_; }
  std::vector<int> const& buckets() const { return buckets_; }
  int max_count() const { return max_count_; }
  double width() const { return max_value_ - min_value_; }

 private:
  int value_to_bucket(double value) const
  {
    ASSERT(min_value_ <= value && value <= max_value_);

    // If number_of_buckets == 3,
    //  |...|...|...|
    //  ^           ^
    // min_value_   max_value_

    return std::min(static_cast<int>(std::floor(buckets_.size() * (value - min_value_) / (max_value_ - min_value_))),
        static_cast<int>(buckets_.size()));
  }
};

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with mpprecision digits precision.
  std::streamsize const old_precision = std::cout.precision(mpprecision);
//  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

#if 0
  double d = 16.0;
  for (int i = 0; i < 1000; ++i)
    d = std::nextafter(d, -100.0);
  for (int i = 0; i < 2000; ++i)
  {
    double dp = std::nextafter(d, 100.0);
    std::cout << d << " += " << ((dp - d) / std::numeric_limits<double>::epsilon()) << " ε" << '\n';
    d = dp;
  }

  return 0;
#endif

  // Handle random seed.
  unsigned long int seed;
  if (argc > 1)
  {
    std::istringstream ss(argv[1]);
    ss >> seed;
  }
  else
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "Seed: " << seed << '\n';

  // Use the mersenne twister engine.
  std::mt19937 engine(seed);

  // Seed the GMP random number generator with the same seed.
  gmp_randstate_t randstate;
  gmp_randinit_mt(randstate);
  gmp_randseed_ui(randstate, seed);

  std::array<mpreal, 4> mpreal_coefficients = { -0.196023, 0.959263, -1.67059, 1.03862 };
  math::CubicPolynomial<mpreal> cubic(mpreal_coefficients);
  std::cout << cubic << '\n';

  // Prepare 1000 buckets to store values in the range [-0.5, 0.5];
  constexpr int number_of_buckets = 1000;
  Buckets buckets(number_of_buckets, -2.5, 2.5);

  constexpr double limit = 13263554;     // 0.5 / sqrt(0.1) * 2^23, where 0.25 < sqrt(0.1) < 0.5 so that floor(log2(limit)) = 23.
  cairowindow::Range const random_x_range{0.999 * limit, limit};
  std::ostringstream formula;

  constexpr int loop_size = 10000000;
  for (int n = 0; n < loop_size; ++n)
  {
    mpreal x = random_mpreal(random_x_range.min(), random_x_range.max(), randstate);

    // d1
    mpreal x2 = utils::square(x);
    // ~ (0.5 / sqrt(0.1) * 2^23)^2 =
    //    0.5^2 / 0.1 * 2^46 =
    //    0.5 / (0.1 * 4) * 0.5 * 4 * 2^46
    //    0.5 / 0.4 * 2^47
    //    when floor(log2(limit^2)) = 47 because 0.25 < 0.4 < 0.5.
    double d1 = x2.toDouble();
    // d2
    mpreal c2_x = cubic[2] * x;
    double d2 = c2_x.toDouble();
    // d3
    mpreal c3_x2 = cubic[3] * x2;
    mpreal c3_d1 = cubic[3] * mpreal{d1};
    double d3 = c3_d1.toDouble();
    // d4
    mpreal c2_x_plus_c3_x2 = c2_x + c3_x2;
    mpreal d2_plus_d3 = mpreal{d2} + mpreal{d3};
    double d4 = d2_plus_d3.toDouble();
    // d5
    mpreal c1_plus_c2_x_plus_c3_x2 = cubic[1] + c2_x_plus_c3_x2;
    mpreal c1_plus_d4 = cubic[1] + mpreal{d4};
    double d5 = c1_plus_d4.toDouble();
    // d6
    mpreal c1_x_plus_c2_x2_plus_c3_x3 = c1_plus_c2_x_plus_c3_x2 * x;
    mpreal d5_x = mpreal{d5} * x;
    double d6 = d5_x.toDouble();
    // f_a
    mpreal c0_plus_c1_x_plus_c2_x2_plus_c3_x3 = cubic[0] + c1_x_plus_c2_x2_plus_c3_x3;
    mpreal c0_plus_d6 = cubic[0] + mpreal{d6};
    double f_a = c0_plus_d6.toDouble();
#define F_A 1
#if D1
    mpreal eta = (mpreal{d1} - x2) / abs(x2);  // Hence we expect this to have a PDF uniform between ±0.4.
    if (n == 0)
      formula << "((double)d1 - x^2) / x^2";
#elif D2
    mpreal eta = (mpreal{d2} - c2_x) / abs(c2_x);
    if (n == 0)
      formula << "((double)d2 - c2 x) / |c2 x|";
#elif D3
    mpreal eta = (mpreal{d3} - c3_x2) / abs(c3_x2);
    if (n == 0)
      formula << "((double)d3 - c3 x^2) / |c3 x^2|";
#elif D4
    mpreal eta = (mpreal{d4} - c2_x_plus_c3_x2) / abs(c2_x_plus_c3_x2);
    if (n == 0)
      formula << "((double)d4 - (c2 x + c3 x^2)) / |c2 x + c3 x^2|";
#elif D5
    mpreal eta = (mpreal{d5} - c1_plus_c2_x_plus_c3_x2) / abs(c1_plus_c2_x_plus_c3_x2);
    if (n== 0)
      formula << "((double)d5 - (c1 + c2 x + c3 x^2) / |c1 + c2 x + c3 x^2|";
#elif D6
    mpreal eta = (mpreal{d6} - c1_x_plus_c2_x2_plus_c3_x3) / abs(c1_x_plus_c2_x2_plus_c3_x3);
    if (n == 0)
      formula << "((double)d6 - (c1 x + c2 x^2 + c3 x^3) / |c1 x + c2 x^2 + c3 x^3|";
#elif F_A
    mpreal eta = (mpreal{f_a} - c0_plus_c1_x_plus_c2_x2_plus_c3_x3) / abs(c0_plus_c1_x_plus_c2_x2_plus_c3_x3);
    if (n == 0)
      formula << "((double)f_a - (c0 + c1 x + c2 x^2 + c3 x^3) / |c0 + c1 x + c2 x^2 + c3 x^3|";
#endif

    buckets.insert(eta_to_machine_epsilons(eta.toDouble()));
  }

  std::cout << "Max count: " << buckets.max_count() << '\n';

  using namespace cairowindow;

  std::ostringstream title;
  title << "η Probability Density Function " << formula.str() << " with x in " << random_x_range;
  double const count_scale = static_cast<double>(number_of_buckets) / (loop_size * buckets.width());
  Range const eta_div_epsilon_range{buckets.min_value(), buckets.max_value()};
  Range const pdf_range{0, count_scale * buckets.max_count()};
  QuickGraph qg(title.str(), "η/ε", "P(η | x)", eta_div_epsilon_range, pdf_range);

  draw::LineStyle red_line_style({.line_color = color::red, .line_width = 2.0});

  std::vector<int> const& buckets_vector = buckets.buckets();
  for (int bucket = 0; bucket != buckets_vector.size(); ++bucket)
  {
    int count = buckets_vector[bucket];
    double eta_x = static_cast<double>(bucket) / buckets_vector.size() * buckets.width() + buckets.min_value();
    qg.add_line(cairowindow::LinePiece{{eta_x, 0}, {eta_x, count_scale * count}});
  }

  // Draw the single-η distribution.
  qg.add_function([&](double z){
    return std::abs(z) < 0.25 ? 1.5 : 0.125 / utils::square(z) - 0.5;
  }, red_line_style);

  qg.wait_for_keypress();
}
