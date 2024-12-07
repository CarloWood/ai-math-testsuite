#include "sys.h"
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
  std::vector<math::CubicPolynomial<mpreal>> cubics;
  std::uniform_real_distribution<double> dist(-10.0, 10.0);
  for (int i = 0; i < number_of_cubics; ++i)
    cubics.emplace_back(dist(engine), dist(engine), dist(engine), dist(engine));
  Dout(dc::notice, "Done (" << cubics.size() << " cubics).");

  for (auto&& cubic : cubics)
    Dout(dc::notice, std::setprecision(mpprecision) << cubic);
}
