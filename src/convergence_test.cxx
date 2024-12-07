#include "sys.h"
#include "math/CubicPolynomial.h"
#include "mpreal/mpreal.h"

using mpreal = mpfr::mpreal;

constexpr mp_prec_t precision_in_bits = 256;
constexpr double log10_of_two = 0.30102999566;
constexpr int mpprecision = precision_in_bits * log10_of_two;

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  math::CubicPolynomial<mpreal> P(1.1, 2.2, 3.1415, -0.01);

  Dout(dc::notice, "P(0.33) = " << std::setprecision(mpprecision) << P(0.33));
}
