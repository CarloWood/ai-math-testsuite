#include "sys.h"
#include <iostream>
#include "utils/square.h"
#include "mpreal/mpreal.h"
#include <array>
#include <limits>
#include "debug.h"

using mpreal = mpfr::mpreal;

constexpr mp_prec_t precision_in_bits = 256;
constexpr double log10_of_two = 0.30102999566;
constexpr int mpprecision = precision_in_bits * log10_of_two;

mpreal evaluate_cubic(mpreal const& x, std::array<mpreal, 4> const& c, std::array<mpreal, 10> const& e, bool derivative = false)
{
  if (derivative)
  {
    // c₁ + (2c₂ + 3c₃x) x
    mpreal exp2 = (1 + e[9]) * (3 * c[3] * x);
    mpreal exp1 = (1 + e[8]) * (2 * c[2] + exp2);
    mpreal exp0 = (1 + e[7]) * (exp1 * x);
    return (1 + e[6]) * (c[1] + exp0);
  }

  // c₀ + (c₁ + (c₂ + c₃x) x) x
  mpreal exp4 = (1 + e[5]) * (c[3] * x);
  mpreal exp3 = (1 + e[4]) * (c[2] + exp4);
  mpreal exp2 = (1 + e[3]) * (exp3 * x);
  mpreal exp1 = (1 + e[2]) * (c[1] + exp2);
  mpreal exp0 = (1 + e[1]) * (exp1 * x);
  return (1 + e[0]) * (c[0] + exp0);
}

mpreal random_between(mpreal const& a, mpreal const& b, gmp_randstate_t& state)
{
  mpreal r;
  mpfr_urandomb(r.mpfr_ptr(), state);
  return a + (b - a) * r;
}

#include <cmath>
#include <utility>
#include <limits>
#include <cstdint>

std::pair<int, int> decompose_float(double value)
{
  // Handle zero case.
  if (std::abs(value) < 0.00001)
    return {0, 0};

  // Handle negative numbers.
  int sign = (value < 0) ? -1 : 1;
  value = std::abs(value);

  // Find k by counting decimal places.
  int k = 0;
  double scaled = value;

  // Scale up until we get close to an integer.
  while (scaled < 0.9 && k < 5)
  {
    scaled *= 10;
    ++k;
  }

  // Round to nearest integer to handle floating point imprecision.
  int n = static_cast<int>(std::round(scaled));

  // Apply the sign to n.
  n *= sign;

  return {n, k};
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Initialize random state.
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, time(NULL));

  // Print everything with mpprecision digits precision.
  std::streamsize const old_precision = std::cout.precision(mpprecision);
  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  // Define a recognizable root.
  mpreal r0("0.1");

  // Define the coefficients of a cubic, such that r is a root.
  // First define two other (real) root (I don't think it matter that they are picked to be real).
  mpreal r1 = random_between(1.0, 2.0, state);
  mpreal r2 = random_between(-1.0, 0.0, state);
  mpreal c3 = random_between(-0.5, 2.0, state);

  //                          -c₃⋅r₀⋅r₁⋅r₂     + c₃⋅(r₀⋅r₁ + r₀⋅r₂ + r₁⋅r₂)⋅x     - c₃⋅(r₀ + r₁ + r₂)⋅x² + c₃⋅x³
  std::array<mpreal, 4> c = { -c3 * r0 * r1 * r2, c3 * (r0 * r1 + r0 * r2 + r1 * r2), -c3 * (r0 + r1 + r2),    c3 };

  // Initialize all epsilons to be zero.
  std::array<mpreal, 10> e = {};
  // Then this should be zero.
  std::cout << "f(0.1) = " << evaluate_cubic(r0, c, e) << std::endl;

  // Define a xi.
  mpreal xi("1e-8");
  mpreal x_n = r0 + xi;

  // Define all epsilons.
  for (int j = 0; j < e.size(); ++j)
    e[j] = std::numeric_limits<double>::epsilon();

  // We want to automate finding each coefficient in the general formula:
  //
  //   fₐ(x_n) = f(x_n) + Σ (Σ (Σ αᵢⱼₖ⋅εⱼ)⋅cᵢ)⋅x_n^k
  //                      k  i  j

  std::array<std::array<std::array<int, c.size()>, e.size()>, c.size()> alpha;        // alpha[i][j][k]

  mpreal delta("1e-24");
  constexpr bool do_derivative = false;
  for (int i = 0; i < c.size(); ++i)
  {
    for (int j = 0; j < e.size(); ++j)
    {
      // Let calculate
      // ∂(∂f/∂cᵢ)/∂εⱼ ≈ [(f(cᵢ + Δ, εⱼ + Δ, ...) - f(cᵢ, εⱼ + Δ, ...))/Δ - (f(cᵢ + Δ, εⱼ, ...) - f(cᵢ, εⱼ, ...))/Δ] / Δ
      //               = (f(cᵢ + Δ, εⱼ + Δ, ...) - f(cᵢ, εⱼ + Δ, ...) - f(cᵢ + Δ, εⱼ, ...) + f(cᵢ, εⱼ, ...)) / Δ²

      mpreal fa_xn = evaluate_cubic(x_n, c, e, do_derivative);

      e[j] += delta;
      mpreal fa_xn_dej = evaluate_cubic(x_n, c, e, do_derivative);

      c[i] += delta;
      mpreal fa_xn_dci_dej = evaluate_cubic(x_n, c, e, do_derivative);

      e[j] = std::numeric_limits<double>::epsilon();
      mpreal fa_xn_dci = evaluate_cubic(x_n, c, e, do_derivative);

      c[i] -= delta;

      // Calculate the derivative ∂(∂f/∂cᵢ)/∂εⱼ.
      mpreal ddfdcidej = (fa_xn_dci_dej - fa_xn_dej - fa_xn_dci + fa_xn) / utils::square(delta);
      double dd = ddfdcidej.toDouble();

//      std::cout << "i = " << i << ", j = " << j << ", ∂(∂f/∂cᵢ)/∂εⱼ = " << dd << std::endl;

      // Make sure we don't leave anything in alpha undefined.
      for (int k = 0; k < c.size(); ++k)
        alpha[i][j][k] = 0;

      // Extract alpha.
      // This partial derivative now has the form αᵢⱼₖ⋅x_n^k, where x_n ≈ 0.1 and αᵢⱼₖ is expected to be a small integer (less than 9).
      // For example:
      //   dd == -0.02  --> k = 2, αᵢⱼ₂ = -2  because -2 * 0.1^2 = 0.02
      //   dd ==  0.003 --> k = 3, αᵢⱼ₂ = 3   because  3 * 0.1^3 = 0.003
      //
      auto [a, k] = decompose_float(dd);

      ASSERT(0 <= k && k < c.size());

      if (a != 0)
        alpha[i][j][k] = a;
    }
  }

  //   fₐ(x_n) = f(x_n) + Σ (Σ (Σ αᵢⱼₖ⋅εⱼ)⋅cᵢ)⋅x_n^k
  //                      k  i  j
  std::cout << "fₐ(x_n) = f(x_n)";
  for (int k = 0; k < c.size(); ++k)
  {
    for (int i = 0; i < c.size(); ++i)
    {
      bool need_parens = false;
      std::ostringstream oss;
      bool saw_epsilon = false;
      for (int j = 0; j < e.size(); ++j)
      {
        if (alpha[i][j][k])
        {
          int a = alpha[i][j][k];
          if (saw_epsilon)
          {
            if (a > 0)
              oss << " + ";
            else
            {
              oss << " - ";
              a = -a;
            }
            need_parens = true;
            if (a == 1)
              oss << "e" << j;
            else
              oss << a << "⋅e" << j;
          }
          else
          {
            if (std::abs(a) == 1)
            {
              if (a == 1)
                oss << "e" << j;
              else
                oss << "-e" << j;
            }
            else
              oss << a << "⋅e" << j;
          }
          saw_epsilon = true;
        }
      }
      if (!saw_epsilon)
        continue;
      std::cout << " + ";
      if (need_parens)
        std::cout << "(";
      std::cout << oss.str();
      if (need_parens)
        std::cout << ")";
      std::cout << "⋅c" << i;
      if (k > 0)
      {
        std::cout << "⋅x_n";
        if (k > 1)
          std::cout << "^" << k;
      }
    }
  }
  std::cout << std::endl;

  // e0⋅(c0 + c1⋅x_n + c2⋅x_n^2 + c3⋅x_n^3) + (e1 + e2)⋅(c1⋅x_n + c2⋅x_n^2 + c3⋅x_n^3) + (e3 + e4)⋅(c2⋅x_n^2 + c3⋅x_n^3) + e5⋅c3⋅x_n^3
}
