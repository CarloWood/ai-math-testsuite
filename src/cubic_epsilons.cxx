#include "sys.h"
#include <iostream>
#include "utils/square.h"
#include "utils/popcount.h"
#include "math/subsuper_string.h"
#include "mpreal/mpreal.h"
#include <cstdint>
#include <array>
#include <limits>
#include <utility>
#include "debug.h"

using mpreal = mpfr::mpreal;

constexpr mp_prec_t precision_in_bits = 256;
constexpr double log10_of_two = 0.30102999566;
constexpr int mpprecision = precision_in_bits * log10_of_two;

constexpr int number_of_coefficients = 4;
using Coefficients = std::array<mpreal, number_of_coefficients>;
using Epsilons = std::array<mpreal, 11>;

#ifdef CWDEBUG
void print(mpreal val)
{
  std::cout << val << std::endl;
}
#endif

struct CoefficientsSet
{
  using Index = int;
  using Mask = uint32_t;

  static constexpr Mask index2mask(Index i1) { return static_cast<Mask>(1) << i1; }

  Mask mask_;

  // Default construct an empty set.
  CoefficientsSet() : mask_(0) { }

  CoefficientsSet(Index i) : mask_(index2mask(i)) { }
  CoefficientsSet(Index i0, Index i1) : mask_(index2mask(i0)|index2mask(i1)) { }

  int count() const
  {
    return utils::popcount(mask_);
  }

  bool empty() const
  {
    return mask_ == 0;
  }
};

class VariablesSet : public CoefficientsSet
{
  Index epsilon_index_;                 // 0 ... Epsilons::size()-1

  // Variabbles used during generation of the permutations.
  mutable Coefficients coefficients_;
  mutable Epsilons epsilons_;
  mutable int empty_sign_;              // The sign to be used when current_mask_ has no bits set.

  mutable uint32_t unused_mask_;
  mutable uint32_t current_mask_{0};
  static constexpr uint32_t epsilon_mask_{1 << number_of_coefficients};

 public:
  static mpreal const delta;

 public:
  VariablesSet(CoefficientsSet const& coefficients_set, Index epsilon_index = -1) :
    CoefficientsSet(coefficients_set), epsilon_index_(epsilon_index) { }

  int size() const
  {
    return CoefficientsSet::count() + (epsilon_index_ == -1 ? 0 : 1);
  }

  bool empty() const
  {
    return epsilon_index_ == -1 && CoefficientsSet::empty() == 0;
  }

  void begin_permutations(Coefficients const& c, Epsilons const& e) const;
  bool next_permutation(Coefficients const& c, Epsilons const& e) const;

  int sign() const
  {
    return utils::popcount(current_mask_) % 2 == 1 ? -empty_sign_ : empty_sign_;
  }

  Coefficients const& c() const
  {
    return coefficients_;
  }

  Epsilons const& e() const
  {
    return epsilons_;
  }
};

//static
mpreal const VariablesSet::delta("1e-24");

enum What {
  Cubic,
  Derivative,
  Newton
};

mpreal f(mpreal const& x, Coefficients const& c, Epsilons const& e, What what)
{
  DoutEntering(dc::notice, "f(" << std::setprecision(mpprecision) << x << ", " << c << ", " << e << ", what)");
  if (what == Cubic)
  {
    // c₀ + (c₁ + (c₂ + c₃x) x) x
    mpreal exp4 = (1 + e[5]) * (c[3] * x);
    mpreal exp3 = (1 + e[4]) * (c[2] + exp4);
    mpreal exp2 = (1 + e[3]) * (exp3 * x);
    mpreal exp1 = (1 + e[2]) * (c[1] + exp2);
    mpreal exp0 = (1 + e[1]) * (exp1 * x);
    return (1 + e[0]) * (c[0] + exp0);
  }
  else if (what == Derivative)
  {
    // c₁ + (2c₂ + 3c₃x) x
    mpreal exp2 = (1 + e[9]) * (3 * c[3] * x);
    mpreal exp1 = (1 + e[8]) * (2 * c[2] + exp2);
    mpreal exp0 = (1 + e[7]) * (exp1 * x);
    return (1 + e[6]) * (c[1] + exp0);
  }

  // x - f(x)/f'(x)
  mpreal exp1 = (1 + e[10]) * (f(x, c, e, Cubic) / f(x, c, e, Derivative));
  return x - exp1;
}

void VariablesSet::begin_permutations(Coefficients const& coefficients, Epsilons const& epsilons) const
{
  coefficients_ = coefficients;
  epsilons_ = epsilons;
  empty_sign_ = size() % 2 == 1 ? -1 : 1;

  unused_mask_ = ~mask_;
  if (epsilon_index_ != -1)
    unused_mask_ &= ~epsilon_mask_;
}

bool VariablesSet::next_permutation(Coefficients const& c, Epsilons const& e) const
{
  for (int l = 0;; ++l)
  {
    if ((current_mask_ & epsilon_mask_) != 0)
    {
      if (l == 0)
        epsilons_[epsilon_index_] = e[epsilon_index_];
      else
        epsilons_[epsilon_index_] += delta;
    }
    for (int i = 0; i < c.size(); ++i)
    {
      if ((current_mask_ & index2mask(i)) != 0)
      {
        if (l == 0)
          coefficients_[i] = c[i];
        else
          coefficients_[i] += delta;
      }
    }
    if (l == 1)
      break;
    current_mask_ |= unused_mask_;
    ++current_mask_;
    current_mask_ &= ~unused_mask_;
  }
  return current_mask_;
}

mpreal derivative_of_f(mpreal const& x, Coefficients const& c, Epsilons const& e, What what, VariablesSet const& vs)
{
  if (vs.empty())
    return f(x, c, e, what);

  mpreal result = 0;
  vs.begin_permutations(c, e);
  do
  {
    result += vs.sign() * f(x, vs.c(), vs.e(), what);
  }
  while (vs.next_permutation(c, e));

  return result / pow(VariablesSet::delta, vs.size());
}

mpreal random_between(mpreal const& a, mpreal const& b, gmp_randstate_t& state)
{
  mpreal r;
  mpfr_urandomb(r.mpfr_ptr(), state);
  return a + (b - a) * r;
}

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
  mpreal r1("1.5");// = random_between(1.0, 2.0, state);
  mpreal r2("-0.5");// = random_between(-1.0, 0.0, state);
  mpreal c3 = random_between(-0.5, 2.0, state);

  Coefficients c = { -c3 * r0 * r1 * r2, c3 * (r0 * r1 + r0 * r2 + r1 * r2), -c3 * (r0 + r1 + r2), c3 };

  // Initialize all epsilons to be zero.
  Epsilons e = {};
  // Then this should be zero.
  std::cout << "f(0.1) = " << f(r0, c, e, Cubic) << std::endl;

  // Define a xi.
  mpreal xi("1e-8");
  mpreal x_n = r0 + xi;

  // Calculate the derivative at x_n.
  mpreal dfx = c[1] + (2 * c[2] + 3 * c[3] * x_n) * x_n;

  // Define all epsilons.
  for (int j = 0; j < e.size(); ++j)
    e[j] = 0; //std::numeric_limits<double>::epsilon();

  What what = Derivative;

  // If what == Cubic or Derivative, we want to automate finding each coefficient in the general formula:
  //
  //   fₐ(x_n) = f(x_n) + Σ (Σ (Σ αᵢⱼₖ⋅εⱼ)⋅cᵢ)⋅x_n^k
  //                      k  i  j
  std::array<std::array<std::array<int, c.size()>, e.size()>, c.size()> alpha;        // alpha[i][j][k]

  // The code below finds that the Approximations (that is, the calculated values including floating point round off errors) are
  //
  // Cubic:
  //   fₐ(x_n) = f(x_n) + Ε₀(x_n)
  // Derivative:
  //   f'ₐ(x_n) = f'(x_n) + E₁(x_n)
  //
  // where
  //   Ε₀(x_n) = ε₀⋅c₀ + (ε₀ + ε₁ + ε₂)⋅  c₁⋅x_n + (ε₀ + ε₁ + ε₂ + ε₃ + ε₄)⋅  c₂⋅x_n² + (ε₀ + ε₁ + ε₂ + ε₃ + ε₄ + ε₅)⋅c₃⋅x_n³
  //   E₁(x_n) = ε₆⋅c₁ + (ε₆ + ε₇ + ε₈)⋅2⋅c₂⋅x_n +      (ε₆ + ε₇ + ε₈ + ε₉)⋅3⋅c₃⋅x_n²
  //
  // If what == Newton we're calculating
  //
  //   x_{n+1} = x_n - f(x_n) / f'(x_n)
  //
  // Once x_n is close enough to a root, the result becomes much less accurate than it could have been:
  // the terms with ε's become dominant. Therefore we might as well stop after a single iteration and investigate
  // the difference g(x_n) = x_{n+1} - x_n.
  //
  //   gₐ(x_n) = (1 + ε₁₀) -fₐ(x_n) / f'ₐ(x_n) = (1 + ε₁₀) -(f(x_n) + Ε₀(x_n)) / (f'(x_n) + Ε₁(x_n)) =
  //
  // Multiply with (f'(x_n) - Ε₁(x_n)) / (f'(x_n) - Ε₁(x_n)):
  //
  //           = (1 + ε₁₀) -(f(x_n) + Ε₀(x_n)) (f'(x_n) - Ε₁(x_n)) / (f'(x_n)² - Ε₁(x_n)²)
  //
  // Neglect any higher order epsilon factors:
  //
  //   gₐ(x_n) = -(1 + ε₁₀) f(x_n) / f'(x_n) + (f(x_n) Ε₁(x_n) - f'(x_n) Ε₀(x_n)) / f'(x_n)²
  //

  for (int i = 0; i < c.size(); ++i)
  {
    for (int j = 0; j < e.size(); ++j)
    {
      // Make sure we don't leave anything in alpha undefined.
      for (int k = 0; k < c.size(); ++k)
        alpha[i][j][k] = 0;

      // Let calculate
      // ∂(∂f/∂cᵢ)/∂εⱼ ≈ [(f(cᵢ + Δ, εⱼ + Δ, ...) - f(cᵢ, εⱼ + Δ, ...))/Δ - (f(cᵢ + Δ, εⱼ, ...) - f(cᵢ, εⱼ, ...))/Δ] / Δ
      //               = (f(cᵢ + Δ, εⱼ + Δ, ...) - f(cᵢ, εⱼ + Δ, ...) - f(cᵢ + Δ, εⱼ, ...) + f(cᵢ, εⱼ, ...)) / Δ²

#if 0
      static mpreal const delta("1e-24");
      mpreal fa_xn = f(x_n, c, e, what);

      e[j] += delta;
      mpreal fa_xn_dej = f(x_n, c, e, what);

      c[i] += delta;
      mpreal fa_xn_dci_dej = f(x_n, c, e, what);

      e[j] = 0; //std::numeric_limits<double>::epsilon();
      mpreal fa_xn_dci = f(x_n, c, e, what);

      c[i] -= delta;

      mpreal ddfdcidej = (fa_xn_dci_dej - fa_xn_dej - fa_xn_dci + fa_xn) / utils::square(delta);
#else
      // Calculate the derivative ∂²f/∂cᵢ∂εⱼ.
      mpreal ddfdcidej = derivative_of_f(x_n, c, e, what, {CoefficientsSet{i}, j});
#endif
      double dd = ddfdcidej.toDouble();

      std::cout << "i = " << i << ", j = " << j << ", ∂(∂f/∂cᵢ)/∂εⱼ = " << (dd * dfx * dfx/ (-c[2] - 2 * c[3] * x_n)) << std::endl;

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
            if (a != 1)
              oss << a << "⋅";
          }
          else
          {
            if (std::abs(a) == 1)
            {
              if (a == -1)
                oss << "-";
            }
            else
              oss << a << "⋅";
          }
          oss << "ε" << math::to_subscript(j);
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
      std::cout << "⋅c" << math::to_subscript(i);
      if (k > 0)
      {
        std::cout << "⋅x_n";
        if (k > 1)
          std::cout << math::to_superscript(k);
      }
    }
  }
  std::cout << std::endl;
}
