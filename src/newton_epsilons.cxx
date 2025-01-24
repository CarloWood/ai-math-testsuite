#include "sys.h"
#include "utils/square.h"
#include "utils/macros.h"
#include "utils/almost_equal.h"
#include "utils/EnumIterator.h"
#include <utility>
#include "debug.h"

#define FOREACH_VARIABLE(X) \
  X(c0) \
  X(c1) \
  X(c2) \
  X(c3) \
  X(x)

#include "FloatingPointRoundOffError.h"

// For gdb.
void print(symbolic::Expression const& expression)
{
  std::cout << expression << std::endl;
}

void print(floating_point_round_off_error::Expression const& expression)
{
  std::cout << expression << std::endl;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using namespace floating_point_round_off_error;
  using FA = floating_point_round_off_error::Expression;

  ExpressionManager expression_manager;

  FOREACH_VARIABLE(DECLARE_VARIABLE);   // Declare all variables (c0, c1, c2, c3, x).

  // Construct the approximation cubic f_a(x_n).
  auto fa = c0 + (c1 + (c2 + c3 * x) * x) * x;

  // Construct the derivative dfa(x_n).
  auto three_c3x = 3 * c3 * x;
  auto dfa = c1 + (2 * c2 + three_c3x) * x;

  std::cout << "dfa = " << dfa << std::endl;

  // And ga.
  FA ga = x - fa / dfa;

  Dout(dc::notice, "fa = " << fa);
  Dout(dc::notice, "dfa = " << dfa);
  Dout(dc::notice, "ga = " << ga);

  x.symbol() = 0.1;
  for (symbolic::Symbol const* epsilon : expression_manager.epsilons())
    *epsilon = 0.0;

  // df^2 =                                   c₁²  + 4c₁c₂x  + 6c₁c₃x² + 4c₂²x² + 12c₂c₃x³ + 9c₃²x⁴

  // gₐ * df^2 =
  //               -c₀c₁ - 2c₀c₂x - 3c₀c₃x²        +  c₁c₂x² + 2c₁c₃x³ + 2c₂²x³ +  7c₂c₃x⁴ + 6c₃²x⁵ +
  //           ε₀ (                                          -  c₁c₃x³          -  2c₂c₃x⁴ - 3c₃²x⁵) +
  //     (ε₁ + ε₂)(                                -  c₁c₂x² -  c₁c₃x³ - 2c₂²x³ -  5c₂c₃x⁴ - 3c₃²x⁵) +
  //     (ε₃ + ε₄)(                          -c₁²x - 3c₁c₂x² - 4c₁c₃x³ - 2c₂²x³ -  5c₂c₃x⁴ - 3c₃²x⁵) +
  //           ε₅ (-c₀c₁ - 2c₀c₂x - 3c₀c₃x² - c₁²x - 3c₁c₂x² - 4c₁c₃x³ - 2c₂²x³ -  5c₂c₃x⁴ - 3c₃²x⁵) +
  //           ε₆ (                 3c₀c₃x²                  + 3c₁c₃x³          +  3c₂c₃x⁴ + 3c₃²x⁵)
  //     (ε₇ + ε₈)(        2c₀c₂x + 3c₀c₃x²        + 2c₁c₂x² + 3c₁c₃x³ + 2c₂²x³ +  5c₂c₃x⁴ + 3c₃²x⁵) +
  //           ε₉ ( c₀c₁ + 2c₀c₂x + 3c₀c₃x² + c₁²x + 3c₁c₂x² + 4c₁c₃x³ + 2c₂²x³ +  5c₂c₃x⁴ + 3c₃²x⁵) +

  // The index of the first epsilon in the above expression per row.
  std::array<size_t, 7> epsilon_index = { 0, 1, 3, 5, 6, 7, 9 };

  //------------------------------------------------------------------------------------------------
  // Verify the above expression.

  // Construct the accurate cubic f(x_n).
  auto& f = fa.zero_order();
  std::cout << "f(x) = " << f << std::endl;

  // And its derivative.
  auto& df = f.derivative(x.symbol());
  std::cout << "df(x) = " << df << std::endl;

  auto ga_df2 = ga.reset_denominator();
  std::cout << "gₐ * df² = " << ga_df2 << std::endl;

  auto& symbolic_ga_df2 = ga_df2.expand();

  // The factor of coefficient pairs.
  std::array<std::array<int, 8>, 10> factors = {{
    { 0,  0,  0,  0,  0,  0,  0,  0},   // c₀²
    {-1,  0,  0,  0, -1,  0,  0,  1},   // c₀c₁
    {-2,  0,  0,  0, -2,  0,  2,  2},   // c₀c₂
    {-3,  0,  0,  0, -3,  3,  3,  3},   // c₀c₃
    { 0,  0,  0, -1, -1,  0,  0,  1},   // c₁²
    { 1,  0, -1, -3, -3,  0,  2,  3},   // c₁c₂
    { 2, -1, -1, -4, -4,  3,  3,  4},   // c₁c₃
    { 2,  0, -2, -2, -2,  0,  2,  2},   // c₂²
    { 7, -2, -5, -5, -5,  3,  5,  5},   // c₂c₃
    { 6, -3, -3, -3, -3,  3,  3,  3}    // c₃²
  }};

  // Run over all cᵢcⱼ.
  int col = 0;
  for (auto i = c0.index_; i <= c3.index_; ++i)
  {
    for (auto j = i; j <= c3.index_; ++j)
    {
      int dc = i == j ? 2 : 1;  // Derivative correction.
      double estimate = symbolic_ga_df2
        .derivative(expression_manager.get_variable(i).symbol())
        .derivative(expression_manager.get_variable(j).symbol()).evaluate();
      ASSERT(utils::almost_equal(estimate, dc * factors[col][0] * std::pow(0.1, i.get_value() + j.get_value() - 1), 1e-6) ||
          (factors[col][0] == 0 && std::abs(estimate) < 1e-9));
      ++col;
    }
  }

  // The top-row is the zero-order part.
  auto& ga_df2_no_epsilons = ga_df2.zero_order();

  // Get the rest; the part with epsilons.
  auto& ga_df2_epsilon = symbolic_ga_df2 - ga_df2_no_epsilons;
  col = 0;
  for (auto i = c0.index_; i <= c3.index_; ++i)
  {
    for (auto j = i; j <= c3.index_; ++j)
    {
      int dc = i == j ? 2 : 1;  // Derivative correction.
      auto& dg_dcidcj = ga_df2_epsilon
        .derivative(expression_manager.get_variable(i).symbol())
        .derivative(expression_manager.get_variable(j).symbol());
      // Run over the table rows (k = 0 is the second row).
      for (int k = 0; k <= 6; ++k)
      {
        FA::epsilon_index_type eps{epsilon_index[k]};
        *expression_manager.epsilons()[eps] = 1.0;
        double estimate = dg_dcidcj.evaluate();
        ASSERT(utils::almost_equal(estimate, dc * factors[col][k + 1] * std::pow(0.1, i.get_value() + j.get_value() - 1), 1e-6) ||
            (factors[col][k + 1] == 0 && std::abs(estimate) < 1e-9));
        *expression_manager.epsilons()[eps] = 0.0;
      }
      ++col;
    }
  }

  std::cout << "Newtons's method:" << std::endl;
  std::cout << "gₐ = " << symbolic::UseUtf8{1} << ga << std::endl;
  std::cout << "Δgₐ = ±√(" << symbolic::UseUtf8{1} << ga.error_squared() << ')' << std::endl;

  //------------------------------------------------------------------------------------------------
  // Halley's method.
  FA x_n_plus_1 = x - (fa * dfa) / (utils::square(dfa) - fa * (c2 + three_c3x));

  std::cout << "Halley's method:" << std::endl;
  std::cout << "xₙ₊₁ = " << symbolic::UseUtf8{1} << x_n_plus_1 << std::endl;
  std::cout << "Δxₙ₊₁ = ±√(" << symbolic::UseUtf8{1} << x_n_plus_1.error_squared() << ')' << std::endl;
}
