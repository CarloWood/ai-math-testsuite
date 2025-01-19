#include "sys.h"
#include "FloatingPointRoundOffError.h"
#include "utils/macros.h"
#include "utils/almost_equal.h"
#include "utils/EnumIterator.h"
#include <utility>
#include "debug.h"

struct VariablesCategory { };
using VariablesIndex = utils::ArrayIndex<VariablesCategory>;

constexpr int number_of_epsilons = 21;

#define FOREACH_EPSILON(X) \
  X(e0) \
  X(e1) \
  X(e2) \
  X(e3) \
  X(e4) \
  X(e5) \
  X(e6) \
  X(e7) \
  X(e8) \
  X(e9) \
  X(e10) \
  X(e11) \
  X(e12) \
  X(e13) \
  X(e14) \
  X(e15) \
  X(e16) \
  X(e17) \
  X(e18) \
  X(e19) \
  X(e20)

#define FOREACH_VARIABLE(X) \
  X(c0) \
  X(c1) \
  X(c2) \
  X(c3) \
  X(x)

#define DECLARE_VARIABLE_ENUM(name) name,

enum class VariablesEnum
{
  FOREACH_VARIABLE(DECLARE_VARIABLE_ENUM)
  FOREACH_EPSILON(DECLARE_VARIABLE_ENUM)
  count
};

VariablesIndex epsilon_index_to_variables_index(FloatingPointRoundOffError<number_of_epsilons>::epsilon_index_type i)
{
  return VariablesIndex{static_cast<int>(VariablesEnum::e0) + i.get_value()};
}

FloatingPointRoundOffError<number_of_epsilons>::epsilon_index_type variables_index_to_epsilon_index(VariablesIndex i)
{
  return FloatingPointRoundOffError<number_of_epsilons>::epsilon_index_type{i.get_value() - static_cast<int>(VariablesEnum::e0)};
}

constexpr size_t number_of_variables = static_cast<size_t>(VariablesEnum::count);

utils::Array<symbolic::Symbol const*, number_of_variables, VariablesIndex> variables;

struct Variable
{
  VariablesIndex index;
  char const* name;
};

#define DECLARE_VARIABLE(name) Variable name{VariablesIndex{static_cast<int>(VariablesEnum::name)}, #name};

namespace Variables {
  FOREACH_VARIABLE(DECLARE_VARIABLE)
  FOREACH_EPSILON(DECLARE_VARIABLE)
} // namespace Variables

#define ADD_CASE_RETURN(name) AI_CASE_RETURN(name);

std::string to_string(VariablesEnum variable)
{
  using enum VariablesEnum;
  switch (variable)
  {
    FOREACH_VARIABLE(ADD_CASE_RETURN)
    FOREACH_EPSILON(ADD_CASE_RETURN)
    case count:   // Suppress compiler warning about count.
      ASSERT(false);
  }
  AI_NEVER_REACHED
}

#define DECLARE_SYMBOL(name) symbolic::Symbol const& name(*variables[Variables::name.index]);

// For gdb.
void print(symbolic::Expression const& expression)
{
  std::cout << expression << std::endl;
}

void print(FloatingPointRoundOffError<number_of_epsilons> const& fa)
{
  std::cout << fa << std::endl;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using FA = FloatingPointRoundOffError<number_of_epsilons>;

  for (VariablesIndex index = variables.ibegin(); index != variables.iend(); ++index)
    variables[index] = &symbolic::Symbol::realize(to_string(static_cast<VariablesEnum>(index.get_value())));

  FOREACH_VARIABLE(DECLARE_SYMBOL);

  using namespace Variables;

  utils::Array<symbolic::Symbol const*, number_of_epsilons, FA::epsilon_index_type> epsilons;
  for (FA::epsilon_index_type i = epsilons.ibegin(); i != epsilons.iend(); ++i)
    epsilons[i] = variables[epsilon_index_to_variables_index(i)];

  // Construct the approximation cubic f_a(x_n).
  FA c3x{c3 * x, variables_index_to_epsilon_index(e5.index)};
  FA c2_plus_c3x{c2 + c3x, variables_index_to_epsilon_index(e4.index)};
  FA c2_plus_c3x_x{c2_plus_c3x * x, variables_index_to_epsilon_index(e3.index)};
  FA c1_plus_c2_plus_c3x_x{c1 + c2_plus_c3x_x, variables_index_to_epsilon_index(e2.index)};
  std::cout << "c1_plus_c2_plus_c3x_x = " << c1_plus_c2_plus_c3x_x << std::endl;
  FA c1_plus_c2_plus_c3x_x_x{c1_plus_c2_plus_c3x_x * x, variables_index_to_epsilon_index(e1.index)};
  FA fa{c0 + c1_plus_c2_plus_c3x_x_x, variables_index_to_epsilon_index(e0.index)};

  // Construct the derivative dfa(x_n).
  FA three_c3x{3 * c3 * x, variables_index_to_epsilon_index(e9.index)};
  FA two_c2_plus_three_c3x{2 * c2 + three_c3x, variables_index_to_epsilon_index(e8.index)};
  FA two_c2_plus_three_c3x_x{two_c2_plus_three_c3x * x, variables_index_to_epsilon_index(e7.index)};
  FA dfa{c1 + two_c2_plus_three_c3x_x, variables_index_to_epsilon_index(e6.index)};

  // And ga.
  FA ga = x - fa / dfa;

  Dout(dc::notice, "fa = " << fa);
  Dout(dc::notice, "dfa = " << dfa);
  Dout(dc::notice, "ga = " << ga);

  x = 0.1;
  for (FA::epsilon_index_type i = epsilons.ibegin(); i != epsilons.iend(); ++i)
    *variables[epsilon_index_to_variables_index(i)] = 0.0;

  // df^2 =                                   c₁²  + 4c₁c₂x  + 6c₁c₃x² + 4c₂²x² + 12c₂c₃x³ + 9c₃²x⁴

  // gₐ * df^2 =
  //               -c₀c₁ - 2c₀c₂x - 3c₀c₃x²        +  c₁c₂x² + 2c₁c₃x³ + 2c₂²x³ +  7c₂c₃x⁴ + 6c₃²x⁵ +
  //           ε₀ (-c₀c₁ - 2c₀c₂x - 3c₀c₃x² - c₁²x - 3c₁c₂x² - 4c₁c₃x³ - 2c₂²x³ -  5c₂c₃x⁴ - 3c₃²x⁵) +
  //     (ε₁ + ε₂)(                          -c₁²x - 3c₁c₂x² - 4c₁c₃x³ - 2c₂²x³ -  5c₂c₃x⁴ - 3c₃²x⁵) +
  //     (ε₃ + ε₄)(                                -  c₁c₂x² -  c₁c₃x³ - 2c₂²x³ -  5c₂c₃x⁴ - 3c₃²x⁵) +
  //           ε₅ (                                          -  c₁c₃x³          -  2c₂c₃x⁴ - 3c₃²x⁵) +
  //           ε₆ ( c₀c₁ + 2c₀c₂x + 3c₀c₃x² + c₁²x + 3c₁c₂x² + 4c₁c₃x³ + 2c₂²x³ +  5c₂c₃x⁴ + 3c₃²x⁵) +
  //     (ε₇ + ε₈)(        2c₀c₂x + 3c₀c₃x²        + 2c₁c₂x² + 3c₁c₃x³ + 2c₂²x³ +  5c₂c₃x⁴ + 3c₃²x⁵) +
  //           ε₉ (                 3c₀c₃x²                  + 3c₁c₃x³          +  3c₂c₃x⁴ + 3c₃²x⁵)

  std::array<int, 7> epsilon_index = { 0, 1, 3, 5, 6, 7, 9 };

  //------------------------------------------------------------------------------------------------
  // Verify the above expression.

  // Construct the accurate cubic f(x_n).
  auto& f = c0 + c1 * x + c2 * (x^2) + c3 * (x^3);
  // And its derivative.
  auto& df = f.derivative(x);

  auto& ga_df2 = (ga * (df^2)).expand(epsilons);

  // Define the top-row.
  auto& ga_df2_no_epsilons = -c0 * c1 - 2 * c0 * c2 * x - 3 * c0 * c3 * (x^2) +
    c1 * c2 * (x^2) + 2 * c1 * c3 * (x^3) + 2 * (c2^2) * (x^3) + 7 * c2 * c3 * (x^4) + 6 * (c3^2) * (x^5);

  // The factor of coefficient pairs.
  std::array<std::array<int, 8>, 10> factors = {{
    {0, 0, 0, 0, 0, 0, 0},              // c₀²
    {-1, -1, 0, 0, 0, 1, 0, 0},         // c₀c₁
    {-2, -2, 0, 0, 0, 2, 2, 0},         // c₀c₂
    {-3, -3, 0, 0, 0, 3, 3, 3},         // c₀c₃
    {0, -1, -1, 0, 0, 1, 0, 0},         // c₁²
    {1, -3, -3, -1, 0, 3, 2, 0},        // c₁c₂
    {2, -4, -4, -1, -1, 4, 3, 3},       // c₁c₃
    {2, -2, -2, -2, 0, 2, 2, 0},        // c₂²
    {7, -5, -5, -5, -2, 5, 5, 3},       // c₂c₃
    {6, -3, -3, -3, -3, 3, 3, 3}        // c₃²
  }};

  // Run over all cᵢcⱼ.
  int col = 0;
  for (auto i = Variables::c0.index; i <= Variables::c3.index; ++i)
  {
    for (auto j = i; j <= Variables::c3.index; ++j)
    {
      int dc = i == j ? 2 : 1;  // Derivative correction.
      double estimate = ga_df2.derivative(*variables[i]).derivative(*variables[j]).evaluate();
      ASSERT(utils::almost_equal(estimate, dc * factors[col][0] * std::pow(0.1, i.get_value() + j.get_value() - 1), 1e-6) ||
          (factors[col][0] == 0 && std::abs(estimate) < 1e-9));
      ++col;
    }
  }
  // Get the part with epsilons.
  auto& ga_df2_epsilon = ga_df2 - ga_df2_no_epsilons;
  col = 0;
  for (auto i = Variables::c0.index; i <= Variables::c3.index; ++i)
  {
    for (auto j = i; j <= Variables::c3.index; ++j)
    {
      int dc = i == j ? 2 : 1;  // Derivative correction.
      auto& dg_dcidcj = ga_df2_epsilon.derivative(*variables[i]).derivative(*variables[j]);
      // Run over the table rows (k = 0 is the second row).
      for (int k = 0; k <= 6; ++k)
      {
        auto eps = epsilon_index_to_variables_index(FA::epsilon_index_type{epsilon_index[k]});
        *variables[eps] = 1.0;
        double estimate = dg_dcidcj.evaluate();
        ASSERT(utils::almost_equal(estimate, dc * factors[col][k + 1] * std::pow(0.1, i.get_value() + j.get_value() - 1), 1e-6) ||
            (factors[col][k + 1] == 0 && std::abs(estimate) < 1e-9));
        *variables[eps] = 0.0;
      }
      ++col;
    }
  }

  FA fa_div_dfa{fa / dfa, variables_index_to_epsilon_index(e10.index)};
  FA ga2{x - fa_div_dfa, variables_index_to_epsilon_index(e11.index)};

  std::cout << "ga2 = " << ga2 << std::endl;
  std::cout << "error = ±(" << ga2.error() << ')' << std::endl;

  //------------------------------------------------------------------------------------------------
  // Halley's method.

  // Q_x = f(x_n);
  // dQ_x = df(x_n);
  // half_ddQ_x = half_second_derivative(x_n);
  // step = Q_x * dQ_x / (utils::square(dQ_x) - Q_x * half_ddQ_x);
  // x_n_plus_1 = x_n - step;

  // Construct half times the second derivative: hddfa(x_n).
  FA hddfa{c2 + three_c3x, variables_index_to_epsilon_index(e10.index)};
  FA square_dQ_x{dfa * dfa, variables_index_to_epsilon_index(e11.index)};
  FA Q_x_dQ_x{fa * dfa, variables_index_to_epsilon_index(e12.index)};
  FA Q_x_half_ddQ_x{fa * hddfa, variables_index_to_epsilon_index(e13.index)};
  FA square_dQ_x_minus_Q_x_half_ddQ_x{square_dQ_x - Q_x_half_ddQ_x, variables_index_to_epsilon_index(e14.index)};
  FA step{Q_x_dQ_x / square_dQ_x_minus_Q_x_half_ddQ_x, variables_index_to_epsilon_index(e15.index)};
  FA x_n_plus_1{x - step, variables_index_to_epsilon_index(e16.index)};

  std::cout << "step = " << step << std::endl;
  std::cout << "step.error = ±(" << step.error() << ')' << std::endl;
}
