#include "sys.h"
#include "FloatingPointRoundOffError.h"
#include "utils/macros.h"
#include "utils/almost_equal.h"
#include "utils/EnumIterator.h"
#include <utility>
#include "debug.h"

constexpr int number_of_epsilons = 21;

#define FOREACH_VARIABLE(X) \
  X(c0) \
  X(c1) \
  X(c2) \
  X(c3) \
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
  X(e20) \
  X(x)

#define DECLARE_VARIABLE(name) name,
#define ADD_CASE_RETURN(name) AI_CASE_RETURN(name);
#define DECLARE_SYMBOL(name) symbolic::Symbol const& name(*variables[to_index(Variables::name)]);

enum class Variables
{
  FOREACH_VARIABLE(DECLARE_VARIABLE)

  count
};

constexpr size_t number_of_variables = static_cast<size_t>(Variables::count);

std::string to_string(Variables variable)
{
  using enum Variables;
  switch (variable)
  {
    FOREACH_VARIABLE(ADD_CASE_RETURN)
    case count:   // Suppress compiler warning about count.
      ASSERT(false);
  }
  AI_NEVER_REACHED
}

std::array<symbolic::Symbol const*, number_of_variables> variables;

int to_index(Variables variable)
{
  return static_cast<int>(variable);
}

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

  for (Variables variable : utils::EnumIterator<Variables, Variables::c0, Variables::x>())
//  for (int v = 0; v < number_of_variables; ++v)
  {
//    Variables variable = static_cast<Variables>(v);
    variables[variable] = &symbolic::Symbol::realize(to_string(variable));
  }

  FOREACH_VARIABLE(DECLARE_SYMBOL);

  FA::epsilon_index_type e0i(0);
  FA::epsilon_index_type e1i(1);
  FA::epsilon_index_type e2i(2);
  FA::epsilon_index_type e3i(3);
  FA::epsilon_index_type e4i(4);
  FA::epsilon_index_type e5i(5);
  FA::epsilon_index_type e6i(6);
  FA::epsilon_index_type e7i(7);
  FA::epsilon_index_type e8i(8);
  FA::epsilon_index_type e9i(9);
  FA::epsilon_index_type e10i(10);
  FA::epsilon_index_type e11i(11);
  FA::epsilon_index_type e12i(12);
  FA::epsilon_index_type e13i(13);
  FA::epsilon_index_type e14i(14);
  FA::epsilon_index_type e15i(15);
  FA::epsilon_index_type e16i(16);
  FA::epsilon_index_type e17i(17);
  FA::epsilon_index_type e18i(18);
  FA::epsilon_index_type e19i(19);
  FA::epsilon_index_type e20i(20);

  utils::Array<symbolic::Symbol const*, number_of_epsilons, FA::epsilon_index_type> epsilons;
  for (FA::epsilon_index_type i = epsilons.ibegin(); i != epsilons.iend(); ++i)
    epsilons[i] = variables[to_index(Variables::e0) + i.get_value()];

  // Construct the approximation cubic f_a(x_n).
  FA c3x{c3 * x, e5i};
  FA c2_plus_c3x{c2 + c3x, e4i};
  FA c2_plus_c3x_x{c2_plus_c3x * x, e3i};
  FA c1_plus_c2_plus_c3x_x{c1 + c2_plus_c3x_x, e2i};
  std::cout << "c1_plus_c2_plus_c3x_x = " << c1_plus_c2_plus_c3x_x << std::endl;
  FA c1_plus_c2_plus_c3x_x_x{c1_plus_c2_plus_c3x_x * x, e1i};
  FA fa{c0 + c1_plus_c2_plus_c3x_x_x, e0i};

  // Construct the derivative dfa(x_n).
  FA three_c3x{3 * c3 * x, e9i};
  FA two_c2_plus_three_c3x{2 * c2 + three_c3x, e8i};
  FA two_c2_plus_three_c3x_x{two_c2_plus_three_c3x * x, e7i};
  FA dfa{c1 + two_c2_plus_three_c3x_x, e6i};

  // And ga.
  FA ga = x - fa / dfa;

  Dout(dc::notice, "fa = " << fa);
  Dout(dc::notice, "dfa = " << dfa);
  Dout(dc::notice, "ga = " << ga);

  x = 0.1;
  for (FA::epsilon_index_type i = epsilons.ibegin(); i != epsilons.iend(); ++i)
    *variables[to_index(Variables::e0) + i.get_value()] = 0.0;

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
  for (int i = 0; i <= 3; ++i)
  {
    for (int j = i; j <= 3; ++j)
    {
      int dc = i == j ? 2 : 1;  // Derivative correction.
      double estimate = ga_df2.derivative(*variables[i]).derivative(*variables[j]).evaluate();
      ASSERT(utils::almost_equal(estimate, dc * factors[col][0] * std::pow(0.1, i + j - 1), 1e-6) ||
          (factors[col][0] == 0 && std::abs(estimate) < 1e-9));
      ++col;
    }
  }
  // Get the part with epsilons.
  auto& ga_df2_epsilon = ga_df2 - ga_df2_no_epsilons;
  col = 0;
  for (int i = 0; i <= 3; ++i)
  {
    for (int j = i; j <= 3; ++j)
    {
      int dc = i == j ? 2 : 1;  // Derivative correction.
      auto& dg_dcidcj = ga_df2_epsilon.derivative(*variables[i]).derivative(*variables[j]);
      // Run over the table rows (k = 0 is the second row).
      for (int k = 0; k <= 6; ++k)
      {
        int eps = epsilon_index[k];
        *variables[eps + 4] = 1.0;
        double estimate = dg_dcidcj.evaluate();
        ASSERT(utils::almost_equal(estimate, dc * factors[col][k + 1] * std::pow(0.1, i + j - 1), 1e-6) ||
            (factors[col][k + 1] == 0 && std::abs(estimate) < 1e-9));
        *variables[eps + 4] = 0.0;
      }
      ++col;
    }
  }

  FA fa_div_dfa{fa / dfa, e10i};
  FA ga2{x - fa_div_dfa, e11i};

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
  FA hddfa{c2 + three_c3x, e10i};
  FA square_dQ_x{dfa * dfa, e11i};
  FA Q_x_dQ_x{fa * dfa, e12i};
  FA Q_x_half_ddQ_x{fa * hddfa, e13i};
  FA square_dQ_x_minus_Q_x_half_ddQ_x{square_dQ_x - Q_x_half_ddQ_x, e14i};
  FA step{Q_x_dQ_x / square_dQ_x_minus_Q_x_half_ddQ_x, e15i};
  FA x_n_plus_1{x - step, e16i};

  std::cout << "step = " << step << std::endl;
  std::cout << "step.error = ±(" << step.error() << ')' << std::endl;
}
