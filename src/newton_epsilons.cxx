#include "sys.h"
#include "utils/macros.h"
#include "utils/has_print_on.h"
#include "utils/almost_equal.h"
#include "cairowindow/symbolic/symbolic.h"
#include <array>
#include <utility>
#include "debug.h"

using utils::has_print_on::operator<<;

class FA
{
 private:
  symbolic::Expression const& zero_order_;              // Term without ε's.
  symbolic::Expression const& first_order_;             // Term with single ε factor(s).
  symbolic::Expression const& common_denominator_;      // Denominator without ε's.

 public:
  // Construct a FirstOrder object where the zero order is an integer.
  FA(int zero_order, symbolic::Expression const& first_order) :
    zero_order_(symbolic::Constant::realize(zero_order)), first_order_(first_order), common_denominator_(symbolic::Constant::s_cached_one) { }

  // Construct a FirstOrder object from two expressions.
  FA(symbolic::Expression const& zero_order, symbolic::Expression const& first_order,
      symbolic::Expression const& common_denominator = symbolic::Constant::s_cached_one) :
    zero_order_(zero_order), first_order_(first_order), common_denominator_(common_denominator) { }

  // Multiply with an expression that does NOT contain ε's.
  FA operator*(symbolic::Expression const& arg) const
  {
    return {zero_order_ * arg, first_order_ * arg, common_denominator_};
  }

  symbolic::Expression const& expand() const
  {
    return (zero_order_ + first_order_) / common_denominator_;
  }

  // Add an expression that does NOT contain ε's.
  FA operator+(symbolic::Expression const& arg) const
  {
    return {zero_order_ + arg * common_denominator_, first_order_, common_denominator_};
  }

  // Subtract an expression that does NOT contain ε's.
  FA operator-(symbolic::Expression const& arg) const
  {
    return {zero_order_ - arg * common_denominator_, first_order_, common_denominator_};
  }

  // Multiply with an expression that can contain ε's.
  // (z1 + f1)/d1 * (z2 + f2)/d2 = (z1 z2, z1 f2 + f1 z2)/(d1 d2)     (neglecting f1 * f2).
  FA operator*(FA const& fa) const
  {
    return {zero_order_ * fa.zero_order_, zero_order_ * fa.first_order_ + first_order_ * fa.zero_order_,
      common_denominator_ * fa.common_denominator_};
  }

  // Divide by an expression that can contain ε's.
  // ((z1 + f1)/d1) / ((z2 + f2)/d2) = ((z1 + f1) / (z2 + f2))/(d1 d2) =
  //                                 = ((z1 + f1)(z2 - f2) / ((z2 + f2)(z2 - f2)))/(d1 d2) =
  //                                 = (z1 z2 - z1 f2 + f1 z2) / (z2^2 d1 d2)
  FA operator/(FA const& fa) const
  {
    return {zero_order_ * fa.zero_order_, first_order_ * fa.zero_order_ - zero_order_ * fa.first_order_,
      (fa.zero_order_^2) * common_denominator_ * fa.common_denominator_};
  }

  // Add an expression that can contain ε's.
  // (z1 + f1)/d1 + (z2 + f2)/d2 = (z1 d2 + f1 d2)/(d1 d2) + (z2 d1 + f2 d1)/(d1 d2) =
  //                             = ((z1 d2 + z2 d1) + (f1 d2 + f2 d1))/(d1 d2)
  FA operator+(FA const& fa) const
  {
    return {zero_order_ * fa.common_denominator_ + fa.zero_order_ * common_denominator_,
      first_order_ * fa.common_denominator_ + fa.first_order_ * common_denominator_,
      common_denominator_ * fa.common_denominator_};
  }

  // Negation.
  FA operator-() const
  {
    return {-zero_order_, -first_order_, common_denominator_};
  }

  // The operations are commutative.
  friend FA operator*(symbolic::Expression const& arg, FA const& fa)
  {
    return fa * arg;
  }

  friend FA operator+(symbolic::Expression const& arg, FA const& fa)
  {
    return fa + arg;
  }

  friend FA operator-(symbolic::Expression const& arg, FA const& fa)
  {
    return -fa + arg;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = false;
    bool need_plus = false;
    if (!symbolic::Constant::is_zero(zero_order_))
    {
      if (!symbolic::Constant::is_zero(first_order_))
      {
        need_parens = !symbolic::Constant::is_one(common_denominator_);
        os << '(';
      }
      os << zero_order_;
      need_plus = true;
    }
    if (!symbolic::Constant::is_zero(first_order_))
    {
      if (need_plus)
        os << " + ";
      os << "\e[33m" << first_order_ << "\e[0m";
    }
    if (need_parens)
      os << ')';
    if (!symbolic::Constant::is_one(common_denominator_))
    {
      os << " / ";
      need_parens = symbolic::needs_parens(common_denominator_.precedence(), symbolic::Precedence::ratio, symbolic::after);
      if (need_parens)
        os << '(';
      os << common_denominator_;
      if (need_parens)
        os << ')';
    }
  }

  double evaluate() const
  {
    return (zero_order_.evaluate() + first_order_.evaluate()) / common_denominator_.evaluate();
  }
};

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

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  for (int v = 0; v < number_of_variables; ++v)
  {
    Variables variable = static_cast<Variables>(v);
    variables[v] = &symbolic::Symbol::realize(to_string(variable));
  }

  FOREACH_VARIABLE(DECLARE_SYMBOL);

  // Construct the accurate cubic f(x_n).
  auto& f = c0 + c1 * x + c2 * (x^2) + c3 * (x^3);
  // And its derivative.
  auto& df = f.derivative(x);
  // And the accurate g(x_n).
  auto& g = x - f / df;
  Dout(dc::notice, "f(x) = " << f);
  Dout(dc::notice, "f'(x) = " << df);
  Dout(dc::notice, "g(x) = x - f(x)/f'(x) = " << g);

  // Construct the 1 + ε factors.
  FA e0p1{1, e0};       // 1 + ε₀
  FA e1p1{1, e1};       // 1 + ε₁
  FA e2p1{1, e2};       // ⋮
  FA e3p1{1, e3};
  FA e4p1{1, e4};
  FA e5p1{1, e5};
  FA e6p1{1, e6};
  FA e7p1{1, e7};
  FA e8p1{1, e8};
  FA e9p1{1, e9};

  // Construct the approximation cubic f_a(x_n).
  FA c3x = e5p1 * c3 * x;
  FA c2_plus_c3x = e4p1 * (c2 + c3x);
  FA c2_plus_c3x_x = e3p1 * c2_plus_c3x * x;
  FA c1_plus_c2_plus_c3x_x = e2p1 * (c1 + c2_plus_c3x_x);
  FA c1_plus_c2_plus_c3x_x_x = e1p1 * c1_plus_c2_plus_c3x_x * x;
  FA fa = e0p1 * (c0 + c1_plus_c2_plus_c3x_x_x);

  // Construct the derivative dfa(x_n).
  FA three_c3x = e9p1 * (3 * c3 * x);
  FA two_c2_plus_three_c3x = e8p1 * (2 * c2 + three_c3x);
  FA two_c2_plus_three_c3x_x = e7p1 * two_c2_plus_three_c3x * x;
  FA dfa = e6p1 * (c1 + two_c2_plus_three_c3x_x);

  // And ga.
  FA ga = x - fa / dfa;

  Dout(dc::notice, "fa = " << fa);
  Dout(dc::notice, "dfa = " << dfa);
  Dout(dc::notice, "ga = " << ga);

  x = 0.1;
  e0 = 0;
  e1 = 0;
  e2 = 0;
  e3 = 0;
  e4 = 0;
  e5 = 0;
  e6 = 0;
  e7 = 0;
  e8 = 0;
  e9 = 0;

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

  // Verify the above expression.
  auto& ga_df2 = (ga * (df^2)).expand();

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
}

  // g * df^2 =
  //             (-c₀c₁        + (-2c₀  + c₁x + 2c₂x²)c₂x + (-3c₀ + 2c₁x + 7c₂x² + 6c₃x³))c₃x² +
  //         -ε₀ ( c₀c₁ + c₁²x + ( 2c₀ + 3c₁x + 2c₂x²)c₂x + ( 3c₀ + 4c₁x + 5c₂x² + 3c₃x³))c₃x² +
  //  -(ε₁ + ε₂) (        c₁²x   (       3c₁x + 2c₂x²)c₂x + (       4c₁x + 5c₂x² + 3c₃x³))c₃x² +
  //  -(ε₃ + ε₄) (               (        c₁x + 2c₂x²)c₂x + (        c₁x + 5c₂x² + 3c₃x³))c₃x² +
  //         -ε₅ (                                          (        c₁x + 2c₂x² + 3c₃x³))c₃x² +
  //          ε₆ ( c₀c₁ + c₁²x + ( 2c₀ + 3c₁x + 2c₂x²)c₂x + ( 3c₀ + 4c₁x + 5c₂x² + 3c₃x³))c₃x² +
  //   (ε₇ + ε₈) (               ( 2c₀ + 2c₁x + 2c₂x²)c₂x + ( 3c₀ + 3c₁x + 5c₂x² + 3c₃x³))c₃x² +
  //          ε₉ (                                          ( 3c₀ + 3c₁x + 3c₂x² + 3c₃x³))c₃x²

