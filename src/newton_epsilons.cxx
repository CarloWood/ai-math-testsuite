#include "sys.h"
#include "utils/macros.h"
#include "utils/has_print_on.h"
#include "utils/almost_equal.h"
#include "utils/Array.h"
#include "cairowindow/symbolic/symbolic.h"
#include <array>
#include <utility>
#include <algorithm>
#include "debug.h"

using utils::has_print_on::operator<<;

struct EpsilonIndexCategory { };
using epsilon_index_type = utils::ArrayIndex<EpsilonIndexCategory>;

constexpr int number_of_epsilons = 21;

class FA
{
 private:
  symbolic::Expression const& zero_order_;                                                              // Term without ε's.
  utils::Array<symbolic::Expression const*, number_of_epsilons, epsilon_index_type> first_order_{};     // Terms with single ε factor(s).
  symbolic::Expression const& common_denominator_;                                                      // Denominator without ε's.

 public:
  // Construct a FirstOrder object where the zero order is an integer.
  FA(int zero_order) :
    zero_order_(symbolic::Constant::realize(zero_order)), common_denominator_(symbolic::Constant::s_cached_one) { }

  // Construct a FirstOrder object from an expressions that is to be multiplied by (1 + ε).
  FA(symbolic::Expression const& expression, epsilon_index_type epsilon_index,
      symbolic::Expression const& common_denominator = symbolic::Constant::s_cached_one) :
    zero_order_(expression), common_denominator_(common_denominator)
  {
    first_order_[epsilon_index] = &expression;
  }

  // Construct a FirstOrder object from an FA object that is to be multiplied by (1 + ε).
  FA(FA const& fa, epsilon_index_type epsilon_index) :
    zero_order_(fa.zero_order_), first_order_(fa.first_order_), common_denominator_(fa.common_denominator_)
  {
    // Don't use an already used epsilon_index.
    ASSERT(first_order_[epsilon_index] == nullptr);
    first_order_[epsilon_index] = &fa.zero_order_;
  }

 private:
  FA(symbolic::Expression const& zero_order, symbolic::Expression const& common_denominator) :
    zero_order_(zero_order), common_denominator_(common_denominator) { }

  FA(symbolic::Expression const& zero_order,
      utils::Array<symbolic::Expression const*, number_of_epsilons, epsilon_index_type> const& first_order,
      symbolic::Expression const& common_denominator) :
    zero_order_(zero_order), first_order_(first_order), common_denominator_(common_denominator) { }

 public:
  // Multiply with an expression that does NOT contain ε's.
  FA operator*(symbolic::Expression const& arg) const
  {
    FA result{zero_order_ * arg, common_denominator_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(*first_order_[i] * arg);
    return result;
  }

  symbolic::Expression const& expand(utils::Array<symbolic::Symbol const*, number_of_epsilons, epsilon_index_type> const& epsilon) const
  {
    symbolic::Expression const* result = &zero_order_;
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result = &(*result + *first_order_[i] * *epsilon[i]);
    return *result / common_denominator_;
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
    FA result{zero_order_ * fa.zero_order_, common_denominator_ * fa.common_denominator_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(*first_order_[i] * fa.zero_order_);
    for (epsilon_index_type i = fa.first_order_.ibegin(); i != fa.first_order_.iend(); ++i)
      if (fa.first_order_[i] != nullptr)
      {
        if (result.first_order_[i] == nullptr)
          result.first_order_[i] = &(zero_order_ * *fa.first_order_[i]);
        else
          result.first_order_[i] = &(*result.first_order_[i] + zero_order_ * *fa.first_order_[i]);
      }
    return result;
  }

  // Divide by an expression that can contain ε's.
  // ((z1 + f1)/d1) / ((z2 + f2)/d2) = ((z1 + f1) / (z2 + f2))/(d1 / d2) =
  //                                 = ((z1 + f1)(z2 - f2) / ((z2 + f2)(z2 - f2)))/(d1 / d2) =
  //                                 = (z1 z2 - z1 f2 + f1 z2) / (z2^2 d1 / d2) =
  //                                 = (z1 z2 d2 - z1 f2 d2 + f1 z2 d2) / (z2^2 d1)
  // Note, here
  //   z1 = zero_order_
  //   f1 = first_order_
  //   d1 = common_denominator_
  //   z2 = fa.zero_order_
  //   f2 = fa.first_order_
  //   d2 = fa.common_denominator_
  //
  // If f2 only depends on a single ε, then it will be equal to z2 * ε2, so that the result will be
  //  (z1 z2 d2 - z1 z2⋅ε2 d2 + f1 z2 d2) / (z2^2 d1) = (z1 d2 - z1 d2⋅ε2 + f1 d2) / (z2 d1)
  //
  FA operator/(FA const& fa) const
  {
    // Check if fa.first_order_ only contains a single non-zero pointer.
    bool f2_is_single_epsilon = std::ranges::count_if(fa.first_order_, [](auto p) { return p != nullptr; }) == 1;
    if (f2_is_single_epsilon)
    {
      FA result{zero_order_ * fa.common_denominator_, fa.zero_order_ * common_denominator_};
      for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
        if (first_order_[i] != nullptr)
          result.first_order_[i] = &(*first_order_[i] * fa.common_denominator_);
      for (epsilon_index_type i = fa.first_order_.ibegin(); i != fa.first_order_.iend(); ++i)
        if (fa.first_order_[i] != nullptr)
        {
          ASSERT(fa.first_order_[i]->equals(fa.zero_order_));
          if (result.first_order_[i] == nullptr)
            result.first_order_[i] = &(-zero_order_ * fa.common_denominator_);
          else
            result.first_order_[i] = &(*result.first_order_[i] - zero_order_ * fa.common_denominator_);
        }
      return result;
    }
    else
    {
      FA result{zero_order_ * fa.zero_order_ * fa.common_denominator_, (fa.zero_order_^2) * common_denominator_};
      for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
        if (first_order_[i] != nullptr)
          result.first_order_[i] = &(*first_order_[i] * fa.zero_order_ * fa.common_denominator_);
      for (epsilon_index_type i = fa.first_order_.ibegin(); i != fa.first_order_.iend(); ++i)
        if (fa.first_order_[i] != nullptr)
        {
          if (result.first_order_[i] == nullptr)
            result.first_order_[i] = &(-zero_order_ * *fa.first_order_[i] * fa.common_denominator_);
          else
            result.first_order_[i] = &(*result.first_order_[i] - zero_order_ * *fa.first_order_[i] * fa.common_denominator_);
        }
      return result;
    }
  }

  // Add an expression that can contain ε's.
  // (z1 + f1)/d1 + (z2 + f2)/d2 = (z1 d2 + f1 d2)/(d1 d2) + (z2 d1 + f2 d1)/(d1 d2) =
  //                             = ((z1 d2 + z2 d1) + (f1 d2 + f2 d1))/(d1 d2)
  FA operator+(FA const& fa) const
  {
    FA result{zero_order_ * fa.common_denominator_ + fa.zero_order_ * common_denominator_,
      common_denominator_ * fa.common_denominator_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(*first_order_[i] * fa.common_denominator_);
    for (epsilon_index_type i = fa.first_order_.ibegin(); i != fa.first_order_.iend(); ++i)
      if (fa.first_order_[i] != nullptr)
      {
        if (result.first_order_[i] == nullptr)
          result.first_order_[i] = &(*fa.first_order_[i] * common_denominator_);
        else
          result.first_order_[i] = &(*result.first_order_[i] + *fa.first_order_[i] * common_denominator_);
      }

    return result;
  }

  FA operator-(FA const& fa) const
  {
    return *this + (-fa);
  }

  // The operations are commutative.
  friend FA operator*(symbolic::Expression const& arg, FA const& fa)
  {
    return fa * arg;
  }

  // Negation.
  FA operator-() const
  {
    FA result{-zero_order_, common_denominator_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(-*first_order_[i]);
    return result;
  }

  friend FA operator+(symbolic::Expression const& arg, FA const& fa)
  {
    return fa + arg;
  }

  friend FA operator-(symbolic::Expression const& arg, FA const& fa)
  {
    return -fa + arg;
  }

  symbolic::Expression const& zero_order() const { return zero_order_ / common_denominator_; }

  void print_on(std::ostream& os) const
  {
    static utils::Array<char const*, number_of_epsilons, epsilon_index_type> subscript = {
      "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉", "₁₀", "₁₁", "₁₂", "₁₃", "₁₄", "₁₅", "₁₆", "₁₇", "₁₈", "₁₉", "₂₀"
    };
    bool need_parens = false;
    bool need_plus = false;
    bool first_order_is_zero = true;
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
    {
      if (first_order_[i] != nullptr && !symbolic::Constant::is_zero(*first_order_[i]))
      {
        first_order_is_zero = false;
        break;
      }
    }
    if (!symbolic::Constant::is_zero(zero_order_))
    {
      if (!first_order_is_zero)
      {
        need_parens = !symbolic::Constant::is_one(common_denominator_);
        if (need_parens)
          os << '(';
      }
      os << zero_order_;
      need_plus = true;
    }
    if (!first_order_is_zero)
    {
      char const* prefix = need_plus ? " + " : "";
      os << "\e[33m";
      for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      {
        if (first_order_[i] != nullptr && !symbolic::Constant::is_zero(*first_order_[i]))
          os << prefix << '(' << *first_order_[i] << ")⋅ε" << subscript[i];
        prefix = " + ";
      }
      os << "\e[0m";
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

  symbolic::Expression const& error() const
  {
    symbolic::Expression const* sum_of_squared_first_order_elements = &symbolic::Constant::s_cached_zero;
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        sum_of_squared_first_order_elements = &(*sum_of_squared_first_order_elements + ((*first_order_[i])^2));
    return *sum_of_squared_first_order_elements / (common_denominator_^2);
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

void print(FA const& fa)
{
  std::cout << fa << std::endl;
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

  epsilon_index_type e0i(0);
  epsilon_index_type e1i(1);
  epsilon_index_type e2i(2);
  epsilon_index_type e3i(3);
  epsilon_index_type e4i(4);
  epsilon_index_type e5i(5);
  epsilon_index_type e6i(6);
  epsilon_index_type e7i(7);
  epsilon_index_type e8i(8);
  epsilon_index_type e9i(9);
  epsilon_index_type e10i(10);
  epsilon_index_type e11i(11);
  epsilon_index_type e12i(12);
  epsilon_index_type e13i(13);
  epsilon_index_type e14i(14);
  epsilon_index_type e15i(15);
  epsilon_index_type e16i(16);
  epsilon_index_type e17i(17);
  epsilon_index_type e18i(18);
  epsilon_index_type e19i(19);
  epsilon_index_type e20i(20);

  utils::Array<symbolic::Symbol const*, number_of_epsilons, epsilon_index_type> epsilons;
  for (epsilon_index_type i = epsilons.ibegin(); i != epsilons.iend(); ++i)
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
  for (epsilon_index_type i = epsilons.ibegin(); i != epsilons.iend(); ++i)
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

  // Verify the above expression.
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
