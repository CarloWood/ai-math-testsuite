#pragma once

#include "utils/has_print_on.h"
#include "utils/Array.h"
#include "utils/Vector.h"
#include "utils/to_digits_string.h"
#include "cairowindow/symbolic/symbolic.h"
#include <algorithm>

#ifndef FOREACH_VARIABLE
// Define FOREACH_VARIABLE to a list of variables that are used in the program.
//
// For example if you have variables x, y, and z, then you should define FOREACH_VARIABLE as:
//
// #define FOREACH_VARIABLE(X) \
//  X(x) \
//  X(y) \
//  X(z)
//
#error "FOREACH_VARIABLE is not defined"
#endif // FOREACH_VARIABLE

namespace floating_point_round_off_error {

struct VariablesCategory { };
using VariablesIndex = utils::ArrayIndex<VariablesCategory>;

// class VariablesEnum is used to convert a name to an index.
#define DECLARE_VARIABLE_ENUM(name) name,

enum class VariablesEnum
{
  FOREACH_VARIABLE(DECLARE_VARIABLE_ENUM)
  count
};

constexpr size_t number_of_variables = static_cast<size_t>(VariablesEnum::count);

class ExpressionManager;

struct Variable
{
  VariablesIndex const index_;
  char const* const name_;
  symbolic::Symbol const* const symbol_;
  ExpressionManager* manager_;

  Variable(VariablesEnum index, char const* name, symbolic::Symbol const* symbol, ExpressionManager& manager) :
    index_(static_cast<int>(index)), name_(name), symbol_(symbol), manager_(&manager) { }

  // Accessors.
  symbolic::Symbol const& symbol() const { return *symbol_; }
  ExpressionManager* manager() const { return manager_; }

  void print_on(std::ostream& os) const
  {
    os << name_ << " (index:" << index_ << ", symbol:" << *symbol_ << ')';
  }
};

// This is used in the constructor of ExpressionManager to create all Variable objects.
#define CREATE_VARIABLE(name) Variable{VariablesEnum::name, #name, &symbolic::Symbol::realize(#name), *this},

// Declare expression_manager before using FOREACH_VARIABLE(DECLARE_VARIABLE).
#define DECLARE_VARIABLE(name) Variable const& name = expression_manager.get_variable(VariablesEnum::name);

#define ADD_CASE_RETURN(name) AI_CASE_RETURN(name);

std::string to_string(VariablesEnum variable)
{
  using enum VariablesEnum;
  switch (variable)
  {
    FOREACH_VARIABLE(ADD_CASE_RETURN)
    case count:   // Suppress compiler warning about count.
      ASSERT(false);
  }
  AI_NEVER_REACHED
}

using utils::has_print_on::operator<<;

struct EpsilonIndexCategory { };

class ExpressionManager
{
 public:
  using epsilon_index_type = utils::VectorIndex<EpsilonIndexCategory>;

 private:
  utils::Array<Variable, number_of_variables, VariablesIndex> const variables_;
  epsilon_index_type next_epsilon_index_{0};                                                    // The next ε index that is still unused.
  using epsilons_type = utils::Vector<symbolic::Symbol const*, epsilon_index_type>;
  epsilons_type epsilons_;                                                                      // The ε Symbols per index.

 public:
  ExpressionManager() : variables_{
      FOREACH_VARIABLE(CREATE_VARIABLE) // Create all Variable objects for c0, c1, c2, c3 and x.
    }
  {
  }

  epsilon_index_type next_epsilon_index()
  {
    epsilons_.push_back(&symbolic::Symbol::realize("ε" + utils::to_subscript_string(next_epsilon_index_.get_value())));
    return next_epsilon_index_++;
  }

  epsilons_type const& epsilons() const { return epsilons_; }

  static VariablesIndex to_index(VariablesEnum variable)
  {
    return VariablesIndex{static_cast<int>(variable)};
  }

  Variable const& get_variable(VariablesIndex index) const
  {
    return variables_[index];
  }

  Variable const& get_variable(VariablesEnum variable) const
  {
    return get_variable(to_index(variable));
  }

  VariablesIndex get_index(char const* name) const
  {
    for (Variable const& variable : variables_)
      if (strcmp(variable.name_, name) == 0)
        return variable.index_;
    AI_NEVER_REACHED
  }

  void print_on(std::ostream& os) const
  {
    os << "{Variables: ";
    char const* separator = "";
    for (Variable const& variable : variables_)
    {
      os << separator << variable;
      separator = ", ";
    }
  }
};

// Forward declare friend functions.
class Expression;
Expression operator+(Variable const& var1, Variable const& var2);
Expression operator-(Variable const& var1, Variable const& var2);
Expression operator*(Variable const& var1, Variable const& var2);
Expression operator/(Variable const& var1, Variable const& var2);
Expression operator+(Variable const& var, Expression const& fa);
Expression operator-(Variable const& var, Expression const& fa);
Expression operator*(Variable const& var, Expression const& fa);
Expression operator/(Variable const& var, Expression const& fa);
Expression operator+(symbolic::Expression const& arg, Expression const& fa);
Expression operator-(symbolic::Expression const& arg, Expression const& fa);
Expression operator*(symbolic::Expression const& arg, Expression const& fa);
Expression operator/(symbolic::Expression const& arg, Expression const& fa);

inline symbolic::Expression const& operator*(int value, Variable const& var)
{
  return value * var.symbol();
}

Expression operator*(symbolic::Expression const& expression, Variable const& var);

class Expression
{
 public:
  using epsilon_index_type = utils::VectorIndex<EpsilonIndexCategory>;

 private:
  ExpressionManager* manager_{nullptr};
  symbolic::Expression const& zero_order_;                                                      // Term without ε's.
  using first_order_type = utils::Vector<symbolic::Expression const*, epsilon_index_type>;
  first_order_type first_order_{};                                                              // Terms with single ε factor(s).
  symbolic::Expression const& common_denominator_;                                              // Denominator without ε's.

 public:
  // Construct an Expression object where the zero order is an integer.
  Expression(int zero_order) :
    zero_order_(symbolic::Constant::realize(zero_order)), common_denominator_(symbolic::Constant::s_cached_one) { }

 private:
  // Used by reset_denominator
  Expression(symbolic::Expression const& zero_order, first_order_type const& first_order, ExpressionManager* manager) :
    manager_(manager), zero_order_(zero_order), first_order_(first_order), common_denominator_(symbolic::Constant::s_cached_one)
  {
  }

  // Construct an Expression object from an expression that is to be multiplied by (1 + ε).
  // The result is therefore zero_order = first_order_[ε] = expression.
  Expression(symbolic::Expression const& expression, ExpressionManager* manager) :
    manager_(manager), zero_order_(expression), common_denominator_(symbolic::Constant::s_cached_one)
  {
    epsilon_index_type epsilon_index = manager->next_epsilon_index();
    grow_first_order(epsilon_index);
    first_order_[epsilon_index] = &expression;
  }

  Expression(symbolic::Expression const& zero_order, symbolic::Expression const& common_denominator, ExpressionManager* manager) :
    manager_(manager), zero_order_(zero_order), common_denominator_(common_denominator)
  {
    epsilon_index_type epsilon_index = manager->next_epsilon_index();
    grow_first_order(epsilon_index);
    first_order_[epsilon_index] = &zero_order_;
  }

  Expression(symbolic::Expression const& zero_order, first_order_type const& first_order,
      symbolic::Expression const& common_denominator, ExpressionManager* manager) :
    manager_(manager), zero_order_(zero_order), first_order_(first_order), common_denominator_(common_denominator)
  {
    epsilon_index_type epsilon_index = manager->next_epsilon_index();
    grow_first_order(epsilon_index);
    first_order_[epsilon_index] = &zero_order_;
  }

 public:
  // Multiply with whatever is in the denominator.
  Expression reset_denominator() const
  {
    return {zero_order_, first_order_, symbolic::Constant::s_cached_one, manager_};
  }

  // Add two variables (i.e. without ε's).
  friend Expression operator+(Variable const& var1, Variable const& var2)
  {
    // Can only add variables that have the same manager.
    ASSERT(var1.manager() == var2.manager());
    return {var1.symbol() + var2.symbol(), var1.manager()};
  }
  // Subtract two variables (i.e. without ε's).
  friend Expression operator-(Variable const& var1, Variable const& var2)
  {
    // Can only add variables that have the same manager.
    ASSERT(var1.manager() == var2.manager());
    return {var1.symbol() - var2.symbol(), var1.manager()};
  }
  // Multiply two variables (i.e. without ε's).
  friend Expression operator*(Variable const& var1, Variable const& var2)
  {
    // Can only multiply variables that have the same manager.
    ASSERT(var1.manager() == var2.manager());
    return {var1.symbol() * var2.symbol(), var1.manager()};
  }
  // Divide two variables (i.e. without ε's).
  friend Expression operator/(Variable const& var1, Variable const& var2)
  {
    // Can only multiply variables that have the same manager.
    ASSERT(var1.manager() == var2.manager());
    return {var1.symbol() / var2.symbol(), var1.manager()};
  }

  Expression operator+(Variable const& var) const
  {
    // Can only add variables that have the same manager.
    ASSERT(var.manager() == manager_);
    return *this + var.symbol();
  }

  Expression operator-(Variable const& var) const
  {
    // Can only subtract variables that have the same manager.
    ASSERT(var.manager() == manager_);
    return *this - var.symbol();
  }

  Expression operator*(Variable const& var) const
  {
    // Can only multiply variables that have the same manager.
    ASSERT(var.manager() == manager_);
    return *this * var.symbol();
  }

  Expression operator/(Variable const& var) const
  {
    // Can only divide variables that have the same manager.
    ASSERT(var.manager() == manager_);
    return *this / var.symbol();
  }

  friend Expression operator+(Variable const& var, Expression const& fa)
  {
    return fa + var.symbol();
  }

  friend Expression operator-(Variable const& var, Expression const& fa)
  {
    return -fa + var.symbol();
  }

  friend Expression operator*(Variable const& var, Expression const& fa)
  {
    return fa * var.symbol();
  }

  friend Expression operator/(Variable const& var, Expression const& fa)
  {
    return var.symbol() / fa;
  }

  friend Expression operator*(symbolic::Expression const& expression, Variable const& var)
  {
    return {expression * var.symbol(), var.manager()};
  }

  // Add an expression that does NOT contain ε's.
  Expression operator+(symbolic::Expression const& arg) const
  {
    return {zero_order_ + arg * common_denominator_, first_order_, common_denominator_, manager_};
  }

  // Subtract an expression that does NOT contain ε's.
  Expression operator-(symbolic::Expression const& arg) const
  {
    return {zero_order_ - arg * common_denominator_, first_order_, common_denominator_, manager_};
  }

  // Multiply with an expression that does NOT contain ε's.
  Expression operator*(symbolic::Expression const& arg) const
  {
    Expression result{zero_order_ * arg, common_denominator_, manager_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(*first_order_[i] * arg);
    return result;
  }

  // Divide by an expression that does NOT contain ε's.
  Expression operator/(symbolic::Expression const& arg) const
  {
    return {zero_order_, first_order_, common_denominator_ * arg, manager_};
  }

  friend Expression operator+(symbolic::Expression const& arg, Expression const& fa)
  {
    return fa + arg;
  }

  friend Expression operator-(symbolic::Expression const& arg, Expression const& fa)
  {
    return -fa + arg;
  }

  // The operations are commutative.
  friend Expression operator*(symbolic::Expression const& arg, Expression const& fa)
  {
    return fa * arg;
  }

  friend Expression operator/(symbolic::Expression const& arg, Expression const& fa)
  {
    bool f2_is_single_epsilon = std::ranges::count_if(fa.first_order_, [](auto p) { return p != nullptr; }) == 1;
    if (f2_is_single_epsilon)
    {
      Expression result{arg * fa.common_denominator_, fa.zero_order_, fa.manager_};
      for (epsilon_index_type i = fa.first_order_.ibegin(); i != fa.first_order_.iend(); ++i)
        if (fa.first_order_[i] != nullptr)
          result.first_order_[i] = &(-arg * fa.common_denominator_);
      return result;
    }
    else
    {
      Expression result{arg * fa.zero_order_ * fa.common_denominator_, (fa.zero_order_^2), fa.manager_};
      for (epsilon_index_type i = fa.first_order_.ibegin(); i != fa.first_order_.iend(); ++i)
        if (fa.first_order_[i] != nullptr)
          result.first_order_[i] = &(-arg * *fa.first_order_[i] * fa.common_denominator_);
      return result;
    }
  }

  symbolic::Expression const& expand() const
  {
    symbolic::Expression const* result = &zero_order_;
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result = &(*result + *first_order_[i] * *manager_->epsilons()[i]);
    return *result / common_denominator_;
  }

  // Add an expression that can contain ε's.
  // (z1 + f1)/d1 + (z2 + f2)/d2 = (z1 d2 + f1 d2)/(d1 d2) + (z2 d1 + f2 d1)/(d1 d2) =
  //                             = ((z1 d2 + z2 d1) + (f1 d2 + f2 d1))/(d1 d2)
  Expression operator+(Expression const& fa) const
  {
    Expression result{zero_order_ * fa.common_denominator_ + fa.zero_order_ * common_denominator_,
      common_denominator_ * fa.common_denominator_, manager_};
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

  Expression operator-(Expression const& fa) const
  {
    return *this + (-fa);
  }

  // Multiply with an expression that can contain ε's.
  // (z1 + f1)/d1 * (z2 + f2)/d2 = (z1 z2, z1 f2 + f1 z2)/(d1 d2)     (neglecting f1 * f2).
  Expression operator*(Expression const& fa) const
  {
    Expression result{zero_order_ * fa.zero_order_, common_denominator_ * fa.common_denominator_, manager_};
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
  //
  // If f2 only depends on a single ε, then it will be equal to z2 * ε2, so that the result will be
  //  (z1 z2 d2 - z1 z2⋅ε2 d2 + f1 z2 d2) / (z2^2 d1) = (z1 d2 - z1 d2⋅ε2 + f1 d2) / (z2 d1)
  //
  Expression operator/(Expression const& fa) const
  {
    auto& z1 = zero_order_;
    auto& f1 = first_order_;
    auto& d1 = common_denominator_;
    auto& z2 = fa.zero_order_;
    auto& f2 = fa.first_order_;
    auto& d2 = fa.common_denominator_;

    // Check if f2 only contains a single non-zero pointer.
    bool f2_is_single_epsilon = std::ranges::count_if(f2, [](auto p) { return p != nullptr; }) == 1;
    if (!f2_is_single_epsilon)
    {
      // Construct ((z1 z2 d2) / (z2^2 d1))⋅(1 + ε_next).
      Expression result{z1 * z2 * d2, (z2^2) * d1, manager_};
      // Add (f1 z2 d2) / (z2^2 d1), by looping over all ε's of f1.
      for (epsilon_index_type i = f1.ibegin(); i != f1.iend(); ++i)
        if (f1[i] != nullptr)
          result.first_order_[i] = &(*f1[i] * z2 * d2);
      // Subtract (z1 f2 d2) / (z2^2 d1), by looping over all ε's of f2.
      for (epsilon_index_type i = f2.ibegin(); i != f2.iend(); ++i)
        if (f2[i] != nullptr)
        {
          if (result.first_order_[i] == nullptr)
            result.first_order_[i] = &(-z1 * *f2[i] * d2);
          else
            result.first_order_[i] = &(*result.first_order_[i] - z1 * *f2[i] * d2);
        }
      return result;
    }
    else
    {
      // Construct ((z1 d2) / (z2 d1))⋅(1 + ε_next).
      Expression result{z1 * d2, z2 * d1, manager_};
      // Add (f1 d2) / (z2 d1), by looping over all ε's of f1.
      for (epsilon_index_type i = f1.ibegin(); i != f1.iend(); ++i)
        if (f1[i] != nullptr)
          result.first_order_[i] = &(*f1[i] * d2);
      // Subtract (z1 d2⋅ε) / (z2 d1) for that single ε of f2 for which we known is equal to z2.
      for (epsilon_index_type i = f2.ibegin(); i != f2.iend(); ++i)
        if (f2[i] != nullptr)
        {
          ASSERT(f2[i]->equals(z2));
          if (result.first_order_[i] == nullptr)
            result.first_order_[i] = &(-z1 * d2);
          else
            result.first_order_[i] = &(*result.first_order_[i] - z1 * d2);
        }
      return result;
    }
  }

  // Negation.
  Expression operator-() const
  {
    Expression result{-zero_order_, common_denominator_, manager_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(-*first_order_[i]);
    return result;
  }

  symbolic::Expression const& zero_order() const { return zero_order_ / common_denominator_; }

  symbolic::Expression const& first_order(epsilon_index_type epsilon_index) const
  {
    return first_order_[epsilon_index] == nullptr ? symbolic::Constant::s_cached_zero : *first_order_[epsilon_index] / common_denominator_;
  }

  void print_on(std::ostream& os) const
  {
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
          os << prefix << '(' << *first_order_[i] << ")⋅ε" << utils::to_subscript_string(i.get_value());
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

  symbolic::Expression const& error_squared() const
  {
    symbolic::Expression const* sum_of_squared_first_order_elements = &symbolic::Constant::s_cached_zero;
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        sum_of_squared_first_order_elements = &(*sum_of_squared_first_order_elements + ((*first_order_[i])^2));
    return *sum_of_squared_first_order_elements / (common_denominator_^2);
  }

 private:
  void grow_first_order(epsilon_index_type epsilon_index)
  {
    if (first_order_.size() <= epsilon_index.get_value())
      first_order_.resize(epsilon_index.get_value() + 1);
  }
};

} // namespace floating_point_round_off_error
