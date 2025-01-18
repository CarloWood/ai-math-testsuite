#pragma once

#include "utils/has_print_on.h"
#include "utils/Array.h"
#include "cairowindow/symbolic/symbolic.h"
#include <algorithm>

using utils::has_print_on::operator<<;

template<int number_of_epsilons>
class FloatingPointRoundOffError
{
 public:
  using epsilon_index_type = utils::ArrayIndex<FloatingPointRoundOffError>;

 private:
  symbolic::Expression const& zero_order_;                                                              // Term without ε's.
  utils::Array<symbolic::Expression const*, number_of_epsilons, epsilon_index_type> first_order_{};     // Terms with single ε factor(s).
  symbolic::Expression const& common_denominator_;                                                      // Denominator without ε's.

 public:
  // Construct a FirstOrder object where the zero order is an integer.
  FloatingPointRoundOffError(int zero_order) :
    zero_order_(symbolic::Constant::realize(zero_order)), common_denominator_(symbolic::Constant::s_cached_one) { }

  // Construct a FirstOrder object from an expressions that is to be multiplied by (1 + ε).
  FloatingPointRoundOffError(symbolic::Expression const& expression, epsilon_index_type epsilon_index,
      symbolic::Expression const& common_denominator = symbolic::Constant::s_cached_one) :
    zero_order_(expression), common_denominator_(common_denominator)
  {
    first_order_[epsilon_index] = &expression;
  }

  // Construct a FirstOrder object from an FloatingPointRoundOffError object that is to be multiplied by (1 + ε).
  FloatingPointRoundOffError(FloatingPointRoundOffError const& fa, epsilon_index_type epsilon_index) :
    zero_order_(fa.zero_order_), first_order_(fa.first_order_), common_denominator_(fa.common_denominator_)
  {
    // Don't use an already used epsilon_index.
    ASSERT(first_order_[epsilon_index] == nullptr);
    first_order_[epsilon_index] = &fa.zero_order_;
  }

 private:
  FloatingPointRoundOffError(symbolic::Expression const& zero_order, symbolic::Expression const& common_denominator) :
    zero_order_(zero_order), common_denominator_(common_denominator) { }

  FloatingPointRoundOffError(symbolic::Expression const& zero_order,
      utils::Array<symbolic::Expression const*, number_of_epsilons, epsilon_index_type> const& first_order,
      symbolic::Expression const& common_denominator) :
    zero_order_(zero_order), first_order_(first_order), common_denominator_(common_denominator) { }

 public:
  // Multiply with an expression that does NOT contain ε's.
  FloatingPointRoundOffError operator*(symbolic::Expression const& arg) const
  {
    FloatingPointRoundOffError result{zero_order_ * arg, common_denominator_};
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
  FloatingPointRoundOffError operator+(symbolic::Expression const& arg) const
  {
    return {zero_order_ + arg * common_denominator_, first_order_, common_denominator_};
  }

  // Subtract an expression that does NOT contain ε's.
  FloatingPointRoundOffError operator-(symbolic::Expression const& arg) const
  {
    return {zero_order_ - arg * common_denominator_, first_order_, common_denominator_};
  }

  // Multiply with an expression that can contain ε's.
  // (z1 + f1)/d1 * (z2 + f2)/d2 = (z1 z2, z1 f2 + f1 z2)/(d1 d2)     (neglecting f1 * f2).
  FloatingPointRoundOffError operator*(FloatingPointRoundOffError const& fa) const
  {
    FloatingPointRoundOffError result{zero_order_ * fa.zero_order_, common_denominator_ * fa.common_denominator_};
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
  FloatingPointRoundOffError operator/(FloatingPointRoundOffError const& fa) const
  {
    // Check if fa.first_order_ only contains a single non-zero pointer.
    bool f2_is_single_epsilon = std::ranges::count_if(fa.first_order_, [](auto p) { return p != nullptr; }) == 1;
    if (f2_is_single_epsilon)
    {
      FloatingPointRoundOffError result{zero_order_ * fa.common_denominator_, fa.zero_order_ * common_denominator_};
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
      FloatingPointRoundOffError result{zero_order_ * fa.zero_order_ * fa.common_denominator_, (fa.zero_order_^2) * common_denominator_};
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
  FloatingPointRoundOffError operator+(FloatingPointRoundOffError const& fa) const
  {
    FloatingPointRoundOffError result{zero_order_ * fa.common_denominator_ + fa.zero_order_ * common_denominator_,
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

  FloatingPointRoundOffError operator-(FloatingPointRoundOffError const& fa) const
  {
    return *this + (-fa);
  }

  // The operations are commutative.
  friend FloatingPointRoundOffError operator*(symbolic::Expression const& arg, FloatingPointRoundOffError const& fa)
  {
    return fa * arg;
  }

  // Negation.
  FloatingPointRoundOffError operator-() const
  {
    FloatingPointRoundOffError result{-zero_order_, common_denominator_};
    for (epsilon_index_type i = first_order_.ibegin(); i != first_order_.iend(); ++i)
      if (first_order_[i] != nullptr)
        result.first_order_[i] = &(-*first_order_[i]);
    return result;
  }

  friend FloatingPointRoundOffError operator+(symbolic::Expression const& arg, FloatingPointRoundOffError const& fa)
  {
    return fa + arg;
  }

  friend FloatingPointRoundOffError operator-(symbolic::Expression const& arg, FloatingPointRoundOffError const& fa)
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
