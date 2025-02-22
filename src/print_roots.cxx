#include "sys.h"
#include "IEEE754.h"
#include "cairowindow/QuickGraph.h"
#include "math/CubicPolynomial.h"
#include "math/bracket_zero.h"
#include "mpreal/mpreal.h"
#include "utils/square.h"
#include "utils/split.h"
#include <array>
#include <string>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <concepts>
#include "utils/debug_ostream_operators.h"
#include "debug.h"

#define FOREACH_VARIABLE(X) \
  X(c0) \
  X(c1) \
  X(c2) \
  X(c3) \
  X(x)

#include "FloatingPointRoundOffError.h"

using namespace mpfr;

constexpr mp_prec_t precision_in_bits = 256;
constexpr double log10_of_two = 0.30102999566;
constexpr int mpprecision = precision_in_bits * log10_of_two;
double const rel_err_threshold = std::pow(10.0, -0.5 * mpprecision);

// For gdb.
void print(mpreal const& val)
{
  std::cout << val << std::endl;
}

void print(symbolic::Expression const& expression)
{
  std::cout << expression << std::endl;
}

void print(floating_point_round_off_error::Expression const& expression)
{
  std::cout << expression << std::endl;
}

// Update the concept definition to correctly check the streamability.
template<typename T>
concept StreamExtractable = requires (std::stringstream& ss, T& value) {
  { ss >> value } -> std::same_as<std::istream&>;
};

// Function to parse a string in the desired form and return an array of coefficients.
template<typename T>
std::array<T, 4> parse_string(std::string const& input)
{
  std::regex term_regex(R"(([+-])?\s*((0[.]|[1-9](\d*[.])?)\d*(e[+-]?\d+)?)(\s*(([a-z])(\^([0-9]+))?)?)?)");
  auto terms_begin = std::sregex_iterator(input.begin(), input.end(), term_regex);
  auto terms_end = std::sregex_iterator();

  std::array<T, 4> coefficients{};

  std::string symbol;
  bool first = true;
  for (std::sregex_iterator i = terms_begin; i != terms_end; ++i, first = false)
  {
    std::smatch const& match = *i;
    double sign = match[1].length() == 1 && match[1].str()[0] == '-' ? -1.0 : 1.0;
    if (sign == 1.0 && !first && (match[1].length() == 0 || match[1].str()[0] != '+'))
      THROW_ALERT("Expected leading + or - before \"[TERM]\".", AIArgs("[TERM]", match.str()));

    T value;
    if constexpr (StreamExtractable<T>)
    {
      std::stringstream ss;
      ss << match[2];
      ss >> value;
    }
    else
      value = T{match[2].str()};
    value *= sign;
    if (symbol.empty())
      symbol = match[8].str();
    else if (match[8].length() > 0)
    {
      if (symbol != match[8].str())
        THROW_ALERT("Mismatching symbols [PREV] != [CURRENT].", AIArgs("[PREV]", symbol)("[CURRENT]", match[8]));
    }
    int exponent = symbol.empty() ? 0 : (match[10].length() == 0 ? 1 : match[10].str()[0] - '0');

    if (0 <= exponent && exponent < coefficients.size())
    {
      if (coefficients[exponent] != 0.0)
        THROW_ALERT("Duplicated term with exponent [EXPONENT].", AIArgs("[EXPONENT]", exponent));
      coefficients[exponent] = value;
    }
  }

  return coefficients;
}

mpreal exact_root(double guess, std::array<mpreal, 4> coefficients)
{
  mpreal u = guess;

  auto evaluate = [&](mpreal u) -> mpreal { return coefficients[0] + (coefficients[1] + (coefficients[2] + coefficients[3] * u) * u) * u; };
  auto derivative = [&](mpreal u) -> mpreal { return coefficients[1] + (2 * coefficients[2] + 3 * coefficients[3] * u) * u; };
  auto half_second_derivative = [&](mpreal u) -> mpreal { return coefficients[2] + 3 * coefficients[3] * u; };

  mpreal prev_u = 0;
  mpreal step = 0;
  do
  {
    prev_u = u;
    mpreal Q_u = evaluate(u);
    mpreal dQ_u = derivative(u);
    mpreal half_ddQ_u = half_second_derivative(u);
    step = Q_u * dQ_u / (utils::square(dQ_u) - Q_u * half_ddQ_u);
    u -= step;
  }
  while (abs((prev_u - u) / u) > rel_err_threshold);

  return u;
}

mpreal next_offset(mpreal root, double offset, std::array<mpreal, 4> coefficients)
{
  DoutEntering(dc::notice, std::setprecision(60) << "next_offset(" << root << ", " << offset << ", coefficients)");
  mpreal offset_n{offset};
  mpreal x_n = root * (1 + offset_n);

  Dout(dc::notice, std::setprecision(60) << "offset_n = " << offset_n);
  Dout(dc::notice, std::setprecision(60) << "x_n = " << x_n);

  auto evaluate = [&](mpreal u) -> mpreal { return coefficients[0] + (coefficients[1] + (coefficients[2] + coefficients[3] * u) * u) * u; };
  auto derivative = [&](mpreal u) -> mpreal { return coefficients[1] + (2 * coefficients[2] + 3 * coefficients[3] * u) * u; };
  auto half_second_derivative = [&](mpreal u) -> mpreal { return coefficients[2] + 3 * coefficients[3] * u; };

  mpreal Q_x = evaluate(x_n);
  mpreal dQ_x = derivative(x_n);
  mpreal half_ddQ_x = half_second_derivative(x_n);
  mpreal step = Q_x * dQ_x / (utils::square(dQ_x) - Q_x * half_ddQ_x);
  mpreal x_n_plus_1 = x_n - step;

  Dout(dc::notice, std::setprecision(60) << "x_n_plus_1 = " << x_n_plus_1 << "; returning " << (x_n_plus_1 / root - 1));

  return x_n_plus_1 / root - 1;
}

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with mpprecision digits precision.
  std::streamsize const old_precision = std::cout.precision(mpprecision);
//  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  // Check if there is exactly one argument provided.
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <argument>" << std::endl;
    return 1;
  }

  // Store the argument in a std::string.
  std::string input = argv[1];

  // Print the string.
  std::cout << "The input argument is: \"" << input << "\"." << std::endl;

  using namespace floating_point_round_off_error;
  using FA = floating_point_round_off_error::Expression;

  ExpressionManager expression_manager;

  FOREACH_VARIABLE(DECLARE_VARIABLE);   // Declare all variables (c0, c1, c2, c3, x).

  // Construct the approximation cubic f_a(x_n).
  auto fa = c0 + (c1 + (c2 + c3 * x) * x) * x;

  // Construct the derivative dfa(x_n).
  auto three_c3x = 3 * c3 * x;
  auto dfa = c1 + (2 * c2 + three_c3x) * x;

  // Prepare the formula required to calculate the floating point round off errors of Halley's method.
  FA x_n_plus_1 = x - (fa * dfa) / (utils::square(dfa) - fa * (c2 + three_c3x));

  // Get the expression for xₙ₊₁.
  auto& expression_x_n_plus_1 = x_n_plus_1.expand();

  std::cout << "Halley's method:" << std::endl;
//  std::cout << "xₙ₊₁ = " << symbolic::UseUtf8{1} << x_n_plus_1 << std::endl;
//  std::cout << "Δxₙ₊₁ = ±√(" << symbolic::UseUtf8{1} << x_n_plus_1.error_squared() << ')' << std::endl;

  // Keep a reference to the symbol x handy.
  symbolic::Symbol const& symbol_x = expression_manager.get_variable(x.index_).symbol();
  // The error function.
  auto& error_squared = x_n_plus_1.error_squared();

  Dout(dc::notice, "error_squared = " << error_squared);

  symbolic::Sum const& sum = static_cast<symbolic::Sum const&>(error_squared);
  Dout(dc::notice, "sum = " << sum.arg1());
  symbolic::Product const& product = static_cast<symbolic::Product const&>(sum.arg1());
  Dout(dc::notice, "product = " << product.arg2());
  symbolic::Product const& product2 = static_cast<symbolic::Product const&>(product.arg2());
  Dout(dc::notice, "product2 = " << product2.arg1());
  Dout(dc::notice, "is power: " << product2.arg1().is_power());
  symbolic::Power const& power = static_cast<symbolic::Power const&>(product2.arg1());
  Dout(dc::notice, "power = " << power.arg1());
  auto& error_squared_times = error_squared * power.arg1();
  Dout(dc::notice, "error_squared_times = " << error_squared_times);

#if 0
  symbolic::Symbol const& symbol_r = symbolic::Symbol::realize("r");
  symbolic::Symbol const& symbol_h = symbolic::Symbol::realize("h");
  auto& r_plus_h = symbol_r + symbol_h;
  auto& error_squared_subs_x = error_squared.subs(symbol_x, r_plus_h);

  Dout(dc::notice, "error_squared_subs_x = " << error_squared_subs_x);
#endif

  try
  {
    auto double_coefficients = parse_string<double>(input);
    math::CubicPolynomial cubic(double_coefficients);
    std::cout << std::setprecision(18) << cubic << '\n';

    // Set the coefficients in the error formula.
    for (auto i = c0.index_; i <= c3.index_; ++i)
      expression_manager.get_variable(i).symbol() = cubic[i.get_value()];       // Assumes c0.index_ == 0, etc.

    double c0 = cubic[0];
    double c1 = cubic[1];
    double c2 = cubic[2];
    double c3 = cubic[3];
    {
      double d = c2 / c3;
      double D = 1.0 - 4.0 * c0 * c2 / utils::square(c1);
      double e = -0.5 * c1 / c2 * (1.0 - std::sqrt(D));
      double epsilon = e / d;
      bool special_case1 = D >= 0.0 && std::abs(epsilon) < 0.1;
      Dout(dc::notice, "special_case1: D = " << D << "; epsilon = " << epsilon << " -->" << (special_case1 ? "" : " NOT") << " special case 1.");
    }
    {
      double c2_div_c3 = c2 / c3;
      double c1_div_c3 = c1 / c3;
      double D = utils::square(c2_div_c3) - 3.0 * c1_div_c3;

      if (D >= 0.0)
      {
        double e_minus = (-c2_div_c3 - std::sqrt(D)) / 3.0;
        double e_plus = (-c2_div_c3 + std::sqrt(D)) / 3.0;
        double eval_e_minus = cubic(e_minus);
        double eval_e_plus = cubic(e_plus);
        double extreme = std::abs(eval_e_minus) < std::abs(eval_e_plus) ? e_minus : e_plus;

        double qc2 = 0.5 * cubic.second_derivative(extreme);
        math::QuadraticPolynomial qp(cubic(extreme), cubic.derivative(extreme), qc2);
        std::array<double, 2> roots;
        int n = qp.get_roots(roots);
        if (n > 1)
        {
          ASSERT(n == 2);
          Dout(dc::notice|continued_cf, "special_case2: qp = " << qp << " where root*c3/qc2 is " << (roots[0] * c3 / qc2) <<
              " and " << (roots[1] * c3 / qc2) << " -->");
          bool special_case2 = false;
          if (std::abs(roots[0] * c3 / qc2) < 0.1 && std::abs(roots[1] * c3 / qc2) < 0.1)
            special_case2 = true;
          Dout(dc::finish, (special_case2 ? "" : " NOT") << " special case 2.");
        }
        else
          Dout(dc::notice, "special_case2: qp = " << qp << " which has " << n << " roots! --> NOT special case 2.");
      }
      else
        Dout(dc::notice, "special_case2: D = " << D << " --> NOT special case 2.");
    }

    std::ostringstream title;
    title << cubic;
    double const y_max_error = 2.5;     // log10(max machine epsilons of relative error)
    cairowindow::Range const x_range{0.5, 1.0};
    cairowindow::Range const y_range{-0.01, 0.01};
    double const error_scale = y_range.max() / y_max_error;
    cairowindow::QuickGraph qg(title.str(), "x", "y", x_range, y_range);
    qg.add_function([&](double x){ return cubic(x); }, cairowindow::color::teal);
    qg.add_function([&](double x){
        symbol_x = x;
        return std::log10(/*numeric_limits<double>::epsilon() * */
            0.866675 * std::sqrt(error_squared.evaluate()) / std::abs(expression_x_n_plus_1.evaluate())) * error_scale;
     }, cairowindow::color::red);
    qg.add_function([&](double x){
        double result = -std::log10(std::abs(utils::square(cubic.derivative(x)) - 0.5 * cubic.second_derivative(x) * cubic(x))) * 0.002;
        return result;
        }, cairowindow::color::blue);

    std::array<mpreal, 4> coefficients = parse_string<mpreal>(input);
    for (int i = 0; i < coefficients.size(); ++i)
      std::cout << "i=" << i << " : " << coefficients[i] << std::endl;

    std::array<double, 3> roots;
    int n = cubic.get_roots(roots);
    std::cout << "This cubic has " << n << " root" << (n != 1 ? "s" : "") << ":\n";
    mpreal root;
    double derivative;
    for (int i = 0; i < n; ++i)
    {
      mpreal root_i = exact_root(roots[i], coefficients);
      double derivative_i = cubic.derivative(root_i.toDouble());
      std::cout << "  " << root_i << " (as double: " << root_i.toDouble() << "; with derivative " << derivative_i << ")" << '\n';
      if (i == 0 || std::abs(derivative_i) < std::abs(derivative))
      {
        root = root_i;
        derivative = derivative_i;
      }
    }

    using namespace cairowindow;
    using Window = cairowindow::Window;

    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});

    double z = std::log10(std::abs(0.866675 / cubic.derivative(root.toDouble()))) * error_scale;
    Dout(dc::notice, "z = " << z);
    qg.add_line({{0.0, z}, Direction::right}, solid_line_style({.line_color = color::blue}));

    // Draw a vertical line where this root is.
    qg.add_line({{root.toDouble(), 0.0}, Direction::up}, solid_line_style({.line_color = color::lime}));

    // Draw the derivative.
    qg.add_function([&](double x){ return cubic.derivative(x); });

    // Create a window.
    Window window("log/log graph " + title.str(), 1200, 900);
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      Debug(NAMESPACE_DEBUG::init_thread("event_loop"));
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Log/log of cubic", {},
        "log₁₀(x - root)", {},
        "log₁₀(y)", {});
//    plot.set_xrange({-17.0, -12.0});
//    plot.set_yrange({-17.0, -12.0});
    plot.set_xrange({-16.0, -12.0});
    plot.set_yrange({-16.0, -12.0});
    plot.add_to(background_layer, false);

    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});

    auto slider_xn = plot.create_slider(second_layer, {1028, 283, 7, 400}, -2.0, plot.xrange().min(), plot.xrange().max());
    auto slider_xn_label = plot.create_text(second_layer, slider_style, Pixel{1028, 683}, "xₙ");

    while(1)
    {
      window.set_send_expose_events(false);

      // Draw the cubic.
      plot::BezierFitter plot_curve1;
      plot_curve1.solve(
          [&](double log10_offset) -> cairowindow::Point {
            double offset = std::pow(10.0, log10_offset);
            double x = root.toDouble() * (1.0 + offset);
            double y = cubic(x);
            return {log10_offset, std::log10(std::abs(y))};
          },
          plot.viewport());
      plot.add_bezier_fitter(second_layer, solid_line_style, plot_curve1);
      plot::BezierFitter plot_curve2;
      plot_curve2.solve(
          [&](double log10_offset) -> cairowindow::Point {
            double offset = std::pow(10.0, log10_offset);
            double x = root.toDouble() * (1.0 - offset);
            double y = cubic(x);
            return {log10_offset, std::log10(std::abs(y))};
          },
          plot.viewport());
      plot.add_bezier_fitter(second_layer, solid_line_style({.line_color = color::red}), plot_curve2);

      // Pick some value for x_n.
      double offset_n = std::pow(10.0, slider_xn.value());
      double x_n = root.toDouble() * (1.0 + offset_n);

      symbol_x = x_n;
      double relative_error =
        0.866675 * std::numeric_limits<double>::epsilon() * std::sqrt(error_squared.evaluate()) / std::abs(expression_x_n_plus_1.evaluate());
      plot::Line plot_error{{std::log10(relative_error), 0.0}, Direction::up};
      plot.add_line(second_layer, solid_line_style({.line_color = color::darkred, .line_width = 2.0}), plot_error);

      double offset_n_plus_1 = next_offset(root, offset_n, coefficients).toDouble();
      double real_root = root.toDouble();
      // offset_n_plus_1 = x_minus_x_n - root.toDouble() -->
      // x_minus_x_n = offset_n_plus_1 + root.toDouble()

      // Draw a vertical line at x_n.
      plot::Line plot_vertical_line{{slider_xn.value(), 0.0}, Direction::up};
      plot.add_line(second_layer, solid_line_style, plot_vertical_line);

      plot::Line plot_initial_root_found{{std::log10(std::abs(x_n / real_root - 1.0)), 0.0}, Direction::up};
      plot.add_line(second_layer, solid_line_style({.line_color = color::lime}), plot_initial_root_found);

      double x_n_plus_1 = root.toDouble() * (1.0 + offset_n_plus_1);
      plot::Line plot_root_found{{std::log10(std::abs(offset_n_plus_1)), 0.0}, Direction::up};
      plot.add_line(second_layer, solid_line_style({.line_color = color::red}), plot_root_found);

      // Draw the Halley "approach" curve (cuts the x-axis at x_{n+1}).
      plot::BezierFitter plot_curve3;
      plot_curve3.solve(
          [&](double log10_offset) -> cairowindow::Point {
            double offset = std::pow(10.0, log10_offset);
            double x = root.toDouble() * (1.0 + offset);                // x   = root * (1 + offset)
                                                                        // x_n = root * (1 + offset_n)
            double x_minus_x_n = root.toDouble() * (offset - offset_n); // x - x_n = root * (offset - offset_n)
            // Halley:
            double gn = cubic(x_n) + utils::square(cubic.derivative(x_n)) * x_minus_x_n /
                (cubic.derivative(x_n) - cubic.half_second_derivative(x_n) * x_minus_x_n);
            // Newton-Raphson:
            //double gn = cubic(x_n) + cubic.derivative(x_n) * x_minus_x_n;
            return {log10_offset, std::log10(std::abs(gn))};
          },
          plot.viewport());
      plot.add_bezier_fitter(second_layer, solid_line_style({.line_color = color::blue}), plot_curve3);
      plot::BezierFitter plot_curve4;
      plot_curve4.solve(
          [&](double log10_offset) -> cairowindow::Point {
            double offset = -std::pow(10.0, log10_offset);
            double x = root.toDouble() * (1.0 + offset);                // x   = root * (1 + offset)
                                                                        // x_n = root * (1 + offset_n)
            double x_minus_x_n = root.toDouble() * (offset - offset_n); // x - x_n = root * (offset - offset_n)
            // Halley:
            double gn = cubic(x_n) + utils::square(cubic.derivative(x_n)) * x_minus_x_n /
                (cubic.derivative(x_n) - cubic.half_second_derivative(x_n) * x_minus_x_n);
            // Newton-Raphson:
            //double gn = cubic(x_n) + cubic.derivative(x_n) * x_minus_x_n;
            return {log10_offset, std::log10(std::abs(gn))};
          },
          plot.viewport());
      plot.add_bezier_fitter(second_layer, solid_line_style({.line_color = color::purple}), plot_curve4);

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.
    }
    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }
}
