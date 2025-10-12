#include "sys.h"
#include "math/Direction.h"
#include "math/Line.h"
#include "math/LinePiece.h"
#include "math/Point.h"
#include "math/Vector.h"

#include <array>
#include <cmath>
#include <concepts>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "debug.h"

namespace {

constexpr double kTolerance = 1e-9;

void require(bool condition, std::string const& test_description, std::string const& message)
{
  if (!condition)
    throw std::runtime_error(test_description + message);
}

void require_near(double actual, double expected, double tolerance, std::string const& test_description, std::string const& message)
{
  if (std::isnan(expected))
  {
    require(std::isnan(actual), test_description, message + ": expected NaN");
    return;
  }

  if (std::isinf(expected))
  {
    require(std::isinf(actual) && std::signbit(actual) == std::signbit(expected), test_description, message + ": expected infinity");
    return;
  }

  if (std::abs(actual - expected) > tolerance)
  {
    std::ostringstream oss;
    oss << test_description << " " << message << ": expected " << expected << ", got " << actual;
    throw std::runtime_error(oss.str());
  }
}

void require_equal(double actual, double expected, std::string const& test_description, std::string const& message)
{
  require_near(actual, expected, 0.0, test_description, message);
}

// Traits that expose dimension and scalar type.
template<typename> struct object_traits;

template<int N, typename T>
struct object_traits<math::Point<N, T>> {
  static constexpr std::size_t dim = N;
  using scalar_type = T;
};

template<int N, typename T>
struct object_traits<math::Direction<N, T>> {
  static constexpr std::size_t dim = N;
  using scalar_type = T;
};

template<int N, typename T>
struct object_traits<math::Vector<N, T>> {
  static constexpr std::size_t dim = N;
  using scalar_type = T;
};

template<typename Obj>
concept ConceptObjectDouble =
  requires { object_traits<Obj>::dim; typename object_traits<Obj>::scalar_type; } &&
  std::same_as<typename object_traits<Obj>::scalar_type, double>;

template<ConceptObjectDouble Obj>
void require_components_near(Obj const& object,
                             std::array<double, object_traits<Obj>::dim> const& expected,
                             double tolerance,
                             std::string const& test_description,
                             std::string const& context)
{
  constexpr int N = static_cast<int>(object_traits<Obj>::dim);
  for (int i = 0; i < N; ++i)
  {
    std::ostringstream label;
    label << context << " component " << i;
    require_near(object[i], expected[i], tolerance, test_description, label.str());
  }
}

struct ConvertibleToDouble
{
  double value_;
  operator double() const { return value_; }
};

template<int N>
void test_point(std::string test_description)
{
  using Point = math::Point<N>;
  using Direction = math::Direction<N>;
  using Vector = math::Vector<N>;

  double const x = 3.9851;
  double const y = -2.7795;
  double const z = 8.3333;

  ConvertibleToDouble const cd_x{x};

  //---------------------------------------------------------------------------
  // Constructor and accessors.
  Point p1;
  {
    if constexpr (N == 2)
    {
      Point tp1{cd_x, y};
      p1 = tp1;
    }
    else
    {
      Point tp1{cd_x, y, z};
      p1 = tp1;
    }
  }
  Point const& cp = p1;

  auto& e = p1.eigen();
  // Check that the type of `e` is eigen_type&.
  static_assert(std::same_as<decltype(e), std::remove_cvref_t<decltype(e)>&>);
  // Check that the type of `ce` is eigen_type const&.
  auto& ce = cp.eigen();
  static_assert(std::same_as<decltype(ce), std::remove_reference_t<decltype(ce)> const&>);

  require_equal(e[0], x, test_description, "eigen x");
  require_equal(ce[1], y, test_description, "eigen y");
  if constexpr (N == 3)
    require_equal(ce[2], z, test_description, "eigen z");

  require_equal(p1[0], x, test_description, "operator[] x");
  require_equal(p1[1], y, test_description, "operator[] y");
  if constexpr (N == 3)
    require_equal(p1[2], z, test_description, "operator[] z");

  require_equal(cp.x(), x, test_description, "x() const");
  require_equal(cp.y(), y, test_description, "y() const");
  if constexpr (N == 3)
    require_equal(cp.z(), z, test_description, "z() const");

  double const x1 = -0.0010982;
  double const y1 = 100.001;
  double const z1 = -2.25;

  //---------------------------------------------------------------------------
  // Copy constructor.
  Point const p2(cp);

  double const x2 = x;
  double const y2 = y;
  double const z2 = z;

  // Construct a Point that is the sum of the coordinates of p1 and p2.

  //---------------------------------------------------------------------------
  // Setters
  p1.x() = x1;
  p1.y() = y1;
  if constexpr (N == 3)
    p1.z() = z1;

  if constexpr (N == 2)
    require_components_near(p1, {x1, y1}, kTolerance, test_description, "setters");
  else
    require_components_near(p1, {x1, y1, z1}, kTolerance, test_description, "setters");

  //---------------------------------------------------------------------------
  // operators

  // Use p2 as a direction.
  Direction const d2{p2};
  double dx = p2.x();
  double dy = p2.y();
  double dz = N == 3 ? p2[2] : 0.0;
  double const len = std::sqrt(dx * dx + dy * dy + dz * dz);
  dx /= len;
  dy /= len;
  dz /= len;

  // Add Direction to a Point.
  Point translated = p1 + d2;
  if constexpr (N == 2)
    require_components_near(translated, {x1 + dx, y1 + dy}, kTolerance, test_description, "operator+(direction)");
  else
    require_components_near(translated, {x1 + dx, y1 + dy, z1 + dz}, kTolerance, test_description, "operator+(direction)");

  // Use p2 as a Vector.
  Vector const shift{p2};

  translated = p1 + shift;
  if constexpr (N == 2)
    require_components_near(translated, {x1 + x2, y1 + y2}, kTolerance, test_description, "operator+(vector)");
  else
    require_components_near(translated, {x1 + x2, y1 + y2, z1 + z2}, kTolerance, test_description, "operator+(vector)");

  translated = p1 - shift;
  if constexpr (N == 2)
    require_components_near(translated, {x1 - x2, y1 - y2}, kTolerance, test_description, "operator-(vector)");
  else
    require_components_near(translated, {x1 - x2, y1 - y2, z1 - z2}, kTolerance, test_description, "operator-(vector)");

  Point mutable_point = p1;
  mutable_point += d2;
  if constexpr (N == 2)
    require_components_near(mutable_point, {x1 + dx, y1 + dy}, kTolerance, test_description, "operator+= direction");
  else
    require_components_near(mutable_point, {x1 + dx, y1 + dy, z1 + dz}, kTolerance, test_description, "operator+= direction");
  mutable_point += shift;
  if constexpr (N == 2)
    require_components_near(mutable_point, {x1 + dx + x2, y1 + dy + y2}, kTolerance, test_description, "operator+= vector");
  else
    require_components_near(mutable_point, {x1 + dx + x2, y1 + dy + y2, z1 + dz + z2}, kTolerance, test_description, "operator+= vector");
  mutable_point -= shift;
  if constexpr (N == 2)
    require_components_near(mutable_point, {x1 + dx, y1 + dy}, kTolerance, test_description, "operator-= vector");
  else
    require_components_near(mutable_point, {x1 + dx, y1 + dy, z1 + dz}, kTolerance, test_description, "operator-= vector");

  Vector difference = p1 - p2;
  if constexpr (N == 2)
    require_components_near(difference, {x1 - x2, y1 - y2}, kTolerance, test_description, "difference operator");
  else
    require_components_near(difference, {x1 - x2, y1 - y2, z1 - z2}, kTolerance, test_description, "difference operator");

  require(!(p1 != p1), test_description, "equality");
  require(p1 != p2, test_description, "inequality");
}

void test_point_2d()
{
  test_point<2>("Point<2>");
}

void test_point_3d()
{
  test_point<3>("Point<3>");
}

void test_vector_2d()
{
  char const* test_description = "Vector<2>";

  using Vector = math::Vector<2>;
  using Direction = math::Direction<2>;
  using Point = math::Point<2>;
  using LinePiece = math::LinePiece<2>;

  Vector v{3.0, 4.0};
  require_components_near(v, {3.0, 4.0}, kTolerance, test_description, "initializer");

  Direction direction{Point{0.0, 0.0}, Point{1.0, 1.0}};
  Vector vdir{direction, std::sqrt(2.0)};
  require_components_near(vdir, {1.0, 1.0}, kTolerance, test_description, "from direction");

  Vector from_points{Point{1.0, 2.0}, Point{4.0, 6.0}};
  require_components_near(from_points, {3.0, 4.0}, kTolerance, test_description, "from two points");

  Vector from_point{Point{3.0, 4.0}};
  require_components_near(from_point, {3.0, 4.0}, kTolerance, test_description, "from single point");

  LinePiece segment{Point{0.0, 0.0}, Point{-3.0, -4.0}};
  Vector from_segment{segment};
  require_components_near(from_segment, {-3.0, -4.0}, kTolerance, test_description, "from line piece");

  Vector::eigen_type eigen;
  eigen << -6.0, 8.0;
  Vector from_eigen{eigen};
  require_components_near(from_eigen, {-6.0, 8.0}, kTolerance, test_description, "from eigen");

  v[0] = 5.0;
  v[1] = -1.0;
  require_components_near(v, {5.0, -1.0}, kTolerance, test_description, "operator[] write");
  require_near(static_cast<Vector const&>(v)[0], 5.0, kTolerance, test_description, "operator[] const");

  require_near(v.x(), 5.0, kTolerance, test_description, "x()");
  require_near(v.y(), -1.0, kTolerance, test_description, "y()");

  require_near(v.dot(Vector{2.0, 3.0}), 7.0, kTolerance, test_description, "dot");
  require_near(v.cross(Vector{2.0, 3.0}), 17.0, kTolerance, test_description, "cross scalar");

  Direction dir = v.direction();
  require_near(dir.x(), v.x() / v.norm(), kTolerance, test_description, "direction x");
  require_near(dir.y(), v.y() / v.norm(), kTolerance, test_description, "direction y");

  require_near(v.norm(), std::sqrt(26.0), kTolerance, test_description, "norm");
  require_near(v.norm_squared(), 26.0, kTolerance, test_description, "norm_squared");

  Point as_point = v.as_point();
  require_components_near(as_point, {5.0, -1.0}, kTolerance, test_description, "as_point");

  Vector nan_vector{std::numeric_limits<double>::quiet_NaN(), 0.0};
  require(nan_vector.isnan(), test_description, "isnan");
  Vector finite_vector{1.0, 1.0};
  require(finite_vector.isfinite(), test_description, "isfinite true");
  Vector infinite_vector{std::numeric_limits<double>::infinity(), 1.0};
  require(!infinite_vector.isfinite(), test_description, "isfinite false");

  require_components_near(v.rotate_90_degrees(), {1.0, 5.0}, kTolerance, test_description, "rotate 90");
  require_components_near(v.rotate_180_degrees(), {-5.0, 1.0}, kTolerance, test_description, "rotate 180");
  require_components_near(v.rotate_270_degrees(), {-1.0, -5.0}, kTolerance, test_description, "rotate 270");

  Vector w{1.0, 2.0};
  Vector u = v;
  u += w;
  require_components_near(u, {6.0, 1.0}, kTolerance, test_description, "operator+=");
  u -= w;
  require_components_near(u, {5.0, -1.0}, kTolerance, test_description, "operator-=");
  u *= 2.0;
  require_components_near(u, {10.0, -2.0}, kTolerance, test_description, "operator*=");
  u /= 4.0;
  require_components_near(u, {2.5, -0.5}, kTolerance, test_description, "operator/=");

  Vector divided = v / 2.0;
  require_components_near(divided, {2.5, -0.5}, kTolerance, test_description, "operator/(scalar)");

  require_components_near(v + w, {6.0, 1.0}, kTolerance, test_description, "operator+");
  require_components_near(v - w, {4.0, -3.0}, kTolerance, test_description, "operator-");
  require_components_near(2.0 * v, {10.0, -2.0}, kTolerance, test_description, "scalar multiplication");

  Point base{1.0, 1.0};
  require_components_near(base + v, {6.0, 0.0}, kTolerance, test_description, "Point + Vector");
  require_components_near(base - v, {-4.0, 2.0}, kTolerance, test_description, "Point - Vector");
}

void test_vector_3d()
{
  char const* test_description = "Vector<3>";

  using Vector = math::Vector<3>;
  using Direction = math::Direction<3>;
  using Point = math::Point<3>;
  using LinePiece = math::LinePiece<3>;

  Vector v{1.0, 2.0, 3.0};
  require_components_near(v, {1.0, 2.0, 3.0}, kTolerance, test_description, "initializer");

  Vector from_points{Point{1.0, 2.0, 3.0}, Point{4.0, 6.0, 9.0}};
  require_components_near(from_points, {3.0, 4.0, 6.0}, kTolerance, test_description, "from two points");

  Vector from_point{Point{3.0, 4.0, 5.0}};
  require_components_near(from_point, {3.0, 4.0, 5.0}, kTolerance, test_description, "from single point");

  Direction direction{Point{0.0, 0.0, 0.0}, Point{0.0, 3.0, 4.0}};
  Vector from_direction{direction, 5.0};
  require_components_near(from_direction, {0.0, 3.0, 4.0}, kTolerance, test_description, "from direction");

  LinePiece segment{Point{0.0, 0.0, 0.0}, Point{-1.0, -2.0, -3.0}};
  Vector from_segment{segment};
  require_components_near(from_segment, {-1.0, -2.0, -3.0}, kTolerance, test_description, "from line piece");

  Vector::eigen_type eigen;
  eigen << 2.0, -1.0, 4.0;
  Vector from_eigen{eigen};
  require_components_near(from_eigen, {2.0, -1.0, 4.0}, kTolerance, test_description, "from eigen");

  require_near(v.dot(Vector{4.0, 5.0, 6.0}), 32.0, kTolerance, test_description, "dot");

  Vector cross = v.cross(Vector{-1.0, 0.0, 2.0});
  require_components_near(cross, {4.0, -5.0, 2.0}, kTolerance, test_description, "cross product");

  Direction dir = v.direction();
  double norm = v.norm();
  require_near(dir.x(), v.x() / norm, kTolerance, test_description, "direction x");
  require_near(dir.y(), v.y() / norm, kTolerance, test_description, "direction y");
  require_near(dir.z(), v.z() / norm, kTolerance, test_description, "direction z");

  require_near(v.norm(), std::sqrt(14.0), kTolerance, test_description, "norm");
  require_near(v.norm_squared(), 14.0, kTolerance, test_description, "norm_squared");

  Point as_point = v.as_point();
  require_components_near(as_point, {1.0, 2.0, 3.0}, kTolerance, test_description, "as_point");

  Vector nan_vector{0.0, std::numeric_limits<double>::quiet_NaN(), 0.0};
  require(nan_vector.isnan(), test_description, "isnan");

  Vector infinite_vector{0.0, std::numeric_limits<double>::infinity(), 0.0};
  require(!infinite_vector.isfinite(), test_description, "isfinite false");

  Vector u = v;
  u += Vector{1.0, 0.0, -1.0};
  require_components_near(u, {2.0, 2.0, 2.0}, kTolerance, test_description, "operator+=");
  u -= Vector{1.0, 1.0, 1.0};
  require_components_near(u, {1.0, 1.0, 1.0}, kTolerance, test_description, "operator-=");
  u *= 3.0;
  require_components_near(u, {3.0, 3.0, 3.0}, kTolerance, test_description, "operator*=");
  u /= 3.0;
  require_components_near(u, {1.0, 1.0, 1.0}, kTolerance, test_description, "operator/=");

  Vector divided = v / 2.0;
  require_components_near(divided, {0.5, 1.0, 1.5}, kTolerance, test_description, "operator/(scalar)");

  require_components_near(v + Vector{1.0, 1.0, 1.0}, {2.0, 3.0, 4.0}, kTolerance, test_description, "operator+");
  require_components_near(v - Vector{1.0, 1.0, 1.0}, {0.0, 1.0, 2.0}, kTolerance, test_description, "operator-");

  Point base{1.0, 2.0, 3.0};
  require_components_near(base + v, {2.0, 4.0, 6.0}, kTolerance, test_description, "Point<3> + Vector<3>");
  require_components_near(base - v, {0.0, 0.0, 0.0}, kTolerance, test_description, "Point<3> - Vector<3>");
}

void test_direction_2d()
{
  char const* test_description = "Direction<2>";

  using Direction = math::Direction<2>;
  using Line = math::Line<2>;
  using LinePiece = math::LinePiece<2>;
  using Point = math::Point<2>;

  Direction angle{std::numbers::pi_v<double> / 4.0};
  double expected = std::sqrt(0.5);
  require_near(angle.x(), expected, kTolerance, test_description, "from angle x");
  require_near(angle.y(), expected, kTolerance, test_description, "from angle y");

  Direction from_points{Point{0.0, 0.0}, Point{0.0, 5.0}};
  require_components_near(from_points, {0.0, 1.0}, kTolerance, test_description, "from points");

  Direction from_point_only{Point{0.0, -3.0}};
  require_components_near(from_point_only, {0.0, -1.0}, kTolerance, test_description, "from single point");

  LinePiece segment{Point{0.0, 0.0}, Point{1.0, 0.0}};
  Direction from_segment{segment};
  require_components_near(from_segment, {1.0, 0.0}, kTolerance, test_description, "from segment");

  Line line{Point{0.0, 0.0}, angle};
  Direction from_line{line};
  require_components_near(from_line, {angle.x(), angle.y()}, kTolerance, test_description, "from line");

  require_near(angle.dot(Direction::right), angle.x(), kTolerance, test_description, "dot");

  require_near(from_points.as_angle(), std::numbers::pi_v<double> / 2.0, kTolerance, test_description, "as_angle");

  Direction normal = angle.normal();
  require_components_near(normal, {-expected, expected}, kTolerance, test_description, "normal");

  Direction inverse = angle.inverse();
  require_components_near(inverse, {-angle.x(), -angle.y()}, kTolerance, test_description, "inverse");

  Direction normal_inverse = angle.normal_inverse();
  require_components_near(normal_inverse, {expected, -expected}, kTolerance, test_description, "normal_inverse");

  require_components_near(Direction::up, {0.0, 1.0}, kTolerance, test_description, "up");
  require_components_near(Direction::down, {0.0, -1.0}, kTolerance, test_description, "down");
  require_components_near(Direction::left, {-1.0, 0.0}, kTolerance, test_description, "left");
  require_components_near(Direction::right, {1.0, 0.0}, kTolerance, test_description, "right");
}

void test_direction_3d()
{
  char const* test_description = "Direction<3>";

  using Direction = math::Direction<3>;
  using Line = math::Line<3>;
  using Point = math::Point<3>;

  Direction from_points{Point{0.0, 0.0, 0.0}, Point{1.0, 2.0, 2.0}};
  double norm = std::sqrt(9.0);
  require_components_near(from_points, {1.0 / norm, 2.0 / norm, 2.0 / norm}, kTolerance, test_description, "from points");

  Direction from_point_only{Point{0.0, 0.0, -3.0}};
  require_components_near(from_point_only, {0.0, 0.0, -1.0}, kTolerance, test_description, "from origin to point");

  Line line{Point{1.0, 1.0, 1.0}, from_points};
  Direction from_line{line};
  require_components_near(from_line, {1.0 / norm, 2.0 / norm, 2.0 / norm}, kTolerance, test_description, "from line");

  require_near(from_points.dot(from_line), 1.0, kTolerance, test_description, "dot");

  require_near(from_points.as_angle(), std::atan2(2.0 / norm, 1.0 / norm), kTolerance, test_description, "as_angle");
}

void test_line_piece_2d()
{
  char const* test_description = "LinePiece<2>";

  using LinePiece = math::LinePiece<2>;
  using Point = math::Point<2>;

  Point a{0.0, 0.0};
  Point b{3.0, 4.0};
  LinePiece segment{a, b};

  require_components_near(segment.from(), {0.0, 0.0}, kTolerance, test_description, "from");
  require_components_near(segment.to(), {3.0, 4.0}, kTolerance, test_description, "to");
  require_near(segment.norm(), 5.0, kTolerance, test_description, "length");

  auto direction = segment.direction();
  require_components_near(direction, {0.6, 0.8}, kTolerance, test_description, "direction");
}

void test_line_piece_3d()
{
  char const* test_description = "LinePiece<3>";

  using LinePiece = math::LinePiece<3>;
  using Point = math::Point<3>;

  Point a{0.0, 0.0, 0.0};
  Point b{1.0, 2.0, 2.0};
  LinePiece segment{a, b};

  require_near(segment.norm(), 3.0, kTolerance, test_description, "length");

  auto direction = segment.direction();
  require_components_near(direction, {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}, kTolerance, test_description, "direction");
}

void test_line_2d()
{
  char const* test_description = "Line<2>";

  using Direction = math::Direction<2>;
  using Line = math::Line<2>;
  using Point = math::Point<2>;

  Direction slope{Point{0.0, 0.0}, Point{1.0, 2.0}};
  Point point{1.0, 1.0};
  Line line{point, slope};

  require_components_near(line.point(), {1.0, 1.0}, kTolerance, test_description, "point");
  require_components_near(line.direction(), {slope.x(), slope.y()}, kTolerance, test_description, "direction");

  Direction const& as_direction = line;
  require(&as_direction == &line.direction(), test_description, "conversion to direction");

  Line other{Point{0.0, 2.0}, Direction::down};
  Point intersection = line.intersection_with(other);
  require_components_near(intersection, {0.0, -1.0}, kTolerance, test_description, "intersection");
}

void test_line_3d()
{
  char const* test_description = "Line<3>";

  using Direction = math::Direction<3>;
  using Line = math::Line<3>;
  using Point = math::Point<3>;

  Direction direction{Point{0.0, 0.0, 0.0}, Point{0.0, 0.0, 1.0}};
  Point point{1.0, 2.0, 3.0};
  Line line{point, direction};

  require_components_near(line.point(), {1.0, 2.0, 3.0}, kTolerance, test_description, "point");
  require_components_near(line.direction(), {0.0, 0.0, 1.0}, kTolerance, test_description, "direction");

  Direction const& as_direction = line;
  require(&as_direction == &line.direction(), test_description, "conversion to direction");
}

} // namespace

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  try
  {
    test_point_2d();
    test_point_3d();
    test_vector_2d();
    test_vector_3d();
    test_direction_2d();
    test_direction_3d();
    test_line_piece_2d();
    test_line_piece_3d();
    test_line_2d();
    test_line_3d();
  }
  catch (std::exception const& error)
  {
    std::cerr << "math_geometry_tests failure: " << error.what() << '\n';
    return 1;
  }

  std::cout << "math_geometry_tests passed" << std::endl;
  return 0;
}

