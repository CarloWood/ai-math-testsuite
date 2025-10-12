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
constexpr double kExact = 0.0;

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

template<typename Object>
concept ConceptObjectDouble =
  requires { object_traits<Object>::dim; typename object_traits<Object>::scalar_type; } &&
  std::same_as<typename object_traits<Object>::scalar_type, double>;

template<ConceptObjectDouble Object>
void require_components_near(Object const& object,
                             std::initializer_list<double> expected_list,
                             double tolerance,
                             std::string const& test_description,
                             std::string const& context)
{
  constexpr int N = static_cast<int>(object_traits<Object>::dim);
  ASSERT(expected_list.size() >= N);
  double const* expected = expected_list.begin();
  for (int i = 0; i < N; ++i, ++expected)
  {
    std::ostringstream label;
    label << context << " component " << i;
    require_near(object[i], *expected, tolerance, test_description, label.str());
  }
}

template<ConceptObjectDouble Object>
void require_components_equal(Object const& object,
                              std::initializer_list<double> expected_list,
                              std::string const& test_description,
                              std::string const& context)
{
  require_components_near(object, expected_list, kExact, test_description, context);
}

template <int I, class... Ts>
decltype(auto) get(Ts&&... ts)
{
  return std::get<I>(std::forward_as_tuple(ts...));
}

template<ConceptObjectDouble Object, typename... Args>
auto make(Args&&... args)
{
  constexpr int N = static_cast<int>(object_traits<Object>::dim);
  if constexpr (N == 2)
    return Object{get<0>(args...), get<1>(args...)};
  else
    return Object{get<0>(args...), get<1>(args...), get<2>(args...)};
};

template<int N>
auto make_array(double x, double y, double z)
{
  std::array<double, N> result{};
  result[0] = x;
  result[1] = y;
  if constexpr (N == 3)
    result[2] = z;
  return result;
};

template<int N>
auto make_norm(double x, double y, double z)
{
  double squared = x * x + y * y;
  if constexpr (N == 3)
    squared += z * z;
  return std::sqrt(squared);
};

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
  Point p1 = make<Point>(cd_x, y, z);
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

  require_components_equal(p1, {x1, y1, z1}, test_description, "setters");

  //---------------------------------------------------------------------------
  // operators

  // Use p2 as a direction.
  Direction const d2{p2};
  double dx = p2.x();
  double dy = p2.y();
  double dz = N == 3 ? p2[2] : 0.0;
  double const len = make_norm<N>(dx, dy, dz);
  dx /= len;
  dy /= len;
  dz /= len;

  // Add Direction to a Point.
  Point translated = p1 + d2;
  require_components_near(translated, {x1 + dx, y1 + dy, z1 + dz}, kTolerance, test_description, "operator+(direction)");

  // Use p2 as a Vector.
  Vector const shift{p2};

  translated = p1 + shift;
  require_components_near(translated, {x1 + x2, y1 + y2, z1 + z2}, kTolerance, test_description, "operator+(vector)");

  translated = p1 - shift;
  require_components_near(translated, {x1 - x2, y1 - y2, z1 - z2}, kTolerance, test_description, "operator-(vector)");

  Point mutable_point = p1;
  mutable_point += d2;
  require_components_near(mutable_point, {x1 + dx, y1 + dy, z1 + dz}, kTolerance, test_description, "operator+= direction");
  mutable_point += shift;
  require_components_near(mutable_point, {x1 + dx + x2, y1 + dy + y2, z1 + dz + z2}, kTolerance, test_description, "operator+= vector");
  mutable_point -= shift;
  require_components_near(mutable_point, {x1 + dx, y1 + dy, z1 + dz}, kTolerance, test_description, "operator-= vector");

  Vector difference = p1 - p2;
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

template<int N>
void test_vector(std::string test_description)
{
  using Vector = math::Vector<N>;
  using Direction = math::Direction<N>;
  using Point = math::Point<N>;
  using LinePiece = math::LinePiece<N>;

  double const x1 = 1.4829374;
  double const y1 = -2.7138465;
  double const z1 = 0.9271354;
  double const x2 = -3.1847265;
  double const y2 = 4.2716384;
  double const z2 = -1.6372945;
  double const x3 = 0.5172834;
  double const y3 = -1.2947382;
  double const z3 = 2.8437591;

  double const sx = -0.9473625;
  double const sy = 1.5738295;
  double const sz = -2.1937453;
  double const ex = 0.6831725;
  double const ey = -1.4928372;
  double const ez = 1.0384726;

  double const eigen_x = -0.4718293;
  double const eigen_y = 2.6173845;
  double const eigen_z = -1.9473825;
  double const scalar = 1.7325;

  Point const start_point = make<Point>(sx, sy, sz);
  Point const end_point = make<Point>(ex, ey, ez);

  double const dx = ex - sx;
  double const dy = ey - sy;
  double const dz = ez - sz;
  double const delta_norm = make_norm<N>(dx, dy, dz);

  // Constructor from component values.
  Vector v = make<Vector>(x1, y1, z1);
  require_components_near(v, {x1, y1, z1}, kTolerance, test_description, "component constructor");

  // Constructor from a Direction and length.
  double const direction_length = delta_norm * 0.642713;
  Direction const direction{start_point, end_point};
  Vector const from_direction{direction, direction_length};
  require_components_near(from_direction,
      {dx * (direction_length / delta_norm),
       dy * (direction_length / delta_norm),
       dz * (direction_length / delta_norm)},
      kTolerance,
      test_description,
      "constructor(direction, length)");

  // Constructor from two Point objects.
  Vector const from_points{start_point, end_point};
  require_components_near(from_points, {dx, dy, dz}, kTolerance, test_description, "constructor(point, point)");

  // Constructor from a single Point.
  Vector const from_point{end_point};
  require_components_near(from_point, {ex, ey, ez}, kTolerance, test_description, "constructor(point)");

  // Constructor from a LinePiece.
  LinePiece const segment{start_point, end_point};
  Vector const from_segment{segment};
  require_components_near(from_segment, {dx, dy, dz}, kTolerance, test_description, "constructor(line_piece)");

  // Constructor from eigen_type.
  typename Vector::eigen_type eigen_input;
  if constexpr (N == 2)
    eigen_input << eigen_x, eigen_y;
  else
    eigen_input << eigen_x, eigen_y, eigen_z;
  Vector const from_eigen{eigen_input};
  require_components_near(from_eigen, {eigen_x, eigen_y, eigen_z}, kTolerance, test_description, "constructor(eigen)");

  // Mutable element access via operator[].
  double const new_x = -4.6138275;
  double const new_y = 1.9375421;
  double const new_z = -0.3827564;
  v[0] = new_x;
  v[1] = new_y;
  if constexpr (N == 3)
    v[2] = new_z;
  require_components_near(v, {new_x, new_y, new_z}, kTolerance, test_description, "operator[] write");

  // Const element access via operator[].
  Vector const& const_ref = v;
  require_near(const_ref[0], new_x, kTolerance, test_description, "operator[] const x");
  require_near(const_ref[1], new_y, kTolerance, test_description, "operator[] const y");
  if constexpr (N == 3)
    require_near(const_ref[2], new_z, kTolerance, test_description, "operator[] const z");

  // Coordinate accessors.
  require_near(v.x(), new_x, kTolerance, test_description, "x()");
  require_near(v.y(), new_y, kTolerance, test_description, "y()");
  if constexpr (N == 3)
    require_near(v.z(), new_z, kTolerance, test_description, "z()");

  Vector const w = make<Vector>(x2, y2, z2);

  // Dot product.
  double const expected_dot = new_x * x2 + new_y * y2 + (N == 3 ? new_z * z2 : 0.0);
  require_near(v.dot(w), expected_dot, kTolerance, test_description, "dot");

  // Cross product (scalar for 2D, Vector for 3D).
  if constexpr (N == 2)
  {
    double const expected_cross = new_x * y2 - new_y * x2;
    require_near(v.cross(w), expected_cross, kTolerance, test_description, "cross");
  }
  else
  {
    require_components_near(v.cross(w),
        {new_y * z2 - new_z * y2,
         new_z * x2 - new_x * z2,
         new_x * y2 - new_y * x2},
        kTolerance,
        test_description,
        "cross");
  }

  // Direction conversion.
  double const norm = make_norm<N>(new_x, new_y, new_z);
  Direction const dir = v.direction();
  require_near(dir.x(), new_x / norm, kTolerance, test_description, "direction x");
  require_near(dir.y(), new_y / norm, kTolerance, test_description, "direction y");
  if constexpr (N == 3)
    require_near(dir.z(), new_z / norm, kTolerance, test_description, "direction z");

  // Norm calculations.
  double const norm_squared = new_x * new_x + new_y * new_y + (N == 3 ? new_z * new_z : 0.0);
  require_near(v.norm(), norm, kTolerance, test_description, "norm");
  require_near(v.norm_squared(), norm_squared, kTolerance, test_description, "norm_squared");

  // Conversion to Point.
  Point const as_point = v.as_point();
  require_components_near(as_point, {new_x, new_y, new_z}, kTolerance, test_description, "as_point");

  // NaN and infinity checks.
  Vector const nan_vector = make<Vector>(std::numeric_limits<double>::quiet_NaN(), y3, z3);
  require(nan_vector.isnan(), test_description, "isnan");
  Vector const finite_vector = make<Vector>(x3, y3, z3);
  require(finite_vector.isfinite(), test_description, "isfinite true");
  Vector const infinite_vector = make<Vector>(std::numeric_limits<double>::infinity(), y3, z3);
  require(!infinite_vector.isfinite(), test_description, "isfinite false");

  // Rotations (2D only).
  if constexpr (N == 2)
  {
    require_components_near(v.rotate_90_degrees(), {-new_y, new_x}, kTolerance, test_description, "rotate_90_degrees");
    require_components_near(v.rotate_180_degrees(), {-new_x, -new_y}, kTolerance, test_description, "rotate_180_degrees");
    require_components_near(v.rotate_270_degrees(), {new_y, -new_x}, kTolerance, test_description, "rotate_270_degrees");
  }

  // Compound assignment operators.
  Vector u = v;
  u += w;
  require_components_near(u, {new_x + x2, new_y + y2, new_z + z2}, kTolerance, test_description, "operator+=");
  u -= w;
  require_components_near(u, {new_x, new_y, new_z}, kTolerance, test_description, "operator-=");
  u *= scalar;
  require_components_near(u, {new_x * scalar, new_y * scalar, new_z * scalar}, kTolerance, test_description, "operator*=");
  u /= scalar;
  require_components_near(u, {new_x, new_y, new_z}, kTolerance, test_description, "operator/=");

  // Binary arithmetic operators.
  Vector const divided = v / scalar;
  require_components_near(divided, {new_x / scalar, new_y / scalar, new_z / scalar}, kTolerance, test_description, "operator/(scalar)");
  require_components_near(v + w, {new_x + x2, new_y + y2, new_z + z2}, kTolerance, test_description, "operator+");
  require_components_near(v - w, {new_x - x2, new_y - y2, new_z - z2}, kTolerance, test_description, "operator-");
  require_components_near(scalar * v, {scalar * new_x, scalar * new_y, scalar * new_z}, kTolerance, test_description, "scalar multiplication");

  // Point and Vector interactions.
  Point const base = make<Point>(x3, y3, z3);
  require_components_near(base + v, {x3 + new_x, y3 + new_y, z3 + new_z}, kTolerance, test_description, "point + vector");
  require_components_near(base - v, {x3 - new_x, y3 - new_y, z3 - new_z}, kTolerance, test_description, "point - vector");
}

void test_vector_2d()
{
  test_vector<2>("Vector<2>");
}

void test_vector_3d()
{
  test_vector<3>("Vector<3>");
}


template<int N>
void test_direction(std::string const& test_description)
{
  using Direction = math::Direction<N>;
  using Line = math::Line<N>;
  using LinePiece = math::LinePiece<N>;
  using Point = math::Point<N>;

  double const angle_value = 1.2743812934672198;
  double const x1 = -0.338451927341;
  double const y1 = 2.19453123741;
  double const z1 = -1.83452987412;
  double const x2 = 1.18234752341;
  double const y2 = -0.71325894732;
  double const z2 = 0.95234712853;
  double const single_x = -0.67123984512;
  double const single_y = 0.44231589743;
  double const single_z = -1.25123987431;
  double const segment_from_x = -0.52734987123;
  double const segment_from_y = 0.38421598273;
  double const segment_from_z = -0.91254378214;
  double const segment_to_x = 0.69321734981;
  double const segment_to_y = -0.18243798231;
  double const segment_to_z = 0.41235678294;
  double const line_point_x = 0.21843987123;
  double const line_point_y = -1.84352987341;
  double const line_point_z = 0.74231578293;
  double const zero = 0.0;
  double const one = 1.0;
  double const minus_one = -1.0;

  //---------------------------------------------------------------------------
  // Constructor from two points.
  Point const point_a = make<Point>(x1, y1, z1);
  Point const point_b = make<Point>(x2, y2, z2);
  Direction const from_points{point_a, point_b};

  double const diff_x = x2 - x1;
  double const diff_y = y2 - y1;
  double const diff_z = z2 - z1;
  double const diff_norm = make_norm<N>(diff_x, diff_y, diff_z);
  require_components_near(from_points,
                          {diff_x / diff_norm, diff_y / diff_norm, diff_z / diff_norm},
                          kTolerance,
                          test_description,
                          "constructor from two points");

  //---------------------------------------------------------------------------
  // Constructor from a single point.
  Point const point_only = make<Point>(single_x, single_y, single_z);
  Direction const from_point_only{point_only};
  double const point_norm = make_norm<N>(single_x, single_y, single_z);
  require_components_near(from_point_only,
                          {single_x / point_norm, single_y / point_norm, single_z / point_norm},
                          kTolerance,
                          test_description,
                          "constructor from point");

  //---------------------------------------------------------------------------
  // Constructor from a line piece.
  Point const segment_from = make<Point>(segment_from_x, segment_from_y, segment_from_z);
  Point const segment_to = make<Point>(segment_to_x, segment_to_y, segment_to_z);
  LinePiece const segment{segment_from, segment_to};
  Direction const from_segment{segment};
  double const segment_diff_x = segment_to_x - segment_from_x;
  double const segment_diff_y = segment_to_y - segment_from_y;
  double const segment_diff_z = segment_to_z - segment_from_z;
  double const segment_norm = make_norm<N>(segment_diff_x, segment_diff_y, segment_diff_z);
  require_components_near(from_segment,
                          {segment_diff_x / segment_norm,
                           segment_diff_y / segment_norm,
                           segment_diff_z / segment_norm},
                          kTolerance,
                          test_description,
                          "constructor from line piece");

  //---------------------------------------------------------------------------
  // Constructor from a line.
  Point const base_point = make<Point>(line_point_x, line_point_y, line_point_z);
  Line const line{base_point, from_points};
  Direction const from_line{line};
  require_components_near(from_line,
                          {diff_x / diff_norm, diff_y / diff_norm, diff_z / diff_norm},
                          kTolerance,
                          test_description,
                          "constructor from line");

  //---------------------------------------------------------------------------
  // Dot product.
  double dot_expected =
    (diff_x / diff_norm) * (diff_x / diff_norm) +
    (diff_y / diff_norm) * (diff_y / diff_norm);
  if constexpr (N == 3)
    dot_expected += (diff_z / diff_norm) * (diff_z / diff_norm);
  require_near(from_points.dot(from_line), dot_expected, kTolerance, test_description, "dot product");

  //---------------------------------------------------------------------------
  // Angle conversion.
  double const expected_angle = std::atan2(diff_y / diff_norm, diff_x / diff_norm);
  require_near(from_points.as_angle(), expected_angle, kTolerance, test_description, "as_angle");

  if constexpr (N == 2)
  {
    //-------------------------------------------------------------------------
    // Constructor from angle and accessors.
    Direction const from_angle{angle_value};
    double const cos_angle = std::cos(angle_value);
    double const sin_angle = std::sin(angle_value);
    require_components_near(from_angle,
                            {cos_angle, sin_angle, zero},
                            kTolerance,
                            test_description,
                            "constructor from angle");

    //-------------------------------------------------------------------------
    // normal().
    Direction const normal = from_angle.normal();
    require_components_near(normal,
                            {-sin_angle, cos_angle, zero},
                            kTolerance,
                            test_description,
                            "normal()");

    //-------------------------------------------------------------------------
    // inverse().
    Direction const inverse = from_angle.inverse();
    require_components_near(inverse,
                            {-cos_angle, -sin_angle, zero},
                            kTolerance,
                            test_description,
                            "inverse()");

    //-------------------------------------------------------------------------
    // normal_inverse().
    Direction const normal_inverse = from_angle.normal_inverse();
    require_components_near(normal_inverse,
                            {sin_angle, -cos_angle, zero},
                            kTolerance,
                            test_description,
                            "normal_inverse()");

    //-------------------------------------------------------------------------
    // dot().
    Direction const& right = Direction::right;
    double const dot_with_right = cos_angle * right[0] + sin_angle * right[1];
    require_near(from_angle.dot(right), dot_with_right, kTolerance, test_description, "dot with right");

    //-------------------------------------------------------------------------
    // Axis aligned constants.
    require_components_near(Direction::up, {zero, one, zero}, kTolerance, test_description, "up");
    require_components_near(Direction::down, {zero, minus_one, zero}, kTolerance, test_description, "down");
    require_components_near(Direction::left, {minus_one, zero, zero}, kTolerance, test_description, "left");
    require_components_near(Direction::right, {one, zero, zero}, kTolerance, test_description, "right");
  }
}

template<int N>
void test_line_piece(std::string const& test_description)
{
  using LinePiece = math::LinePiece<N>;
  using Point = math::Point<N>;

  double const from_x = -0.84123984512;
  double const from_y = 0.31245789341;
  double const from_z = -1.2841293475;
  double const to_x = 1.27453782341;
  double const to_y = -0.52431987412;
  double const to_z = 0.71345872341;

  //---------------------------------------------------------------------------
  // Constructor and from()/to() accessors.
  Point const from_point = make<Point>(from_x, from_y, from_z);
  Point const to_point = make<Point>(to_x, to_y, to_z);
  LinePiece const segment{from_point, to_point};
  require_components_near(segment.from(), {from_x, from_y, from_z}, kTolerance, test_description, "from()");
  require_components_near(segment.to(), {to_x, to_y, to_z}, kTolerance, test_description, "to()");

  //---------------------------------------------------------------------------
  // norm().
  double const delta_x = to_x - from_x;
  double const delta_y = to_y - from_y;
  double const delta_z = to_z - from_z;
  double const expected_norm = make_norm<N>(delta_x, delta_y, delta_z);
  require_near(segment.norm(), expected_norm, kTolerance, test_description, "norm()");

  //---------------------------------------------------------------------------
  // direction().
  require_components_near(segment.direction(),
                          {delta_x / expected_norm, delta_y / expected_norm, delta_z / expected_norm},
                          kTolerance,
                          test_description,
                          "direction()");
}

template<int N>
void test_line(std::string const& test_description)
{
  using Direction = math::Direction<N>;
  using Line = math::Line<N>;
  using Point = math::Point<N>;

  double const intersection_x = -0.41235897231;
  double const intersection_y = 0.91234782941;
  double const intersection_z = -1.12458723914;
  double const dir1_dx = 0.68123984723;
  double const dir1_dy = -0.43218973412;
  double const dir1_dz = 0.29341872341;
  double const offset1 = -0.76234987123;
  double const dir2_dx = -0.38412987342;
  double const dir2_dy = -0.59124873214;
  double const dir2_dz = 0.71235897431;
  double const offset2 = 0.54873912873;

  double const line_point_x = intersection_x + offset1 * dir1_dx;
  double const line_point_y = intersection_y + offset1 * dir1_dy;
  double const line_point_z = intersection_z + offset1 * dir1_dz;
  double const line_point_next_x = line_point_x + dir1_dx;
  double const line_point_next_y = line_point_y + dir1_dy;
  double const line_point_next_z = line_point_z + dir1_dz;

  Point const base_point = make<Point>(line_point_x, line_point_y, line_point_z);
  Point const forward_point = make<Point>(line_point_next_x, line_point_next_y, line_point_next_z);
  Direction const direction{base_point, forward_point};
  Line const line{base_point, direction};

  //---------------------------------------------------------------------------
  // point().
  require_components_near(line.point(), {line_point_x, line_point_y, line_point_z}, kTolerance, test_description, "point()");

  //---------------------------------------------------------------------------
  // direction().
  double const dir1_norm = make_norm<N>(dir1_dx, dir1_dy, dir1_dz);
  require_components_near(line.direction(),
                          {dir1_dx / dir1_norm, dir1_dy / dir1_norm, dir1_dz / dir1_norm},
                          kTolerance,
                          test_description,
                          "direction()");

  //---------------------------------------------------------------------------
  // Conversion to Direction.
  Direction const& as_direction = line;
  require(&as_direction == &line.direction(), test_description, "conversion to direction");

  if constexpr (N == 2)
  {
    //-------------------------------------------------------------------------
    // intersection_with().
    double const other_point_x = intersection_x + offset2 * dir2_dx;
    double const other_point_y = intersection_y + offset2 * dir2_dy;
    double const other_point_z = intersection_z + offset2 * dir2_dz;
    double const other_point_next_x = other_point_x + dir2_dx;
    double const other_point_next_y = other_point_y + dir2_dy;
    double const other_point_next_z = other_point_z + dir2_dz;

    Point const other_point = make<Point>(other_point_x, other_point_y, other_point_z);
    Point const other_forward = make<Point>(other_point_next_x, other_point_next_y, other_point_next_z);
    Direction const other_direction{other_point, other_forward};
    Line const other_line{other_point, other_direction};

    double const delta_x = other_point_x - line_point_x;
    double const delta_y = other_point_y - line_point_y;
    double const determinant = dir1_dx * dir2_dy - dir1_dy * dir2_dx;
    double const t1 = (delta_x * dir2_dy - delta_y * dir2_dx) / determinant;
    double const expected_x = line_point_x + t1 * dir1_dx;
    double const expected_y = line_point_y + t1 * dir1_dy;
    double const expected_z = line_point_z + t1 * dir1_dz;

    Point const intersection = line.intersection_with(other_line);
    require_components_near(intersection, {expected_x, expected_y, expected_z}, kTolerance, test_description, "intersection_with()");
  }
}

void test_direction_2d()
{
  test_direction<2>("Direction<2>");
}

void test_direction_3d()
{
  test_direction<3>("Direction<3>");
}


void test_line_piece_2d()
{
  test_line_piece<2>("LinePiece<2>");
}

void test_line_piece_3d()
{
  test_line_piece<3>("LinePiece<3>");
}

void test_line_2d()
{
  test_line<2>("Line<2>");
}

void test_line_3d()
{
  test_line<3>("Line<3>");
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

