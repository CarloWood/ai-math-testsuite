#include "sys.h"
#include "math/Direction.h"
#include "math/Line.h"
#include "math/LinePiece.h"
#include "math/Point.h"
#include "math/Vector.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include "debug.h"

namespace {

constexpr double kTolerance = 1e-9;

void require(bool condition, std::string const& message)
{
  if (!condition)
    throw std::runtime_error(message);
}

template<typename T>
void require_near(T actual, T expected, T tolerance, std::string const& message)
{
  if (std::isnan(expected))
  {
    require(std::isnan(actual), message + ": expected NaN");
    return;
  }

  if (std::isinf(expected))
  {
    require(std::isinf(actual) && std::signbit(actual) == std::signbit(expected), message + ": expected infinity");
    return;
  }

  if (std::abs(actual - expected) > tolerance)
  {
    std::ostringstream oss;
    oss << message << ": expected " << expected << ", got " << actual;
    throw std::runtime_error(oss.str());
  }
}

template<int N, typename T, typename Object>
void require_components_near(Object const& object, std::array<T, N> const& expected, T tolerance, std::string const& context)
{
  for (int i = 0; i < N; ++i)
  {
    std::ostringstream label;
    label << context << " component " << i;
    require_near(object[i], expected[i], tolerance, label.str());
  }
}

void test_point_2d()
{
  using math::Direction;
  using math::Point;
  using math::Vector;

  Point<2, double> p{1.0, 2.0};
  require_near(p.eigen()(0), 1.0, kTolerance, "Point<2> eigen x");
  require_near(p.eigen()(1), 2.0, kTolerance, "Point<2> eigen y");

  require_near(p[0], 1.0, kTolerance, "Point<2> operator[] x");
  require_near(p[1], 2.0, kTolerance, "Point<2> operator[] y");

  Point<2, double> const& cp = p;
  require_near(cp.x(), 1.0, kTolerance, "Point<2> x() const");
  require_near(cp.y(), 2.0, kTolerance, "Point<2> y() const");

  p.x() = 3.0;
  p.y() = 4.0;
  require_components_near<2>(p, {3.0, 4.0}, kTolerance, "Point<2> setters");

  Direction<2, double> dx{Point<2, double>{0.0, 0.0}, Point<2, double>{1.0, 0.0}};

  Point<2, double> translated = p + dx;
  require_components_near<2>(translated, {4.0, 4.0}, kTolerance, "Point<2> operator+(direction)");

  translated = p + Vector<2, double>{-1.0, 2.0};
  require_components_near<2>(translated, {2.0, 6.0}, kTolerance, "Point<2> operator+(vector)");

  translated = p - Vector<2, double>{-1.0, 2.0};
  require_components_near<2>(translated, {4.0, 2.0}, kTolerance, "Point<2> operator-(vector)");

  Point<2, double> mutable_point = p;
  mutable_point += dx;
  require_components_near<2>(mutable_point, {4.0, 4.0}, kTolerance, "Point<2> operator+= direction");
  mutable_point += Vector<2, double>{1.0, -1.0};
  require_components_near<2>(mutable_point, {5.0, 3.0}, kTolerance, "Point<2> operator+= vector");
  mutable_point -= Vector<2, double>{2.0, 2.0};
  require_components_near<2>(mutable_point, {3.0, 1.0}, kTolerance, "Point<2> operator-= vector");

  Point<2, double> q{3.0, 1.0};
  Vector<2, double> difference = q - p;
  require_components_near<2>(difference, {0.0, -3.0}, kTolerance, "Point<2> difference operator");

  require(!(q != q), "Point<2> equality");
  require(p != q, "Point<2> inequality");
}

void test_point_3d()
{
  using math::Direction;
  using math::Point;
  using math::Vector;

  Point<3, double> p{1.0, 2.0, 3.0};
  require_components_near<3>(p, {1.0, 2.0, 3.0}, kTolerance, "Point<3> initial values");

  p.eigen()(0) = 4.0;
  p.eigen()(1) = 5.0;
  p.eigen()(2) = 6.0;
  require_components_near<3>(p, {4.0, 5.0, 6.0}, kTolerance, "Point<3> eigen write access");

  require_near(static_cast<Point<3, double> const&>(p).x(), 4.0, kTolerance, "Point<3> x() const");
  require_near(static_cast<Point<3, double> const&>(p).y(), 5.0, kTolerance, "Point<3> y() const");
  require_near(static_cast<Point<3, double> const&>(p).z(), 6.0, kTolerance, "Point<3> z() const");

  Direction<3, double> dz{Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 1.0}};
  Vector<3, double> shift{1.0, -2.0, 0.5};

  Point<3, double> translated = p + dz;
  require_components_near<3>(translated, {4.0, 5.0, 7.0}, kTolerance, "Point<3> operator+(direction)");

  translated = p + shift;
  require_components_near<3>(translated, {5.0, 3.0, 6.5}, kTolerance, "Point<3> operator+(vector)");

  translated = p - shift;
  require_components_near<3>(translated, {3.0, 7.0, 5.5}, kTolerance, "Point<3> operator-(vector)");

  Point<3, double> mutable_point = p;
  mutable_point += dz;
  require_components_near<3>(mutable_point, {4.0, 5.0, 7.0}, kTolerance, "Point<3> operator+= direction");
  mutable_point += shift;
  require_components_near<3>(mutable_point, {5.0, 3.0, 7.5}, kTolerance, "Point<3> operator+= vector");
  mutable_point -= shift;
  require_components_near<3>(mutable_point, {4.0, 5.0, 7.0}, kTolerance, "Point<3> operator-= vector");

  Point<3, double> q{4.0, 5.0, 7.0};
  Vector<3, double> difference = q - p;
  require_components_near<3>(difference, {0.0, 0.0, 1.0}, kTolerance, "Point<3> difference operator");

  require(!(q != q), "Point<3> equality");
  require(p != q, "Point<3> inequality");
}

void test_vector_2d()
{
  using math::Direction;
  using math::LinePiece;
  using math::Point;
  using math::Vector;

  Vector<2, double> v{3.0, 4.0};
  require_components_near<2>(v, {3.0, 4.0}, kTolerance, "Vector<2> initializer");

  Direction<2, double> direction{Point<2, double>{0.0, 0.0}, Point<2, double>{1.0, 1.0}};
  Vector<2, double> vdir{direction, std::sqrt(2.0)};
  require_components_near<2>(vdir, {1.0, 1.0}, kTolerance, "Vector<2> from direction");

  Vector<2, double> from_points{Point<2, double>{1.0, 2.0}, Point<2, double>{4.0, 6.0}};
  require_components_near<2>(from_points, {3.0, 4.0}, kTolerance, "Vector<2> from two points");

  Vector<2, double> from_point{Point<2, double>{3.0, 4.0}};
  require_components_near<2>(from_point, {3.0, 4.0}, kTolerance, "Vector<2> from single point");

  LinePiece<2, double> segment{Point<2, double>{0.0, 0.0}, Point<2, double>{-3.0, -4.0}};
  Vector<2, double> from_segment{segment};
  require_components_near<2>(from_segment, {-3.0, -4.0}, kTolerance, "Vector<2> from line piece");

  Vector<2, double>::eigen_type eigen;
  eigen << -6.0, 8.0;
  Vector<2, double> from_eigen{eigen};
  require_components_near<2>(from_eigen, {-6.0, 8.0}, kTolerance, "Vector<2> from eigen");

  v[0] = 5.0;
  v[1] = -1.0;
  require_components_near<2>(v, {5.0, -1.0}, kTolerance, "Vector<2> operator[] write");
  require_near(static_cast<Vector<2, double> const&>(v)[0], 5.0, kTolerance, "Vector<2> operator[] const");

  require_near(v.x(), 5.0, kTolerance, "Vector<2> x()");
  require_near(v.y(), -1.0, kTolerance, "Vector<2> y()");

  require_near(v.dot(Vector<2, double>{2.0, 3.0}), 7.0, kTolerance, "Vector<2> dot");
  require_near(v.cross(Vector<2, double>{2.0, 3.0}), 17.0, kTolerance, "Vector<2> cross scalar");

  Direction<2, double> dir = v.direction();
  require_near(dir.x(), v.x() / v.norm(), kTolerance, "Vector<2> direction x");
  require_near(dir.y(), v.y() / v.norm(), kTolerance, "Vector<2> direction y");

  require_near(v.norm(), std::sqrt(26.0), kTolerance, "Vector<2> norm");
  require_near(v.norm_squared(), 26.0, kTolerance, "Vector<2> norm_squared");

  Point<2, double> as_point = v.as_point();
  require_components_near<2>(as_point, {5.0, -1.0}, kTolerance, "Vector<2> as_point");

  Vector<2, double> nan_vector{std::numeric_limits<double>::quiet_NaN(), 0.0};
  require(nan_vector.isnan(), "Vector<2> isnan");
  Vector<2, double> finite_vector{1.0, 1.0};
  require(finite_vector.isfinite(), "Vector<2> isfinite true");
  Vector<2, double> infinite_vector{std::numeric_limits<double>::infinity(), 1.0};
  require(!infinite_vector.isfinite(), "Vector<2> isfinite false");

  require_components_near<2>(v.rotate_90_degrees(), {1.0, 5.0}, kTolerance, "Vector<2> rotate 90");
  require_components_near<2>(v.rotate_180_degrees(), {-5.0, 1.0}, kTolerance, "Vector<2> rotate 180");
  require_components_near<2>(v.rotate_270_degrees(), {-1.0, -5.0}, kTolerance, "Vector<2> rotate 270");

  Vector<2, double> w{1.0, 2.0};
  Vector<2, double> u = v;
  u += w;
  require_components_near<2>(u, {6.0, 1.0}, kTolerance, "Vector<2> operator+=");
  u -= w;
  require_components_near<2>(u, {5.0, -1.0}, kTolerance, "Vector<2> operator-=");
  u *= 2.0;
  require_components_near<2>(u, {10.0, -2.0}, kTolerance, "Vector<2> operator*=");
  u /= 4.0;
  require_components_near<2>(u, {2.5, -0.5}, kTolerance, "Vector<2> operator/=");

  Vector<2, double> divided = v / 2.0;
  require_components_near<2>(divided, {2.5, -0.5}, kTolerance, "Vector<2> operator/(scalar)");

  require_components_near<2>(v + w, {6.0, 1.0}, kTolerance, "Vector<2> operator+");
  require_components_near<2>(v - w, {4.0, -3.0}, kTolerance, "Vector<2> operator-");
  require_components_near<2>(2.0 * v, {10.0, -2.0}, kTolerance, "Vector<2> scalar multiplication");

  Point<2, double> base{1.0, 1.0};
  require_components_near<2>(base + v, {6.0, 0.0}, kTolerance, "Point + Vector");
  require_components_near<2>(base - v, {-4.0, 2.0}, kTolerance, "Point - Vector");
}

void test_vector_3d()
{
  using math::Direction;
  using math::LinePiece;
  using math::Point;
  using math::Vector;

  Vector<3, double> v{1.0, 2.0, 3.0};
  require_components_near<3>(v, {1.0, 2.0, 3.0}, kTolerance, "Vector<3> initializer");

  Vector<3, double> from_points{Point<3, double>{1.0, 2.0, 3.0}, Point<3, double>{4.0, 6.0, 9.0}};
  require_components_near<3>(from_points, {3.0, 4.0, 6.0}, kTolerance, "Vector<3> from two points");

  Vector<3, double> from_point{Point<3, double>{3.0, 4.0, 5.0}};
  require_components_near<3>(from_point, {3.0, 4.0, 5.0}, kTolerance, "Vector<3> from single point");

  Direction<3, double> direction{Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 3.0, 4.0}};
  Vector<3, double> from_direction{direction, 5.0};
  require_components_near<3>(from_direction, {0.0, 3.0, 4.0}, kTolerance, "Vector<3> from direction");

  LinePiece<3, double> segment{Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{-1.0, -2.0, -3.0}};
  Vector<3, double> from_segment{segment};
  require_components_near<3>(from_segment, {-1.0, -2.0, -3.0}, kTolerance, "Vector<3> from line piece");

  Vector<3, double>::eigen_type eigen;
  eigen << 2.0, -1.0, 4.0;
  Vector<3, double> from_eigen{eigen};
  require_components_near<3>(from_eigen, {2.0, -1.0, 4.0}, kTolerance, "Vector<3> from eigen");

  require_near(v.dot(Vector<3, double>{4.0, 5.0, 6.0}), 32.0, kTolerance, "Vector<3> dot");

  Dout(dc::notice, "v = " << v);
  Vector<3, double> cross = v.cross(Vector<3, double>{-1.0, 0.0, 2.0});
  require_components_near<3>(cross, {4.0, -5.0, 2.0}, kTolerance, "Vector<3> cross product");

  Direction<3, double> dir = v.direction();
  double norm = v.norm();
  require_near(dir.x(), v.x() / norm, kTolerance, "Vector<3> direction x");
  require_near(dir.y(), v.y() / norm, kTolerance, "Vector<3> direction y");
  require_near(dir.z(), v.z() / norm, kTolerance, "Vector<3> direction z");

  require_near(v.norm(), std::sqrt(14.0), kTolerance, "Vector<3> norm");
  require_near(v.norm_squared(), 14.0, kTolerance, "Vector<3> norm_squared");

  Point<3, double> as_point = v.as_point();
  require_components_near<3>(as_point, {1.0, 2.0, 3.0}, kTolerance, "Vector<3> as_point");

  Vector<3, double> nan_vector{0.0, std::numeric_limits<double>::quiet_NaN(), 0.0};
  require(nan_vector.isnan(), "Vector<3> isnan");

  Vector<3, double> infinite_vector{0.0, std::numeric_limits<double>::infinity(), 0.0};
  require(!infinite_vector.isfinite(), "Vector<3> isfinite false");

  Vector<3, double> u = v;
  u += Vector<3, double>{1.0, 0.0, -1.0};
  require_components_near<3>(u, {2.0, 2.0, 2.0}, kTolerance, "Vector<3> operator+=");
  u -= Vector<3, double>{1.0, 1.0, 1.0};
  require_components_near<3>(u, {1.0, 1.0, 1.0}, kTolerance, "Vector<3> operator-=");
  u *= 3.0;
  require_components_near<3>(u, {3.0, 3.0, 3.0}, kTolerance, "Vector<3> operator*=");
  u /= 3.0;
  require_components_near<3>(u, {1.0, 1.0, 1.0}, kTolerance, "Vector<3> operator/=");

  Vector<3, double> divided = v / 2.0;
  require_components_near<3>(divided, {0.5, 1.0, 1.5}, kTolerance, "Vector<3> operator/(scalar)");

  require_components_near<3>(v + Vector<3, double>{1.0, 1.0, 1.0}, {2.0, 3.0, 4.0}, kTolerance, "Vector<3> operator+");
  require_components_near<3>(v - Vector<3, double>{1.0, 1.0, 1.0}, {0.0, 1.0, 2.0}, kTolerance, "Vector<3> operator-");

  Point<3, double> base{1.0, 2.0, 3.0};
  require_components_near<3>(base + v, {2.0, 4.0, 6.0}, kTolerance, "Point<3> + Vector<3>");
  require_components_near<3>(base - v, {0.0, 0.0, 0.0}, kTolerance, "Point<3> - Vector<3>");
}

void test_direction_2d()
{
  using math::Direction;
  using math::Line;
  using math::LinePiece;
  using math::Point;

  Direction<2, double> angle{std::numbers::pi_v<double> / 4.0};
  double expected = std::sqrt(0.5);
  require_near(angle.x(), expected, kTolerance, "Direction<2> from angle x");
  require_near(angle.y(), expected, kTolerance, "Direction<2> from angle y");

  Direction<2, double> from_points{Point<2, double>{0.0, 0.0}, Point<2, double>{0.0, 5.0}};
  require_components_near<2>(from_points, {0.0, 1.0}, kTolerance, "Direction<2> from points");

  Direction<2, double> from_point_only{Point<2, double>{0.0, -3.0}};
  require_components_near<2>(from_point_only, {0.0, -1.0}, kTolerance, "Direction<2> from single point");

  LinePiece<2, double> segment{Point<2, double>{0.0, 0.0}, Point<2, double>{1.0, 0.0}};
  Direction<2, double> from_segment{segment};
  require_components_near<2>(from_segment, {1.0, 0.0}, kTolerance, "Direction<2> from segment");

  Line<2, double> line{Point<2, double>{0.0, 0.0}, angle};
  Direction<2, double> from_line{line};
  require_components_near<2>(from_line, {angle.x(), angle.y()}, kTolerance, "Direction<2> from line");

  require_near(angle.dot(Direction<2, double>::right), angle.x(), kTolerance, "Direction<2> dot");

  require_near(from_points.as_angle(), std::numbers::pi_v<double> / 2.0, kTolerance, "Direction<2> as_angle");

  Direction<2, double> normal = angle.normal();
  require_components_near<2>(normal, {-expected, expected}, kTolerance, "Direction<2> normal");

  Direction<2, double> inverse = angle.inverse();
  require_components_near<2>(inverse, {-angle.x(), -angle.y()}, kTolerance, "Direction<2> inverse");

  Direction<2, double> normal_inverse = angle.normal_inverse();
  require_components_near<2>(normal_inverse, {expected, -expected}, kTolerance, "Direction<2> normal_inverse");

  require_components_near<2>(Direction<2, double>::up, {0.0, 1.0}, kTolerance, "Direction<2> up");
  require_components_near<2>(Direction<2, double>::down, {0.0, -1.0}, kTolerance, "Direction<2> down");
  require_components_near<2>(Direction<2, double>::left, {-1.0, 0.0}, kTolerance, "Direction<2> left");
  require_components_near<2>(Direction<2, double>::right, {1.0, 0.0}, kTolerance, "Direction<2> right");
}

void test_direction_3d()
{
  using math::Direction;
  using math::Line;
  using math::Point;

  Direction<3, double> from_points{Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{1.0, 2.0, 2.0}};
  double norm = std::sqrt(9.0);
  require_components_near<3>(from_points, {1.0 / norm, 2.0 / norm, 2.0 / norm}, kTolerance, "Direction<3> from points");

  Direction<3, double> from_point_only{Point<3, double>{0.0, 0.0, -3.0}};
  require_components_near<3>(from_point_only, {0.0, 0.0, -1.0}, kTolerance, "Direction<3> from origin to point");

  Line<3, double> line{Point<3, double>{1.0, 1.0, 1.0}, from_points};
  Direction<3, double> from_line{line};
  require_components_near<3>(from_line, {1.0 / norm, 2.0 / norm, 2.0 / norm}, kTolerance, "Direction<3> from line");

  require_near(from_points.dot(from_line), 1.0, kTolerance, "Direction<3> dot");

  require_near(from_points.as_angle(), std::atan2(2.0 / norm, 1.0 / norm), kTolerance, "Direction<3> as_angle");
}

void test_line_piece_2d()
{
  using math::LinePiece;
  using math::Point;

  Point<2, double> a{0.0, 0.0};
  Point<2, double> b{3.0, 4.0};
  LinePiece<2, double> segment{a, b};

  require_components_near<2>(segment.from(), {0.0, 0.0}, kTolerance, "LinePiece<2> from");
  require_components_near<2>(segment.to(), {3.0, 4.0}, kTolerance, "LinePiece<2> to");
  require_near(segment.norm(), 5.0, kTolerance, "LinePiece<2> length");

  auto direction = segment.direction();
  require_components_near<2>(direction, {0.6, 0.8}, kTolerance, "LinePiece<2> direction");
}

void test_line_piece_3d()
{
  using math::LinePiece;
  using math::Point;

  Point<3, double> a{0.0, 0.0, 0.0};
  Point<3, double> b{1.0, 2.0, 2.0};
  LinePiece<3, double> segment{a, b};

  require_near(segment.norm(), 3.0, kTolerance, "LinePiece<3> length");

  auto direction = segment.direction();
  require_components_near<3>(direction, {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}, kTolerance, "LinePiece<3> direction");
}

void test_line_2d()
{
  using math::Direction;
  using math::Line;
  using math::Point;

  Direction<2, double> slope{Point<2, double>{0.0, 0.0}, Point<2, double>{1.0, 2.0}};
  Point<2, double> point{1.0, 1.0};
  Line<2, double> line{point, slope};

  require_components_near<2>(line.point(), {1.0, 1.0}, kTolerance, "Line<2> point");
  require_components_near<2>(line.direction(), {slope.x(), slope.y()}, kTolerance, "Line<2> direction");

  Direction<2, double> const& as_direction = line;
  require(&as_direction == &line.direction(), "Line<2> conversion to direction");

  Line<2, double> other{Point<2, double>{0.0, 2.0}, Direction<2, double>::down};
  Point<2, double> intersection = line.intersection_with(other);
  require_components_near<2>(intersection, {0.0, 2.0}, kTolerance, "Line<2> intersection");
}

void test_line_3d()
{
  using math::Direction;
  using math::Line;
  using math::Point;

  Direction<3, double> direction{Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 1.0}};
  Point<3, double> point{1.0, 2.0, 3.0};
  Line<3, double> line{point, direction};

  require_components_near<3>(line.point(), {1.0, 2.0, 3.0}, kTolerance, "Line<3> point");
  require_components_near<3>(line.direction(), {0.0, 0.0, 1.0}, kTolerance, "Line<3> direction");

  Direction<3, double> const& as_direction = line;
  require(&as_direction == &line.direction(), "Line<3> conversion to direction");
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

