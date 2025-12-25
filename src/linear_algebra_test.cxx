#include "sys.h"

#include "math/Matrix.h"
#include "math/Point.h"
#include "math/Vector.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "debug.h"

namespace {

constexpr double kTolerance = 1e-12;

void require(bool condition, std::string const& context)
{
  if (!condition)
    throw std::runtime_error(context);
}

void require_near(double actual, double expected, std::string const& context)
{
  if (std::abs(actual - expected) > kTolerance)
  {
    std::ostringstream oss;
    oss << context << ": expected " << expected << ", got " << actual;
    throw std::runtime_error(oss.str());
  }
}

template<int N, typename P>
void require_point_near(P const& p, std::array<double, N> const& expected, std::string const& context)
{
  for (int i = 0; i < N; ++i)
    require_near(p[i], expected[i], context + ", i=" + std::to_string(i));
}

template<int N, typename V>
void require_vector_near(V const& v, std::array<double, N> const& expected, std::string const& context)
{
  for (int i = 0; i < N; ++i)
    require_near(v[i], expected[i], context + ", i=" + std::to_string(i));
}

//------------------------------------------------------------------------------
// Wrapper types to exercise generic *Ops paths.

template<int N, typename T = double>
class WrappedPoint;

template<int N, typename T = double>
class WrappedVector;

template<int N, typename T = double>
class WrappedDirection;

template<int N, int M, typename T = double>
class WrappedMatrix;

template<int N, typename T>
struct WrappedPointTypes
{
  static constexpr int n = N;
  using scalar_type = T;
  using derived_type = WrappedPoint<N, T>;
  using vector_type = WrappedVector<N, T>;
  using direction_type = WrappedDirection<N, T>;
};

template<int N, typename T>
struct WrappedVectorTypes
{
  static constexpr int n = N;
  using scalar_type = T;
  using derived_type = WrappedVector<N, T>;
  using point_type = WrappedPoint<N, T>;
  using direction_type = WrappedDirection<N, T>;
};

template<int N, int M, typename T>
struct WrappedMatrixTypes
{
  static constexpr int rows = N;
  static constexpr int cols = M;
  using scalar_type = T;
  using derived_type = WrappedMatrix<N, M, T>;
  template<int R, int C>
  using matrix_type = WrappedMatrix<R, C, T>;
  template<int D>
  using vector_type = WrappedVector<D, T>;
};

template<int N, typename T>
class WrappedDirection
{
 public:
  using raw_type = math::Direction<N, T>;

 private:
  raw_type raw_;

 public:
  WrappedDirection() = default;
  explicit WrappedDirection(raw_type const& raw) : raw_(raw) { }
  WrappedDirection(WrappedPoint<N, T> const& from, WrappedPoint<N, T> const& to) : raw_(from.raw(), to.raw()) { }

  raw_type& raw() { return raw_; }
  raw_type const& raw() const { return raw_; }
};

template<int N, typename T>
class WrappedVector : public math::VectorOps<WrappedVectorTypes<N, T>>
{
 public:
  using raw_type = math::Vector<N, T>;
  using eigen_type = typename raw_type::eigen_type;

 private:
  raw_type raw_;

 public:
  WrappedVector() = default;
  WrappedVector(T x, T y) requires (N == 2) : raw_(x, y) { }
  WrappedVector(WrappedPoint<N, T> const& from, WrappedPoint<N, T> const& to) : raw_(from.raw(), to.raw()) { }

  WrappedVector(eigen_type const& v) : raw_(v) { }
  explicit WrappedVector(raw_type const& raw) : raw_(raw) { }

  raw_type& raw() { return raw_; }
  raw_type const& raw() const { return raw_; }
};

template<int N, typename T>
class WrappedPoint : public math::PointOps<WrappedPointTypes<N, T>>
{
 public:
  using raw_type = math::Point<N, T>;
  using eigen_type = typename raw_type::eigen_type;

 private:
  raw_type raw_;

 public:
  WrappedPoint() = default;
  WrappedPoint(T x, T y) requires (N == 2) : raw_(x, y) { }

  WrappedPoint(eigen_type const& point) : raw_(point) { }
  explicit WrappedPoint(raw_type const& raw) : raw_(raw) { }

  raw_type& raw() { return raw_; }
  raw_type const& raw() const { return raw_; }
};

template<int N, int M, typename T>
class WrappedMatrix : public math::MatrixOps<WrappedMatrixTypes<N, M, T>>
{
 public:
  using raw_type = math::Matrix<N, M, T>;
  using eigen_type = typename raw_type::eigen_type;

 private:
  raw_type raw_;

 public:
  WrappedMatrix() = default;
  WrappedMatrix(eigen_type const& m) : raw_(m) { }
  explicit WrappedMatrix(raw_type const& raw) : raw_(raw) { }

  raw_type& raw() { return raw_; }
  raw_type const& raw() const { return raw_; }
};

//------------------------------------------------------------------------------
// Tests

void test_point_ops_math()
{
  math::Point<2> const p{1.0, -2.5};
  math::Point<2> const origin{0.0, 0.0};
  math::Point<2> const x1{-0.3, 0.4};
  math::Direction<2> const dx{origin, x1};
  math::Vector<2> const v{3.0, -1.0};

  require_point_near<2>(p + dx, {0.4, -1.7}, "Point + Direction");
  {
    auto q = p;
    q += dx;
    require_point_near<2>(q, {0.4, -1.7}, "Point += Direction");
  }

  require_point_near<2>(p + v, {4.0, -3.5}, "Point + Vector");
  require_point_near<2>(p - v, {-2.0, -1.5}, "Point - Vector");

  {
    auto q = p;
    q += v;
    require_point_near<2>(q, {4.0, -3.5}, "Point += Vector");
    q -= v;
    require_point_near<2>(q, {1.0, -2.5}, "Point -= Vector");
  }

  math::Vector<2> diff = math::Point<2>{4.0, 7.0} - math::Point<2>{1.0, 2.0};
  require_vector_near<2>(diff, {3.0, 5.0}, "Point - Point (free)");

  require(math::Point<2>{1.0, 2.0} != math::Point<2>{1.0, 3.0}, "Point != (free)");
  require(!(math::Point<2>{1.0, 2.0} != math::Point<2>{1.0, 2.0}), "Point != (free) false");
}

void test_point_ops_wrapped()
{
  WrappedPoint<2> const p{1.0, -2.5};
  WrappedPoint<2> const origin{0.0, 0.0};
  WrappedPoint<2> const x1{-0.3, 0.4};
  WrappedDirection<2> const dx{origin, x1};
  WrappedVector<2> const v{3.0, -1.0};

  require_point_near<2>(p + dx, {0.4, -1.7}, "WrappedPoint + WrappedDirection");
  {
    auto q = p;
    q += dx;
    require_point_near<2>(q, {0.4, -1.7}, "WrappedPoint += WrappedDirection");
  }

  require_point_near<2>(p + v, {4.0, -3.5}, "WrappedPoint + WrappedVector");
  require_point_near<2>(p - v, {-2.0, -1.5}, "WrappedPoint - WrappedVector");

  {
    auto q = p;
    q += v;
    require_point_near<2>(q, {4.0, -3.5}, "WrappedPoint += WrappedVector");
    q -= v;
    require_point_near<2>(q, {1.0, -2.5}, "WrappedPoint -= WrappedVector");
  }

  WrappedVector<2> diff = WrappedPoint<2>{4.0, 7.0} - WrappedPoint<2>{1.0, 2.0};
  require_vector_near<2>(diff, {3.0, 5.0}, "WrappedPoint - WrappedPoint (free)");

  require(WrappedPoint<2>{1.0, 2.0} != WrappedPoint<2>{1.0, 3.0}, "WrappedPoint != (free)");
  require(!(WrappedPoint<2>{1.0, 2.0} != WrappedPoint<2>{1.0, 2.0}), "WrappedPoint != (free) false");
}

void test_vector_ops_math()
{
  math::Vector<2> const v{3.0, 4.0};
  require_near(v.norm(), 5.0, "Vector norm");
  require_near(v.norm_squared(), 25.0, "Vector norm_squared");

  math::Vector<2> const w{1.0, 2.0};
  require_near(v.dot(w), 11.0, "Vector dot");
  require_near(v.cross(w), 2.0, "Vector cross (2D)");

  require_vector_near<2>(v.rotate_90_degrees(), {-4.0, 3.0}, "Vector rotate_90_degrees");
  require_vector_near<2>(v.rotate_180_degrees(), {-3.0, -4.0}, "Vector rotate_180_degrees");
  require_vector_near<2>(v.rotate_270_degrees(), {4.0, -3.0}, "Vector rotate_270_degrees");

  {
    auto u = v;
    u += w;
    require_vector_near<2>(u, {4.0, 6.0}, "Vector +=");
    u -= w;
    require_vector_near<2>(u, {3.0, 4.0}, "Vector -=");
    u *= 2.0;
    require_vector_near<2>(u, {6.0, 8.0}, "Vector *=");
    u /= 2.0;
    require_vector_near<2>(u, {3.0, 4.0}, "Vector /=");
  }

  {
    auto u = v;
    u.negate();
    require_vector_near<2>(u, {-3.0, -4.0}, "Vector negate()");
    require_vector_near<2>(-u, {3.0, 4.0}, "Vector unary -");
  }

  require_vector_near<2>(v / 2.0, {1.5, 2.0}, "Vector / scalar");

  // Free functions from Vector.h.
  require_vector_near<2>(2.0 * v, {6.0, 8.0}, "scalar * Vector (free)");
  require_point_near<2>(math::Point<2>{1.0, 2.0} + v, {4.0, 6.0}, "Point + Vector (free)");
  require_point_near<2>(math::Point<2>{1.0, 2.0} - v, {-2.0, -2.0}, "Point - Vector (free)");
  require_vector_near<2>(v + w, {4.0, 6.0}, "Vector + Vector (free)");
  require_vector_near<2>(v - w, {2.0, 2.0}, "Vector - Vector (free)");

  // direction() and as_point().
  auto d = v.direction();
  require_near(d.dot(math::Direction<2>{math::Point<2>{0.0, 0.0}, math::Point<2>{1.0, 0.0}}), 3.0 / 5.0, "Vector direction()");
  require_point_near<2>(v.as_point(), {3.0, 4.0}, "Vector as_point()");

  // normalize()
  {
    auto u = v;
    require(u.normalize(), "Vector normalize returns true");
    require_near(u.norm(), 1.0, "Vector normalize makes unit length");
  }
}

void test_vector_ops_wrapped()
{
  WrappedVector<2> const v{3.0, 4.0};
  require_near(v.norm(), 5.0, "WrappedVector norm");
  require_near(v.norm_squared(), 25.0, "WrappedVector norm_squared");

  WrappedVector<2> const w{1.0, 2.0};
  require_near(v.dot(w), 11.0, "WrappedVector dot");
  require_near(v.cross(w), 2.0, "WrappedVector cross (2D)");

  require_vector_near<2>(v.rotate_90_degrees(), {-4.0, 3.0}, "WrappedVector rotate_90_degrees");

  {
    auto u = v;
    u += w;
    require_vector_near<2>(u, {4.0, 6.0}, "WrappedVector +=");
    u -= w;
    require_vector_near<2>(u, {3.0, 4.0}, "WrappedVector -=");
    u *= 2.0;
    require_vector_near<2>(u, {6.0, 8.0}, "WrappedVector *=");
    u /= 2.0;
    require_vector_near<2>(u, {3.0, 4.0}, "WrappedVector /=");
  }

  {
    auto u = v;
    u.negate();
    require_vector_near<2>(u, {-3.0, -4.0}, "WrappedVector negate()");
    require_vector_near<2>(-u, {3.0, 4.0}, "WrappedVector unary -");
  }

  require_vector_near<2>(v / 2.0, {1.5, 2.0}, "WrappedVector / scalar");

  auto d = v.direction();
  WrappedDirection<2> dx{WrappedPoint<2>{0.0, 0.0}, WrappedPoint<2>{1.0, 0.0}};
  require_near(d.raw().dot(dx.raw()), 3.0 / 5.0, "WrappedVector direction()");

  require_point_near<2>(v.as_point(), {3.0, 4.0}, "WrappedVector as_point()");

  {
    auto u = v;
    require(u.normalize(), "WrappedVector normalize returns true");
    require_near(u.norm(), 1.0, "WrappedVector normalize makes unit length");
  }
}

void test_matrix_ops_math()
{
  math::Matrix<2, 2> const I{1};        // identity
  require_near(I.determinant(), 1.0, "Matrix determinant");

  math::Matrix<2, 2> A;
  A(0, 0) = 1.0; A(0, 1) = 2.0;
  A(1, 0) = 3.0; A(1, 1) = 4.0;

  require_near(A.determinant(), -2.0, "Matrix determinant (non-identity)");

  auto const Ainv = A.inverse();
  auto const should_be_I = A * Ainv;
  require_near(should_be_I(0, 0), 1.0, "Matrix inverse check (0,0)");
  require_near(should_be_I(1, 1), 1.0, "Matrix inverse check (1,1)");

  auto const At = A.transpose();
  require_near(At(0, 1), 3.0, "Matrix transpose (0,1)");
  require_near(At(1, 0), 2.0, "Matrix transpose (1,0)");

  math::Vector<2> const v{5.0, 6.0};
  auto const Av = A * v;
  require_vector_near<2>(Av, {17.0, 39.0}, "Matrix * Vector (member)");

  auto const vA = v * A;
  require_vector_near<2>(vA, {23.0, 34.0}, "Vector * Matrix (free)");

  auto const twoA = 2.0 * A;
  require_near(twoA(1, 1), 8.0, "scalar * Matrix (free)");

  {
    auto B = A;
    B += A;
    require_near(B(1, 1), 8.0, "Matrix +=");
    B -= A;
    require_near(B(1, 1), 4.0, "Matrix -=");
    B *= 2.0;
    require_near(B(1, 1), 8.0, "Matrix *=");
    B /= 2.0;
    require_near(B(1, 1), 4.0, "Matrix /=");
    B.negate();
    require_near(B(1, 1), -4.0, "Matrix negate()");
    require_near((-B)(1, 1), 4.0, "Matrix unary -");
  }

  auto const ApA = A + A;
  require_near(ApA(1, 1), 8.0, "Matrix +");
  auto const AmA = A - A;
  require_near(AmA(0, 0), 0.0, "Matrix -");
  auto const As2 = A * 2.0;
  require_near(As2(0, 1), 4.0, "Matrix * scalar");
  auto const Ad2 = A / 2.0;
  require_near(Ad2(1, 0), 1.5, "Matrix / scalar");
}

void test_matrix_ops_wrapped()
{
  WrappedMatrix<2, 2> A;
  A(0, 0) = 1.0; A(0, 1) = 2.0;
  A(1, 0) = 3.0; A(1, 1) = 4.0;

  require_near(A.determinant(), -2.0, "WrappedMatrix determinant");

  auto At = A.transpose();
  require_near(At(0, 1), 3.0, "WrappedMatrix transpose (0,1)");

  WrappedVector<2> v{5.0, 6.0};
  auto Av = A * v;
  require_vector_near<2>(Av, {17.0, 39.0}, "WrappedMatrix * WrappedVector (generic ops)");

  auto vA = v * A;
  require_vector_near<2>(vA, {23.0, 34.0}, "WrappedVector * WrappedMatrix (free)");

  auto twoA = 2.0 * A;
  require_near(twoA(1, 1), 8.0, "scalar * WrappedMatrix (free)");
}

} // namespace

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  try
  {
    test_point_ops_math();
    test_point_ops_wrapped();
    test_vector_ops_math();
    test_vector_ops_wrapped();
    test_matrix_ops_math();
    test_matrix_ops_wrapped();
  }
  catch (std::exception const& error)
  {
    std::cerr << "linear_algebra_test failed: " << error.what() << '\n';
    return 1;
  }

  std::cout << "linear_algebra_test passed\n";
  return 0;
}
