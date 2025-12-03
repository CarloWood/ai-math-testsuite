#include "sys.h"
#include "math/Hyperblock.h"
#include "math/Hyperplane.h"
#include "math/Vector.h"
#include "utils/almost_equal.h"
#include "utils/RandomNumber.h"
#include "utils/square.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "debug.h"

#include "cwds/debug_ostream_operators.h"
#include "utils/debug_ostream_operators.h"

namespace {

constexpr double kTolerance = 1e-9;

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

template<int N>
using vector_type = math::Vector<N>;

template<int N>
vector_type<N> vector_from_corners(vector_type<N> const& c1,
                                   vector_type<N> const& c2,
                                   std::size_t mask)
{
  vector_type<N> result = c1;
  for (int d = 0; d < N; ++d)
  {
    if (mask & (1UL << d))
      result[d] = c2[d];
  }
  return result;
}

template<int N>
math::Hyperplane<N> make_plane(vector_type<N> const& normal, vector_type<N> const& point)
{
  return {normal, -normal.dot(point)};
}

template<int N>
std::string vector_to_string(vector_type<N> const& v)
{
  std::ostringstream oss;
  oss << '[';
  for (int i = 0; i < N; ++i)
  {
    if (i != 0)
      oss << ',';
    oss << v[i];
  }
  oss << ']';
  return oss.str();
}

template<int N>
void require_vector_near(vector_type<N> const& actual,
                         vector_type<N> const& expected,
                         std::string const& context)
{
  double difference = (actual - expected).norm();
  if (difference > kTolerance)
  {
    std::ostringstream oss;
    oss << context << ": expected " << vector_to_string(expected)
        << ", got " << vector_to_string(actual)
        << ", difference=" << difference;
    throw std::runtime_error(oss.str());
  }
}

template<int N>
void require_same_points(std::vector<vector_type<N>> actual,
                         std::vector<vector_type<N>> expected,
                         std::string const& context)
{
  require(actual.size() == expected.size(),
          context + ": expected " + std::to_string(expected.size()) +
              " intersections, got " + std::to_string(actual.size()));
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    require_vector_near(actual[i], expected[i], context + ", index=" + std::to_string(i));
  }
}

template<typename Scalar>
struct Entry
{
  std::size_t index;
  Scalar      key;
};

template<typename Scalar>
std::ostream& operator<<(std::ostream& os, Entry<Scalar> const& entry)
{
  os << "{index:" << entry.index << ", key:" << entry.key << "}";
  return os;
}

template<int N>
std::vector<vector_type<N>> reference_intersections(vector_type<N> const& c1,
                                                    vector_type<N> const& c2,
                                                    math::Hyperplane<N> const& plane)
{
  std::vector<vector_type<N>> expected;
  for (std::size_t mask = 0; mask < static_cast<std::size_t>(1UL << N); ++mask)
  {
    for (int d = 0; d < N; ++d)
    {
      if (mask & (1UL << d))
        continue;
      vector_type<N> const corner1 = vector_from_corners(c1, c2, mask);
      vector_type<N> const corner2 = vector_from_corners(c1, c2, mask | (1UL << d));
      if (plane.side(corner1) != plane.side(corner2))
        expected.push_back(plane.intersection(corner1, corner2));
    }
  }

  // Order intersection points starting from expected[0],
  // then around the plane, by angle.
  if (expected.size() > 2)
  {
    using Vec    = vector_type<N>;
    using Scalar = typename Vec::scalar_type;

    // Project onto the plane and compute centroid.
    std::vector<Vec> proj(expected.size());
    Vec centroid;
    centroid.eigen().setZero();
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
      proj[i] = plane.project(expected[i]);
      centroid += proj[i];
    }
    centroid /= static_cast<Scalar>(expected.size());

    // Build an orthonormal basis (t1, t2) in the plane.
    Vec n      = plane.normal();
    Scalar nn  = n.norm();
    Vec n_hat  = n / nn;

    // First basis vector: direction from centroid to first point, projected into the plane.
    Vec t1 = proj[0] - centroid;
    // Make sure t1 is orthogonal to n_hat (for safety in higher dimensions).
    t1 -= t1.dot(n_hat) * n_hat;
    Scalar t1n = t1.norm();
    if (t1n == Scalar(0))
      return {};                // Plane must be touching a corner; return nothing.
    t1 /= t1n;

    Vec t2;
    {
      Vec n = plane.normal();
      Scalar n_norm = n.norm();
      if (n_norm == Scalar(0))
        return expected;
      Vec n_hat = n / n_norm;

      if constexpr (N == 3)
      {
        // Use a cross product to fix orientation:
        // t2 = t1 × n_hat gives angle -pi/2 w.r.t. t1 under atan2(dot(v,t2), dot(v,t1)).
        Vec tmp = t1.cross(n_hat);
        Scalar tmp_norm = tmp.norm();
        if (tmp_norm == Scalar(0))
          return expected;          // degenerate, don't reorder
        t2 = tmp / tmp_norm;
      }
      else
      {
        // Generic fallback for N != 3: Gram–Schmidt from a coordinate axis.
        bool found = false;
        for (int k = 0; k < N && !found; ++k)
        {
          Vec v;
          v.eigen().setZero();
          v[k] = Scalar(1);
          v -= v.dot(n_hat) * n_hat;
          v -= v.dot(t1)    * t1;
          Scalar vn = v.norm();
          if (vn != Scalar(0))
          {
            t2 = v / vn;
            found = true;
          }
        }
        if (!found)
          return expected;          // truly degenerate
      }
    }

    std::vector<Entry<Scalar>> entries;
    entries.reserve(expected.size());

    // Compute raw angles for all points.
    for (std::size_t i = 0; i < proj.size(); ++i)
    {
      Vec v = proj[i] - centroid;
      Scalar x = v.dot(t1);
      Scalar y = v.dot(t2);
      Scalar angle = std::atan2(y, x);
      entries.push_back({i, angle});
    }

    // Normalize angles so that entries[0] becomes angle 0,
    // and sort others by increasing angle.
    Scalar two_pi = static_cast<Scalar>(2.0 * std::acos(-1.0));
    Scalar a0     = entries[0].key;

    for (auto &e : entries)
    {
      Scalar delta = e.key - a0;
      while (delta > two_pi)
        delta -= two_pi;
      while (delta <= Scalar(0))
        delta += two_pi;

      // Go CW around the normal.
      e.key = two_pi - delta;
    }

    // Keep the first point fixed; sort the rest by angle.
    std::stable_sort(entries.begin(), entries.end(), [](Entry<Scalar> const& a, Entry<Scalar> const& b) { return a.key < b.key; });

    // Remove near-duplicate angles, keeping the first occurrence.
    auto new_end = std::unique(entries.begin(), entries.end(), [](Entry<Scalar> const& a, Entry<Scalar> const& b) { return utils::almost_equal(a.key, b.key, 10e-6); });
    entries.erase(new_end, entries.end());

    // Build reordered vector; size is identical to expected.size().
    std::vector<Vec> ordered;
    ordered.reserve(expected.size());
    for (auto const& e : entries)
      ordered.push_back(expected[e.index]);

    expected.swap(ordered);
  }

  return expected;
}

void axis_aligned_slice_test()
{
  vector_type<3> const c1{0.0, 0.0, 0.0};
  vector_type<3> const c2{2.0, 3.0, 4.0};
  math::Hyperblock<3> block{c1, c2};
  math::Hyperplane<3> const plane = make_plane<3>({1.0, 0.0, 0.0}, {1.0, 0.0, 0.0});

  std::vector<vector_type<3>> const expected{
      {1.0, 0.0, 0.0},
      {1.0, 3.0, 0.0},
      {1.0, 3.0, 4.0},
      {1.0, 0.0, 4.0}};

  auto const intersections = block.intersection_points(plane);
  require_same_points(intersections, expected, "axis_aligned_slice_test intersections");
}

void diagonal_slice_test()
{
  vector_type<3> const c1{0.0, 0.0, 0.0};
  vector_type<3> const c2{2.0, 2.0, 2.0};
  math::Hyperblock<3> block{c1, c2};
  math::Hyperplane<3> const plane = make_plane<3>({1.0, 1.0, 1.0}, {1.0, 1.0, 1.0});

  std::vector<vector_type<3>> const expected{
      {2.0, 1.0, 0.0},
      {1.0, 2.0, 0.0},
      {0.0, 2.0, 1.0},
      {0.0, 1.0, 2.0},
      {1.0, 0.0, 2.0},
      {2.0, 0.0, 1.0}};

  auto const intersections = block.intersection_points(plane);
  require_same_points(intersections, expected, "diagonal_slice_test intersections");
}

void disjoint_plane_test()
{
  vector_type<3> const c1{0.0, 0.0, 0.0};
  vector_type<3> const c2{2.0, 2.0, 2.0};
  math::Hyperblock<3> block{c1, c2};
  math::Hyperplane<3> const plane = make_plane<3>({1.0, 0.0, 0.0}, {5.0, 0.0, 0.0});

  auto const intersections = block.intersection_points(plane);
  require(intersections.empty(), "disjoint_plane_test intersections should be empty");
}

void vertex_touching_test()
{
  vector_type<3> const c1{0.0, 0.0, 0.0};
  vector_type<3> const c2{1.0, 1.0, 1.0};
  math::Hyperblock<3> block{c1, c2};
  math::Hyperplane<3> const plane = make_plane<3>({1.0, 1.0, 1.0}, {0.0, 0.0, 0.0});

  // An in-plane corner is considered to not cause any intersections if the plane is only touching.
  std::vector<vector_type<3>> const expected{};

  auto const intersections = block.intersection_points(plane);
  require_same_points(intersections, expected, "vertex_touching_test intersections");
}

void general_reference_test()
{
  vector_type<3> const c1{-1.5, 2.0, -0.5};
  vector_type<3> const c2{3.0, 4.5, 2.5};
  math::Hyperblock<3> block{c1, c2};

  struct PlaneSpec
  {
    vector_type<3> normal;
    vector_type<3> point;
  };

  std::vector<PlaneSpec> planes = {
      {{0.75, -1.0, 0.5}, {0.0, 3.0, 0.0}},
      {{-0.5, 0.2, 1.0}, {1.2, 2.2, 1.0}},
      {{1.0, 1.5, 0.25}, {0.5, 3.5, 0.5}},
      {{-0.75, -0.25, 0.6}, {2.5, 2.0, 0.1}}};

  utils::RandomNumber rng(42);
  std::uniform_real_distribution<> dist_0_1(0, 1);

  auto random_point = [&]() -> vector_type<3> {
    vector_type<3> delta = c2 - c1;
    for (int d = 0; d < 3; ++d)
      delta[d] *= rng.generate(dist_0_1);
    return c1 + delta;
  };

  std::uniform_real_distribution<> dist_phi(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<> dist_cos_theta(-1.0, 1.0);

  auto random_normal = [&]() -> vector_type<3> {
    double phi = rng.generate(dist_phi);
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    double cos_theta = rng.generate(dist_cos_theta);
    double sin_theta = std::sqrt(1 - utils::square(cos_theta));
    double x = sin_phi * cos_theta;
    double y = sin_phi * sin_theta;
    double z = cos_phi;
    return {x, y, z};
  };

  for (int i = 0; i < 32; ++i)
    planes.emplace_back(random_normal(), random_point());

  for (std::size_t i = 0; i < planes.size(); ++i)
  {
    math::Hyperplane<3> const plane = make_plane<3>(planes[i].normal, planes[i].point);
    auto const intersections = block.intersection_points(plane);
    auto const expected = reference_intersections<3>(c1, c2, plane);
    require_same_points(intersections, expected, "general_reference_test plane=" + std::to_string(i));
  }
}

} // namespace

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  try
  {
    axis_aligned_slice_test();
    diagonal_slice_test();
    disjoint_plane_test();
    vertex_touching_test();
    general_reference_test();
  }
  catch (std::exception const& error)
  {
    std::cerr << "Hyperblock tests failed: " << error.what() << '\n';
    return 1;
  }

  std::cout << "Hyperblock tests passed" << std::endl;
  return 0;
}

