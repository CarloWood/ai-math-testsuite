#include "sys.h"
#include "math/Hyperblock.h"
#include "math/Hyperplane.h"
#include "math/Vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "debug.h"

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
void sort_vectors(std::vector<vector_type<N>>& vectors)
{
  std::sort(vectors.begin(), vectors.end(), [](vector_type<N> const& lhs, vector_type<N> const& rhs) {
    for (int i = 0; i < N; ++i)
    {
      if (lhs[i] < rhs[i] - kTolerance)
        return true;
      if (lhs[i] > rhs[i] + kTolerance)
        return false;
    }
    return false;
  });
}

template<int N>
void require_same_points(std::vector<vector_type<N>> actual,
                         std::vector<vector_type<N>> expected,
                         std::string const& context)
{
  sort_vectors(actual);
  sort_vectors(expected);
  require(actual.size() == expected.size(),
          context + ": expected " + std::to_string(expected.size()) +
              " intersections, got " + std::to_string(actual.size()));
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    require_vector_near(actual[i], expected[i], context + ", index=" + std::to_string(i));
  }
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
      {1.0, 0.0, 4.0},
      {1.0, 3.0, 0.0},
      {1.0, 3.0, 4.0}};

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
      {2.0, 0.0, 1.0},
      {1.0, 2.0, 0.0},
      {0.0, 2.0, 1.0},
      {1.0, 0.0, 2.0},
      {0.0, 1.0, 2.0}};

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

  std::vector<vector_type<3>> const expected{
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}};

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

  std::vector<PlaneSpec> const planes = {
      {{0.75, -1.0, 0.5}, {0.0, 3.0, 0.0}},
      {{-0.5, 0.2, 1.0}, {1.2, 2.2, 1.0}},
      {{1.0, 1.5, 0.25}, {0.5, 3.5, 0.5}},
      {{-0.75, -0.25, 0.6}, {2.5, 2.0, 0.1}}};

  for (std::size_t i = 0; i < planes.size(); ++i)
  {
    math::Hyperplane<3> const plane = make_plane<3>(planes[i].normal, planes[i].point);
    auto const intersections = block.intersection_points(plane);
    auto const expected = reference_intersections<3>(c1, c2, plane);
    require_same_points(intersections, expected,
                        "general_reference_test plane=" + std::to_string(i));
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

