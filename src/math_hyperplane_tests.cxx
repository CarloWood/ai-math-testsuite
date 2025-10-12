#include "sys.h"
#include "math/Hyperplane.h"
#include "math/Vector.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include "debug.h"

namespace {

constexpr double kTolerance = 1e-9;
constexpr std::array<double, 8> kCoordinates = {
    -7.897572, -3.1415, -2.0, -1.0, 0.0, 1.0, 2.0, 8.947193};
constexpr std::array<double, 5> kSignedDistances = {-100.12, -2.0, -1.0, 0.0, 4.58374};

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
vector_type<N> make_zero()
{
  vector_type<N> result;
  result.eigen().setZero();
  return result;
}

template<int N>
std::vector<vector_type<N>> generate_vectors(bool exclude_zero)
{
  std::vector<vector_type<N>> vectors;
  std::size_t combinations = 1;
  for (int i = 0; i < N; ++i)
    combinations *= kCoordinates.size();
  vectors.reserve(combinations);
  for (std::size_t index = 0; index < combinations; ++index)
  {
    std::size_t remainder = index;
    vector_type<N> candidate = make_zero<N>();
    bool all_zero = true;
    for (int dim = 0; dim < N; ++dim)
    {
      double value = kCoordinates[remainder % kCoordinates.size()];
      remainder /= kCoordinates.size();
      candidate[dim] = value;
      all_zero &= value == 0.0;
    }
    if (!exclude_zero || !all_zero)
      vectors.push_back(candidate);
  }
  return vectors;
}

template<int N>
std::vector<vector_type<N>> select_subset(std::vector<vector_type<N>> const& source, std::size_t maximum)
{
  if (source.size() <= maximum)
    return source;

  std::vector<vector_type<N>> subset;
  subset.reserve(maximum);
  double step = static_cast<double>(source.size()) / static_cast<double>(maximum);
  double cursor = 0.0;
  for (std::size_t i = 0; i < maximum; ++i)
  {
    std::size_t index = static_cast<std::size_t>(cursor);
    if (index >= source.size())
      index = source.size() - 1;
    subset.push_back(source[index]);
    cursor += step;
  }
  return subset;
}

template<int N>
vector_type<N> normalized(vector_type<N> v)
{
  double const norm = v.norm();
  require(norm > 0.0, "Attempted to normalize a zero vector");
  v /= norm;
  return v;
}

template<int N>
std::string description(int normal_index, int point_index, std::string_view details)
{
  std::ostringstream oss;
  oss << "dimension=" << N << ", normal_index=" << normal_index
      << ", point_index=" << point_index << ": " << details;
  return oss.str();
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
        << ", got " << vector_to_string(actual) << ", difference=" << difference;
    throw std::runtime_error(oss.str());
  }
}

template<int N>
std::vector<vector_type<N>> find_tangents(int n, vector_type<N> const& normal, std::vector<vector_type<N>> const& candidates)
{
  ASSERT(n < candidates.size());
  double const normal_norm_squared = normal.norm_squared();
  std::vector<vector_type<N>> tangents;
  for (auto const& candidate : candidates)
  {
    if (candidate.norm() == 0.0)
      continue;
    double const projection_length = candidate.dot(normal) / normal_norm_squared;
    vector_type<N> tangent = candidate - projection_length * normal;
    if (tangent.norm() > 1e-8)
    {
      tangents.push_back(tangent);
      if (tangents.size() == n)
        return tangents;
    }
  }

  ASSERT(false);
  AI_NEVER_REACHED
}

template<int N>
void run_dimension_tests()
{
  Dout(dc::notice, "Running run_dimension_tests<" << N << ">()...");

  auto const normals_all = generate_vectors<N>(true);
  auto const points_all = generate_vectors<N>(false);

  std::size_t kMaxNormals = std::min(normals_all.size(), 1000UL);
  std::size_t kMaxPoints = std::min(points_all.size(), 1000UL);

  auto const normals = select_subset<N>(normals_all, kMaxNormals);
  auto const points = select_subset<N>(points_all, kMaxPoints);

  vector_type<N> const zero = make_zero<N>();

  int normal_index = 0;
  for (auto const& normal : normals)
  {
    double const normal_norm = normal.norm();
    ASSERT(normal_norm > 0.0);

    vector_type<N> const unit_normal = normal / normal_norm;
    double const normal_norm_squared = normal_norm * normal_norm;

    int point_index = 0;
    for (auto const& base_point : points)
    {
      double const b = -normal.dot(base_point);
      math::Hyperplane<N> const plane{normal, b};

      auto const context = description<N>(normal_index, point_index, "base checks");

      require_near(plane.signed_distance(base_point), 0.0, context + " signed_distance(base_point)");
      require_near(plane.height_along_N(base_point), 0.0, context + " height_along_N(base_point)");
      require(math::in_plane == plane.side(base_point), context + " side(base_point)");
      require_vector_near<N>(plane.project(base_point), base_point, context + " project(base_point)");

      vector_type<N> const origin_projection = plane.project(zero);
      vector_type<N> const expected_origin_projection = (-b / normal_norm_squared) * normal;
      require_vector_near<N>(origin_projection,
                             expected_origin_projection,
                             description<N>(normal_index, point_index, "origin projection"));
      require_near(plane.signed_distance(origin_projection), 0.0,
                  description<N>(normal_index, point_index, "signed_distance(origin_projection)"));

      std::vector<vector_type<N>> const tangents = find_tangents<N>(10, normal, points_all);

      for (vector_type<N> const& tangent : tangents)
      {
        vector_type<N> const other_point = base_point + tangent;
        auto const other_context = description<N>(normal_index, point_index, "other point in plane");
        require_near(plane.signed_distance(other_point), 0.0, other_context + " signed_distance");
        require_vector_near<N>(plane.project(other_point), other_point, other_context + " project");

        for (double const distance : kSignedDistances)
        {
          vector_type<N> const shifted = base_point + distance * unit_normal;
          auto const shift_context = description<N>(normal_index, point_index, "shift along normal");

          require_near(plane.signed_distance(shifted), distance,
                       shift_context + " signed_distance");
          require_near(plane.height_along_N(shifted), distance / normal_norm,
                       shift_context + " height_along_N");
          require_vector_near<N>(plane.project(shifted), base_point,
                                  shift_context + " project onto base");

          math::Sign const expected_side =
              (std::abs(distance) < 1e-6) ? math::in_plane : (distance > 0.0 ? math::positive : math::negative);
          require(plane.side(shifted) == expected_side, shift_context + " side");
        }

        vector_type<N> const mixed_point = other_point + 3.25 * unit_normal + 1.5 * tangent;
        auto const mixed_context = description<N>(normal_index, point_index, "mixed offset");
        double const mixed_height = plane.height_along_N(mixed_point);
        vector_type<N> const projected_mixed = plane.project(mixed_point);
        require_vector_near<N>(projected_mixed, other_point + 1.5 * tangent, mixed_context + " project");
        require_near(mixed_height, 3.25 / normal_norm, mixed_context + " height_along_N");
        require_near(plane.signed_distance(mixed_point), 3.25, mixed_context + " signed_distance");
        require_vector_near<N>(projected_mixed + mixed_height * normal, mixed_point,
                                mixed_context + " reconstruct from projection");
        require_near(normal.dot(projected_mixed) + b, 0.0, mixed_context + " projection lies in plane");

        vector_type<N> const far_positive = base_point + 12.0 * unit_normal;
        vector_type<N> const far_negative = base_point - 7.0 * unit_normal;
        auto const intersection_context = description<N>(normal_index, point_index, "intersection");
        vector_type<N> const intersection = plane.intersection(far_negative, far_positive);
        require_vector_near<N>(intersection, base_point, intersection_context + " intersection point");
        require(plane.side(far_positive) == math::positive, intersection_context + " side positive");
        require(plane.side(far_negative) == math::negative, intersection_context + " side negative");
      }

      ++point_index;
    }
    ++normal_index;
  }
}

} // namespace

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  try
  {
    run_dimension_tests<2>();
    run_dimension_tests<3>();

    bool running_locally = std::getenv("LIBCWD_RCFILE_OVERRIDE_NAME");
    if (running_locally)
    {
      run_dimension_tests<4>();
      run_dimension_tests<7>();
    }
  }
  catch (std::exception const& error)
  {
    std::cerr << "Hyperplane tests failed: " << error.what() << '\n';
    return 1;
  }

  std::cout << "Hyperplane tests passed" << std::endl;
  return 0;
}
