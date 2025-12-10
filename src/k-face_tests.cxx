#include "sys.h"
#include "utils/VectorIndex.h"
#include "utils/uint_leastN_t.h"
#include "utils/BitSet.h"
#include "math/binomial.h"
#include "debug.h"

namespace math {
namespace category {

template<int n, int k> struct kFace;

} // namespace category

// n : number of dimensions of the hypercube.
// k : number of dimensions of the encoded k-face.
//
// Note that a corner has k=0, an edge has k=1...
//
template<int n, int k>
class kFaceRank : public utils::VectorIndex<category::kFace<n, k>>
{
 public:
  // A bitset that can hold at least n bits, where each bit represents an axis.
  using axes_type = utils::BitSet<utils::uint_leastN_t<n>>;

  // The number of k-faces of an n-cube.
  static constexpr size_t size = static_cast<size_t>(binomial(n, k)) << (n - k);
  // The number of (k-1)-faces of a given k-face.
  static constexpr size_t number_of_facets = 2 * k;

  std::array<kFaceRank<n, k - 1>, number_of_facets> facet_indexes()
  {
    return {};
  }
};

} // namespace math

// Returns the rank for this choice of axes for the set of all possible ways one can choose k axes out of n.
template<int n, int k>
int rank(typename math::kFaceRank<n, k>::axes_type k_axes)
{
  // There should be one bit set for each of the k axes.
  ASSERT(k_axes.count() == k);

  using namespace utils::bitset;

  // Calculate the rank according to:
  //
  //                         k  dᵢ-1 ⎛n-1-x⎞
  //   rank(d₁, d₂, …, dₖ) = ∑   ∑   ⎜     ⎟
  //                        i=1 x=mᵢ ⎝ k-i ⎠
  //
  // where m₁=0, mᵢ = dᵢ₋₁ + 1 for i>1.
  int sum = 0;

  int mi = 0;
  Index d_i = index_pre_begin;
  for (int i = 1; i <= k; ++i)
  {
    // Advance d_i to the first/next axis index.
    d_i.next_bit_in(k_axes);
    for (int x = mi; x <= d_i() - 1; ++x)
      sum += math::binomial(n - 1 - x, k - i);
    mi = d_i() + 1;
  }
  return sum;
}

#if 0
template<int n, int k>
int rank_old(std::array<int, k> const& k_axes)
{
  DoutEntering(dc::notice, "*** rank_old<" << n << ", " << k << ">(" << k_axes << ")");
  int sum = 0;
  for (int i = 1; i <= k; ++i)
  {
    Dout(dc::notice, "*** i=" << i << "; adding axis " << k_axes[i - 1]);
    int mi = i == 1 ? 0 : k_axes[i - 2] + 1;
    Dout(dc::notice, "*** mi = " << mi << "; running x from " << mi << " up till and including " << (k_axes[i - 1] - 1));
    for (int x = mi; x <= k_axes[i - 1] - 1; ++x)
    {
      Dout(dc::notice, "***  x = " << x);
      sum += math::binomial(n - 1 - x, k - i);
    }
  }
  return sum;
}
#endif

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  //    A9  6   2
  //    ||  |   |
  // 10100101110100
  //    ^^  ^   ^
  // 101  10 110 00
  //     1011011000
  //  ^^^
  //   \
  //    rank(2, 6, 9, 10) = ...
  //

  constexpr int n = 7;
  constexpr int k = 3;
  using k_face_rank = math::kFaceRank<n, k>;
  Dout(dc::notice, "Number of " << k << "-faces of an " << n << "-cube: " << k_face_rank::size);
  Dout(dc::notice, "The number of facets of such a " << k << "-face: " << k_face_rank::number_of_facets);

  k_face_rank kfr;
  Dout(dc::notice, "Indexes of the facets: " << kfr.facet_indexes());

  using namespace utils::bitset;
  IndexPOD const index_end{n};
  for (Index d0 = index_begin; d0 != index_end; ++d0)
    for (Index d1 = d0 + 1; d1 != index_end; ++d1)
      for (Index d2 = d1 + 1; d2 != index_end; ++d2)
      {
        std::array<Index, k> face_axes = { d0, d1, d2 };
        k_face_rank::axes_type k_axes;
        k_axes.reset();
        for (int i = 0; i < k; ++i)
          k_axes.set(face_axes[i]);
        Dout(dc::notice, "rank(" << face_axes << " (" << k_axes << ")) = " << (rank<n, k>(k_axes)));
      }
  Dout(dc::notice, n << " choose " << k << " = " << math::binomial(n, k));
}
