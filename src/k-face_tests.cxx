#include "sys.h"
#include "utils/VectorIndex.h"
#include "utils/uint_leastN_t.h"
#include "utils/BitSet.h"
#include "utils/create_mask.h"
#include "math/Vector.h"
#include "math/binomial.h"
#include "debug.h"

namespace math {
namespace category {

template<int n, int k> struct kFace;

} // namespace category

// n : number of dimensions of the hypercube.
// k : number of dimensions of the encoded k-face.
//
// This encodes a k-face as a bitmask that can be used as an index into a Vector.
// The bitmask exists of a 'rank' for the k axes that the k-face aligns with,
// followed by n-k 'fixed bits' for the other axes:
//
//     <rrr><fffffff>
//      ^^^  ^^^^^^^
//
// The fixed bits (<fffffff>) are those bits from the coordinate mask of each
// corner of the k-face that do not change: the bits that belong the one of the
// k axis (that are thus not fixed) are removed.
//
// The rank bits (<rrr>) encode the lexiographic rank of the k axes that the k-face
// aligns with, as returned by the member function `rank`.
//
// Note that a 0-face is a corner (k=0), a 1-face is an edge (k=1).
//
template<int n, int k>
class kFaceRank : public utils::VectorIndex<category::kFace<n, k>>
{
 public:
  // A bitset that can hold at least n bits, where each bit represents an axis.
  using axes_type = utils::BitSet<utils::uint_leastN_t<n>>;
  static constexpr std::size_t fixed_mask = utils::create_mask<std::size_t, n - k>();

  // The number of k-faces of an n-cube.
  static constexpr size_t size = static_cast<size_t>(binomial(n, k)) << (n - k);
  // The number of (k-1)-faces of a given k-face.
  static constexpr size_t number_of_facets = 2 * k;

  static uint32_t rank(axes_type k_axes);
  static axes_type unrank(uint32_t r);

  // FIXME: shouldn't this be deleted?
  kFaceRank() = default;

  // Construct a corner.
  kFaceRank(axes_type coordinate_mask) requires (k == 0) : utils::VectorIndex<category::kFace<n, k>>(coordinate_mask())
  {
  }

  kFaceRank(axes_type k_axis, kFaceRank<n, 0> corner) : utils::VectorIndex<category::kFace<n, k>>(static_cast<size_t>(rank(k_axis) << (n - k) | corner.remove_bits(k_axis)))
  {
  }

  std::array<kFaceRank<n, k - 1>, number_of_facets> facet_indexes();

 public: // FIXME: should be private
  int remove_bits(axes_type bits) requires (k == 0)
  {
    ASSERT(!this->undefined());

    using mask_type = typename axes_type::mask_type;

    mask_type const mask = static_cast<mask_type>(this->get_value());
    mask_type const b = bits();

    mask_type result = 0;
    int removed = 0;

    // Walk over all n coordinate bits; skip the ones present in 'b' and
    // compact the remaining bits downward.
    for (int i = 0; i < n; ++i)
    {
      mask_type const bit = static_cast<mask_type>(mask_type{1} << i);
      if (b & bit)
      {
        ++removed;
      }
      else if (mask & bit)
      {
        result |= static_cast<mask_type>(mask_type{1} << (i - removed));
      }
    }

    return static_cast<int>(result);
  }
};

template<int n, int k>
std::array<kFaceRank<n, k - 1>, kFaceRank<n, k>::number_of_facets>
kFaceRank<n, k>::facet_indexes()
{
  ASSERT(!this->undefined());

  uint32_t const mask = this->get_value();
  uint32_t const rrr = mask >> (n - k);
  uint32_t const fff = mask & fixed_mask;
  axes_type k_axes = unrank(rrr);

  Dout(dc::notice, "k_axes = " << k_axes);

  return {};
}

// Returns the rank for this choice of axes for the set of all possible ways one can choose k axes out of n.
//static
template<int n, int k>
uint32_t kFaceRank<n, k>::rank(axes_type k_axes)
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
  uint32_t sum = 0;

  int mi = 0;
  Index d_i = index_pre_begin;
  for (int i = 1; i <= k; ++i)
  {
    // Advance d_i to the first/next axis index.
    d_i.next_bit_in(k_axes);
    for (int x = mi; x <= d_i() - 1; ++x)
      sum += binomial(n - 1 - x, k - i);
    mi = d_i() + 1;
  }
  return sum;
}

// Given a rank in [0, C(n,k)) return the corresponding axes BitSet.
//static
template<int n, int k>
typename kFaceRank<n, k>::axes_type kFaceRank<n, k>::unrank(uint32_t r)
{
  using namespace utils::bitset;
  // Index of the least significant bit.
  constexpr IndexPOD ilsb = { 0 };
  // Initialize the result variable.
  axes_type axes{axes_type::zero};

  // Our binomials are never negative (n and k are positive).
  ASSERT(r < (uint32_t)binomial(n, k));

  Index mi = ilsb;
  for (int i = 1; i <= k; ++i)
  {
    Index d = mi;
    // Find smallest d such that the block size falls past r.
    for (;; ++d)
    {
      uint32_t const block = binomial(n - 1 - d(), k - i);
      if (r < block)
        break;
      r -= block;
    }
    axes.set(d);
    mi = d + 1;
  }

  return axes;
}

} // namespace math

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

  math::Vector<n> const origin(-1);
  math::Vector<n> const far_corner(1, 1, 1, 1, 1, 1, 1);
  using CornerIndex = math::kFaceRank<n, 0>;
  k_face_rank::axes_type::pod_type const c1_pod{0b1001101};
  CornerIndex const c1{c1_pod};
  k_face_rank::axes_type::pod_type k_axes{0b0101010};
  k_face_rank kfr{k_axes, c1};

  //101001011
  //  256 + 64 + 8 + 3 =

  Dout(dc::notice, "kfr = " << kfr);
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
        uint32_t r = k_face_rank::rank(k_axes);
        Dout(dc::notice, "rank(" << face_axes << " (" << k_axes << ")) = " << r);
        ASSERT((k_face_rank::unrank(r) == k_axes));
      }
  Dout(dc::notice, n << " choose " << k << " = " << math::binomial(n, k));
}
