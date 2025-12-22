#include "sys.h"
#include "math/Vector.h"
#include "math/kFace.h"

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  constexpr int n = 7;
  constexpr int k = 3;
  using k_face_index = math::kFaceIndex<n, k>;
  Dout(dc::notice, "Number of " << k << "-faces of an " << n << "-cube: " << k_face_index::size);
  Dout(dc::notice, "The number of facets of such a " << k << "-face: " << k_face_index::number_of_facets);

  math::Vector<n> const origin(-1);
  math::Vector<n> const far_corner(1, 1, 1, 1, 1, 1, 1);
  math::CornerIndex<n> const ci{0b1001101};
  k_face_index::axes_type::pod_type k_axes{0b0101010};
  k_face_index kfr({k_axes, ci});

  Dout(dc::notice, "The facets of the k-face " << kfr << " aka " << std::format("{:b}", kfr.get_value()) << " aka " << kfr.as_kface() << " are:");
  auto facet_indexes = kfr.facet_indexes();
  for (auto&& facet : facet_indexes)
    Dout(dc::notice, facet << " aka " << facet.get_value() << " aka " << std::format("{:b}", facet.get_value()) << " = " << facet.as_kface());

  using namespace utils::bitset;
  IndexPOD const index_end{n};
  for (Index d0 = index_begin; d0 != index_end; ++d0)
    for (Index d1 = d0 + 1; d1 != index_end; ++d1)
      for (Index d2 = d1 + 1; d2 != index_end; ++d2)
      {
        std::array<Index, k> face_axes = { d0, d1, d2 };
        k_face_index::axes_type k_axes;
        k_axes.reset();
        for (int i = 0; i < k; ++i)
          k_axes.set(face_axes[i]);
        uint32_t r = k_face_index::rank(k_axes);
        Dout(dc::notice, "rank(" << face_axes << " (" << k_axes << ")) = " << r);
        ASSERT((k_face_index::unrank(r) == k_axes));
      }
  Dout(dc::notice, n << " choose " << k << " = " << math::binomial(n, k));
}
