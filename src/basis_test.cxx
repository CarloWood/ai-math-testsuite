#include "sys.h"
#include "math/Universe.h"
#include "debug.h"

namespace universe {
enum A {};
struct B {};
} // namespace universe

using UniverseA = math::Universe<universe::A, 32>;      // Orthonormal basis: { e1, e2, e3, ... e32 }
using UniverseB = math::Universe<universe::B, 5>;       // Orthonormal basis: { e1, e2, e3, ... e5 }

constexpr int playground7 = 7;                          // Use 7 dimensions for 'playground7'.

template<size_t N>
using Permutation = math::Permutation<N>;

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  //Similitude sim;

  // Construct a 7-dimensional Coordinate-Subspace in UniverseA.
  using playground7_coordinate_subspace_basis = UniverseA::CoordinateSubspace::basis<playground7>;

  Dout(dc::notice, "UniverseA::basis_type = " << debug::type_name_of<UniverseA::basis_type>());
  Dout(dc::notice, "playground7_coordinate_subspace_basis = " << debug::type_name_of<playground7_coordinate_subspace_basis>());
  Dout(dc::notice, "with number of dimensions: " << playground7_coordinate_subspace_basis::n);

  // Construct a basis for playground3.
  auto basis = UniverseB::CoordinateSubspace::from_permutation(Permutation{0, 4, 3}, 1, 2);
  Dout(dc::notice, "basis = " << basis);

  Dout(dc::notice, "Leaving main()");
}
