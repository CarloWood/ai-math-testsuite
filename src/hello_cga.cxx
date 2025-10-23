#include <iostream>
#include <vsr/space/vsr_cga3D.h>
#include <vsr/space/vsr_cga3D_op.h>     // Gen::rotor / translator / motor, etc.

using namespace vsr::cga;

int main ()
{
  // Make a Euclidean direction and a conformal point.
  Vec const ex{1.0, 0.0, 0.0};
  Pnt const p = Round::null(ex * 2.0);     // Same as PT(2,0,0).  A conformal point.

  // 90° rotation in the xy-plane: bivector e12 times angle/2.
  Biv const B = Biv::xy * (M_PI * 0.25);
  Rot const R = Gen::rotor(B);             // Euclidean rotor from a bivector.

  // Translator by (0, 1, 0).
  Trs const T = Gen::translator(Vec{0.0, 1.0, 0.0});

  // A motor = translation then rotation (compose as you like).
  Mot const M = R * T;

  // Apply motor to a point (spin / sandwich).
  Pnt const q = p.sp(M);                    // Equivalent to M * p * ~M.

  // Extract Euclidean coordinates back from a conformal point.
  Vec const q_xyz = Round::location(q);     // Returns the 3D location.
  std::cout << q_xyz << '\n';
}
