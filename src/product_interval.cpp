#include <iostream>
#include <cassert>
#include <array>
#include <algorithm>
#include <ranges>

//         |                  |                  |                 |
//         |     A_M < 0      | A_m < 0  A_M > 0 |     A_m > 0     |
//         |                  |                  |                 |
// --------+------------------+------------------+-----------------+
//         |                  |                  |                 |
// B_m > 0 |   M = A_M * B_m  |   M = A_M * B_M  |  M = A_M * B_M  |
//         |   m = A_m * B_M  |   m = A_m * B_M  |  m = A_m * B_m  |
// --------+------------------+------------------+-----------------+
//         |                  | \              / |                 |
// B_m < 0 |                  |   \          /   |                 |
//         |                  |     \      /     |                 |
//         |  M = A_m * B_m   |       \  /       |  M = A_M * B_M  |
//         |  m = A_m * B_M   |       /  \       |  m = A_M * B_m  |
// B_M > 0 |                  |     /      \     |                 |
//         |                  |   /          \   |                 |
//         |                  | /              \ |                 |
// --------+------------------+------------------+-----------------+
//         |                  |                  |                 |
// B_M < 0 |   M = A_m * B_m  |   M = A_m * B_m  |  M = A_m * B_M  |
//         |   m = A_M * B_M  |   m = A_M * B_m  |  m = A_M * B_m  |

std::array<int, 4> products(int Am, int AM, int Bm, int BM)
{
  return { Am * Bm, Am * BM, AM * Bm, AM * BM };
}

int main()
{
  for (int A_m = -9; A_m < 9; ++A_m)
    for (int A_M = A_m + 1; A_M <= 9; ++A_M)
      for (int B_m = -9; B_m < 9; ++B_m)
        for (int B_M = B_m + 1; B_M <= 9; ++B_M)
        {
          auto p = products(A_m, A_M, B_m, B_M);
          int m = *std::ranges::min_element(p);
          int M = *std::ranges::max_element(p);

          if (B_m >= 0)         // Also B_M >= 0. The center also works if B_m == 0, but this block is cheaper.
          {
            if (A_m >= 0)       // Also A_M >= 0.
                                        //   BM>=0 ----.
                                        //   Bm>=0 ---.|
                                        //   AM>=0 --.||
            {                           //   Am>=0 -.|||
              assert(m == A_m * B_m);   // 1x       1111
              assert(M == A_M * B_M);   // 3x       1111
            }
            else if (A_M < 0)
            {
              assert(m == A_m * B_M);   // 3x       0011
              assert(M == A_M * B_m);   // 1x       0011
            }
            else // A_m < 0 && A_M > 0.
            {
              assert(m == A_m * B_M);   // 3x       0111
              assert(M == A_M * B_M);   // 3x       0111
            }
          }
          else if (B_M <= 0)    // The center also works if B_M == 0, but this block is cheaper.
          {
            if (A_m >= 0)
                                        //   BM>0  ----.
                                        //   Bm>=0 ---.|
                                        //   AM>=0 --.||
            {                           //   Am>=0 -.|||
              assert(m == A_M * B_m);   // 3x       1100
              assert(M == A_m * B_M);   // 1x       1100
            }
            else if (A_M <= 0)
            {
              assert(m == A_M * B_M);   // 1x       0000
              assert(M == A_m * B_m);   // 3x       0000
            }
            else // A_m < 0 && A_M > 0
            {
              assert(m == A_M * B_m);   // 3x       0100
              assert(M == A_m * B_m);   // 3x       0100
            }
          }
          else // B_m < 0 && B_M > 0
          {
            if (A_m >= 0)       // The center also works if A_m == 0, but this block is cheaper.
                                        //   BM>0  ----.
                                        //   Bm>=0 ---.|
                                        //   AM>=0 --.||
            {                           //   Am>=0 -.|||
              assert(m == A_M * B_m);   // 3x       1101
              assert(M == A_M * B_M);   // 3x       1101
            }
            else if (A_M <= 0)  // The center also works if A_M == 0, but this block is cheaper.
            {
              assert(m == A_m * B_M);   // 3x       0001
              assert(M == A_m * B_m);   // 3x       0001
            }
            else // A_m < 0 && A_M > 0
            {
                                                                //   BM>0  ----.
                                                                //   Bm>=0 ---.|
                                                                //   AM>0  --.||
              // Center.                                        //   Am>=0 -.|||
              assert(m == std::min(A_m * B_M, A_M * B_m));      // 1x       0101
              assert(M == std::max(A_M * B_M, A_m * B_m));      // 1x       0101
            }
          }
        }
}
