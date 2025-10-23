include(AICxxProject)

find_package(Eigen3 3.3...<6.0 REQUIRED NO_MODULE)

find_library(GMP_LIBRARY gmp)
find_library(MPFR_LIBRARY mpfr)
find_library(GLPK_LIBRARY glpk)
find_package(versor REQUIRED CONFIG)

add_executable(bracket_zero_test bracket_zero_test.cxx)
target_link_libraries(bracket_zero_test AICxx::math ${AICXX_OBJECTS_LIST})

add_executable(math_geometry_tests math_geometry_tests.cxx)
target_link_libraries(math_geometry_tests AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(math_hyperplane_tests math_hyperplane_tests.cxx)
target_link_libraries(math_hyperplane_tests AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(math_hyperblock_tests math_hyperblock_tests.cxx)
target_link_libraries(math_hyperblock_tests AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(roots_test roots_test.cxx)
target_link_libraries(roots_test AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(hello_cga hello_cga.cxx)
target_link_libraries(hello_cga versor::vsr)

if (EnableCairoWindow)

add_executable(cubic_roots cubic_roots.cxx)
target_link_libraries(cubic_roots AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST})

add_executable(random_cubics random_cubics.cxx)
target_link_libraries(random_cubics AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(random_cubics2 random_cubics2.cxx)
target_link_libraries(random_cubics2 AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(cubic_roots2 cubic_roots2.cxx)
target_link_libraries(cubic_roots2 AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(small_root small_root.cxx)
target_link_libraries(small_root AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

endif (EnableCairoWindow)

if (MPFR_LIBRARY)
if (GMP_LIBRARY)
if (GLPK_LIBRARY AND EnableCairoWindow)

add_executable(Halley_iterations Halley_iterations.cxx)
target_link_libraries(Halley_iterations AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${GLPK_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen)

add_executable(minus_three_case minus_three_case.cxx)
target_link_libraries(minus_three_case AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${GLPK_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen)

add_executable(plus_three_case plus_three_case.cxx)
target_link_libraries(plus_three_case AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${GLPK_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen)

endif (GLPK_LIBRARY AND EnableCairoWindow)

add_executable(cubic_epsilons cubic_epsilons.cxx)
target_link_libraries(cubic_epsilons AICxx::math ${AICXX_OBJECTS_LIST}
  ${GMP_LIBRARY} ${MPFR_LIBRARY})

if (EnableCairoWindow)

add_executable(eta_correlation eta_correlation.cxx)
target_link_libraries(eta_correlation AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${GMP_LIBRARY} ${MPFR_LIBRARY})

add_executable(newton_epsilons newton_epsilons.cxx)
target_link_libraries(newton_epsilons AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${GMP_LIBRARY} ${MPFR_LIBRARY})

endif (EnableCairoWindow)

endif (GMP_LIBRARY)

add_executable(IEEE754_test IEEE754_test.cxx)
target_link_libraries(IEEE754_test ${AICXX_OBJECTS_LIST}
  ${MPFR_LIBRARY})

if (EnableCairoWindow)

add_executable(convergence_test convergence_test.cxx)
target_link_libraries(convergence_test AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${MPFR_LIBRARY})

add_executable(epsilon_sums epsilon_sums.cxx)
target_link_libraries(epsilon_sums AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${MPFR_LIBRARY})

add_executable(print_roots print_roots.cxx)
target_link_libraries(print_roots AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST}
  ${MPFR_LIBRARY} Eigen3::Eigen)

endif (EnableCairoWindow)

endif (MPFR_LIBRARY)
