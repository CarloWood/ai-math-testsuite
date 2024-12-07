include(AICxxProject)

find_library(GMP_LIBRARY gmp)
find_library(MPFR_LIBRARY mpfr)
find_library(GLPK_LIBRARY glpk REQUIRED)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(bracket_zero_test bracket_zero_test.cxx)
target_link_libraries(bracket_zero_test AICxx::math ${AICXX_OBJECTS_LIST})

add_executable(roots_test roots_test.cxx)
target_link_libraries(roots_test AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(Halley_iterations Halley_iterations.cxx)
target_link_libraries(Halley_iterations AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} ${GMP_LIBRARY} ${MPFR_LIBRARY}      ${GLPK_LIBRARY} Eigen3::Eigen)

add_executable(random_cubics random_cubics.cxx)
target_link_libraries(random_cubics AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(random_cubics2 random_cubics2.cxx)
target_link_libraries(random_cubics2 AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(print_roots print_roots.cxx)
target_link_libraries(print_roots AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} ${MPFR_LIBRARY} Eigen3::Eigen)

add_executable(cubic_roots cubic_roots.cxx)
target_link_libraries(cubic_roots AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST})

add_executable(cubic_roots2 cubic_roots2.cxx)
target_link_libraries(cubic_roots2 AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(IEEE754_test IEEE754_test.cxx)
target_link_libraries(IEEE754_test ${AICXX_OBJECTS_LIST} ${MPFR_LIBRARY})

if (GMP_LIBRARY AND MPFR_LIBRARY AND GLPK_LIBRARY)
add_executable(minus_three_case minus_three_case.cxx)
target_link_libraries(minus_three_case AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} ${GMP_LIBRARY} ${MPFR_LIBRARY}       ${GLPK_LIBRARY} Eigen3::Eigen)

add_executable(plus_three_case plus_three_case.cxx)
target_link_libraries(plus_three_case AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} ${GMP_LIBRARY} ${MPFR_LIBRARY}        ${GLPK_LIBRARY} Eigen3::Eigen)
endif ()

add_executable(small_root small_root.cxx)
target_link_libraries(small_root AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(convergence_test convergence_test.cxx)
target_link_libraries(convergence_test AICxx::cairowindow AICxx::symbolic AICxx::math ${AICXX_OBJECTS_LIST} ${MPFR_LIBRARY})
