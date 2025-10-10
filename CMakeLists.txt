cmake_minimum_required(VERSION 3.14...3.26)

project(ai-math-testsuite
  LANGUAGES CXX
  DESCRIPTION "Testsuite for git submodule math.")

# Compile all subdirectories with the same standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use gitache for libcwd.
set(GITACHE_PACKAGES libcwd_r)
include(cwm4/cmake/StableGitache)

# This project uses aicxx modules.
include(cwm4/cmake/AICxxProject)

# Exit if someone tries to contaminates the source directory with an in-source build.
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "Please do out-of-source builds.\nCleanup: \"rm -rf CMake CMakeCache.txt CMakeFiles/\"")
endif ()

find_package(Boost CONFIG COMPONENTS filesystem)

include(AICxxSubmodules)

add_subdirectory(src)
