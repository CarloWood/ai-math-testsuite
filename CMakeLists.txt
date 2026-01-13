cmake_minimum_required(VERSION 3.14...3.26)

project(ai-math-testsuite
  LANGUAGES CXX
  DESCRIPTION "Testsuite for git submodule math.")

# Compile all subdirectories with the same standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# If a local gitache directory is present then use that rather than downloading one.
if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/gitache/gateway.cmake)
  # This will disable the use of the GIT_REPOSITORY/GIT_TAG below, and disable the
  # FetchContent- download and update step. Instead, use the gitache directory as-is.
  set(FETCHCONTENT_SOURCE_DIR_GITACHE "${CMAKE_CURRENT_LIST_DIR}/gitache" CACHE INTERNAL "" FORCE)
endif ()

# Use gitache for libcwd.
set(GITACHE_PACKAGES libcwd_r)
FetchContent_Declare(
  gitache
  GIT_REPOSITORY "https://github.com/CarloWood/gitache.git"
  GIT_TAG        "master"
)
FetchContent_MakeAvailable(gitache)

# This project uses aicxx modules.
include(cwm4/cmake/AICxxProject NO_POLICY_SCOPE)

# Exit if someone tries to contaminates the source directory with an in-source build.
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "Please do out-of-source builds.\nCleanup: \"rm -rf CMake CMakeCache.txt CMakeFiles/\"")
endif ()

find_package(Boost CONFIG COMPONENTS filesystem)

include(AICxxSubmodules)

add_subdirectory(src)
# We use utils/to_string.h
add_subdirectory(enchantum)
