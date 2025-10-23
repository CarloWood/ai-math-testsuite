# Every package in GITACHE_PACKAGES must have a config file.
set(GITACHE_CONFIGS_DIR "${CMAKE_SOURCE_DIR}/cmake/gitache-configs" CACHE PATH "Directory containing configs for gitache packages.")
message(STATUS "Reading package configurations from \"${GITACHE_CONFIGS_DIR}/\"...")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cwm4/cmake")

include(color_vars)
include(gitache_config)
include(ExternalProject)

# Create necessary directories.
set(_package_directory "${CMAKE_CURRENT_BINARY_DIR}/packages")
file(MAKE_DIRECTORY ${_package_directory})

# Generate a unique ID for file locking.
string(TIMESTAMP _current_time "%Y-%m-%dT%H:%M")
string(RANDOM LENGTH 8 _random_string)
set(gitache_package_LOCK_ID "${CMAKE_SOURCE_DIR} - configured ${_current_time} - ${_random_string}")

function(lock_directory package_root)
  message(DEBUG "Locking directory \"${package_root}\".")
  execute_process(
    COMMAND "${GITACHE_CORE_SOURCE_DIR}/lock.sh" ${package_root} ${gitache_package_LOCK_ID}
    WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    RESULT_VARIABLE _exit_code
  )
  if(_exit_code)
    message(FATAL_ERROR "Could not lock directory.")
  endif()
endfunction()

function(unlock_directory package_root)
  message(DEBUG "Unlocking directory \"${package_root}\".")
  execute_process(
    COMMAND "${GITACHE_CORE_SOURCE_DIR}/unlock.sh" ${package_root} ${gitache_package_LOCK_ID}
    WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    RESULT_VARIABLE _exit_code
  )
  if(_exit_code)
    message(FATAL_ERROR "Could not unlock directory.")
  endif()
endfunction()

set(gitache_log_level)
if(DEFINED CACHE{CMAKE_MESSAGE_LOG_LEVEL})
  set(gitache_log_level "-DCMAKE_MESSAGE_LOG_LEVEL=${CMAKE_MESSAGE_LOG_LEVEL}")
endif()

set(package_independent_seed "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
foreach (gitache_package ${GITACHE_PACKAGES})
  if(NOT EXISTS "${GITACHE_CONFIGS_DIR}/${gitache_package}.cmake")
    message(FATAL_ERROR
      " No configuration file found for package ${gitache_package}.\n"
      " Please add a file \"${GITACHE_CONFIGS_DIR}/${gitache_package}.cmake\"."
    )
  endif()
  # Load the user specified configuration for this package.
  set(GIT_TAG)
  set(gitache_package_CMAKE_CONFIG "Release")
  include("${GITACHE_CONFIGS_DIR}/${gitache_package}.cmake") # Sets GIT_TAG if any is specified.
  if (NOT "${cmake_arguments}" STREQUAL "" AND NOT "${configure_arguments}" STREQUAL "")
    message(FATAL_ERROR
      " ${GITACHE_CONFIGS_DIR}/${gitache_package}.cmake: can only specify one of CMAKE_ARGS or CONFIGURE_ARGS."
    )
  endif()
  # Calculate a hash that determines everything that might have an influence on the build result.
  set(gitache_package_HASH_CONTENT "${package_independent_seed}|${arguments_to_FetchContent_Declare}|${bootstrap_command}|${cmake_arguments}${configure_arguments}")
  string(SHA256 ${gitache_package}_config_hash "${gitache_package_HASH_CONTENT}")

  Dout("${gitache_package}_config_hash = \"${${gitache_package}_config_hash}\"")
  Dout("${gitache_package}: package_independent_seed = \"${package_independent_seed}\".")
  Dout("${gitache_package}: arguments_to_FetchContent_Declare = \"${arguments_to_FetchContent_Declare}\".")
  Dout("${gitache_package}: bootstrap_command = \"${bootstrap_command}\".")
  Dout("${gitache_package}: cmake_arguments = \"${cmake_arguments}\".")
  Dout("${gitache_package}: configure_arguments = \"${configure_arguments}\".")

  set(gitache_package_ROOT "${GITACHE_ROOT}/${gitache_package}")
  set(gitache_package_NAME "gitache_package_${gitache_package}")
  set(gitache_package_INSTALL_PREFIX "${gitache_package_ROOT}/${${gitache_package}_config_hash}")
  set(gitache_package_BOOTSTRAP_COMMAND ${bootstrap_command})
  set(gitache_package_CMAKE_ARGS "${cmake_arguments} -DCMAKE_INSTALL_PREFIX=${gitache_package_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${gitache_package_CMAKE_CONFIG}")
  set(gitache_package_CONFIGURE_ARGS ${configure_arguments})
  list(APPEND gitache_package_CONFIGURE_ARGS "--prefix=${gitache_package_INSTALL_PREFIX}")

  set(_gitupdate_stamp_file "${gitache_package_ROOT}/src/gitache_package_${gitache_package}-stamp/gitache_package_${gitache_package}-gitupdate")
  set(_lock_stamp_file "${gitache_package_ROOT}/src/gitache_package_${gitache_package}-stamp/gitache_package_${gitache_package}-lock")

  file(MAKE_DIRECTORY ${gitache_package_ROOT})
  lock_directory(${gitache_package_ROOT})
  include("${CMAKE_CURRENT_LIST_DIR}/package.cmake")
  unlock_directory(${gitache_package_ROOT})

  if(ERROR_MESSAGE)
    break()
  endif()

  # The find_package will be done in the parent scope.
  set(${gitache_package}_ROOT "${gitache_package_INSTALL_PREFIX}" PARENT_SCOPE)
  # Also set this variable in the current scope, so that subsequent packages can find packages we already configured.
  set(${gitache_package}_ROOT "${gitache_package_INSTALL_PREFIX}")
  Dout("${gitache_package}: ${gitache_package}_ROOT = \"${${gitache_package}_ROOT}\".")

  # Remove possible old entries from the cache.
  if(DEFINED CACHE{${gitache_package}_DIR} AND NOT ${${gitache_package}_DIR} STREQUAL "${${gitache_package}_ROOT}/lib/cmake/${gitache_package}")
    message(NOTICE ">> ${Red}Removing old cache value ${gitache_package}_DIR (\"${${gitache_package}_DIR}\")!${ColourReset}")
    unset(${gitache_package}_DIR CACHE)
    unset(${gitache_package}_VERSION CACHE)
  endif()

endforeach ()
