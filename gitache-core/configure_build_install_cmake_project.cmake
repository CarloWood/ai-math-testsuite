# The following variables need to be defined before including this file:
#
# SOURCE_DIR                    - The source directory of the package.
# BINARY_DIR                    - The binary directory of the package.
# CMAKE_CONFIG                  - The build config (for multi-target generators).
# CMAKE_ARGS                    - The arguments that should be passed to cmake during configuration.
#                                 This already includes -DCMAKE_BUILD_TYPE=${CMAKE_CONFIG}, for single-target generators.
# INSTALL_PREFIX                - The prefix used for installation of the package.
# CMAKE_MESSAGE_LOG_LEVEL       - An optional log-level for the cmake child processes.
# GITACHE_CORE_SOURCE_DIR       - The directory containing gitache-core.
# PACKAGE_NAME                  - Used for message output (short name of the gitache package).
# HASH_CONTENT_ES               - The string that was used to calculate the SHA256 that is part of the install prefix.
#                                 But with semicolons replaces with <-:-:->.

list(PREPEND CMAKE_MODULE_PATH "${GITACHE_CORE_SOURCE_DIR}/utils")
include(debug_support)  # For Dout.
include(color_vars)
include(ProcessorCount)

#include(${CMAKE_SOURCE_DIR}/../cwm4/cmake/dump_cmake_variables.cmake)
#dump_cmake_variables(.* RelWith)

# Convert HASH_CONTENT back.
string(REPLACE "<-:-:->" ";" HASH_CONTENT "${HASH_CONTENT_ES}")

# Process log-level, again.
set(_gitache_log_level)
set(_where NONE)
if(DEFINED CACHE{CMAKE_MESSAGE_LOG_LEVEL})
  message(STATUS "CMAKE_MESSAGE_LOG_LEVEL = \"${CMAKE_MESSAGE_LOG_LEVEL}\".")
  set(_gitache_log_level "-DCMAKE_MESSAGE_LOG_LEVEL=${CMAKE_MESSAGE_LOG_LEVEL}")
  if(${CMAKE_MESSAGE_LOG_LEVEL} STREQUAL "DEBUG")
    set(_where "STDOUT")
  endif()
endif()

# Determine the number of cores we have, again.
ProcessorCount(_cpus)
if(_cpus EQUAL 0)
  set(_cpus 4)
endif()

# Configure step.
message("${BoldCyan}Running configure step for '${PACKAGE_NAME}' [${HASH_CONTENT}].${ColourReset}")
separate_arguments(_cmake_args UNIX_COMMAND ${CMAKE_ARGS})
execute_process(
  COMMAND
    ${CMAKE_COMMAND}
      -S ${SOURCE_DIR}
      -B ${BINARY_DIR}
      ${_gitache_log_level}
      ${_cmake_args}
  COMMAND_ECHO ${_where}
  RESULT_VARIABLE _exit_code
)

if(_exit_code)
  message(FATAL_ERROR "Failed to configure cmake project at \"${SOURCE_DIR}\".")
endif()

# Build step.
message("${BoldCyan}Running build step for '${PACKAGE_NAME}' [${HASH_CONTENT}].${ColourReset}")
execute_process(
  COMMAND
    ${CMAKE_COMMAND}
    --build ${BINARY_DIR}
    --config ${CMAKE_CONFIG}
    --parallel ${_cpus}
  COMMAND_ECHO ${_where}
  RESULT_VARIABLE _exit_code
)

if(_exit_code)
  message(FATAL_ERROR "Failed to build cmake project at \"${BINARY_DIR}\".")
endif()

# Install step.
message("${BoldCyan}Running install step for '${PACKAGE_NAME}' [${HASH_CONTENT}].${ColourReset}")
execute_process(
  COMMAND
    ${CMAKE_COMMAND}
      --install ${BINARY_DIR}
      --config ${CMAKE_CONFIG}
      --prefix ${INSTALL_PREFIX}
  COMMAND_ECHO ${_where}
  RESULT_VARIABLE _exit_code
)

if(_exit_code)
  message(FATAL_ERROR "Failed to install cmake project at \"${BINARY_DIR}\".")
endif()
