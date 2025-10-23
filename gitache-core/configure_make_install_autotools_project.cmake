# The following variables need to be defined before including this file:
#
# SOURCE_DIR                    - The source directory of the package.
# BINARY_DIR                    - The binary directory of the package.
# CMAKE_MESSAGE_LOG_LEVEL       - An optional log-level for the cmake child processes.
# GITACHE_CORE_SOURCE_DIR       - The directory containing gitache-core.
# PACKAGE_NAME                  - Used for message output (short name of the gitache package).
# HASH_CONTENT_ES               - The string that was used to calculate the SHA256 that is part of the install prefix.
#                                 But with semicolons replaced with <-:-:->.
# CONFIGURE_ARGS                - User defined arguments that should be passed to configure.
#                                 This already includes the --prefix=${INSTALL_PREFIX}.

list(PREPEND CMAKE_MODULE_PATH "${GITACHE_CORE_SOURCE_DIR}/utils")
include(color_vars)
include(debug_support)  # For Dout.
include(ProcessorCount)

#include(${CMAKE_SOURCE_DIR}/../cwm4/cmake/dump_cmake_variables.cmake)
#dump_cmake_variables(.* RelWith)

# Convert HASH_CONTENT back.
string(REPLACE "<-:-:->" ";" HASH_CONTENT "${HASH_CONTENT_ES}")

# Process log-level, again.
set(_where NONE)
if(DEFINED CACHE{CMAKE_MESSAGE_LOG_LEVEL})
  message(STATUS "CMAKE_MESSAGE_LOG_LEVEL = \"${CMAKE_MESSAGE_LOG_LEVEL}\".")
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
execute_process(
  COMMAND
    ${SOURCE_DIR}/configure ${CONFIGURE_ARGS}
  COMMAND_ECHO ${_where}
  WORKING_DIRECTORY ${BINARY_DIR}
  RESULT_VARIABLE _exit_code
)

if(_exit_code)
  message(FATAL_ERROR "Failed to configure autotools project at \"${SOURCE_DIR}\".")
endif()

# Build step.
message("${BoldCyan}Running build step for '${PACKAGE_NAME}' [${HASH_CONTENT}].${ColourReset}")
execute_process(
  COMMAND
    make -j ${_cpus}
  COMMAND_ECHO ${_where}
  WORKING_DIRECTORY ${BINARY_DIR}
  RESULT_VARIABLE _exit_code
)

if(_exit_code)
  message(FATAL_ERROR "Failed to make autotools project at \"${BINARY_DIR}\".")
endif()

# Install step.
message("${BoldCyan}Running install step for '${PACKAGE_NAME}' [${HASH_CONTENT}].${ColourReset}")
execute_process(
  COMMAND
    make install
  COMMAND_ECHO ${_where}
  WORKING_DIRECTORY ${BINARY_DIR}
  RESULT_VARIABLE _exit_code
)

if(_exit_code)
  message(FATAL_ERROR "Failed to install autotools project at \"${BINARY_DIR}\".")
endif()
