# The following variables need to be defined before including this file:
#
# GITACHE_PACKAGES                      - Space separated list of packages required by the main project.
# gitache_where                         - NONE or STDOUT, depending on log-level.
# GITACHE_CORE_SOURCE_DIR               - The directory containing gitache-core.

# Every package in GITACHE_PACKAGES must have a config file in _gitache_default_config_dir.
set(_gitache_default_config_dir "${CMAKE_SOURCE_DIR}/cmake/gitache-configs" CACHE PATH "Default directory containing configs for gitache packages.")
set(GITACHE_CONFIG_DIRS "${_gitache_default_config_dir}" CACHE STRING "Semicolon-separated list of directories containing configs for gitache packages.")

set(package_independent_seed "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
Dout("package_independent_seed = \"${package_independent_seed}\"")

# Keep track of (other) registered configuration directories and processed packages.
set_property(GLOBAL PROPERTY GITACHE_CONFIG_DIRS "${GITACHE_CONFIG_DIRS}")
set_property(GLOBAL PROPERTY GITACHE_CORE_MAIN_DIR "${CMAKE_CURRENT_LIST_DIR}")
set_property(GLOBAL PROPERTY GITACHE_COMMAND_ECHO "${gitache_where}")
set_property(GLOBAL PROPERTY GITACHE_PROCESSED_PACKAGES "")
set_property(GLOBAL PROPERTY GITACHE_PACKAGE_INDEPENDENT_SEED "${package_independent_seed}")
set_property(GLOBAL PROPERTY GITACHE_CORE_SOURCE_DIR "${GITACHE_CORE_SOURCE_DIR}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cwm4/cmake")

include(color_vars)
include(gitache_config)
include(ExternalProject)

# Create necessary directories.
set(_package_directory "${CMAKE_CURRENT_BINARY_DIR}/packages")
file(MAKE_DIRECTORY ${_package_directory})

function(lock_directory package_root)
  message(DEBUG "Locking directory \"${package_root}\".")
  get_property(_core_dir GLOBAL PROPERTY GITACHE_CORE_SOURCE_DIR)
  # Generate a unique ID for file locking.
  string(TIMESTAMP _current_time "%Y-%m-%dT%H:%M")
  string(RANDOM LENGTH 8 _random_string)
  set(_lock_id "${CMAKE_CURRENT_SOURCE_DIR} - configured ${_current_time} - ${_random_string}")
  set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY GITACHE_LOCK_ID "${_lock_id}")
  execute_process(
    COMMAND "${_core_dir}/lock.sh" ${package_root} ${_lock_id}
    WORKING_DIRECTORY ${_core_dir}
    RESULT_VARIABLE _exit_code
  )
  if(_exit_code)
    message(FATAL_ERROR "Could not lock directory.")
  endif()
endfunction()

function(unlock_directory package_root)
  message(DEBUG "Unlocking directory \"${package_root}\".")
  get_property(_core_dir GLOBAL PROPERTY GITACHE_CORE_SOURCE_DIR)
  get_property(_lock_id DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY GITACHE_LOCK_ID)
  execute_process(
    COMMAND "${_core_dir}/unlock.sh" ${package_root} ${_lock_id}
    WORKING_DIRECTORY ${_core_dir}
    RESULT_VARIABLE _exit_code
  )
  if(_exit_code)
    message(FATAL_ERROR "Could not unlock directory.")
  endif()
endfunction()

function(_gitache_get_registered_config_dirs out_var)
  get_property(_dirs GLOBAL PROPERTY GITACHE_CONFIG_DIRS)
  if(NOT _dirs)
    set(_dirs)
  endif()
  set(${out_var} "${_dirs}" PARENT_SCOPE)
endfunction()

function(_gitache_register_config_dir absolute_dir)
  get_property(_dirs GLOBAL PROPERTY GITACHE_CONFIG_DIRS)
  if(_dirs)
    list(FIND _dirs "${absolute_dir}" _index)
  else()
    set(_index -1)
  endif()
  if(_index EQUAL -1)
    if(_dirs)
      list(APPEND _dirs "${absolute_dir}")
    else()
      set(_dirs "${absolute_dir}")
    endif()
    set_property(GLOBAL PROPERTY GITACHE_CONFIG_DIRS "${_dirs}")
  endif()
endfunction()

function(gitache_register_config_dir dir)
  get_filename_component(_absolute "${dir}" ABSOLUTE)
  if(NOT IS_DIRECTORY "${_absolute}")
    message(FATAL_ERROR "gitache_register_config_dir: \"${dir}\" is not a directory.")
  endif()

  _gitache_register_config_dir(${_absolute})
  Dout("Registered gitache config dir \"${_absolute}\".")
endfunction()

function(_gitache_locate_package_config package out_var)
  # If a CMakeLists.txt file calls `gitache_require_packages` directly, CMAKE_CURRENT_LIST_DIR will be the directory where this CMakeLists.txt resides.
  _gitache_register_config_dir(${CMAKE_CURRENT_LIST_DIR}/cmake/gitache-configs)
  get_property(_dirs GLOBAL PROPERTY GITACHE_CONFIG_DIRS)
  foreach(_dir ${_dirs})
    set(_candidate "${_dir}/${package}.cmake")
    if(EXISTS "${_candidate}")
      set(${out_var} "${_candidate}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${out_var} "" PARENT_SCOPE)
endfunction()

set(gitache_log_level)
if(DEFINED CACHE{CMAKE_MESSAGE_LOG_LEVEL})
  set(gitache_log_level "-DCMAKE_MESSAGE_LOG_LEVEL=${CMAKE_MESSAGE_LOG_LEVEL}")
endif()

function(gitache_process_package package)
  if("${package}" STREQUAL "")
    message(FATAL_ERROR "gitache_process_package called with an empty package name.")
  endif()

  get_property(_processed GLOBAL PROPERTY GITACHE_PROCESSED_PACKAGES)
  if(_processed)
    list(FIND _processed "${package}" _already_index)
  else()
    set(_already_index -1)
  endif()
  if(NOT _already_index EQUAL -1)
    message(DEBUG "gitache package \"${package}\" already processed; skipping.")
    return()
  endif()

  set(gitache_package "${package}")

  _gitache_locate_package_config("${package}" _config_file)
  if("${_config_file}" STREQUAL "")
    get_property(_search_dirs GLOBAL PROPERTY GITACHE_CONFIG_DIRS)
    if(_search_dirs)
      string(JOIN "\n  - " _formatted_dirs ${_search_dirs})
      set(_search_dirs_message "\n  - ${_formatted_dirs}")
    else()
      set(_search_dirs_message " <no directories registered>")
    endif()
    message(FATAL_ERROR
      " No configuration file found for package ${package}.\n"
      " Searched:${_search_dirs_message}.\n"
      " Please add a file \"<config_dir>/${package}.cmake\" and/or add 'CONFIG_DIRS \"<config_dir>\"' to `gitache_require_packages` (or call `gitache_register_config_dir(\"<config_dir>\") first)`.")
  else()
    Dout(${gitache_package}: "using config file: \"${_config_file}\"")
  endif()

  set(GIT_TAG)
  set(gitache_package_CMAKE_CONFIG "Release")
  include("${_config_file}")
  if(NOT "${cmake_arguments}" STREQUAL "" AND NOT "${configure_arguments}" STREQUAL "")
    message(FATAL_ERROR
      " ${_config_file}: can only specify one of CMAKE_ARGS or CONFIGURE_ARGS."
    )
  endif()

  get_property(_package_seed GLOBAL PROPERTY GITACHE_PACKAGE_INDEPENDENT_SEED)
  if(NOT _package_seed)
    set(_package_seed "${package_independent_seed}")
  endif()

  set(gitache_package_HASH_CONTENT "${_package_seed}|${arguments_to_FetchContent_Declare}|${bootstrap_command}|${cmake_arguments}${configure_arguments}")
  string(SHA256 ${gitache_package}_config_hash "${gitache_package_HASH_CONTENT}")

  Dout("${gitache_package}_config_hash = \"${${gitache_package}_config_hash}\"")
  Dout("${gitache_package}: package_independent_seed = \"${_package_seed}\".")
  Dout("${gitache_package}: arguments_to_FetchContent_Declare: ${arguments_to_FetchContent_Declare}")
  Dout("${gitache_package}: bootstrap_command: ${bootstrap_command}")
  Dout("${gitache_package}: cmake_arguments: ${cmake_arguments}")
  Dout("${gitache_package}: configure_arguments: ${configure_arguments}")

  set(gitache_package_ROOT "${GITACHE_ROOT}/${gitache_package}")
  set(gitache_package_NAME "gitache_package_${gitache_package}")
  set(gitache_package_INSTALL_PREFIX "${gitache_package_ROOT}/${${gitache_package}_config_hash}")
  set(gitache_package_BOOTSTRAP_COMMAND ${bootstrap_command})
  set(gitache_package_CMAKE_ARGS "${cmake_arguments} -DCMAKE_INSTALL_PREFIX=${gitache_package_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${gitache_package_CMAKE_CONFIG}")
  set(gitache_package_CONFIGURE_ARGS ${configure_arguments})
  list(APPEND gitache_package_CONFIGURE_ARGS "--prefix=${gitache_package_INSTALL_PREFIX}")

  set(_gitupdate_stamp_file "${gitache_package_ROOT}/src/gitache_package_${gitache_package}-stamp/gitache_package_${gitache_package}-gitupdate")
  set(_lock_stamp_file "${gitache_package_ROOT}/src/gitache_package_${gitache_package}-stamp/gitache_package_${gitache_package}-lock")

  get_property(_main_dir GLOBAL PROPERTY GITACHE_CORE_MAIN_DIR)
  if(NOT _main_dir)
    message(FATAL_ERROR "gitache main directory is unknown.")
  endif()

  get_property(_command_echo GLOBAL PROPERTY GITACHE_COMMAND_ECHO)
  if(DEFINED _command_echo AND NOT "${_command_echo}" STREQUAL "")
    set(gitache_where "${_command_echo}")
  else()
    set(gitache_where NONE)
  endif()

  file(MAKE_DIRECTORY ${gitache_package_ROOT})
  lock_directory(${gitache_package_ROOT})
  # This sets ERROR_MESSAGE upon a fatal error.
  include("${_main_dir}/package.cmake")
  unlock_directory(${gitache_package_ROOT})
  if(ERROR_MESSAGE)
    message(FATAL_ERROR "${ERROR_MESSAGE}")
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

  get_property(_processed GLOBAL PROPERTY GITACHE_PROCESSED_PACKAGES)
  if(_processed)
    list(APPEND _processed "${package}")
  else()
    set(_processed "${package}")
  endif()
  set_property(GLOBAL PROPERTY GITACHE_PROCESSED_PACKAGES "${_processed}")
endfunction()

function(gitache_require_packages)
  message(STATUS "*** CMAKE_SOURCE_DIR = \"${CMAKE_SOURCE_DIR}\"")
  message(STATUS "*** CMAKE_CURRENT_SOURCE_DIR = \"${CMAKE_CURRENT_SOURCE_DIR}\"")
  message(STATUS "*** CMAKE_CURRENT_LIST_DIR = \"${CMAKE_CURRENT_LIST_DIR}\"")
  message(STATUS "*** PROJECT_SOURCE_DIR = \"${PROJECT_SOURCE_DIR}\"")

  # This function can be called from a submodules CMakeLists.txt;
  # we have to load all required global properties that were set by the super project!
  get_property(_core_dir GLOBAL PROPERTY GITACHE_CORE_SOURCE_DIR)
  set(GITACHE_CORE_SOURCE_DIR ${_core_dir})

  set(options)
  set(oneValueArgs)
  set(multiValueArgs CONFIG_DIRS)
  cmake_parse_arguments(gitache_require "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  foreach(_dir ${gitache_require_CONFIG_DIRS})
    gitache_register_config_dir("${_dir}")
  endforeach()

  set(_packages ${gitache_require_UNPARSED_ARGUMENTS})
  if(NOT _packages)
    return()
  endif()

  set(_processed_in_call)
  foreach(_package ${_packages})
    gitache_process_package("${_package}")
    list(APPEND _processed_in_call "${_package}")
  endforeach()

  foreach(_package ${_processed_in_call})
    if(DEFINED ${_package}_ROOT)
      set(${_package}_ROOT "${${_package}_ROOT}" PARENT_SCOPE)
    endif()
  endforeach()
endfunction()

foreach(_config_dir IN LISTS GITACHE_CONFIG_DIRS)
  gitache_register_config_dir("${_config_dir}")
endforeach()
_gitache_get_registered_config_dirs(_gitache_effective_dirs)
if(_gitache_effective_dirs)
  string(JOIN "\", \"" _gitache_dir_message ${_gitache_effective_dirs})
  message(STATUS "Reading package configurations from \"${_gitache_dir_message}\".")
else()
  message(STATUS "No gitache configuration directories registered.")
endif()

if(GITACHE_PACKAGES)
  gitache_require_packages(${GITACHE_PACKAGES})
endif()
