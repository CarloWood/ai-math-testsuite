# The following variables need to be defined before including this file:
#
# gitache_package                       - The short name of the package.
# gitache_package_NAME                  - The name of the package as used by FetchContent.
# gitache_package_ROOT                  - The '_deps' base directory that FetchContent should use.
# gitache_package_CMAKE_CONFIG          - The build config (for multi-target generators).
# gitache_package_CMAKE_ARGS            - If the package is a cmake project, arguments that should be passed to cmake during configuration.
#                                         This already includes -DCMAKE_BUILD_TYPE=${gitache_package_CMAKE_CONFIG}, for single-target.
# gitache_package_CONFIGURE_ARGS        - If the package is an autotools project, arguments that should be passed to configure.
#                                         This already includes --prefix=${gitache_package_INSTALL_PREFIX}.
# gitache_package_INSTALL_PREFIX        - The prefix used for installation of the package.
# arguments_to_FetchContent_Declare     - Space separated list of arguments that need to be passed to FetchContent_Declare.
# gitache_package_LOCK_ID               - An identifying string for this install (this is written to the DONE file).
# gitache_log_level                     - An optional -DCMAKE_MESSAGE_LOG_LEVEL=* argument for a cmake child process.
# gitache_where                         - NONE or STDOUT, depending on log-level.
# GITACHE_CORE_SOURCE_DIR               - The directory containing gitache-core.
# gitache_package_HASH_CONTENT          - The string over which the hash is calculated (this is written to the DONE file). Contains semi-colons!
# gitache_package_BOOTSTRAP_COMMAND     - A user defined command to run before configuration of a the package.

if (NOT DEFINED GITACHE_CORE_SOURCE_DIR OR "${GITACHE_CORE_SOURCE_DIR}" STREQUAL "")
  #  message(FATAL_ERROR "gitache-core/package.cmake: GITACHE_CORE_SOURCE_DIR is not set!")
  set(ERROR_MESSAGE "gitache-core/package.cmake:22: GITACHE_CORE_SOURCE_DIR is not set!")
  return()
endif ()

# This is an output variable.
set(ERROR_MESSAGE False)

# Use $GITACHE_ROOT/package_name as base dir (this replaces the usual ${CMAKE_CURRENT_BINARY_DIR}/_deps).
set(FETCHCONTENT_BASE_DIR "${gitache_package_ROOT}")

# Following how FetchContent works, the source code of a packages is put into
set(_source_dir ${FETCHCONTENT_BASE_DIR}/${gitache_package_NAME}-src)

# And the build directory then is prepared at
set(_build_dir ${FETCHCONTENT_BASE_DIR}/${gitache_package_NAME}-build)

# The existence of this file marks successful installation.
set(_done_file "${gitache_package_INSTALL_PREFIX}/DONE")

# We already have FETCHCONTENT_BASE_DIR locked before including this file.
# Hence if this file exists then we're done.
if(EXISTS ${_done_file})
  message(STATUS "${gitache_package_NAME} is already installed.")
else()
  set(FETCHCONTENT_QUIET OFF)

  Dout("Calling FetchContent_Declare(${gitache_package_NAME} SOURCE_SUBDIR \"doesnotexist\" ${arguments_to_FetchContent_Declare})")
  separate_arguments(_args UNIX_COMMAND ${arguments_to_FetchContent_Declare})
  FetchContent_Declare(
    ${gitache_package_NAME}
    SOURCE_SUBDIR "doesnotexist"
    # E.g. GIT_TAG and GIT_REPOSITORY. See gitache_config.cmake for a full list.
    ${_args}
  )

  string(TOLOWER ${gitache_package_NAME} gitache_package_NAME_lc)
  FetchContent_GetProperties(${gitache_package_NAME})

  if(NOT ${gitache_package_NAME_lc}_POPULATED)
    # Unset CMAKE_GENERATOR while running FetchContent_MakeAvailable, so that it
    # uses the (platform) default generator (you can still pass a -G to its
    # CMAKE_ARGS of course.
    set(CMAKE_GENERATOR_STORE ${CMAKE_GENERATOR})
    unset(CMAKE_GENERATOR CACHE)

    # FetchContent_MakeAvailable might do nothing even if the source tree is empty,
    # because "stamp" files have been stored in the build directory of the super project.
    if(NOT EXISTS "${FETCHCONTENT_BASE_DIR}/${gitache_package_NAME_lc}-src/.git")
      set(_fc_stamp_dir "${CMAKE_BINARY_DIR}/CMakeFiles/fc-stamp/${gitache_package_NAME_lc}")
      set(_fc_tmp_dir   "${CMAKE_BINARY_DIR}/CMakeFiles/fc-tmp/${gitache_package_NAME_lc}")
      if(EXISTS "${_fc_stamp_dir}" OR EXISTS "${_fc_tmp_dir}")
        message(STATUS "Warning: ${FETCHCONTENT_BASE_DIR}/${gitache_package_NAME_lc}-src/.git does not exist. Removing stamp files.")
        file(REMOVE "${_fc_stamp_dir}/download.stamp"
                    "${_fc_stamp_dir}/update.stamp"
                    "${_fc_stamp_dir}/patch.stamp"
                    "${_fc_stamp_dir}/gitache_package_versor-gitclone-lastrun.txt")
        file(REMOVE_RECURSE "${_fc_tmp_dir}")
      endif()
    endif()

    # Populate the source and binary directories of this package.
    Dout("Calling FetchContent_MakeAvailable(${gitache_package_NAME})")
    FetchContent_MakeAvailable(${gitache_package_NAME})
    set(CMAKE_GENERATOR "${CMAKE_GENERATOR_STORE}" CACHE INTERNAL "The projects generator")
    set(FETCHCONTENT_QUIET ON)

    if(NOT "${gitache_package_BOOTSTRAP_COMMAND}" STREQUAL "")
      # Bootstrap step.
      message("${BoldCyan}Running bootstrap step for '${gitache_package}' [${gitache_package_HASH_CONTENT}].${ColourReset}")
      message(STATUS "gitache_package_BOOTSTRAP_COMMAND = \"${gitache_package_BOOTSTRAP_COMMAND}\".")
      set(_exit_code 2)
      execute_process(
        COMMAND
          ${gitache_package_BOOTSTRAP_COMMAND}
        COMMAND_ECHO ${gitache_where}
        WORKING_DIRECTORY ${${gitache_package_NAME_lc}_SOURCE_DIR}
        RESULT_VARIABLE _exit_code
      )
      if(_exit_code)
        if(_exit_code EQUAL 2)
          set(ERROR_MESSAGE "Fatal error: execute_process() did not run. See CMake Error above.")
        else()
          set(ERROR_MESSAGE "Failed to bootstrap autotools project at \"${${gitache_package_NAME_lc}_SOURCE_DIR}\".")
        endif()
        return()
      endif()
    endif()
    string(REPLACE ";" "<-:-:->" _hash_content_encoded_semicolon "${gitache_package_HASH_CONTENT}")
    if(EXISTS ${${gitache_package_NAME_lc}_SOURCE_DIR}/CMakeLists.txt)
      Dout("Attempting to configure/build/install \"${gitache_package_NAME}\" as cmake project; running:")
      # Start a separate process to configure, build and install this cmake package.
      execute_process(
        COMMAND
          ${CMAKE_COMMAND} ${gitache_log_level}
            -DCMAKE_ARGS=${gitache_package_CMAKE_ARGS}
            -DCMAKE_CONFIG=${gitache_package_CMAKE_CONFIG}
            -DPACKAGE_NAME=${gitache_package}
            -DGITACHE_CORE_SOURCE_DIR=${GITACHE_CORE_SOURCE_DIR}
            -DSOURCE_DIR=${${gitache_package_NAME_lc}_SOURCE_DIR}
            -DBINARY_DIR=${${gitache_package_NAME_lc}_BINARY_DIR}
            -DINSTALL_PREFIX=${gitache_package_INSTALL_PREFIX}
            -DHASH_CONTENT_ES=${_hash_content_encoded_semicolon}
            -P "${CMAKE_CURRENT_LIST_DIR}/configure_build_install_cmake_project.cmake"
        COMMAND_ECHO ${gitache_where}
        RESULT_VARIABLE
          _exit_code
      )
      if(_exit_code) # A cmake script always returns either 0 (success) or 1 (failure).
        set(ERROR_MESSAGE "Failed to config/build/install cmake package \"${gitache_package_NAME}\".")
        return()
      endif()
    elseif(EXISTS ${${gitache_package_NAME_lc}_SOURCE_DIR}/Makefile.am AND
           EXISTS ${${gitache_package_NAME_lc}_SOURCE_DIR}/configure.ac)
      Dout("Attempting to configure/build/install \"${gitache_package_NAME}\" as autotools project; running:")
      # Start a separate process to configure, make and install this autotools package.
      string(REPLACE ";" "\\;" escaped_CONFIGURE_ARGS "${gitache_package_CONFIGURE_ARGS}")
      execute_process(
        COMMAND
          ${CMAKE_COMMAND} ${gitache_log_level}
            -DCONFIGURE_ARGS=${escaped_CONFIGURE_ARGS}
            -DPACKAGE_NAME=${gitache_package}
            -DGITACHE_CORE_SOURCE_DIR=${GITACHE_CORE_SOURCE_DIR}
            -DSOURCE_DIR=${${gitache_package_NAME_lc}_SOURCE_DIR}
            -DBINARY_DIR=${${gitache_package_NAME_lc}_BINARY_DIR}
            -DHASH_CONTENT_ES=${_hash_content_encoded_semicolon}
            -P "${CMAKE_CURRENT_LIST_DIR}/configure_make_install_autotools_project.cmake"
        COMMAND_ECHO ${gitache_where}
        RESULT_VARIABLE
          _exit_code
      )
      if(_exit_code)
        set(ERROR_MESSAGE "Failed to configure/make/install autotools package \"${gitache_package_NAME}\".")
        return()
      endif()
    else()
      set(ERROR_MESSAGE "Don't know how to build gitache package \"${gitache_package_NAME}\".")
      return()
    endif()
  else()
    set(FETCHCONTENT_QUIET ON)
  endif()

  # The above only has to be done once.
  Dout("Creating ${_done_file}")
  file(WRITE ${_done_file} "${gitache_package_HASH_CONTENT} - ${gitache_package_LOCK_ID}\n")
endif()

# Restore default values.
unset(FETCHCONTENT_BASE_DIR)
