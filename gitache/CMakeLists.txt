cmake_minimum_required(VERSION 3.18...4.1.1)

message(DEBUG "DEBUG: Entering ${CMAKE_CURRENT_LIST_FILE}")

string(ASCII 27 _escape)
set(_red "${_escape}[31m")
set(_reset "${_escape}[m")

# If the environment variable GITACHE_ROOT is not set then gitache is disabled.
if (NOT DEFINED ENV{GITACHE_ROOT})
  message(STATUS "Environment variable GITACHE_ROOT is not set: ${_red}gitache disabled${_reset}.")
  return()
endif ()
set(GITACHE_ROOT $ENV{GITACHE_ROOT} CACHE INTERNAL "")

# Sanity check.
if (NOT GITACHE_PACKAGES)
  message(WARNING "${_red}GITACHE_PACKAGES is not set or empty${_reset}. Bailing out. ${_red}gitache disabled${_reset}.")
  return()
endif ()

project(gitache
  LANGUAGES NONE
  DESCRIPTION "The gateway to gitache-core."
)

# Each gitache sha1 (commit) is directly linked to a certain gitache-core commit, so that if you do not
# change gitache version (sha1) you also won't get a different gitache-core.
if (DEFINED ENV{GITACHE_CORE_SHA1})
  # This feature is intended for the developer of gitache.
  message(STATUS "Using environment variable GITACHE_CORE_SHA1 as gitache-core sha1.")
  set(GITACHE_CORE_SHA1 "$ENV{GITACHE_CORE_SHA1}")
else ()
  # The master branch of gitache follows the master branch of gitache-core.
  # Use the branch 'stable' of gitache to get a fixed (stable) SHA1 here.
  set(GITACHE_CORE_SHA1 "master")
endif ()
set(GITACHE_CORE_GIT_REPOSITORY "https://github.com/CarloWood/gitache-core.git")
set(GITACHE_CORE_DIR "${GITACHE_ROOT}/gitache-core" CACHE INTERNAL "")

function(_gitache_lock_core_directory)
  message(STATUS "Locking directory \"${GITACHE_CORE_DIR}\"...")
  file(LOCK ${GITACHE_CORE_DIR} DIRECTORY
    GUARD PROCESS
    RESULT_VARIABLE _error_result
    TIMEOUT 10
  )
  if (_error_result)
    message(FATAL_ERROR "  ${_error_result}.")
  endif ()
endfunction()

function(_gitache_unlock_core_directory)
  message(STATUS "Unlocking directory \"${GITACHE_CORE_DIR}\".")
  file(LOCK ${GITACHE_CORE_DIR} DIRECTORY RELEASE)
endfunction()

function(fetch_gitache_core)
  # The following lock is obtained while cloning and/or updating the gitache-core repository.
  _gitache_lock_core_directory()
  message(STATUS "  success.\n** Populating \"${GITACHE_CORE_DIR}\" with ${GITACHE_CORE_SHA1}.")
  FetchContent_Populate(gitache-core
    QUIET
    GIT_REPOSITORY "${GITACHE_CORE_GIT_REPOSITORY}"
    GIT_TAG        "${GITACHE_CORE_SHA1}"
    SUBBUILD_DIR   "${GITACHE_CORE_DIR}/subbuild"
    SOURCE_DIR     "${GITACHE_CORE_DIR}/source"
  )
  _gitache_unlock_core_directory()
endfunction()

message(STATUS "----------------------------------------------------\n** Populating gitache:")
message(STATUS "GITACHE_ROOT = \"${GITACHE_ROOT}\"")
# This seems informative enough to warrant it always being printed:
# Use -DCMAKE_MESSAGE_LOG_LEVEL=DEBUG for maximum output.
message(STATUS "CMAKE_MESSAGE_LOG_LEVEL = \"${CMAKE_MESSAGE_LOG_LEVEL}\"")
if (NOT IS_ABSOLUTE ${GITACHE_ROOT})
  message(FATAL_ERROR " GITACHE_ROOT must be an absolute path.")
endif ()
if (NOT EXISTS ${GITACHE_ROOT})
  message(FATAL_ERROR " ${GITACHE_ROOT}: no such file or directory.")
endif ()
if (NOT IS_DIRECTORY ${GITACHE_ROOT})
  message(FATAL_ERROR " ${GITACHE_ROOT}: not a directory.")
endif ()
# We're not testing if the directory is writable because cmake isn't supporting that
# (other than creating a file). We'll just proceed and error out when it isn't.
if (NOT EXISTS ${GITACHE_CORE_DIR})
  message(STATUS "Creating \"${GITACHE_CORE_DIR}\".")
  file(MAKE_DIRECTORY ${GITACHE_CORE_DIR})
  if (NOT EXISTS ${GITACHE_CORE_DIR})
    message(FATAL_ERROR " Please make GITACHE_ROOT writable for current user.")
  endif ()
  # Fetch gitache-core for the first time.
  fetch_gitache_core()
endif ()

# Should we use a submodule that is already added to the current source directory?
if (EXISTS ${CMAKE_SOURCE_DIR}/gitache-core/main.cmake)
  set(gitache_core_is_local TRUE)
else ()
  set(gitache_core_is_local FALSE)
endif ()

if (gitache_core_is_local)
  message(STATUS "Found local gitache-core submodule at \"${CMAKE_SOURCE_DIR}/gitache-core\". Using that.")
  # Set the gitache-core source directory.
  set(GITACHE_CORE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/gitache-core")
  include("${GITACHE_CORE_SOURCE_DIR}/updater.cmake")
else ()
  # Set the gitache-core source directory.
  set(GITACHE_CORE_SOURCE_DIR "${GITACHE_CORE_DIR}/source")
  # Pass control to gitache-core: updater.cmake, when successful, ends with including main.cmake.
  # If "updater.cmake" is reloaded, it sets gitache_need_include to TRUE and we do the include twice.
  set(gitache_need_include TRUE)
  while (gitache_need_include)
    set(gitache_need_include FALSE)
    file(COPY "${GITACHE_CORE_DIR}/source/updater.cmake" DESTINATION "${GITACHE_CORE_DIR}")
    include("${GITACHE_CORE_DIR}/updater.cmake") # Sets gitache_need_include to TRUE if source/updater.cmake got replaced.
    # Note: here we could REMOVE the copied file ${GITACHE_CORE_DIR}/updater.cmake again, but
    # because COPY checks the timestamp on source and destination and skips the copy when
    # those are equal, we can use the copy as a cache.
    if (gitache_need_include)
      message(DEBUG "DEBUG: gitache-core was updated: reloading \"${GITACHE_CORE_DIR}/source/updater.cmake\"!")
    endif ()
  endwhile ()
endif ()

# Propagate the package_ROOT's to the parent.
foreach(_package ${GITACHE_PACKAGES})
  set(${_package}_ROOT "${${_package}_ROOT}" PARENT_SCOPE)
endforeach()
