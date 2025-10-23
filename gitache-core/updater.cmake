# Copyright (c) 2020, Carlo Wood (carlo@alinoe.com)
# All rights reserved.

# This module makes sure the gitache-core source with sha1
# of ${GITACHE_CORE_SHA1} is available.
#
# Input:
#    gitache_core_is_local              True if this is a local submodule.
#    GITACHE_CORE_SHA1                  The SHA1 of the required commit of gitache-core.
#
# Output:
#    CMAKE_MODULE_PATH                  This file prepends the necessary module load paths.
#    git_where                          Depending on CMAKE_MESSAGE_LOG_LEVEL, set to either NONE or STDOUT.
#    git_executable                     Set to the absolute path to the git executable.
#    gitache_need_include               If GITACHE_CORE_SHA1 wasn't checked out yet, then
#                                       if gitache_core_is_local is FALSE, this file will
#                                       do the (fetch and) checkout and set
#                                       gitache_need_include to TRUE. Then return to reload.
#                                       Otherwise, when gitache_core_is_local is TRUE, that
#                                       is an error.
#
# If GITACHE_CORE_SHA1 is already checked out, then GITACHE_CORE_SOURCE_DIR
# is locked and control passed to main.cmake.

message(DEBUG "DEBUG: Entering `${CMAKE_CURRENT_LIST_FILE}`")

# Make sure that GITACHE_CORE_SHA1 is checked out.
# This is done by executing git commands directly; so find git.

# Add utils subdirectory to CMAKE_MODULE_PATH.
list(PREPEND CMAKE_MODULE_PATH "${GITACHE_CORE_SOURCE_DIR}/utils")

# Include utility modules and find git executable.
include(gitache_get_git_executable)
include(debug_support)
gitache_get_git_executable(git_executable)
Dout("git_executable = \"${git_executable}\".")

# Show COMMANDs if log-level is DEBUG.
set(gitache_where NONE)
if(DEFINED CACHE{CMAKE_MESSAGE_LOG_LEVEL} AND "${CMAKE_MESSAGE_LOG_LEVEL}" STREQUAL "DEBUG")
  set(gitache_where "STDOUT")
endif()

# Stop other processes from changing the SHA1.
lock_core_directory()

# Get the SHA1 that is checked out right now.
execute_process(COMMAND ${git_executable} rev-parse HEAD
  COMMAND_ECHO ${gitache_where}
  WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
  OUTPUT_VARIABLE head_sha1
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# In case the passed "GITACHE_CORE_SHA1" isn't a SHA1, try to do something sane.
set(_fetch_done false)
set(_commit_sha1)
if(NOT GITACHE_CORE_SHA1 MATCHES
    "^[0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F]\
[0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F]\
[0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F]\
[0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F][0-9a-fA-F]$")
  # Fetch upstream.
  execute_process(COMMAND ${git_executable} fetch --tags
    COMMAND_ECHO ${gitache_where}
    WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    RESULT_VARIABLE _exit_code
  )
  if(NOT _exit_code)
    set(_fetch_done true)
  endif()
  # Is it a tag?
  execute_process(COMMAND ${git_executable} show-ref --hash --verify refs/tags/${GITACHE_CORE_SHA1}
    COMMAND_ECHO ${gitache_where}
    WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    RESULT_VARIABLE _exit_code
    OUTPUT_VARIABLE _commit_sha1
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(_exit_code)
    # Is it a branch?
    execute_process(COMMAND ${git_executable} show-ref --hash --verify refs/remotes/origin/${GITACHE_CORE_SHA1}
      COMMAND_ECHO ${gitache_where}
      WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
      RESULT_VARIABLE _exit_code
      OUTPUT_VARIABLE _commit_sha1
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
  endif()
endif()
if(NOT _commit_sha1)
  set(_fetched FALSE)
  while (TRUE)
    # Is it anything that refers to an existing commit?
    execute_process(COMMAND ${git_executable} rev-parse --verify "${GITACHE_CORE_SHA1}^{commit}"
      COMMAND_ECHO ${gitache_where}
      WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
      RESULT_VARIABLE _exit_code
      OUTPUT_VARIABLE _commit_sha1
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (NOT _exit_code)
      break()
    endif ()
    if (_fetched)
      message(FATAL_ERROR "The environment variable GITACHE_CORE_SHA1 is set to \"${GITACHE_CORE_SHA1}\", which does not exist in the gitache-core repository.")
    endif ()
    # Fetch once, then loop to retry the rev-parse.
    execute_process(
      COMMAND ${git_executable} fetch
      COMMAND_ECHO ${gitache_where}
      WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
      RESULT_VARIABLE _fetch_exit
    )
    if (_fetch_exit)
      message(FATAL_ERROR "git fetch failed with exit code ${_fetch_exit}.")
    endif ()
    set(_fetched TRUE)
  endwhile ()
endif()

# If the right SHA1 is not already checked out,
if(NOT head_sha1 STREQUAL _commit_sha1)
  Dout("head_sha1 = \"${head_sha1}\" != _commit_sha1 = \"${_commit_sha1}\" (= GITACHE_CORE_SHA1 = \"${GITACHE_CORE_SHA1}\").")
  if(gitache_core_is_local)
    if(GITACHE_CORE_SHA1 STREQUAL "")
      set(_fatal_message "Local ${PROJECT_NAME} detected.")
    else()
      set(_fatal_message "The local submodule has checked out ${head_sha1}, but ${GITACHE_CORE_SHA1} is requested.")
    endif()
    message(FATAL_ERROR
      " ${_fatal_message} Please set the environment variable GITACHE_CORE_SHA1 to the sha1 that is checked out before calling cmake:\n"
      " \n     export GITACHE_CORE_SHA1=$(git -C ${GITACHE_CORE_SOURCE_DIR} rev-parse HEAD)\n")
  endif()
  # check if the SHA1 is in the local repository.
  execute_process(COMMAND ${git_executable} cat-file -e "${GITACHE_CORE_SHA1}^{commit}"
    COMMAND_ECHO ${gitache_where}
    WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    RESULT_VARIABLE _exit_code
    ERROR_QUIET
  )
  if(_exit_code AND NOT _fetch_done)
    # That SHA1 is not known yet. Fetch it from upstream.
    execute_process(COMMAND ${git_executable} fetch
      COMMAND_ECHO ${gitache_where}
      WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    )
  endif()
  set(gitache_need_include TRUE)
  # Now checkout the needed SHA1.
  execute_process(COMMAND ${git_executable} checkout ${_commit_sha1}
    COMMAND_ECHO ${gitache_where}
    WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
    RESULT_VARIABLE _exit_code
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(_exit_code)
    message(FATAL_ERROR " Failed to checkout ${_commit_sha1} of gitache-core!")
  endif()
  # This file was changed. Reload it!
  unlock_core_directory()
  return()
elseif(_commit_sha1 STREQUAL GITACHE_CORE_SHA1)
  message(STATUS "Gitache-core is already at ${GITACHE_CORE_SHA1}.")
else()
  message(STATUS "Gitache-core is already at \"${GITACHE_CORE_SHA1}\" (${_commit_sha1}).")
endif()

# Now, with ${GITACHE_CORE_SOURCE_DIR} process locked, start the real thing.
list(PREPEND CMAKE_MODULE_PATH "${GITACHE_CORE_SOURCE_DIR}")
set(ERROR_MESSAGE False)
include(main)   # Also uses gitache_where.

# We're finished with gitache-core.
unlock_core_directory()

if(ERROR_MESSAGE)
  # This happens when ERROR_MESSAGE was set in package.cmake.
  message(FATAL_ERROR ${ERROR_MESSAGE})
endif()
