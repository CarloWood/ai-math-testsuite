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
_gitache_lock_core_directory()
set(CORE_DIRECTORY_LOCKED True)

while (TRUE) # So we can break out of it after setting ERROR_MESSAGE.
  set(ERROR_MESSAGE False)

# This might not even be a git repository; in that case there is nothing to update.
if (NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/.git")
  Dout("Not a git repository: skipping the update of gitache-core.")
else (NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/.git")
  if (NOT gitache_core_is_local)
    set(ERROR_MESSAGE "Can't find a .git directory in ${CMAKE_CURRENT_LIST_DIR}.")
    break()
  endif ()

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
        set(ERROR_MESSAGE "The environment variable GITACHE_CORE_SHA1 is set to \"${GITACHE_CORE_SHA1}\", which does not exist in the gitache-core repository.")
        break()
      endif ()
      # Fetch once, then loop to retry the rev-parse.
      execute_process(
        COMMAND ${git_executable} fetch
        COMMAND_ECHO ${gitache_where}
        WORKING_DIRECTORY ${GITACHE_CORE_SOURCE_DIR}
        RESULT_VARIABLE _fetch_exit
      )
      if (_fetch_exit)
        set(ERROR_MESSAGE "git fetch failed with exit code ${_fetch_exit}.")
        break()
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
      set(ERROR_MESSAGE
        " ${_fatal_message} Please set the environment variable GITACHE_CORE_SHA1 to the sha1 that is checked out before calling cmake:\n"
        " \n     export GITACHE_CORE_SHA1=$(git -C ${GITACHE_CORE_SOURCE_DIR} rev-parse HEAD)\n")
      break()
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
      set(ERROR_MESSAGE " Failed to checkout ${_commit_sha1} of gitache-core!")
      break()
    endif()
    # This file was changed. Reload it!
    _gitache_unlock_core_directory()
    return()
  elseif(_commit_sha1 STREQUAL GITACHE_CORE_SHA1)
    message(STATUS "Gitache-core is already at ${GITACHE_CORE_SHA1}.")
  else()
    message(STATUS "Gitache-core is already at \"${GITACHE_CORE_SHA1}\" (${_commit_sha1}).")
  endif()

endif (NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/.git")

  # No error occured. Just continue running (leave this fake "loop").
  break()
endwhile()

if (NOT ERROR_MESSAGE)
  # Now, with ${GITACHE_CORE_SOURCE_DIR} process locked, start the real thing.
  list(PREPEND CMAKE_MODULE_PATH "${GITACHE_CORE_SOURCE_DIR}")
  include(main)   # Also uses gitache_where.

function(gitache_register_config_dir)
  _gitache_lock_core_directory()
  set(CORE_DIRECTORY_LOCKED True)
  gitache_register_config_dir_locked(${ARGV})
  _gitache_unlock_core_directory()
  set(CORE_DIRECTORY_LOCKED False)

  if (ERROR_MESSAGE)
    message(FATAL_ERROR ${ERROR_MESSAGE})
  endif ()
endfunction()

function(gitache_require_packages)
  # This function can be called from a submodules CMakeLists.txt;
  # we have to load all required global properties that were set by the super project!
  get_property(_core_dir GLOBAL PROPERTY GITACHE_CORE_SOURCE_DIR)
  set(GITACHE_CORE_SOURCE_DIR ${_core_dir})

  set(options)
  set(oneValueArgs)
  set(multiValueArgs CONFIG_DIRS)
  cmake_parse_arguments(gitache_require "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(_packages ${gitache_require_UNPARSED_ARGUMENTS})

  _gitache_lock_core_directory()
  set(CORE_DIRECTORY_LOCKED True)

  foreach(_dir ${gitache_require_CONFIG_DIRS})
    _gitache_register_config_dir_locked("${_dir}")
    if (ERROR_MESSAGE)
      set(ERROR_MESSAGE "${ERROR_MESSAGE}" PARENT_SCOPE)
      break()
    endif ()
  endforeach()

  if (NOT ERROR_MESSAGE)
    _gitache_require_packages_locked(${_packages})
  endif ()

  _gitache_unlock_core_directory()
  set(CORE_DIRECTORY_LOCKED False)

  if (ERROR_MESSAGE)
    message(FATAL_ERROR ${ERROR_MESSAGE})
  endif ()

  # Propagate the package_ROOT's to the parent.
  foreach(_package ${_packages})
    set(${_package}_ROOT "${${_package}_ROOT}" PARENT_SCOPE)
  endforeach()
endfunction()

endif (NOT ERROR_MESSAGE)

# We're finished with gitache-core.
_gitache_unlock_core_directory()
set(CORE_DIRECTORY_LOCKED False)

if (ERROR_MESSAGE)
  message(FATAL_ERROR ${ERROR_MESSAGE})
endif ()
