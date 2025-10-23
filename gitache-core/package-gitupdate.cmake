message(DEBUG "Running package-gitupdate.cmake with \
GIT_EXE = \"${GIT_EXE}\"; \
GIT_TAG = \"${GIT_TAG}\"; \
PACKAGE_ROOT = \"${PACKAGE_ROOT}\"; \
STAMP_FILE = \"${STAMP_FILE}\"; \
PACKAGE_NAME = \"${PACKAGE_NAME}\"; \
PWD = $ENV{PWD}.")

execute_process(
  COMMAND "${GIT_EXE}" rev-list --max-count=1 HEAD
  RESULT_VARIABLE _exit_code
  OUTPUT_VARIABLE head_sha
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
if(_exit_code)
  message(FATAL_ERROR "${PACKAGE_NAME}: failed to get the hash for HEAD.")
endif()

execute_process(
  COMMAND "${GIT_EXE}" show-ref ${GIT_TAG}
  OUTPUT_VARIABLE show_ref_output
  )
# If a remote ref is asked for, which can possibly move around,
# we must always do a fetch and checkout.
if("${show_ref_output}" MATCHES "remotes")
  set(is_remote_ref 1)
else()
  set(is_remote_ref 0)
endif()

# Tag is in the form <remote>/<tag> (i.e. origin/master) we must strip
# the remote from the tag.
if("${show_ref_output}" MATCHES "refs/remotes/${GIT_TAG}")
  string(REGEX MATCH "^([^/]+)/(.+)$" _unused "${GIT_TAG}")
  set(git_remote "${CMAKE_MATCH_1}")
  set(git_tag "${CMAKE_MATCH_2}")
else()
  set(git_remote "origin")
  set(git_tag "${GIT_TAG}")
endif()

# This will fail if the tag does not exist (it probably has not been fetched
# yet).
execute_process(
  COMMAND "${GIT_EXE}" rev-list --max-count=1 ${GIT_TAG}
  RESULT_VARIABLE _exit_code
  OUTPUT_VARIABLE tag_sha
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# Is the hash checkout out that we want?
if(_exit_code OR NOT tag_sha STREQUAL head_sha)
  message(STATUS "Need rebuild ${PACKAGE_NAME}: requested commit is not checked out.")
  set(_need_touch_stamp_file TRUE)
else()
  set(_need_touch_stamp_file FALSE)
endif()
if(is_remote_ref OR _need_touch_stamp_file)
  execute_process(
    COMMAND "${GIT_EXE}" fetch
    RESULT_VARIABLE _exit_code
    )
  if(_exit_code)
    message(FATAL_ERROR "${PACKAGE_NAME}: failed to fetch repository '$ENV{PWD}'.")
  endif()

  if(is_remote_ref)
    # Check if stash is needed
    execute_process(
      COMMAND "${GIT_EXE}" status --porcelain
      RESULT_VARIABLE _exit_code
      OUTPUT_VARIABLE repo_status
      )
    if(_exit_code)
      message(FATAL_ERROR "${PACKAGE_NAME}: failed to get the status.")
    endif()
    string(LENGTH "${repo_status}" need_stash)

    # If not in clean state, stash changes in order to be able to perform git pull --rebase
    if(need_stash)
      message(STATUS "Need rebuild ${PACKAGE_NAME}: the source tree has been editted manually!")
      set(_need_touch_stamp_file TRUE)  # Once editted, we don't know if more changes were made (without calculating a hash of the whole source tree anyway).
      execute_process(
        COMMAND "${GIT_EXE}" stash save --all;--quiet
        RESULT_VARIABLE _exit_code
        )
      if(_exit_code)
        message(FATAL_ERROR "${PACKAGE_NAME}: failed to stash changes.")
      endif()
    endif()

    # Pull changes from the remote branch
    execute_process(
      COMMAND "${GIT_EXE}" rebase ${git_remote}/${git_tag}
      RESULT_VARIABLE _exit_code
      )
    if(_exit_code)
      # Rebase failed: Restore previous state.
      execute_process(
        COMMAND "${GIT_EXE}" rebase --abort
      )
      if(need_stash)
        execute_process(
          COMMAND "${GIT_EXE}" stash pop --index --quiet
          )
      endif()
      message(FATAL_ERROR "\nFailed to rebase in '$ENV{PWD}'.\nYou will have to resolve the conflicts manually.")
    endif()

    execute_process(
      COMMAND "${GIT_EXE}" rev-list --max-count=1 HEAD
      RESULT_VARIABLE _exit_code
      OUTPUT_VARIABLE new_head_sha
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    if(_exit_code)
      message(FATAL_ERROR "${PACKAGE_NAME}: failed to get the hash for HEAD after rebase.")
    endif()
    if(NOT new_head_sha STREQUAL head_sha)
      message(STATUS "Need rebuild ${PACKAGE_NAME}: HEAD changed from ${head_sha} to ${new_head_sha}.")
      set(_need_touch_stamp_file TRUE)
    endif()

    if(need_stash)
      execute_process(
        COMMAND "${GIT_EXE}" stash pop --index --quiet
        RESULT_VARIABLE _exit_code
        )
      if(_exit_code)
        # Stash pop --index failed: Try again dropping the index
        execute_process(
          COMMAND "${GIT_EXE}" reset --hard --quiet
          RESULT_VARIABLE _exit_code
          )
        execute_process(
          COMMAND "${GIT_EXE}" stash pop --quiet
          RESULT_VARIABLE _exit_code
          )
        if(_exit_code)
          # Stash pop failed: Restore previous state.
          execute_process(
            COMMAND "${GIT_EXE}" reset --hard --quiet ${head_sha}
          )
          execute_process(
            COMMAND "${GIT_EXE}" stash pop --index --quiet
          )
          message(FATAL_ERROR "\n${PACKAGE_NAME}: failed to unstash changes in: '$ENV{PWD}'.\nYou will have to resolve the conflicts manually.")
        endif()
      endif()
    endif()
  else()
    execute_process(
      COMMAND "${GIT_EXE}" checkout ${GIT_TAG}
      RESULT_VARIABLE _exit_code
      )
    if(_exit_code)
      message(FATAL_ERROR "${PACKAGE_NAME}: failed to checkout tag: '${GIT_TAG}'.")
    endif()
  endif()

  set(init_submodules TRUE)
  if(init_submodules)
    execute_process(
      COMMAND "${GIT_EXE}" submodule update --recursive --init 
      RESULT_VARIABLE _exit_code
      )
  endif()
  if(_exit_code)
    message(FATAL_ERROR "${PACKAGE_NAME}: failed to update submodules in: '$ENV{PWD}'.")
  endif()
endif()

if(_need_touch_stamp_file OR NOT EXISTS ${STAMP_FILE})
  file(TOUCH ${STAMP_FILE})
  file(COPY ${STAMP_FILE} DESTINATION ${PACKAGE_ROOT})
endif()
