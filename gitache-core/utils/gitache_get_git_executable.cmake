function(gitache_get_git_executable git_executable_varname)
  if(NOT DEFINED ENV{GITACHE_GIT_EXECUTABLE})
    find_package(Git REQUIRED)
    set(git_executable ${GIT_EXECUTABLE})
    set(git_version ${GIT_VERSION_STRING})
  else()
    set(git_executable $ENV{GITACHE_GIT_EXECUTABLE})
    execute_process(COMMAND ${git_executable} --version
        RESULT_VARIABLE _exit_code
        OUTPUT_VARIABLE output
        ERROR_VARIABLE output
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
    )
    if(_exit_code EQUAL 0)
      if(output MATCHES "^git version [0-9.]+$")
        string(REPLACE "git version " "" git_version "${output}")
      else()
        message(FATAL_ERROR "`${git_executable} --version` gives unexpected output: \"${output}\".")
      endif()
    else()
      message(FATAL_ERROR "Failed to run `${git_executable} --version`: ${_exit_code} ${output}")
    endif()
  endif()

  message(DEBUG "Using git executable: ${git_executable}")
  if(git_version STREQUAL "")
    message(FATAL_ERROR "Failed to determine version of git.")
  endif()

  # 2.13 is the first version that moved to a hardened SHA-1 implementation
  # by default, which isn’t vulnerable to the SHAttered attack.
  # See https://git-scm.com/docs/hash-function-transition/
  if(git_version VERSION_LESS "2.13.0")
    message(FATAL_ERROR "Git version 2.13.0 or higher is required.")
  endif()

  set(${git_executable_varname} ${git_executable} PARENT_SCOPE)
endfunction()
