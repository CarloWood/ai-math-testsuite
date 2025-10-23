## Gitache package configuration files

Each package listed in `GITACHE_PACKAGES` must come with a config
file in `cmake/gitache-configs` in the root of the main project.

For example if `myproject/CMakeLists.txt` contains the line

    set(GITACHE_PACKAGES foobar blah)

then the files `myproject/cmake/gitache-configs/foobar.cmake` and
`myproject/cmake/gitache-configs/blah.cmake` must exist.

Their contents should be, for example,

    gitache_config(
      GIT_REPOSITORY
        "https://github.com/CarloWood/libcwd.git"
      GIT_TAG
        "master"
      CMAKE_ARGS
        "-DEnableLibcwdAlloc:BOOL=OFF -DEnableLibcwdLocation:BOOL=ON"
    )

Again, you are encouraged to use a SHA1 for the `GIT_TAG` for reasons
of reproduciblity and stability, but a tag or branch work too. A branch
will NOT result in a `git fetch` every time the project is (re)configured
(see below); a tag or sha1 only if that tag or sha1 doesn't exist (yet).

Currently supported are:

    GIT_REPOSITORY     - Git URL to clone / fetch.
    GIT_TAG            - The thing to checkout with git.
    BOOTSTRAP_COMMAND  - Command to run before configuration.
    CMAKE_ARGS         - Arguments to pass to cmake (cmake projects only).
    CONFIGURE_ARGS     - Arguments to pass to configure (autotools projects only).

Both, `CMAKE_ARGS` and `CONFIGURE_ARGS` are intended to pass package configurations
only. Do not attempt to pass an install prefix.

See https://github.com/CarloWood/ai-statefultask-testsuite/tree/master/cmake/gitache-configs
for example configurations.

Also, see https://github.com/CarloWood/ai-statefultask-testsuite/blob/master/CMakeLists.txt
for an example of a project that uses gitache.

## Using a branch instead of SHA1

In order to 'git update' a `GIT_TAG` branch (like "master" in the example above),
you need to remove the 'DONE' file from the ROOT of the gitache project
(the value of the ROOT is printed during configuration when cmake DEBUG
output is turned on (add `-DCMAKE_MESSAGE_LOG_LEVEL=DEBUG` to cmake when
configuring));

For example, the above would print something like

    -- DEBUG: libcwd_r: libcwd_r_ROOT = "/opt/gitache/libcwd_r/1289adbf8a106a19b815fdf539814aaf5ffdf1eac652290a38af5677e08f34e8".

then, in order to pull the `master` again, remove the file
`/opt/gitache/libcwd_r/1289adbf8a106a19b815fdf539814aaf5ffdf1eac652290a38af5677e08f34e8/DONE`.

## Gitache developers

If you want to experiment with making changes to gitache itself,
and want to support both (a gitache submodule being present, or not) you
can use the first version given under `Basic usage` (see https://github.com/CarloWood/gitache)
but add the following immediately after the `include(FetchContent)` line:

    # If a local gitache submodule is present then use that rather than downloading one.
    if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/gitache/.git)
      # This will disable the use of the GIT_REPOSITORY/GIT_TAG below, and disable the
      # FetchContent- download and update step. Instead, use the gitache submodule as-is.
      set(FETCHCONTENT_SOURCE_DIR_GITACHE "${CMAKE_CURRENT_LIST_DIR}/gitache" CACHE INTERNAL "" FORCE)
    endif ()

This more verbose version allows you to clone gitache in the root
of your project, without adding it to the project as a submodule, and
make changes to it as you see fit, to experiment, while users of the
project don't need a local gitache directory.

Cloning a local gitache-core in the root of the project is supported
by default. This too will disable any downloading of gitache-core and
use the local gitache-core directory as-is.
