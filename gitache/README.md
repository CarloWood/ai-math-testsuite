# Usage

## Basic usage

Add to the `CMakeLists.txt` of the main project, after its `project()` line (requires cmake 3.14 or later):

    include(FetchContent)  
    
    set(GITACHE_PACKAGES libcwd_r)  # List of requested packages (example).
    FetchContent_Declare(  
      gitache  
      GIT_REPOSITORY "https://github.com/CarloWood/gitache.git"  
      GIT_TAG        "stable"  
    )
    FetchContent_MakeAvailable(gitache)  

Where, if you're paranoid, you should replace `stable` with its
git commit sha1 (from the `stable` branch, the `master` branch might
be unstable). Both, `stable` and `master` will cause a `git fetch`
every time the project is configured to see if the project was updated,
so for a large project please just use the sha1 of stable (or master,
if that works at that moment).

`GITACHE_PACKAGES` must be set to the list of packages that need
to be downloaded, configured and compiled (if not already in the cache).

Alternatively you can add gitache as a submodule to the root
of your project and replace all of the above with,

    set(GITACHE_PACKAGES libcwd_r)  # List of requested packages (example).
    include(gitache/gateway.cmake)

That will not update gitache automatically but simply use
whatever you have checked out.

Each package listed in `GITACHE_PACKAGES` must have a configuration
file in one of the registered gitache configuration directories.  By
default gitache looks in the cmake cache variable `GITACHE_CONFIGS_DIR`,
which in turn defaults to `${CMAKE_SOURCE_DIR}/cmake/gitache-configs`,
but additional directories can be registered.

A submodule can call `gitache_require_packages` to request (extra)
packages:

```
gitache_require_packages(foobar foobaz
  # optional; defaults to cmake/gitache-configs
  CONFIG_DIRS
    cmake/extra-configs1
    cmake/extra-configs2
)
```

Note that `CONFIG_DIRS` must be the last argument; the package names must
given first. The paths specified are relative the directory that the
current `CMakeLists.txt` file resides in.

See https://github.com/CarloWood/gitache-core for more information.

## Requirements

The minimum cmake version required is 3.14.

Furthermore, the user must have the environment variable
`GITACHE_ROOT` set to a writable directory, or gitache will do
nothing except print the message:

    -- Environment variable GITACHE_ROOT is not set: gitache disabled.   

Finally, the user must know what they are doing: warn them not
to set `LD_LIBRARY_PATH` pointing to a library that is also installed
by gitache. This wouldn't work since the project is linked against
the version compiled by gitache, so at runtime `-rpath` must be honored
and (unfortunately?) `LD_LIBRARY_PATH` has a higher precedence.

As an alternative you could add the following line to your project
(assuming you're creating an executable and your linker honors that flag):

    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")

after which the `-rpath` will get precedence so the value of `LD_LIBRARY_PATH`
in the users environment becomes unimportant.

## Gitache developers

See the README.md of https://github.com/CarloWood/gitache-core
