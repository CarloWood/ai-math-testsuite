# Allowing git submodules to request gitache packages

## Background

The top-level project currently drives gitache by setting `GITACHE_PACKAGES` before including `cwm4/cmake/StableGitache`, which in turn fetches gitache and executes `gitache/CMakeLists.txt`.
That entry point stops immediately when `GITACHE_PACKAGES` is empty and hardcodes a single `GITACHE_CONFIGS_DIR` rooted at the super-project (`${CMAKE_SOURCE_DIR}/cmake/gitache-configs`).
Inside `gitache-core/main.cmake` the variable is dereferenced once in a `foreach` loop; each package is resolved exactly one time while `main.cmake` is executing.
As soon as that call returns there is no exported API that would let later code (for example a git submodule) contribute additional package names or configuration directories.

When we move the `versor` configuration file into the `math` submodule the current design therefore fails: the submodule cannot register its config search path, nor can it ask gitache to fetch a package after the initial pass.

## Goals

We want a git submodule to be able to:

1. Publish additional configuration directories (for example `math/cmake/gitache-configs`).
2. Request one or more gitache packages after it is added to the build (for example `gitache` should fetch and expose `versor`).
3. Reuse the existing locking/build/install logic so packages are still configured exactly once and exposed through `<package>_ROOT` variables.

## Suggested gitache changes

1. **Turn the one-shot package loop into a reusable function.**  Factor the body of the `foreach (gitache_package ${GITACHE_PACKAGES})` loop in `gitache-core/main.cmake` into a helper such as `gitache_process_package(package_name)`.  The helper should continue to read the package configuration, compute the hash, lock the installation directory and call `package.cmake` exactly as it does now.  Track processed packages in a global property (for example `GITACHE_PROCESSED_PACKAGES`) so repeated requests are ignored.
2. **Maintain a search path list for configuration files.**  Replace the single `GITACHE_CONFIGS_DIR` cache path with a list (for backwards compatibility, initialize the list with the current default).  Provide a function like `gitache_register_config_dir(<path>)` that appends an absolute directory to the list while normalizing duplicates.  `gitache_process_package` should scan this list until it finds `<path>/<package>.cmake`.
3. **Export a public registration API.**  After defining the helpers, add a function `gitache_require_packages(<pkg> [...])` that (a) optionally accepts a `CONFIG_DIR` keyword for ad-hoc directories, (b) calls `gitache_register_config_dir` when necessary, and (c) invokes `gitache_process_package` for each package.  Call this function once from `main.cmake` with the packages the super project provided so existing projects keep working.
4. **Preserve thread-safety and logging.**  The helper should keep using the existing `lock_directory`, `unlock_directory`, `gitache_log_level` and debug helpers so parallel configurations still serialize package installation and the diagnostics remain unchanged.
5. **Update documentation.**  Extend `gitache/README.md` and `gitache-core/README.md` with a short section that shows how a submodule can call `gitache_register_config_dir` and `gitache_require_packages`.

## Suggested changes in `math`

Once gitache exports the new helpers, the `math` submodule can own the dependency with the following pattern:

1. Move `cmake/gitache-configs/versor.cmake` into `math/cmake/gitache-configs/versor.cmake`.
2. At the top of `math/CMakeLists.txt`, after including `AICxxProject`, call
   ```cmake
   gitache_register_config_dir("${CMAKE_CURRENT_LIST_DIR}/cmake/gitache-configs")
   gitache_require_packages(versor)
   ```
   so the submodule asks gitache for `versor` when it is added to the build.
3. Add `versor::vsr` to the public link interface of `AICxx::math` so consumers inherit the dependency automatically:
   ```cmake
   target_link_libraries(math_ObjLib
     PUBLIC
       Eigen3::Eigen
       versor::vsr
     PRIVATE
       AICxx::cwds
       AICxx::utils)
   ```

With these changes the super project no longer needs to list `versor` in its `GITACHE_PACKAGES`, while any project that pulls in the `math` submodule automatically receives the versor include directories and link libraries through `AICxx::math`.
