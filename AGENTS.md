## Build instruction

Paths in this file that start with a '/' are relative to the project root.
The full path to the project root is stored in $REPOBASE - prioritize using $REPOBASE over the full path.

### Configuration Quick Start

- First check CODEX_DIRECTORY; if empty, treat cwd as repo root and avoid redundant cd.
- Prefer path relative to the cwd, otherwise prioritize $REPOBASE and $BUILDDIR in commands instead of absolute paths.
- Run `./autogen.sh`, which lives in the repository root, and follow the configure and build instructions verbatim.
- If the user says `build` he means the `cmake --build` line (only). Otherwise he'll say `configure`.
  If unsure if the project was already configured at all; just assume it is and try to build first (if the user asks you to build the project).

There is no need to inspect autogen.sh, just configure and build by running what the output of `autogen.sh` tells you.

## Helper shell functions

- `findsymbol <symbol> [--kinds=<kinds-list>] [--scope=<scope>] [--subpath=<sub-path>] [--prefix] [--help]`
       --kinds : a comma separated list of kinds: `c|class, s|struct, u|union, f|func|method, m|member|field, v|var|variable, t|typedef|using, g|Enum, e|enum, n|namespace, d|macro|define`
       --subpath : filters on a contiguous subsequence of pathname components in the output location.
       --scope : filters on scopes that begin with given substring.
       --prefix : also match symbols that begin with `<symbol>`.
  Paths are abbreviated by being either relative to the cwd or the path environment variables (`REPOBASE`, `BUILDDIR`, `CODEX_WORKSPACE`).
  `findsymbol` uses `$BUILDDIR/tags` which can be refreshed by running `make ctags` (only do that if it starts failing).

- Always use `findsymbol` to locate the definition of a class, (member)function, enum etc. Fallback to `rg` if it doesn't find what you are looking for.
- Use `rg` only if you need to find *every* occurrence including comments, usage, forward declarations etc.

For example, to find all locations where a class (or struct) `Vector` is defined in the project, you will use: `findsymbol Vector --kinds=c,s`.
The reason for this tool is to limit the generated output; following up with an `rg` command is strongly discouraged, as that would defeat the whole purpose of `findsymbol`.

## What to do in case of an error during configuration/building

- If any error occurs, at any stage (running `./autogen.sh`, running cmake for configuration,
or during running `cmake --build ...`) then STOP immediately and do NOT try to work around the error.

## Project Structure

This project contains several git submodules. From all existing (so called) "aicxx" submodules,
`ai-math-testsuite` is using `cwm4`, `cwds`, `utils`, `threadsafe`, `math` and `cairowindow`.

A lot of these submodules use the other submodules; for example every submodule uses cwds
and utils (those even use eachother).

### Overview of subdirectories and submodules

- `/math`: the git submodule under test by this project (ai-math-testsuite).
- `/src`: various tests for the `math` submodule.
  Most tests are written by Carlo Wood and not really unit tests.
  Their requirement is to compile without errors and run without asserting.
  `/src/math_geometry_tests.cxx` is the actual testsuite that OpenAI Codex should be working on.
  This file should contain unit tests that thoroughly test the classes (in namespace `math`):
  Point, Vector, Direction, Line and LinePiece.

### Remaining subdirectories that are not of interest to AI Agents

- `/cmake`: Contains instructions for gitache on how to download, configure and install libcwd (and any other github repository that might be required by the project).
- `/cwds`: This is a git submodule containing debugging support code for C++ projects in general, on top of libcwd.
- `/utils`: this is a stable git submodule containing various C++ utilities that might be used by the project and other git submodules.
- `/cwm4`: Contains build system support.
- `/threadsafe` and `/cairowindow`: These submodules are only here because some of the tests in `/src` written by Carlo Wood use them.

## Test instructions

### When working on a test in /src

1. Whenever a .cxx file in /src was changed, run the resulting executable to check if it works.
2. If there are errors, then the test failed and the applied changes must be wrong, or you
   discovered a bug in /math. Determine which is the case by checking the test and verifying
   that it is correct, or fix it.
3. If you decide that the test failed because of a bug in /math - then do NOT attempt to fix it
   unless the user explicitly requests to do so.

### When working, by request, on files in /math.

1. After making the change check if the project still builds. If so, run `math_geometry_tests`.
2. If `math_geometry_tests` fails, then make sure this is not an error in the test itself.
3. If `math_geometry_tests.cxx` in `/src` is doing the right thing, then there must be a problem
   with `/math`. If the problem is caused by the last changes just made, then fix them. Otherwise,
   only report back on what fails where you think the bug is.
