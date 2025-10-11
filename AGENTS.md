# Project Agents.md Guide for OpenAI Codex

This Agents.md file provides comprehensive guidance for OpenAI Codex and other AI agents working with this codebase.

## Project Structure for OpenAI Codex Navigation

This project contains several git submodules. From all existing so called "aicxx" submodules,
ai-math-testsuite is using `cwm4`, `cwds`, `utils`, `threadsafe`, `math` and `cairowindow`.

A lot of these submodules use the other submodules; for example every submodule uses cwds
and utils (they even use eachother).

- `/math`: the git submodule under test by this project (ai-math-testsuite).
- `/src`: various tests for the `math` submodule.
  Most tests are written by Carlo Wood and not really unit tests.
  Their requirement is to compile without errors and run without asserting.
  `/src/math_geometry_tests.cxx` is the actual testsuite that OpenAI Codex should be working on.
  This file must contain unit tests that thoroughly test the classes (in namespace math) Point,
  Vector, Direction, Line and LinePiece.
- `/cmake`: not interesting for AI agents. It contains instructions for gitache on how to download,
  configure and install libcwd (and any other github repository that might be required by the project).
- `/cwds`: not interesting for AI agents. This is a git submodule containing debugging support code
  for C++ projects in general, on top of libcwd.
- `/utils`: this is a stable git submodule containing various C++ utilities that might be used by the
  project and other git submodules.
- `/cwm4`: not interesting for AI agents. Contains build system support.
- `/threadsafe` and `/cairowindow`: not interesting for AI agents. These submodules are only here
  because some of the tests in `/src` written by Carlo Wood use them.

# Build instruction

Immediatly after cloning this repository, and/or after something serious changed in the submodules,
it is normally required to run `autogen.sh`; however - this should not be required in the Codex
environment because all this does is checkout and sync the git submodules which is already taken
care of in the setup and maintenance scripts.

However running `autogen.sh` from the root of the project also prints build instructions.
If build instructions need to be changed then please edit /.build-instructions
for example if something is wrong with a path or configure option.

There is no need to inspect autogen.sh and/or /.build-instructions, just
configure by running what the output of `autogen.sh` tells you.

# What to do in case of an error during configuration/building

If any error occurs, at any stage (running autogen.sh, running cmake for configuration, or
during running cmake --build) then stop immediately and do not try to work around the error.

For example, if a library is missing or a header can't be found, then this requires updating
the build environment and only I can do that. You CAN however suggest how to fix it by
editting /.build-instructions and/or /Codex-environment. The latter contains the instructions
related to, for example, installing missing dependencies.
