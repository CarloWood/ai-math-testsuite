# Gitache-core uses SOURCE_SUBDIR in FetchContent_Declare.
cmake_minimum_required(VERSION 3.18)

# To enable DEBUG messsages, pass -DCMAKE_MESSAGE_LOG_LEVEL=DEBUG to cmake on the command line.
message(DEBUG "DEBUG: Entering ${CMAKE_CURRENT_LIST_FILE}")

if (NOT DEFINED GITACHE_PACKAGES)
  message(FATAL_ERROR "GITACHE_PACKAGES must be set to a list of packages before including gitache/gateway.cmake.")
endif ()

include(FetchContent)  

set(FETCHCONTENT_SOURCE_DIR_GITACHE "gitache" CACHE INTERNAL "" FORCE)  
FetchContent_Declare(gitache)
message(DEBUG "DEBUG: gitache/gateway.cmake: calling `FetchContent_MakeAvailable(gitache)`")
FetchContent_MakeAvailable(gitache)  
message(DEBUG "DEBUG: gitache/gateway.cmake: returned from `FetchContent_MakeAvailable(gitache)`")
