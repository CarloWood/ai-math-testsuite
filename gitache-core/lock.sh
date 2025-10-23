#! /bin/bash
# Copyright (c) 2020  Carlo Wood.
#
# Lock directory $1 by first acquiring an flock(1) on
# $1/gitache.flock and while holding that writing $2
# to the file $1/gitache.lock.
#
# The directory is considered locked once the flock
# is obtained or (afterwards) when gitache.lock exists.

# $2 should be an unique but descriptive string.
# If a lock cannot be obtained because gitache.lock
# already exists, its content (the unique string)
# is printed to help the user understand what project
# and/or process is holding the lock, and removing the
# lock file if it turns out to be a stale one.

function fatal_error()
{
  echo "$0: FATAL_ERROR: $1" >&2
  exit 1
}

if [ $# -ne 2 -o ! -d "$1" ]; then
  if [ $# -ne 2 ]; then
    echo "Called with $# arguments." >&2
  elif [ ! -d "$1" ]; then
    echo "No such directory: \"$1\"" >&2
  fi
  echo "Usage: $0 <directory> <unique_key>" >&2
  exit 1
fi

#echo "-- Locking directory \"$1\"."
(
  while true; do
    flock 9 || fatal_error "Couldn't lock \"$1/gitache.flock\"!?!" 

    if [ -e "$1/gitache.lock" ]; then
      KEY=$(cat "$1/gitache.lock")
      test "$KEY" != "$2" || fatal_error "Calling lock recursively ($2)?!"
      echo "** Can't lock \"$1\" with key '$2' because \"$1/gitache.lock\" already exists [$KEY]. Sleeping 1 second..."
      flock -u 9
      sleep 1
      continue
    fi

    echo "$2" > "$1/gitache.lock" || fatal_error "Failed to write to \"$1/gitache.lock\""
    break;
  done

) 9>"$1/gitache.flock"
