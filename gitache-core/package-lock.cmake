message(DEBUG "Running package-lock.cmake with \
LOCK_SCRIPT=\"${LOCK_SCRIPT}\"; \
PACKAGE_ROOT=\"${PACKAGE_ROOT}\"; \
LOCK_ID=\"${LOCK_ID}\"; \
STAMP_FILE=\"${STAMP_FILE}\".")

execute_process(
  COMMAND ${LOCK_SCRIPT} ${PACKAGE_ROOT} ${LOCK_ID}
)

if(NOT EXISTS "${STAMP_FILE}")
  file(TOUCH "${STAMP_FILE}")
endif()
if(NOT EXISTS "${STAMP_FILE}-trigger")
  file(TOUCH "${STAMP_FILE}-trigger")
endif()
