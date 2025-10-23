function(Dout)
  string(JOIN " " _message ${ARGV})
  message(DEBUG "DEBUG: ${_message}")
endfunction()
