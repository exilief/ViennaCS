add_custom_target(ViennaCS_Tests ALL)

file(
  GLOB entries
  LIST_DIRECTORIES true
  "*")

foreach(entry ${entries})
  if(NOT IS_DIRECTORY ${entry})
    continue()
  endif()

  message(STATUS "[ViennaCS] Adding Test ${entry}")
  add_subdirectory(${entry})
endforeach()
