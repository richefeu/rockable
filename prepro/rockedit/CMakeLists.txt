# rockedit needs nothing but a C++17 compiler and the standard library. The
# compiler, the standard, the build type, the output directory and the install
# prefix all come from the top-level CMakeLists.txt of Rockable.

# prepro/CMakeLists.txt exports this; fall back to the sibling directory so that
# this CMakeLists can still be configured on its own.
if(NOT rockable_lang_DIR)
  set(rockable_lang_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common")
endif()

add_executable(rockedit src/main.cpp)

# rockable_lang.hpp lives in prepro/common, next to rockable.lang, because
# confedit reads the very same language description.
target_include_directories(rockedit PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${rockable_lang_DIR}
)

# The top-level CMakeLists.txt gathers every binary in the build folder, where
# rockable.lang would not be found; put a copy next to them.
add_custom_command(TARGET rockedit POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${rockable_lang_DIR}/rockable.lang
          $<TARGET_FILE_DIR:rockedit>/rockable.lang
  COMMENT "Copying rockable.lang next to the executables")

install(TARGETS rockedit)

# rockedit and confedit look for rockable.lang next to the executable (and in
# its parent directories), so ship it alongside them.
install(FILES ${rockable_lang_DIR}/rockable.lang DESTINATION bin)
