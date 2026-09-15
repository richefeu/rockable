set(cpptools_HDRS "${CMAKE_CURRENT_LIST_DIR}/genesis/cpptools")

# Shared description of the Rockable input-file language (rockable_lang.hpp and
# rockable.lang), read by both confedit and rockedit.
set(rockable_lang_DIR "${CMAKE_CURRENT_LIST_DIR}/common")

add_subdirectory(sweepable)
add_subdirectory(confedit)
add_subdirectory(rockedit)
add_subdirectory(converters)
add_subdirectory(genesis)
add_subdirectory(shapesurvey)
