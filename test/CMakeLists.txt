## tests

set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# The references hold for every build type and platform: rockable compares
# conf-files with a relative tolerance (see compareConf in src/Apps/run.cpp).
set(REGRESSION_DIR ${TEST_DIR}/regression)

# add_regression_test(<name> <case> <conf> <result> <reference>)
#
# Runs rockable on input/<case>/<conf> in a folder of its own, BUILD/test/<name>,
# and compares the conf-file <result> it writes with regression/<reference>.
# A setup test first removes the conf-files of a previous run, so that a
# simulation stopping early cannot be compared with a stale result.
function(add_regression_test NAME CASE CONF RESULT REFERENCE)
  set(RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${NAME})
  file(COPY ${TEST_DIR}/input/${CASE} DESTINATION ${RUN_DIR}/input)
  file(WRITE ${RUN_DIR}/redirection.txt "redirection\ninput/${CASE}/\n${CONF}\n")

  add_test(NAME ${NAME}_clean COMMAND rockable -c WORKING_DIRECTORY ${RUN_DIR})
  set_tests_properties(${NAME}_clean PROPERTIES FIXTURES_SETUP ${NAME})

  add_test(NAME ${NAME}
           COMMAND rockable redirection.txt -n ${RESULT} -r ${REGRESSION_DIR}/${REFERENCE}
           WORKING_DIRECTORY ${RUN_DIR})
  set_tests_properties(${NAME} PROPERTIES FIXTURES_REQUIRED ${NAME})
endfunction()

# Unit test of the r-shape signed distance function (skin-mesh, phase 0)
add_executable(test_ShapeSDF ${TEST_DIR}/ShapeSDF/test_ShapeSDF.cpp)
target_link_libraries(test_ShapeSDF PUBLIC Rockable_core)
add_test(NAME ShapeSDF COMMAND test_ShapeSDF)

# Unit test of the r-shape skin mesher (skin-mesh, phase 1)
add_executable(test_ShapeMesh ${TEST_DIR}/ShapeSDF/test_ShapeMesh.cpp)
target_link_libraries(test_ShapeMesh PUBLIC Rockable_core)
add_test(NAME ShapeMesh COMMAND test_ShapeMesh ${TEST_DIR}/input/518_poly/shape.shp)

add_regression_test(2_prismes_gravity
  2_prismes_gravity 2_prismes_gravity.conf conf1 regression_2_prismes_gravity.conf)
add_regression_test(50_prismes_gravity
  50_prismes_gravity 50_prismes_gravity.conf conf1 regression_50_prismes_gravity.conf)
add_regression_test(518_poly_bruteforce_default_noPeriodicity
  518_poly 518_poly_bruteforce_default_noperiod.conf conf3001 regression_518_poly_default_noperiod.conf)
add_regression_test(512_particles_bruteforce_LawDefault_noPeriodicity
  512_particles 512_particles_bruteforce_LawDefault_noPeriodicity.conf conf101
  regression_512_particles_LawDefault_noPeriodicity.conf)
add_regression_test(512_particles_bruteforce_LawAvalanche_noPeriodicity
  512_particles 512_particles_bruteforce_LawAvalanche_noPeriodicity.conf conf101
  regression_512_particles_LawAvalanche_noPeriodicity.conf)

# Two cubes glued by stickBCM and pulled apart until their interface breaks (BCM law).
# The restart from the conf-file at t = 0.006 must reach the same state, so both use the same reference.
add_regression_test(bcm_two_cubes_stick
  bcm_two_cubes bcm_two_cubes_stick.conf conf15 regression_bcm_two_cubes.conf)
add_regression_test(bcm_two_cubes_restart
  bcm_two_cubes bcm_two_cubes_restart.conf conf15 regression_bcm_two_cubes.conf)

if(ROCKABLE_ENABLE_SOFT_PARTICLES)
  # A soft cube squeezed between two platens, elastic and viscoelastic
  add_regression_test(soft_elastic_cube
    soft_cube soft_elastic_cube.conf conf2 regression_soft_elastic_cube.conf)
  add_regression_test(soft_viscous_cube
    soft_cube soft_viscous_cube.conf conf2 regression_soft_viscous_cube.conf)
endif()
