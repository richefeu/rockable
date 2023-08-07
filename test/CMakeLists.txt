
## tests
#enable_testing()
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(rockable_test run_regression_test.cpp)
target_link_libraries(rockable_test Rockable_core )

##add_test(NAME rockable_brute_force COMMAND run ${TEST_DIR}/redirection_test_brute_force.conf)
set(WORK_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(INPUT_512_PARTICLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/512_particles)
set(TMP_COPY_INPUT_512_PARTICLES_DIR ${WORK_DIR}/input/) 
set(INPUT_518_POLY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/518_poly)
set(TMP_COPY_INPUT_518_POLY_DIR ${WORK_DIR}/input/) 
set(INPUT_50_PRISMES_GRAVITY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/50_prismes_gravity)
set(TMP_COPY_50_PRISMES_GRAVITY_DIR ${WORK_DIR}/input/)
set(INPUT_2_PRISMES_GRAVITY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/2_prismes_gravity)
set(TMP_COPY_2_PRISMES_GRAVITY_DIR ${WORK_DIR}/input) 


if(CMAKE_BUILD_TYPE MATCHES "Release")
	message("-- The regression tests used are in regression/release")
	set(REGRESSION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/regression/release)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
	message("-- The regression tests used are in regression/debug")
	set(REGRESSION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/regression/debug)
elseif()
	message("-- Error: this mode is not supported by our non-regresssion tests")
endif()

file(COPY ${INPUT_512_PARTICLES_DIR}  DESTINATION ${TMP_COPY_INPUT_512_PARTICLES_DIR})
file(COPY ${INPUT_518_POLY_DIR}  DESTINATION ${TMP_COPY_INPUT_518_POLY_DIR})
file(COPY ${INPUT_50_PRISMES_GRAVITY_DIR}  DESTINATION ${TMP_COPY_50_PRISMES_GRAVITY_DIR})
file(COPY ${INPUT_2_PRISMES_GRAVITY_DIR}  DESTINATION ${TMP_COPY_2_PRISMES_GRAVITY_DIR})

add_test(NAME 2_prismes_gravity_exec COMMAND rockable ${TEST_DIR}/redirection_2_prismes_gravity.conf)
add_test(NAME 2_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_2_prismes_gravity.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_2_prismes_gravity.conf)
add_test(NAME 50_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_50_prismes_gravity.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_50_prismes_gravity.conf)
add_test(NAME 518_poly_bruteforce_default_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_518_poly_bruteforce_default_noperiod.conf -n ${WORK_DIR}/conf3001 -r ${REGRESSION_DIR}/regression_518_poly_default_noperiod.conf)
add_test(NAME 518_poly_bruteforce_avalanches_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_518_poly_bruteforce_avalanches_noperiod.conf -n ${WORK_DIR}/conf3001 -r ${REGRESSION_DIR}/regression_518_poly_avalanches_noperiod.conf)
add_test(NAME 512_particles_bruteforce_LawDefault_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_bruteforce_LawDefault_noPeriodicity.conf -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawDefault_noPeriodicity.conf)
add_test(NAME 512_particles_bruteforce_LawAvalanche_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_bruteforce_LawAvalanche_noPeriodicity.conf -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawAvalanche_noPeriodicity.conf)


