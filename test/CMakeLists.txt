
## tests
#enable_testing()
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(rockable_test run_regression_test.cpp)
target_link_libraries(rockable_test Rockable_core )

##add_test(NAME rockable_brute_force COMMAND run ${TEST_DIR}/redirection_test_brute_force.conf)
set(WORK_DIR ${CMAKE_CURRENT_BINARY_DIR})
#set(INPUT_512_PARTICLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/512_particles)
#set(TMP_COPY_INPUT_512_PARTICLES_DIR ${WORK_DIR}/input/) 
set(INPUT_518_POLY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/518_poly)
set(TMP_COPY_INPUT_518_POLY_DIR ${WORK_DIR}/input/) 
#set(INPUT_2_PLASTIC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/2_sphs_LudingPlastic)
#set(TMP_COPY_INPUT_2_PLASTIC_DIR ${WORK_DIR}/input/) 
set(INPUT_50_PRISMES_GRAVITY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/50_prismes_gravity)
set(TMP_COPY_50_PRISMES_GRAVITY_DIR ${WORK_DIR}/input/)
#set(INPUT_50_PRISMES_TRIAXIAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/50_prismes_triaxial)
#set(TMP_COPY_50_PRISMES_TRIAXIAL_DIR ${WORK_DIR}/input/)
set(INPUT_2_PRISMES_GRAVITY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/2_prismes_gravity)
set(TMP_COPY_2_PRISMES_GRAVITY_DIR ${WORK_DIR}/input) 
#set(INPUT_500_POLY_DEPOSIT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/500_poly_deposit_in_periodicBox)
#set(TMP_COPY_500_POLY_DEPOSIT_DIR ${WORK_DIR}/input/) 
#set(INPUT_555_SPHERE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/input/555_sphere_periodic)
#set(TMP_COPY_555_SPHERE_DIR ${WORK_DIR}/input/) 

if(CMAKE_BUILD_TYPE MATCHES "Release")
	message("-- The regression tests used are in regression/release")
	set(REGRESSION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/regression/release)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
	message("-- The regression tests used are in regression/debug")
	set(REGRESSION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/regression/debug)
elseif()
	message("-- Error: this mode is not supported by our non-regresssion tests")
endif()

#file(COPY ${INPUT_512_PARTICLES_DIR}  DESTINATION ${TMP_COPY_INPUT_512_PARTICLES_DIR})
file(COPY ${INPUT_518_POLY_DIR}  DESTINATION ${TMP_COPY_INPUT_518_POLY_DIR})
#file(COPY ${INPUT_2_PLASTIC_DIR}  DESTINATION ${TMP_COPY_INPUT_2_PLASTIC_DIR})
file(COPY ${INPUT_50_PRISMES_GRAVITY_DIR}  DESTINATION ${TMP_COPY_50_PRISMES_GRAVITY_DIR})
#file(COPY ${INPUT_50_PRISMES_TRIAXIAL_DIR}  DESTINATION ${TMP_COPY_50_PRISMES_TRIAXIAL_DIR})
file(COPY ${INPUT_2_PRISMES_GRAVITY_DIR}  DESTINATION ${TMP_COPY_2_PRISMES_GRAVITY_DIR})
#file(COPY ${INPUT_500_POLY_DEPOSIT_DIR}  DESTINATION ${TMP_COPY_500_POLY_DEPOSIT_DIR})
#file(COPY ${INPUT_555_SPHERE_DIR}  DESTINATION ${TMP_COPY_555_SPHERE_DIR})


message("test_dir" ${TEST_DIR})
message("work_dir" ${WORK_DIR})
message("regression_dir" ${REGRESSION_DIR})

add_executable(minitest miniTest.cpp)

add_test(NAME toto COMMAND ./minitest)
add_test(NAME 2_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_2_prismes_gravity2.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_2_prismes_gravity.conf)
#if(ROCKABLE_ENABLE_PERIODIC)
#	add_test(NAME 500_poly_deposit_periodic_box COMMAND rockable_test ${TEST_DIR}/redirection_500_poly_deposit_in_periodicBox.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_500_poly_deposit_in_periodicBox.conf)
#	add_test(NAME 555_sphere_periodic COMMAND rockable_test ${TEST_DIR}/redirection_555_sphere_periodic.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_555_sphere_periodic.conf)
#
#	if (OPENMP_FOUND)
#	add_test(NAME 4threads_500_poly_deposit_periodic_box COMMAND rockable_test ${TEST_DIR}/redirection_500_poly_deposit_in_periodicBox.conf -j 4 -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_500_poly_deposit_in_periodicBox.conf)
#	add_test(NAME 4threads_555_sphere_periodic COMMAND rockable_test ${TEST_DIR}/redirection_555_sphere_periodic.conf -j 4 -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_555_sphere_periodic.conf)
#	endif()
#else()
#	# conf101 is created in build/test/512_particles
#	add_test(NAME 518_poly_bruteforce_avalanches_noperiod COMMAND rockable_test ${TEST_DIR}/redirection_test_518_poly_bruteforce_avalanches_noperiod.conf -n ${WORK_DIR}/conf3001 -r ${REGRESSION_DIR}/regression_518_poly_avalanches_noperiod.conf)
#	add_test(NAME 518_poly_linkcells_avalanches_noperiod COMMAND rockable_test ${TEST_DIR}/redirection_test_518_poly_linkcells_avalanches_noperiod.conf -n ${WORK_DIR}/conf3001 -r ${REGRESSION_DIR}/regression_518_poly_avalanches_noperiod.conf)
#	add_test(NAME 518_poly_bruteforce_default_noperiod COMMAND rockable_test ${TEST_DIR}/redirection_test_518_poly_bruteforce_default_noperiod.conf -n ${WORK_DIR}/conf3001 -r ${REGRESSION_DIR}/regression_518_poly_default_noperiod.conf)
#	add_test(NAME 518_poly_linkcells_default_noperiod COMMAND rockable_test ${TEST_DIR}/redirection_test_518_poly_linkcells_default_noperiod.conf -n ${WORK_DIR}/conf3001 -r ${REGRESSION_DIR}/regression_518_poly_default_noperiod.conf)
#
#	add_test(NAME 512_particles_bruteforce_LawAvalanche_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_bruteforce_LawAvalanche_noPeriodicity.conf -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawAvalanche_noPeriodicity.conf)
#	add_test(NAME 512_particles_linkCells_LawAvalanche_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_linkCells_LawAvalanche_noPeriodicity.conf -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawAvalanche_noPeriodicity.conf)
#	add_test(NAME 512_particles_bruteforce_LawHertz_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_bruteforce_LawHertz_noPeriodicity.conf -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawHertz_noPeriodicity.conf)
#	add_test(NAME 512_particles_linkCells_LawHertz_noPeriodicity COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_linkCells_LawHertz_noPeriodicity.conf -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawHertz_noPeriodicity.conf)
#
#
#	add_test(NAME 2_sphs_LudingPlastic_exec COMMAND run ${TEST_DIR}/redirection_2_sphs_LudingPlastic.conf)
#	add_test(NAME 2_sphs_LudingPlastic_reg COMMAND rockable_test ${TEST_DIR}/redirection_2_sphs_LudingPlastic.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_2_plastic_Luding2Phases_noperoid.conf)
#
#	add_test(NAME 50_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_50_prismes_gravity.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_50_prismes_gravity.conf)
#	add_test(NAME 2_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_2_prismes_gravity2.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_2_prismes_gravity.conf)

	#	add_test(NAME 50_prismes_triaxial_reg COMMAND rockable_test ${TEST_DIR}/redirection_50_prismes_triaxial.conf -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_50_prismes_triaxial.conf)
#
	#if (OPENMP_FOUND)
#		add_test(NAME 4threads_512part_bruteforce_Avalanche_noPeriod COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_bruteforce_LawAvalanche_noPeriodicity.conf -j 4 -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawAvalanche_noPeriodicity.conf)
#		add_test(NAME 4threads_512part_linkCells_Avalanche_noPeriod COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_linkCells_LawAvalanche_noPeriodicity.conf -j 4 -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawAvalanche_noPeriodicity.conf)
#		add_test(NAME 4threads_512part_bruteforce_Hertz_noPeriod COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_bruteforce_LawHertz_noPeriodicity.conf -j 4 -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawHertz_noPeriodicity.conf)
#		add_test(NAME 4threads_512part_linkCells_Hertz_noPeriod COMMAND rockable_test ${TEST_DIR}/redirection_test_512_particles_linkCells_LawHertz_noPeriodicity.conf -j 4 -n ${WORK_DIR}/conf101 -r ${REGRESSION_DIR}/regression_512_particles_LawHertz_noPeriodicity.conf)
#		add_test(NAME 4threads_50_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_50_prismes_gravity.conf -j 4 -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_50_prismes_gravity.conf)
#		add_test(NAME 4threads_2_prismes_gravity_reg COMMAND rockable_test ${TEST_DIR}/redirection_2_prismes_gravity.conf -j 4 -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_2_prismes_gravity.conf)
#		add_test(NAME 4threads_50_prismes_triaxial_reg COMMAND rockable_test ${TEST_DIR}/redirection_50_prismes_triaxial.conf -j 4 -n ${WORK_DIR}/conf1 -r ${REGRESSION_DIR}/regression_50_prismes_triaxial.conf)
	#endif()
#endif()
