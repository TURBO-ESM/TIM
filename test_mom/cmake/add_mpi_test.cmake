
function(add_mpi_test TEST_TARGET_NAME TEST_FILE DATA_DIR)

  set(prefix TIM)
  set(optionalValues "")
  set(singleValues MAX_PES)
  set(multiValues "")

  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} "${optionalValues}" "${singleValues}" "${multiValues}" ${ARGN})

  if(DEFINED TIM_MAX_PES)
    if (TIM_MAX_PES NOT MATCHES "^[1-9][0-9]*$")
      message(FATAL_ERROR "MAX_PES expected a valid positive integer; recieved ${TIM_MAX_PES}")
    endif()
  else()
    set(TIM_MAX_PES 4)
  endif()

  add_executable(${TEST_TARGET_NAME} ${TEST_FILE})
  target_sources(${TEST_TARGET_NAME} PRIVATE $<TARGET_OBJECTS:data_driven_shared_main>)
  target_link_libraries(${TEST_TARGET_NAME} PRIVATE common)

  add_test(NAME MPI_${TEST_TARGET_NAME}
         COMMAND ${MPIEXEC_EXECUTABLE} --tag-output ${MPIEXEC_NUMPROC_FLAG} ${TIM_MAX_PES} ./${TEST_TARGET_NAME} --data-dir=${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR})
  set_tests_properties(MPI_${TEST_TARGET_NAME} PROPERTIES PROCESSORS ${TIM_MAX_PES})
endfunction()
