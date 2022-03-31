# Want test modules in the build directory.
#unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(gather_test gather_test.F90)
target_link_libraries(gather_test common)
if(USE_CAF)
  add_test(gather ./gather_test)
else()
  add_test(gather ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./gather_test ${MPIEXEC_POSTFLAGS})
endif()

add_executable(scatter_test scatter_test.F90)
target_link_libraries(scatter_test common)
if(USE_CAF)
  add_test(scatter ./scatter_test)
else()
  add_test(scatter ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./scatter_test ${MPIEXEC_POSTFLAGS})
endif()

add_executable(localize_test localize_test.F90)
target_link_libraries(localize_test common)
if(USE_CAF)
  add_test(localize ./localize_test)
else()
  add_test(localize ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./localize_test ${MPIEXEC_POSTFLAGS})
endif()

add_executable(collate_test collate_test.F90)
target_link_libraries(collate_test common)
if(USE_CAF)
  add_test(collate ./collate_test)
else()
  add_test(collate ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./collate_test ${MPIEXEC_POSTFLAGS})
endif()

add_executable(distribute_test distribute_test.F90)
target_link_libraries(distribute_test common)
if(USE_CAF)
  add_test(distribute ./distribute_test)
else()
  add_test(distribute ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./distribute_test ${MPIEXEC_POSTFLAGS})
endif()

set_tests_properties(localize gather scatter distribute collate
    PROPERTIES PROCESSORS 4)
if(USE_CAF)
  set_tests_properties(localize gather scatter distribute collate
      PROPERTIES ENVIRONMENT "NAGFORTRAN_NUM_IMAGES=4")
endif()
