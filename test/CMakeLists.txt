# Want test modules in the build directory.
#unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(reduce_test reduce_test.F90)
target_link_libraries(reduce_test common MPI::MPI_Fortran)
add_test(reduce ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./reduce_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(reduce PROPERTIES PROCESSORS 4)

add_executable(test-broadcast broadcast_test.F90)
target_link_libraries(test-broadcast common MPI::MPI_Fortran)
add_test(broadcast ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test-broadcast ${MPIEXEC_POSTFLAGS})
set_tests_properties(broadcast PROPERTIES PROCESSORS 4)

add_executable(distribute_test distribute_test.F90)
target_link_libraries(distribute_test common MPI::MPI_Fortran)
add_test(distribute ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./distribute_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(distribute PROPERTIES PROCESSORS 4)

add_executable(collate_test collate_test.F90)
target_link_libraries(collate_test common MPI::MPI_Fortran)
add_test(collate ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./collate_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(collate PROPERTIES PROCESSORS 4)

add_executable(test-gather gather_test.F90)
target_link_libraries(test-gather common MPI::MPI_Fortran)
add_test(gather ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test-gather ${MPIEXEC_POSTFLAGS})
set_tests_properties(gather PROPERTIES PROCESSORS 4)

add_executable(test-scatter scatter_test.F90)
target_link_libraries(test-scatter common MPI::MPI_Fortran)
add_test(scatter ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test-scatter ${MPIEXEC_POSTFLAGS})
set_tests_properties(scatter PROPERTIES PROCESSORS 4)

add_executable(test_localize localize_test.F90)
target_link_libraries(test_localize common MPI::MPI_Fortran)
add_test(localize ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test_localize ${MPIEXEC_POSTFLAGS})
set_tests_properties(localize PROPERTIES PROCESSORS 4)

add_executable(index_map_collate_test index_map_collate_test.F90)
target_link_libraries(index_map_collate_test common MPI::MPI_Fortran)
add_test(index_map_collate ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./index_map_collate_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(index_map_collate PROPERTIES PROCESSORS 4)

add_executable(index_map_distribute_test index_map_distribute_test.F90)
target_link_libraries(index_map_distribute_test common MPI::MPI_Fortran)
add_test(index_map_distribute ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./index_map_distribute_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(index_map_distribute PROPERTIES PROCESSORS 4)
