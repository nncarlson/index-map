set(SRC
  f90_assert.F90
  index_map_type.F90
  index_map_type-gather_impl.F90
  index_map_type-scatter_impl.F90
#  index_map_type-localize_impl.F90
  parallel_communication.F90
  parallel_communication-broadcast_impl.F90
  parallel_communication-distribute_impl.F90
  parallel_communication-collate_impl.F90
  parallel_communication-reduce_impl.F90
)
add_library(common OBJECT ${SRC})
target_link_libraries(common MPI::MPI_Fortran)
target_include_directories(common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
