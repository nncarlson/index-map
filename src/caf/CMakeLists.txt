set(SRC
  f90_assert.F90
  integer_set_type.F90
  integer_map_type.F90
  coarray_collectives.F90
  index_map_type.F90
  index_map_type-gather_offp_impl.F90
  index_map_type-scatter_offp_impl.F90
  index_map_type-distribute_impl.F90
  index_map_type-collate_impl.F90
  index_map_type-localize_impl.F90
)

add_library(common OBJECT ${SRC})
target_include_directories(common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
