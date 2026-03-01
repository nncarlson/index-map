set(FYPP_SRC
  index_map_type.F90.fypp
  index_map_type-collate_impl.F90.fypp
  index_map_type-distribute_impl.F90.fypp
  index_map_type-gather_offp_impl.F90.fypp
  index_map_type-localize_impl.F90.fypp
  index_map_type-scatter_offp_impl.F90.fypp
)

add_fypp_sources(GEN_SRC
  SOURCES ${FYPP_SRC}
  INCLUDES f90_assert.fpp
)

set(SRC
  f90_assert.F90
  integer_set_type.F90
  integer_map_type.F90
  coarray_collectives.F90
  index_map_type-core_impl.F90
  ${GEN_SRC}
)

add_library(common OBJECT ${SRC})
target_include_directories(common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
