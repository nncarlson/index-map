add_custom_target(docs
  COMMAND ${SPHINX_EXECUTABLE} -b html
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
  COMMENT "Building Sphinx documentation"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
