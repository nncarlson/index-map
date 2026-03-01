find_program(SPHINX_BUILD_EXECUTABLE NAMES sphinx-build)
if(NOT SPHINX_BUILD_EXECUTABLE)
  message(FATAL_ERROR
      "Could not find 'sphinx-build'.\n"
      "Install Sphinx and rerun CMake.")
endif()

execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import furo"
  RESULT_VARIABLE FURO_IMPORT_STATUS
  OUTPUT_QUIET
  ERROR_QUIET
)
if(NOT FURO_IMPORT_STATUS EQUAL 0)
  message(FATAL_ERROR
      "Could not import the Python package 'furo'.\n"
      "Install it with 'python3 -m pip install furo' and rerun CMake.")
endif()

set(SPHINX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(SPHINX_DOCTREE_DIR ${CMAKE_CURRENT_BINARY_DIR}/doctrees)
set(SPHINX_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/conf.py
  ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
  ${CMAKE_CURRENT_SOURCE_DIR}/index_map_type.rst
)

add_custom_command(
  OUTPUT ${SPHINX_OUTPUT_DIR}/index.html
  COMMAND ${CMAKE_COMMAND} -E make_directory ${SPHINX_OUTPUT_DIR}
  COMMAND ${SPHINX_BUILD_EXECUTABLE}
          -b html
          -d ${SPHINX_DOCTREE_DIR}
          ${SPHINX_SOURCE_DIR}
          ${SPHINX_OUTPUT_DIR}
  DEPENDS ${SPHINX_SOURCES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)

add_custom_target(doc_html ALL DEPENDS ${SPHINX_OUTPUT_DIR}/index.html)
