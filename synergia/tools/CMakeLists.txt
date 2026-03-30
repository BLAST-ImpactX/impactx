file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/three_bump.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/print_map_o3.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES __init__.py three_bump.py print_map_o3.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/tools)

add_subdirectory(tests)
