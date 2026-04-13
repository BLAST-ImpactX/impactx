file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/matching.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/one_turn_map.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES __init__.py matching.py one_turn_map.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/optics)

add_subdirectory(tests)
