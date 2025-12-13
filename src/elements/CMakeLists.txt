target_sources(lib
  PRIVATE
    Aperture.cpp
    PolygonAperture.cpp
    Programmable.cpp
    Source.cpp
)

add_subdirectory(diagnostics)
add_subdirectory(transformation)
