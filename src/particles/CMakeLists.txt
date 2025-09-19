target_sources(lib
  PRIVATE
    ChargeDeposition.cpp
    CollectLost.cpp
    ImpactXParticleContainer.cpp
    Push.cpp
    SplitEqually.cpp
)

add_subdirectory(spacecharge)
add_subdirectory(transformation)
add_subdirectory(wakefields)
