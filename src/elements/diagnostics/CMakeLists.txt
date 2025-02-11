target_sources(lib
  PRIVATE
    AdditionalProperties.cpp
    BeamMonitor.cpp
)

if(ImpactX_OPENPMD)
    target_sources(lib
      PRIVATE
        openPMD.cpp
    )
endif()
