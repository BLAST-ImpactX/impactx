###############################################################################
# These are the files equivalent to the primary C++ implementation - but here
# we define how they will appear in our Python module (aka Python bindings).
#
target_sources(pyImpactX
  PRIVATE
    CompensatedReal.cpp
    distribution.cpp
    elements.cpp
    flatten_rho.cpp
    ImpactX.cpp
    ImpactXParticleContainer.cpp
    lattice.cpp
    ReferenceParticle.cpp
    transformation.cpp
    WakeConvolution.cpp
)
