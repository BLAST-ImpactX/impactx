from impactx import ImpactX, distribution, elements

sim = ImpactX()

sim.particle_shape = 2

sim.slice_step_diagnostics = True

sim.init_grids()

# Initialize particle beam
kin_energy_MeV = 2000.0
bunch_charge_C = 1e-09
npart = 10000

pc = sim.particle_container()

# Reference particle
ref = pc.ref_particle()
ref.set_charge_qe(-1).set_mass_MeV(0.51099895).set_kin_energy_MeV(kin_energy_MeV)

distr = distribution.Waterbag(
    lambdaX=3.998488477e-05,
    lambdaY=3.998488477e-05,
    lambdaT=0.001,
    lambdaPx=2.662353876e-05,
    lambdaPy=2.662353876e-05,
    lambdaPt=0.002,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)

ns = 25
lattice_configuration = [
    elements.Drift(ds=0.25, nslice=ns, name="drift1"),
    elements.Quad(ds=1.0, k=1.0, nslice=ns, name="quad1"),
    elements.Drift(ds=0.5, nslice=ns, name="drift2"),
    elements.Quad(ds=1.0, k=-1.0, nslice=ns, name="quad2"),
    elements.Drift(ds=0.25, nslice=ns, name="drift3"),
]

sim.lattice.extend(lattice_configuration)
sim.periods = 2

# Simulate
sim.add_particles(bunch_charge_C, distr, npart)
sim.track_particles()

# Clean Shutdown
sim.finalize()
