from impactx import ImpactX, distribution, elements, twiss

sim = ImpactX()

sim.space_charge = True  # deprecated name of "3D"
sim.particle_shape = 2

sim.slice_step_diagnostics = True

sim.init_grids()

# Initialize particle beam
kin_energy_MeV = 10.0
bunch_charge_C = 1e-12
npart = 1000

# Reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

distr = distribution.Waterbag(
    **twiss(
        beta_x=1.0,  # see (1)
        beta_y=2.0,
        beta_t=3.0,
        emitt_x=1e-06,
        emitt_y=2e-06,
        emitt_t=1e-06,
        alpha_x=-0.5,
        alpha_y=0.25,
        alpha_t=0.0,
    )
)
sim.init_envelope(ref, distr, bunch_charge_C)

cell = [
    elements.Drift(name="d", ds=0.2),
    elements.Quad(name="q", ds=0.1, k=2.0),
]
sim.lattice.append(elements.Marker("start"))  # name by position
sim.lattice.append(elements.DipEdge(psi=0.0, rc=1.0, g=0.0, modify_ref_part=False))
sim.lattice.append(elements.Drift(name="pre", ds=0.4))
sim.lattice.extend(cell * 3)

# Simulate
sim.track_envelope()

# Clean Shutdown
sim.finalize()
