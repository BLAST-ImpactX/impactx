import numpy as np
import amrex.space3d as amr

from impactx import Config, ImpactX, elements


KIN_ENERGY_MEV = 400.0
MASS_MEV = 938.27208816
QM_EEV = 1.0 / (MASS_MEV * 1.0e6)


def make_vec(values):
    vec_cls = amr.PODVector_real_arena if Config.have_gpu else amr.PODVector_real_std
    out = vec_cls()
    for vv in values:
        out.push_back(float(vv))
    return out


def run_one(element):
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.verbose = 0
    sim.init_grids()

    pc = sim.particle_container()
    ref = pc.ref_particle()
    ref.set_charge_qe(1.0).set_mass_MeV(MASS_MEV).set_kin_energy_MeV(KIN_ENERGY_MEV)

    # One particle exactly on the reference trajectory
    pc.add_n_particles(
        make_vec([0.0]),  # x
        make_vec([0.0]),  # y
        make_vec([0.0]),  # t
        make_vec([0.0]),  # px
        make_vec([0.0]),  # py
        make_vec([0.0]),  # pt
        QM_EEV,
        bunch_charge=0.0,
    )

    sim.lattice.append(element)
    sim.track_particles()

    df = pc.to_df(local=True)
    x = df["position_x"]
    px = df["momentum_x"]
    y = df["position_y"]
    py = df["momentum_y"]
    t = df["position_t"]
    pt = df["momentum_t"]
    
    print(x, px, y, py, t, pt)

    atol = 1.0e-13
    assert np.allclose(  
        [x, px, y, py, t, pt],
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        atol=atol,
    )

    sim.finalize()

print("Quad")
run_one(elements.Quad(name="q", ds=0.5, k=1.2))

print("\nExactQuad")
run_one(elements.ExactQuad(name="eq",ds=0.5,k=1.2))

print("\nExactMultipole")
run_one(elements.ExactMultipole(name="em",ds=0.5,k_normal=[0,1.2],k_skew=[0.0,0.0]))
