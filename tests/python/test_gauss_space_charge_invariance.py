"""
Test that analytic Gaussian space-charge models center their field
at the bunch centroid (and not on the reference axis).

The test compares the space charge kicks computed when the bunch
is initialized with zero centroid to the space charge kicks
computed when there is a global centroid offset.

This test was contributed by N. Banerjee.

"""
from multiprocessing import Process, Queue # Child processes
import numpy as np # Arrays
from scipy.constants import c, e, m_p # Useful constants
from impactx import ImpactX, elements # ImpactX interface

Mass = m_p*c*c*1e-6/e # MeV, proton
KineticEnergy = 2.5 # MeV
SigmaX = 2.0e-3 # rms sizes of the round beam (m)
SigmaT = 0.1 # rms length (m)
Charge = 1e-10 # Beam charge (C)i
Displacement = 1e-3 # Beam dispacement from reference orbit (m)
NCells = 64 # PIC grid size for comparison
NParticles = 10000 # Number of macro-particles 
KickLength = 0.1 # m, the drift of test
Seed = 20260915 # Random number seed for beam
Models = ["Gauss2p5D", "Gauss3D", "2p5D"] # SC Models to test


# Utilities
def beam(npart:int, displacement:float):
    """
    A round Gaussian beam with displacement from reference axis.

    Parameters
    ----------
    npart: int
        Particles.

    displacement: float
        Offset of the whole beam in x (m).

    Returns
    -------
    x, y, t, px, py, pt: np.ndarray
    """
    rng = np.random.default_rng(Seed)
    x = rng.normal(0.0, SigmaX, npart)
    y = rng.normal(0.0, SigmaX, npart)
    t = rng.normal(0.0, SigmaT, npart)
    px = np.zeros(npart)
    py = np.zeros(npart)
    pt = np.zeros(npart)
    # Remove the sampling offset, so that the displacement is exact
    for a in (x, y, t, px, py):
        a -= a.mean()
    return x + displacement, y, t, px, py, pt


def coordinates(sim):
    """
    Coords x, y, px, py in the order of particle ids.

    Parameters
    ----------
    sim: impactx.ImpactX

    Returns
    -------
    x, y, px, py: np.ndarray
    """
    df = sim.beam.to_df()
    order = np.argsort(df["idcpu"].to_numpy())
    return tuple(df[name].to_numpy()[order] for name in
                 ("position_x", "position_y", "momentum_x", "momentum_y"))


def worker(q: Queue, model:str, displacement:float):
    """
    Run the test case and return results.

    Parameters
    ----------
    q: Queue
          Queue object to transfer field data back into the host process

    model: str
        The space-charge model.

    displacement: float
        Offset of the beam in x (m).

    Returns
    -------
    None
    """
    # Setup impactx run
    sim = ImpactX()
    sim.verbose = 0
    sim.tiny_profiler = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.particle_shape = 2
    assert model in Models, f"Models must be one of {Models}"
    sim.space_charge = model
    if model == "2p5D":
        sim.n_cell = [NCells, NCells, NCells]
        sim.poisson_solver = "fft"
    sim.init_grids()
    # Setup beam
    sim.beam.ref.set_charge_qe(1.0).set_mass_MeV(Mass)\
        .set_kin_energy_MeV(KineticEnergy)
    x, y, t, px, py, pt = beam(NParticles, displacement)
    weights = np.full(NParticles, Charge/e/NParticles)
    sim.beam.add_n_particles(x, y, t, px, py, pt, sim.beam.ref.qm_ratio_SI,
                             w=weights)
    x0, y0, px0, py0 = coordinates(sim) # Properly ordered coords
    # Track beam through drift
    sim.lattice.append(elements.Drift(ds=KickLength, nslice=1))
    sim.track_particles()
    x1, y1, px1, py1 = coordinates(sim)
    # Send back results
    q.put({"x0": x0, "y0": y0, "dpx": px1 - px0, "dpy": py1 - py0})


def _checkmodel(model:str):
    """
    Run a test case with the requested model and print results.

    model: str
        The space-charge model.

    Returns
    -------
    None
    """
    # First run 2 simulations
    sim_data_queue = Queue() # Queue to shuttle data
    data = [] # Empty list to collect data
    for displacement in [0.0, Displacement]: 
        # Simulation process
        sim_process = Process(target=worker,
                              args=(sim_data_queue, model, displacement))
        sim_process.start() # Start the simulation process
        data.append(sim_data_queue.get()) # Import the data through the queue
        sim_process.join() # End the child process
    centered = data[0]
    moved = data[1]
    # Now analyze and print results
    rms = np.sqrt(np.mean(centered["dpx"]**2 + centered["dpy"]**2))
    mean_kick = moved["dpx"].mean()
    difference = np.sqrt(np.mean((moved["dpx"] - centered["dpx"])**2 +
                                 (moved["dpy"] - centered["dpy"])**2))
    print("%-10s %14.4e %14.4f %16.4f" % (model, rms, mean_kick/rms,
                                          difference/rms))
    return difference/rms


def test_models():
    """
    Run a test case with the requested model and print results.
    """

    print(f"Test: one space-charge kick in a {KickLength:.2f} m drift, {NParticles} particles,"
          f" beam of {Charge:.2e} C displaced by {Displacement/SigmaX:.2f} sigma_x")
    print("%-10s %14s %14s %16s" % ("model", "rms kick", "mean kick/rms",
                                    "invariance/rms"))
    atol = 1.0e-6
    for model in Models:
        difference_kick_rms = _checkmodel(model)
        np.testing.assert_allclose([difference_kick_rms], 0.0, atol=atol)
