#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements, twiss

verbose = False

sim = ImpactX()

# set numerical parameters for space charge
sim.max_level = 0
# sim.n_cell = [64, 64, 64]  #use this for high-resolution runs
sim.n_cell = [32, 32, 32]
sim.particle_shape = 2  # B-spline order
sim.space_charge = "3D"
sim.poisson_solver = "fft"
sim.dynamic_size = True
sim.prob_relative = [1.1, 1.1]

# set parameters for IO control
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# npart = 2  # number of macro particles
# npart = 1000000  #use this for high-resolution runs
npart = 10000
ns = 25  # default number of slices per element

# beam reference parameters
kin_energy_MeV = 2.1226695  # reference energy
bunch_charge_C = 2.9824923076852509e-11  # used with space charge
rest_mass_MeV = 939.294308  # particle rest energy
charge_qe = -1.0  # particle charge

# reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(charge_qe).set_mass_MeV(rest_mass_MeV).set_kin_energy_MeV(
    kin_energy_MeV
)

# particle bunch
distr = distribution.Gaussian(
    **twiss(
        beta_x=1.7285172802864157,
        beta_y=1.8679530886750253,
        beta_t=1801.398618485624,
        emitt_x=3.0537240180190796e-06,
        emitt_y=2.9611736825829232e-06,
        emitt_t=5.4974953691649138e-06,
        alpha_x=0.13233244124654187,
        alpha_y=0.039004242136615824,
        alpha_t=-1.2272478360491141,
    )
)

# initialize particle bunch
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# field coefficients
QWR_sin_coefs = [
    0.0000000000000000,
    -0.44552128900111937,
    -0.46341657636313111,
    -0.20231040922531793,
    5.1588444808049428e-003,
    5.8844715958123073e-002,
    3.4824638238541664e-002,
    8.9094492663513401e-003,
    1.1498524643627059e-003,
    1.7073943163752061e-003,
    2.0526237634658114e-003,
    7.4066308552227407e-004,
    -3.1559085909346231e-004,
    -6.0472607331999882e-004,
    -2.9523109071238765e-004,
    -9.2848523466919987e-005,
    3.6596205157362737e-005,
    -1.0615530581730304e-005,
    1.5947093117194282e-005,
    -2.5026267107225051e-005,
    1.4176336245546890e-005,
    -1.4237084172757583e-005,
    2.4230112875582230e-005,
    -8.7542451259873567e-006,
    1.5747147629040920e-005,
    -1.7289474688373437e-005,
    1.1053141257897324e-005,
    -1.3365651893844886e-005,
    1.7088725904341150e-005,
    -6.9045871183576885e-006,
]

QWR_cos_coefs = [
    1.4817546261092218e-006,
    2.8017487205400926e-003,
    5.8369359230740769e-003,
    3.8183683433996862e-003,
    -1.2845097142791251e-004,
    -1.8561682988446904e-003,
    -1.3146193274420626e-003,
    -3.9608853919712028e-004,
    -5.6483925331224372e-005,
    -1.0009076752079404e-004,
    -1.2795452515182437e-004,
    -5.4591022951439028e-005,
    2.5375542690636577e-005,
    4.6438447566049534e-005,
    2.7634072263776899e-005,
    5.6507816244277809e-006,
    -2.2331226970095930e-006,
    -2.0902475061645998e-006,
    -4.0028352235599307e-007,
    -1.7573212882093792e-007,
    -2.7559960050216326e-007,
    -1.2290845811946394e-006,
    -1.8015143666005873e-006,
    -1.8601874063051138e-006,
    -8.7539741192022014e-007,
    -4.7471486960826148e-007,
    -3.3075456642994716e-007,
    -8.8472914991366558e-007,
    -1.4898207046087908e-006,
    -1.8563793804105000e-006,
]

SOL_sin_coefs = [
    0.0000000000000000,
    8.4583007312587222e-009,
    1.5625118966917928e-008,
    1.2544039813496867e-008,
    1.2951204561920115e-008,
    -1.3824319466948509e-009,
    -2.3654450397774321e-008,
    -2.3981556296348572e-008,
    -4.3793988876300594e-008,
    5.3089692633387364e-009,
    6.5186324602062307e-009,
    3.9130100049078465e-008,
    -6.6213967512065985e-010,
    2.9653165256604552e-008,
    -2.7894202503375709e-008,
    -5.8284914868870358e-008,
    -9.8081727628596127e-008,
    -4.8879883252084255e-008,
    6.2507581315845995e-009,
    -4.0900704334490001e-007,
    -7.8168397919797794e-007,
    -6.5719126518537507e-007,
    -2.5097324396483600e-007,
    -3.0195224098861217e-010,
    3.9559326480385903e-008,
    -8.2029146142303944e-008,
    6.2776962295174599e-008,
    1.1066968561013191e-007,
    7.6001015258952975e-008,
    -3.2838215702213347e-008,
]

SOL_cos_coefs = [
    0.86107056619615296,
    0.55971853194956966,
    9.5165555696743537e-002,
    -7.1096582004436276e-002,
    -3.3379302016024598e-002,
    1.1432386673232042e-002,
    1.1410602237751005e-002,
    -6.1644831212256945e-004,
    -3.6726479972148118e-003,
    -5.8336394378343155e-004,
    9.6514088030205669e-004,
    4.8042818016283753e-004,
    -2.5709054504289078e-004,
    -1.8039534282018529e-004,
    4.4409228704099524e-006,
    1.0008305620550287e-004,
    -9.3173909280520884e-006,
    -8.2580999710389624e-006,
    -2.6499259243726453e-005,
    2.1681775087243477e-005,
    -9.1528851643742082e-006,
    1.2202909413342594e-005,
    -1.4992933452893237e-005,
    1.1952168037430294e-005,
    -1.1120027616860068e-005,
    1.1698954757780738e-005,
    -1.1543046335580481e-005,
    1.0876363164911279e-005,
    -1.0739105769905916e-005,
    1.0832271198166904e-005,
]

HWR_sin_coefs = [
    0.0000000000000000,
    0.62446804950371848,
    0.39543633732295935,
    -8.0003693987048971e-002,
    -0.14321597590146887,
    -2.2354907085085453e-002,
    -6.5045696342092103e-003,
    -2.7625438898537705e-002,
    -8.3495199535108845e-003,
    1.0349654439745908e-002,
    5.8109366862368443e-003,
    6.7645812805355431e-005,
    1.2488286275393725e-003,
    1.6152921567933968e-003,
    -2.3731332273832940e-004,
    -7.3735142731588073e-004,
    -1.7385804784192361e-004,
    8.6718934908533787e-007,
    -1.7379313943036016e-004,
    -5.3503333766242609e-005,
    5.5700992587132386e-005,
    5.0154828533339330e-005,
    -1.1599281414194189e-005,
    3.0462696269272361e-005,
    9.1556008232301467e-006,
    4.6620667631519768e-006,
    -1.8897102253325956e-005,
    5.8585163714951205e-006,
    -1.2633708997350154e-005,
    7.5902876048951134e-006,
]

HWR_cos_coefs = [
    -1.6429636351400490e-006,
    -3.9298641775114274e-003,
    -4.9810064959918410e-003,
    1.5128643575186007e-003,
    3.6056898690726424e-003,
    7.0558052307109764e-004,
    2.4417054683947836e-004,
    1.2198738355022332e-003,
    4.1923771228294382e-004,
    -5.8539325032239131e-004,
    -3.6796369613131819e-004,
    -3.0374723325697628e-006,
    -9.6171251282330811e-005,
    -1.3084677545929443e-004,
    1.9291398171614876e-005,
    7.1480486311979036e-005,
    1.5985530655879959e-005,
    1.6333992165253053e-006,
    1.8109445470568947e-005,
    8.0417488665697479e-006,
    -8.6861802619207784e-006,
    -5.0126344126859984e-006,
    -5.8424457905284499e-008,
    -2.7977139079021079e-006,
    -3.0325780533083796e-006,
    9.1565069457351378e-007,
    1.4630386344968022e-006,
    6.7824413083235335e-007,
    5.9435230451675147e-007,
    2.7696933433025983e-007,
]

# lattice segments

# Fragment of MEBT lattice
D1 = elements.Drift(name="D197", ds=0.17414, nslice=ns)
RF1 = elements.RFCavity(
    name="RF4",
    ds=0.24,
    escale=-1.09922948665414467808720822e-3,
    freq=162500000.0,
    phase=43.82865362,
    cos_coefficients=QWR_cos_coefs,
    sin_coefficients=QWR_sin_coefs,
    mapsteps=100,
    nslice=ns,
)
D2 = elements.Drift(name="D204", ds=0.06586, nslice=ns)
Q1 = elements.ChrQuad(name="Q29", ds=0.05, k=11.4, unit=1, nslice=ns)
D3 = elements.Drift(name="D205", ds=0.07, nslice=ns)
Q2 = elements.ChrQuad(name="Q30", ds=0.1, k=-10.47, unit=1, nslice=ns)
D4 = elements.Drift(name="D208", ds=0.070, nslice=ns)
Q3 = elements.ChrQuad(name="Q31", ds=0.05, k=11.4, unit=1, nslice=ns)
D5 = elements.Drift(name="D208", ds=1.0106078, nslice=ns)

mebt_fragment = [D1, RF1, D2, Q1, D3, Q2, D4, Q3, D5]

# Fragment of HWR lattice
D1b = elements.Drift(name="D1b", ds=0.177422, nslice=ns)
SOL1b = elements.SoftSolenoid(
    name="SOL1b",
    ds=0.3,
    bscale=2.3192,
    cos_coefficients=SOL_cos_coefs,
    sin_coefficients=SOL_sin_coefs,
    unit=1,
    mapsteps=100,
    nslice=ns,
)
D2b = elements.Drift(name="D2b", ds=0.0679, nslice=ns)
RF1b = elements.RFCavity(
    name="RF1b",
    ds=0.25,
    escale=-0.007252998279640,
    freq=162500000.0,
    phase=-34.6544361,
    cos_coefficients=HWR_cos_coefs,
    sin_coefficients=HWR_sin_coefs,
    mapsteps=100,
    nslice=ns,
)
D3b = elements.Drift(name="D3b", ds=0.0679, nslice=ns)
SOL2b = elements.SoftSolenoid(
    name="SOL2b",
    ds=0.3,
    bscale=-2.3969,
    cos_coefficients=SOL_cos_coefs,
    sin_coefficients=SOL_sin_coefs,
    unit=1,
    mapsteps=100,
    nslice=ns,
)
D4b = elements.Drift(name="D4b", ds=0.0679, nslice=ns)

hwr_fragment = [D1b, SOL1b, D2b, RF1b, D3b, SOL2b, D4b]

# assign lattice
sim.lattice.append(monitor)
sim.lattice.extend(mebt_fragment)
sim.lattice.extend(hwr_fragment)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
