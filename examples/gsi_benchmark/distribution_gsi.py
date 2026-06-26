from impactx import distribution, twiss

# beam parameters
kin_energy_MeV = 11.4  # reference kinetic energy
bunch_charge_C = 1.0e-9  # used with space charge
beam_current_A = 1.5  # used with 2D space charge
charge_qe = -1.0  # particle charge

# particle distribution
#
# The matched Twiss functions below are obtained from running gsi.madx using MAD-X.
# These values agree with the values reported here:
#    https://web-docs.gsi.de/~giuliano/research_activity/trapping_benchmarking/twiss-sis18,
# except for the horizontal dispersion and its derivative, which are reported as
# (Dx,Dx') = (1.946,-0.325).  Those values differ from the values obtained in MAD-X
# due to the convention of using pt vs delta as the canonical momentum.  Our dispersion
# convention is consistent with that used in MAD-X.
#
distr = distribution.Gaussian(
    **twiss(
        beta_x=12.79711091,
        beta_y=13.4860516,
        beta_t=6.7632723266784568737e6,
        emitt_x=1.953582871e-6,
        emitt_y=1.8537742844e-6,
        emitt_t=0.0100875,
        alpha_x=1.29174698,
        alpha_y=0.4270522414,
        alpha_t=0.0,
        dispersion_x=12.59966194,
        dispersion_px=-2.106017322,
    )
)
