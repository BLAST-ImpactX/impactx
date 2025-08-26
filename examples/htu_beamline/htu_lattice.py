"""
Define the lattice for the HTU accelerator.

References:
- https://doi.org/10.1103/vh62-gz1p
- https://doi.org/10.1117/12.3056776
"""

from scipy.constants import c, e, m_e

from impactx import elements

impactx_available = True


def get_lattice(
    code,
    vs_current_x=[0] * 8,
    vs_current_y=[0] * 8,
    emq_currents=[0.683822195172598, -0.882403221217054, 1.085116293768873],
    chicane_r56=200.0,
    from_element=None,
    to_element=None,
):
    """
    Get the lattice for the HTU accelerator.

    Parameters:
        code (str):
            "impactx", ...
        vs_current_x (list of 8 elements):
            The current in the VISA horizontal steering magnets, in amps.
        vs_current_y (list of 8 elements):
            The current in the VISA vertical steering magnets, in amps.
        emq_currents (list of 3 elements):
            The current in EMQ1H, EMQ2V and EMQ3H, in amps.
        chicane_r56 (float):
            The r56 setting of the chicane, in microns. This is the parameter that is entered in the
            control system, in BELLA operations. It does not necessarily correspond to the actual R56
            of the chicane, esp. if the mean energy of the beam differs from 100 MeV.
        from_element (str):
            The name of the element to start from. If None, start from the beginning.
        to_element (str):
            The name of the element to end at. If None, end at the last element.

    Returns:
        list: The lattice for the HTU accelerator.
    """
    TCPhosphor = screen("TCPhosphor", code)
    ChicaneSlit = screen("ChicaneSlit", code)
    DCPhosphor = screen("DCPhosphor", code)
    Phosphor1 = screen("Phosphor1", code)
    Aline1 = screen(
        "UC_ALineEbeam1",
        code,
    )
    Aline2 = screen(
        "UC_ALineEBeam2",
        code,
    )
    Aline3 = screen(
        "UC_ALineEBeam3",
        code,
    )
    VisaEBeam1 = screen(
        "UC_VisaEBeam1",
        code,
    )
    VisaEBeam2 = screen(
        "UC_VisaEBeam2",
        code,
    )
    VisaEBeam3 = screen(
        "UC_VisaEBeam3",
        code,
    )
    VisaEBeam4 = screen(
        "UC_VisaEBeam4",
        code,
    )
    VisaEBeam5 = screen(
        "UC_VisaEBeam5",
        code,
    )
    VisaEBeam6 = screen(
        "UC_VisaEBeam6",
        code,
    )
    VisaEBeam7 = screen(
        "UC_VisaEBeam7",
        code,
    )
    VisaEBeam8 = screen(
        "UC_VisaEBeam8",
        code,
    )

    # Distance from the plasma source to the first PMQ, when Jetz=6.8
    # Jetz moves the jet Upstream for positive changes.  Therefore, we can
    # get ~6.8mm closer to the PMQ1, and about ~15mm further away.
    SrcToPMQ1 = drift("SrcToPMQ1", 0.052, code)
    PMQ1V = quadrupole("PMQ1V", L=0.02903, bore_radius=0.006, B=-1.242, code=code)
    PMQ2H = quadrupole("PMQ2H", L=0.02890, bore_radius=0.006, B=1.242, code=code)
    PMQ3V = quadrupole("PMQ3V", L=0.016321, bore_radius=0.006, B=-1.107, code=code)

    L1 = drift("L1", 0.029035, code)
    L2 = drift("L2", 0.0473895, code)

    PMQTrip = [PMQ1V, L1, PMQ2H, L2, PMQ3V]

    PMQTripToTCPhos = drift("PMQTripToTCPhos", 0.2158, code)
    TCPhosToChicane = drift("TCPhosToChicane", 0.42, code)

    # Integrated field (converted from G.cm to T.m) and max current from HTU_Kickers_SteeringMagnets.pdf
    # https://drive.google.com/drive/u/0/folders/1r9InjfW7-92OUpZo4xADUm6rVM5u13yt
    S1 = kicker(
        "S1", 0.0, 0.0, max_integrated_field=1970e-6, max_current=5.0, code=code
    )
    S2 = kicker(
        "S2", 0.0, 0.0, max_integrated_field=1970e-6, max_current=5.0, code=code
    )
    S3 = kicker(
        "S3", 0.0, 0.0, max_integrated_field=1970e-6, max_current=5.0, code=code
    )
    S4 = kicker(
        "S4", 0.0, 0.0, max_integrated_field=1970e-6, max_current=5.0, code=code
    )

    # Calculation of the bending field and bend angle are now found in the dipole definitions below.
    BEND1 = dipole("BEND1", 0.175, r56=chicane_r56, bend=1, code=code)
    BEND2 = dipole("BEND2", 0.175, r56=chicane_r56, bend=2, code=code)
    BEND3 = dipole("BEND3", 0.175, r56=chicane_r56, bend=3, code=code)
    BEND4 = dipole("BEND4", 0.175, r56=chicane_r56, bend=4, code=code)

    L12 = drift("L12", 0.125, code)
    L23 = drift("L23", 0.15, code)

    Chicane = [BEND1, L12, BEND2, L23, ChicaneSlit, L23, BEND3, L12, BEND4]

    DriftToDCPhos = drift("DriftToDCPhos", 0.27, code)
    DriftToEMQTrip = drift("DriftToEMQTrip", 0.405, code)
    EMQ1H = quadrupole(
        "EMQ1H", L=0.1408, current=emq_currents[0], design="EMQD-113-394", code=code
    )
    EMQL1 = drift("EMQL1", 0.112735, code)
    EMQ2V = quadrupole(
        "EMQ2V", L=0.28141, current=emq_currents[1], design="EMQD-113-949", code=code
    )
    EMQL2 = drift("EMQL2", 0.112735, code)
    EMQ3H = quadrupole(
        "EMQ3H", L=0.1409, current=emq_currents[2], design="EMQD-113-394", code=code
    )

    EMQTriplet = [EMQ1H, EMQL1, EMQ2V, EMQL2, EMQ3H]

    DriftToPhos1 = drift("DriftToPhos1", 0.084325, code)
    DriftToSpec = drift("DriftToSpec", 0.28, code)
    MagSpec = dipole("MagSpec", 0.4826, angle=0.0, code=code)

    DriftToAline1 = drift("DriftToAline1", 0.4009, code)
    DriftToAline2 = drift("DriftToAline2", 0.3825, code)
    DriftToAline3 = drift("DriftToAline3", 0.4191, code)

    DriftToUndulator = drift("DriftToUndulator", 0.2945, code)
    # TODO: UndulatorAperture = aperture("UndulatorAperture", 0.005, 0.003, code)

    Aline = [
        DriftToAline1,
        Aline1,
        DriftToAline2,
        Aline2,
        DriftToAline3,
        Aline3,
        DriftToUndulator,
    ]

    # Integrated field (converted from G.cm to T.m) is from HTU_Kickers_SteeringMagnets.pdf
    # but scaled by 0.5 to account for lower peak field, per Sam's recommendation
    VS1 = kicker(
        "VS1",
        vs_current_x[0],
        vs_current_y[0],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS2 = kicker(
        "VS2",
        vs_current_x[1],
        vs_current_y[1],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS3 = kicker(
        "VS3",
        vs_current_x[2],
        vs_current_y[2],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS4 = kicker(
        "VS4",
        vs_current_x[3],
        vs_current_y[3],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS5 = kicker(
        "VS5",
        vs_current_x[4],
        vs_current_y[4],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS6 = kicker(
        "VS6",
        vs_current_x[5],
        vs_current_y[5],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS7 = kicker(
        "VS7",
        vs_current_x[6],
        vs_current_y[6],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )
    VS8 = kicker(
        "VS8",
        vs_current_x[7],
        vs_current_y[7],
        max_integrated_field=1970e-6 * 0.5,
        max_current=5.0,
        code=code,
    )

    VD = drift("VD", 0.018, code)

    VQ1 = quadrupole("VQ1", L=0.0504, bore_radius=0.004, B=0.132, code=code)
    VQ2 = quadrupole("VQ2", L=0.0504, bore_radius=0.004, B=-0.132, code=code)
    FODOCell1 = [VQ1, VQ1, VD, VQ2, VQ2, VD]
    UndulatorSegment1 = [
        *FODOCell1,
        VS1,
        VisaEBeam1,
        *FODOCell1,
        *FODOCell1,
        VS2,
        VisaEBeam2,
        *FODOCell1,
    ]

    VQ3 = quadrupole("VQ3", L=0.0504, bore_radius=0.004, B=0.132, code=code)
    VQ4 = quadrupole("VQ4", L=0.0504, bore_radius=0.004, B=-0.132, code=code)
    FODOCell2 = [VQ3, VQ3, VD, VQ4, VQ4, VD]
    UndulatorSegment2 = [
        *FODOCell2,
        VS3,
        VisaEBeam3,
        *FODOCell2,
        *FODOCell2,
        VS4,
        VisaEBeam4,
        *FODOCell2,
    ]

    VQ5 = quadrupole("VQ5", L=0.0504, bore_radius=0.004, B=0.132, code=code)
    VQ6 = quadrupole("VQ6", L=0.0504, bore_radius=0.004, B=-0.132, code=code)
    FODOCell3 = [VQ5, VQ5, VD, VQ6, VQ6, VD]
    UndulatorSegment3 = [
        *FODOCell3,
        VS5,
        VisaEBeam5,
        *FODOCell3,
        *FODOCell3,
        VS6,
        VisaEBeam6,
        *FODOCell3,
    ]

    VQ7 = quadrupole("VQ7", L=0.0504, bore_radius=0.004, B=0.132, code=code)
    VQ8 = quadrupole("VQ8", L=0.0504, bore_radius=0.004, B=-0.132, code=code)
    FODOCell4 = [VQ7, VQ7, VD, VQ8, VQ8, VD]
    UndulatorSegment4 = [
        *FODOCell4,
        VS7,
        VisaEBeam7,
        *FODOCell4,
        *FODOCell4,
        VS8,
        VisaEBeam8,
        *FODOCell4,
    ]

    B0 = [
        SrcToPMQ1,
        *PMQTrip,
        PMQTripToTCPhos,
        TCPhosphor,
        TCPhosToChicane,
        S1,
        *Chicane,
        S2,
        DriftToDCPhos,
        DCPhosphor,
        DriftToEMQTrip,
        *EMQTriplet,
        S3,
        DriftToPhos1,
        Phosphor1,
        DriftToSpec,
        MagSpec,
        S4,
        *Aline,
    ]
    Undulator = (
        UndulatorSegment1 + UndulatorSegment2 + UndulatorSegment3 + UndulatorSegment4
    )
    full_beamline = B0 + Undulator

    # Only return the elements between from_element and to_element
    if from_element is not None:
        element_indices = [
            i for i, elem in enumerate(full_beamline) if elem.name == from_element
        ]
        if len(element_indices) == 0:
            raise ValueError(f"Element {from_element} not found in beamline.")
        elif len(element_indices) > 1:
            raise ValueError(
                f"Multiple elements with name {from_element} found in beamline."
            )
        # Get the first index of the element
        from_index = element_indices[0]
    else:
        from_index = 0
    if to_element is not None:
        element_indices = [
            i for i, elem in enumerate(full_beamline) if elem.name == to_element
        ]
        if len(element_indices) == 0:
            raise ValueError(f"Element {to_element} not found in beamline.")
        elif len(element_indices) > 1:
            raise ValueError(
                f"Multiple elements with name {to_element} found in beamline."
            )
        # Get the first index of the element
        to_index = element_indices[0]
    else:
        to_index = len(full_beamline) - 1
    beamline = full_beamline[from_index : to_index + 1]

    return beamline


# Define a screen element
def screen(name, code):
    """
    Define a screen element.
    """
    return elements.BeamMonitor(name=name, backend="h5")


# Define a quadrupole element
def peakfield_to_Bgradient(bore_radius, B):
    """
    Convert the peak field in T (and corresponding bore radius) to the field gradient in T/m of a quadrupole.

    Parameters:
        bore_radius (float):
            The bore radius in m.
        B (float):
            The peak field in T.
    """
    Bgradient = B / bore_radius
    return Bgradient


def current_to_Bgradient(current, design):
    """
    Convert a current in A to the field gradient in T/m of a quadrupole.

    Parameters:
        current (float):
            The current in A.
        design (str):
            Indicates which EMQ design is being used
    """
    if design == "EMQD-113-394":
        # From "LBM6_03 Calibration" in EMQD-113-394 Testing Report-FINAL.pdf
        Bgradient = 2.9217 * current + 0.0965  # T/m
    elif design == "EMQD-113-949":
        # From "LBM6_03 Calibration" in EMQD-113-949 Testing Report-FINAL.pdf
        Bgradient = 2.9318 * current + 0.0077  # T/m
    else:
        raise ValueError(f"Unsupported design: {design}")
    return Bgradient


def get_rigidity(reference_energy_eV):
    """
    Return the magnetic rigidity associated with reference_energy_eV in
    units of T-m, to be used for field strength normalization.

    Parameters:
        Bfield (float):
            The reference energy in eV.
    """
    gamma = reference_energy_eV * e / (m_e * c**2)
    beta = (1 - 1 / gamma**2) ** 0.5
    rigidity = m_e * gamma * beta * c / e
    return rigidity


def quadrupole(
    name,
    L,
    k1=None,
    current=None,
    design=None,
    bore_radius=None,
    B=None,
    code=None,
    reference_energy_eV=100e6,
):
    """
    Define a quadrupole element.

    Need to provide k1, current, or bore_radius and B.
    """
    if k1 is None and current is None:
        Bgradient = -1.0 * peakfield_to_Bgradient(bore_radius, B)
        unit = 1
    elif k1 is None:
        Bgradient = -1.0 * current_to_Bgradient(current, design)
        unit = 1
    else:
        Bgradient = k1
        unit = 0
    return elements.ChrQuad(name=name, ds=L, k=Bgradient, unit=unit, nslice=1)


# Define a drift element
def drift(name, L, code):
    """
    Define a drift element.
    """
    return elements.ExactDrift(name=name, ds=L, nslice=1)


# Define a kicker element
def current_to_integrated_field(current, max_current, max_integrated_field):
    integrated_field = current * max_integrated_field / max_current
    return integrated_field


def kicker(
    name,
    current_h,
    current_v,
    max_current,
    max_integrated_field,
    code,
    reference_energy_eV=100e6,
):
    """
    Define a kicker element.
    """
    integrated_field_h = current_to_integrated_field(
        current_h, max_current, max_integrated_field
    )
    integrated_field_v = current_to_integrated_field(
        current_v, max_current, max_integrated_field
    )
    return elements.Kicker(
        name=name, xkick=integrated_field_h, ykick=integrated_field_v, unit="T-m"
    )


# Define a dipole element
def chicane_r56_to_field(r56, L, reference_energy_eV):
    """
    Return the magnetic field in T associated with the chicane setting of r56
    defined at the reference energy reference_energy_eV.

    Parameters:
        r56 (float):
            The chicane r56 at nominal energy in microns.
        L (float):
            The bend length in m.
        reference_energy_eV (float):
            The reference energy in eV.
    """

    # Polynomial fit to find the dipole angle as a function of chicane r56 (no longer used)
    # angle_rad = 0.00743627 + 6.32812981e-05*chicane_r56 -4.83395592e-08*chicane_r56**2 + 1.91504783e-11*chicane_r56**3
    # Updated r56 formula valid near r56=0
    angle_rad = 0.001438389904456 * r56**0.5
    # TODO: This angle explicitly assumes that the reference energy of the beam is 100 MeV! We should make
    # sure that this is the case. In particular, for beams with fluctuations in energy, we should keep the
    # energy of the reference particles constant.

    B = -1.0 * get_rigidity(reference_energy_eV) * angle_rad / L
    return B


def Bfield_to_angle(Bfield, L, reference_energy_eV):
    """
    Return the bend angle in radians associated with a magnetic field Bfield
    in units of T, for a bend of length L and specified energy.

    Parameters:
        Bfield (float):
            The value of the magnetic field in T.
        L (float):
            The bend length in m.
        reference_energy_eV (float):
            The reference energy in eV.
    """
    angle = -1.0 * Bfield * L / get_rigidity(reference_energy_eV)
    return angle


def dipole(
    name, L, angle=None, r56=None, bend=None, code=None, reference_energy_eV=100e6
):
    """
    Define a dipole element.
    """
    if angle is None:
        Bfield = chicane_r56_to_field(r56, L, reference_energy_eV=100e6)
        angle = Bfield_to_angle(Bfield, L, reference_energy_eV)
    else:
        Bfield = 0.0
    if bend == 1 or bend == 4:
        Bfield = -Bfield
        angle = -angle
    angle_deg = angle * 180.0 / (3.1415926535898)
    return elements.ExactSbend(name=name, ds=L, phi=angle_deg, B=Bfield, nslice=1)
