#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_

#include <Kokkos_MathematicalConstants.hpp>
#include <string_view>

#ifndef PDG_VERSION
#define PDG_VERSION 2012
#endif

namespace pconstants {
#if PDG_VERSION == 2012
    constexpr std::string_view pdg_year("2012");
    //  J. Beringer et al. (Particle Data Group), PR D86, 010001 (2012) and 2013
    //  partial
    //     update for the 2014 edition (URL: http://pdg.lbl.gov)
    constexpr double mp = 0.938272046;          // Mass of proton [GeV/c^2]
    constexpr double proton_mass = 0.938272046; // Mass of proton [GeV/c^2]
    constexpr double me = 0.510998928e-3;       // Mass of electron [GeV/c^2]
    constexpr double electron_mass =
        0.510998928e-3;                        // Mass of electron [GeV/c^2]
    constexpr double mmu = 0.1056583715;       // Mass of muon [GeV/c^2]
    constexpr double muon_mass = 0.1056583715; // Mass of muon [GeV/c^2]
    constexpr double e = 1.602176565e-19;      // Charge of proton [C]
#else
#if PDG_VERSION == 2010
    //  K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010)
    constexpr std::string_view pdg_year("2010");
    constexpr double mp = 0.938272013;          // Mass of proton [GeV/c^2]
    constexpr double proton_mass = 0.938272013; // Mass of proton [GeV/c^2]
    constexpr double me = 0.510998910e-3;       // Mass of electron [GeV/c^2]
    constexpr double electron_mass =
        0.510998910e-3;                       // Mass of electron [GeV/c^2]
    constexpr double mmu = 0.105658367;       // Mass of muon [GeV/c^2]
    constexpr double muon_mass = 0.105658367; // Mass of muon [GeV/c^2]
    constexpr double e = 1.602176487e-19;     // Charge of proton [C]
#else
#if PDG_VERSION == 2008 // PDG_VERSION2008
    //  C. Amsler et al. (Particle Data Group), Physics Letters B667, 1 (2008)
    constexpr std::string_view pdg_year("2008");
    constexpr double mp = 0.938272013;          // Mass of proton [GeV/c^2]
    constexpr double proton_mass = 0.938272013; // Mass of proton [GeV/c^2]
    constexpr double me = 0.510998910e-3;       // Mass of electron [GeV/c^2]
    constexpr double electron_mass =
        0.510998910e-3;                       // Mass of electron [GeV/c^2]
    constexpr double mmu = 0.105658367;       // Mass of muon [GeV/c^2]
    constexpr double muon_mass = 0.105658367; // Mass of muon [GeV/c^2]
    constexpr double e = 1.602176487e-19;     // Charge of proton [C]
#if PDG_VERSION == 2006 // PDG_VERSION2006
                        //  S. Eidelman et al. (Particle Data Group), Phys.
                        //  Lett. B 592, 1 (2004)
    constexpr std::string_view pdg_year("2006");
    constexpr double mp = 0.938272029;          // Mass of proton [GeV/c^2]
    constexpr double proton_mass = 0.938272029; // Mass of proton [GeV/c^2]
    constexpr double me = 0.510998918e-3;       // Mass of electron [GeV/c^2]
    constexpr double electron_mass =
        0.510998918e-3;                       // Mass of electron [GeV/c^2]
    constexpr double mmu = 0.105658369;       // Mass of muon [GeV/c^2]
    constexpr double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]
    constexpr double e = 1.60217653e-19;      // Charge of proton [C]
#else
#if PDG_VERSION == 2004
                                          //  S. Eidelman et al. (Particle Data
                                          //  Group), Phys.
                                          //  Lett. B 592, 1 (2004)
    constexpr std::string_view pdg_year("2004");
    constexpr double mp = 0.938272029;          // Mass of proton [GeV/c^2]
    constexpr double proton_mass = 0.938272029; // Mass of proton [GeV/c^2]
    constexpr double me = 0.510998918e-3;       // Mass of electron [GeV/c^2]
    constexpr double electron_mass =
        0.510998918e-3;                       // Mass of electron [GeV/c^2]
    constexpr double mmu = 0.105658369;       // Mass of muon [GeV/c^2]
    constexpr double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]
    constexpr double e = 1.60217653e-19;      // Charge of proton [C]
#else
#if PDG_VERSION == -1 // legacy
    constexpr std::string_view pdg_year("legacy");
    constexpr double mp = 0.93827203;          // Mass of proton [GeV/c^2]
    constexpr double proton_mass = 0.93827203; // Mass of proton [GeV/c^2]
    constexpr double me = 0.51099892e-3;       // Mass of electron [GeV/c^2]
    constexpr double electron_mass =
        0.51099892e-3;                        // Mass of electron [GeV/c^2]
    constexpr double mmu = 0.105658369;       // Mass of muon [GeV/c^2]
    constexpr double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]
    constexpr double e = 1.6021892e-19;       // Charge of proton [C]
#else
#error "No selection for physicsl_constants PDG_VERSION set"
#endif                                // legacy
#endif                                // PDG_VERSION 2004
#endif                                // PDG_VERSION 2006
#endif                                // PDG_VERSION 2008
#endif                                // PDG_VERSION 2010
#endif                                // PDG_VERSION 2012
    constexpr double c = 299792458.0; // Speed of light [m/s]

    constexpr double mu0 = 4 * Kokkos::numbers::pi_v<double> *
                           1.0e-7; // Permittivity of free space [H/m]
    constexpr double epsilon0 =
        1.0 / (c * c * mu0); // Permeability of free space [F/m]

    // Classical radius of a particle = e^2/(4 pi epsilon0 m c^2)
    constexpr double re = e / (4 * Kokkos::numbers::pi_v<double> * epsilon0 *
                               me * 1.0e9); // Classical
                                            // radius of electron [m]
    constexpr double rp = e / (4 * Kokkos::numbers::pi_v<double> * epsilon0 *
                               mp * 1.0e9); // Classical
                                            // radius of proton [m]
    constexpr double rmu = e / (4 * Kokkos::numbers::pi_v<double> * epsilon0 *
                                mmu * 1.0e9); // Classical
                                              // radius of muon [m]

    constexpr int proton_charge = 1;      // Charge in units of e
    constexpr int antiproton_charge = -1; // Charge in units of e
    constexpr int electron_charge = -1;   // Charge in units of e
    constexpr int positron_charge = 1;    // Charge in units of e
    constexpr int muon_charge = -1;       // Charge in units of e
    constexpr int antimuon_charge = 1;    // Charge in units of e

    // mass [kg] * c**2 [m^2/s^2] = E [J]
    // E[j]/q[C} = E [eV]
    // conversion of mass from SI units to GeV used
    constexpr double kg_to_GeV = 1.0e-9 * c * c / e;

}

#endif /* PHYSICAL_CONSTANTS_H_ */
