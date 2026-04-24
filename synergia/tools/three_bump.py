#!/usr/bin/env python

import h5py
import numpy as np
from scipy.optimize import least_squares

import synergia


# class that calculates corrector settings to create a 3 kick local orbit bump
class Three_bump:
    ##################################################

    #  lattice you wnat
    #  to create a bump in.
    # start_name is the name of the starting element for the bump
    #    (use a unique marker if necessary)
    # end_name is the name of the ending element for the bump
    #    (use a unique marker if necessary)
    # hcorr_names is a length=3 sequence of names of horizontal corrector
    #    elements that have a "kick=" attribute
    # vcorr_names is a length=3 sequence of names of vertical corrector
    #    elements that have a "kick=" attribute
    # target name is the name of the element at which the bump offset will
    #    be specified
    # verbose = False/True on whether the module is chatty

    def __init__(
        self,
        lattice,
        start_name,
        end_name,
        hcorr_names,
        vcorr_names,
        target_name,
        verbose=0,
    ):

        self.lattice = lattice
        # this is a reference to the actual working lattice
        # I keep the elements separately so I can adjust them at the end when I know
        # what the settings are
        self.lattice_elements = self.lattice.get_elements()
        self.lattice_elements_idx = range(len(self.lattice_elements))
        self.start_name = start_name
        self.end_name = end_name
        self.hcorr_names = hcorr_names
        self.vcorr_names = vcorr_names
        self.target_name = target_name
        self.verbose = verbose

        self._construct()

    def _construct(self):

        # extract the bump region out of the lattice
        self.bump_lattice = synergia.lattice.Lattice("bump")
        self.bump_lattice.set_reference_particle(self.lattice.get_reference_particle())

        if self.verbose:
            print("length of lattice_elements: ", len(self.lattice_elements))
        elem_names = [e.get_name() for e in self.lattice_elements]
        if self.verbose:
            print("elem_names: ", elem_names)

        # find the start and end of the region by element name
        try:
            start_idx = elem_names.index(self.start_name)
            if self.verbose:
                print("start_idx: ", start_idx)
        except:
            raise RuntimeError("Three_bump: start_name: %s not found" % self.start_name)
        try:
            end_idx = elem_names.index(self.end_name)
            if self.verbose:
                print("end_idx: ", end_idx)
        except:
            raise RuntimeError("Three_bump: end_name: %s not found" % self.end_name)

        # start is before end so I can just peel elements in that range
        if start_idx < end_idx:
            for elem in self.lattice_elements[start_idx : end_idx + 1]:
                self.bump_lattice.append(elem)
            self.bump_idx = self.lattice_elements_idx[start_idx : end_idx + 1]
        else:
            for elem in (
                self.lattice_elements[start_idx:] + self.lattice_elements[: end_idx + 1]
            ):
                self.bump_lattice.append(elem)
            self.bump_idx = (
                self.lattice_elements_idx[start_idx:]
                + self.lattice_elements_idx[: end_idx + 1]
            )

        # self.bump_idx[] is the index pointing to the original elements in self.lattice_elements

        if self.verbose > 2:
            print("bump lattice:")
            print(self.bump_lattice)
            print("bump_idx: ", self.bump_idx)

        bump_elements = self.bump_lattice.get_elements()
        bump_enames = [e.get_name() for e in bump_elements]

        # get the corrector elements and keep track of their original positions in the lattice

        # hcorr_names or vcorr_names may be None, but not both of them
        if (self.hcorr_names is None) and (self.vcorr_names is None):
            raise RuntimeError("no corrector element names specified")

        # get horizontal corrector elements
        if self.hcorr_names is None:
            self.hcorr_idx = None
            self.hcorr_elements = None
        else:
            if len(self.hcorr_names) != 3:
                raise RuntimeError("length of hcorr_names needs to be 3")

            self.hcorr_elements = []
            self.hcorr_idx = []
            for i in range(3):
                try:
                    hc_idx = bump_enames.index(self.hcorr_names[i])
                except:
                    raise RuntimeError(
                        "Three_bump: hcorr_name: %s not found" % self.hcorr_names[i]
                    )

                if self.verbose:
                    print("hcorrector element ", i, " at index: ", hc_idx)
                self.hcorr_elements.append(bump_elements[hc_idx])
                self.hcorr_idx.append(self.bump_idx[hc_idx])
            for e in self.hcorr_elements:
                e.set_double_attribute("hkick", 0.0)

        # get vertical corrector elements
        if self.vcorr_names is None:
            self.vcorr_idx = None
            self.vcorr_elements = None
        else:
            if len(self.vcorr_names) != 3:
                raise RuntimeError("length of hcorr_names needs to be 3")

            self.vcorr_elements = []
            self.vcorr_idx = []
            for i in range(3):
                try:
                    vc_idx = bump_enames.index(self.vcorr_names[i])
                except:
                    raise RuntimeError(
                        "Three_bump: vcorr_name: %s not found" % self.vcorr_names[i]
                    )

                self.vcorr_elements.append(bump_elements[vc_idx])
                self.vcorr_idx.append(self.bump_idx[vc_idx])
            for e in self.vcorr_elements:
                e.set_double_attribute("vkick", 0.0)

        try:
            target_idx = bump_enames.index(self.target_name)
        except:
            raise RuntimeError(
                "Three_bump: target_name: %s not found" % self.target_name
            )

        self.target_elem = bump_elements[target_idx]

    ##################################################

    # print out information about the bump settings

    def information(self):
        print(
            "bump_lattice: ",
            len(self.bump_lattice.get_elements()),
            " elements, length: ",
            self.bump_lattice.get_length(),
        )
        print("horizontal correctors: ")
        if self.hcorr_elements is None:
            print("\t No horizontal correctors")
        else:
            for i in range(3):
                print(
                    f"\t{self.hcorr_elements[i].get_name()}, index: {self.hcorr_idx[i]}, hkick = {self.hcorr_elements[i].get_double_attribute('hkick'):.16g}"
                )

        print("vertical correctors: ")
        if self.vcorr_elements is None:
            print("\t No vertical correctors")
        else:
            for i in range(3):
                print(
                    f"\t{self.vcorr_elements[i].get_name()}, index: {self.vcorr_idx[i]}, vkick = {self.vcorr_elements[i].get_double_attribute('vkick'):.16g}"
                )

        print("target element name: ", self.target_name)

    ##################################################

    # Propagate particle at 0,0,0,0,0,0 through the bump section saving
    #     diagnostics and returning the orbit position at the target and
    #     at the end of the line as a tuple of array [6]

    def propagate_zero(self):
        if self.verbose:
            verbosity = 1
        else:
            verbosity = 0
        comm = synergia.utils.Commxx()
        refpart = self.bump_lattice.get_reference_particle()

        sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            refpart, 8, 0.5e11
        )
        bunch = sim.get_bunch(0, 0)
        lp = bunch.get_particles_numpy()
        lp[:, 0:6] = 0.0
        bunch.checkin_particles()

        # register diagnostics at the target
        diag_target = synergia.bunch.Diagnostics_bulk_track("target.h5", 1)
        diag_orbit = synergia.bunch.Diagnostics_bulk_track("orbit.h5", 1)
        sim.reg_diag_at_element(diag_target, self.target_elem)
        sim.reg_diag_per_turn(diag_orbit)

        stepper = synergia.simulation.Independent_stepper_elements(1)
        bump_propagator = synergia.simulation.Propagator(self.bump_lattice, stepper)

        simlog = synergia.utils.parallel_utils.Logger(
            0, synergia.utils.parallel_utils.LoggerV.ERROR
        )
        bump_propagator.propagate(sim, simlog, 1)

        del simlog
        del bump_propagator
        del stepper
        del diag_orbit
        del diag_target
        del lp
        del bunch
        del sim

        h5_target = h5py.File("target.h5", "r")
        h5_orbit = h5py.File("orbit.h5", "r")

        results = np.hstack(
            (
                h5_target.get("track_coords")[0, 0, 0:3:2],
                h5_orbit.get("track_coords")[1, 0, 0:4],
            )
        )

        h5_target.close()
        h5_orbit.close()

        return results

    ##################################################

    # adjust the correctors to achieve an orbit bump through the bump
    # section with the desired_position at the target element.

    #  Returns the
    #     final values for the correctors as array
    #     [hc1, hc2, hc3, vc1, vc2, vc3]
    #     desired position is a sequence of the desired (x,y) position
    #     at the target element

    def set_bump(self, desired_position):

        targets = np.array(
            [desired_position[0], desired_position[1], 0.0, 0.0, 0.0, 0.0]
        )

        ##################################################
        # Set the horizontal and vertical corrector values in preparation
        #    for the propagation of particles.
        # hcorr_values and vcorr_values are a length 3 sequence

        def set_corrector_elements(hcorr_values, vcorr_values):
            if self.hcorr_elements is not None:
                if self.verbose:
                    print("setting hcorr settings to ", hcorr_values)
                if len(self.hcorr_elements) != len(hcorr_values):
                    raise RuntimeError(
                        "set_corrector_elements: len(hcorr_elements) != len(hcorr_values)"
                    )
                for i in range(len(hcorr_values)):
                    self.hcorr_elements[i].set_double_attribute(
                        "hkick", hcorr_values[i]
                    )

            if self.vcorr_elements is not None:
                if self.verbose:
                    print("setting vcorr settings to ", vcorr_values)
                if len(self.vcorr_elements) != len(vcorr_values):
                    raise RuntimeError(
                        "set_corrector_elements: len(vcorr_elements) != len(vcorr_values)"
                    )
                for i in range(len(vcorr_values)):
                    self.vcorr_elements[i].set_double_attribute(
                        "vkick", vcorr_values[i]
                    )

            if self.verbose:
                print("bump_lattice after settings:", self.bump_lattice)

        ##################################################

        # utility function for fitting the bump.  given corrector settings
        # returns the positions and momenta of the 0 particle

        # x are the corrector settings (hc1, hc2, hc3, vc1, vc2, vc3).
        def bump_f(x):
            hc = (x[0], x[1], x[2])
            vc = (x[3], x[4], x[5])
            if self.verbose:
                print("propagating with hc: ", hc, ", vc: ", vc)
            set_corrector_elements(hc, vc)
            mean = self.propagate_zero()
            if self.verbose:
                print("\tresult: ", mean)
            residuals = mean - targets
            return residuals

        ##################################################

        init_guess = np.zeros(6)
        x_scale = np.array([1.0, 1.0, 1.0, 1.0 / 30.5, 1.0, 1.0 / 7.5])
        result = least_squares(
            bump_f,
            init_guess,
            ftol=1.0e-12,
            xtol=1.0e-12,
            gtol=1.0e-12,
            x_scale=x_scale,
            verbose=2,
        )
        if self.verbose:
            print("corrector values: ", result.x)
            print("cost: ", result.cost)
            print("residuals: ", result.fun)

        # copy the settings in the bump_lattice corrector elements
        # back to the original lattice elements
        if self.hcorr_elements is not None:
            for i in range(3):
                elem = self.lattice.get_elements()[self.hcorr_idx[i]]
                elem.set_double_attribute(
                    "hkick", self.hcorr_elements[i].get_double_attribute("hkick")
                )
        if self.vcorr_elements is not None:
            for i in range(3):
                elem = self.lattice.get_elements()[self.vcorr_idx[i]]
                elem.set_double_attribute(
                    "vkick", self.vcorr_elements[i].get_double_attribute("vkick")
                )

        if self.verbose:
            print("lattice")
            print(self.lattice)

        return result.x

    ##################################################
    ##################################################


#  just a little tester for the class


def create_lattice():
    lattice = synergia.lattice.Lattice("foo")
    dr1 = synergia.lattice.Lattice_element("drift", "dr1")
    dr1.set_double_attribute("l", 2.0)
    dr2a = synergia.lattice.Lattice_element("drift", "dr2a")
    dr2a.set_double_attribute("l", 1.5)
    dr2b = synergia.lattice.Lattice_element("drift", "dr2b")
    dr2b.set_double_attribute("l", 1.5)
    dr3 = synergia.lattice.Lattice_element("drift", "dr3")
    dr3.set_double_attribute("l", 3.0)
    dr4 = synergia.lattice.Lattice_element("drift", "dr4")
    dr4.set_double_attribute("l", 2.0)

    hk1 = synergia.lattice.Lattice_element("hkicker", "hk1")
    hk2 = synergia.lattice.Lattice_element("hkicker", "hk2")
    hk3 = synergia.lattice.Lattice_element("hkicker", "hk3")
    vk1 = synergia.lattice.Lattice_element("vkicker", "vk1")
    vk2 = synergia.lattice.Lattice_element("vkicker", "vk2")
    vk3 = synergia.lattice.Lattice_element("vkicker", "vk3")

    mbegin = synergia.lattice.Lattice_element("marker", "mbegin")
    mend = synergia.lattice.Lattice_element("marker", "mend")
    m = synergia.lattice.Lattice_element("marker", "m")

    lattice.append(dr1)
    lattice.append(mbegin)
    lattice.append(hk1)
    lattice.append(vk1)
    lattice.append(dr2a)
    lattice.append(m)
    lattice.append(dr2b)
    lattice.append(hk2)
    lattice.append(vk2)
    lattice.append(dr3)
    lattice.append(hk3)
    lattice.append(vk3)
    lattice.append(mend)
    lattice.append(dr4)

    refpart = synergia.foundation.Reference_particle(
        1, synergia.foundation.pconstants.mp, synergia.foundation.pconstants.mp + 0.4
    )
    lattice.set_reference_particle(refpart)

    return lattice


#####################################################################################################

if __name__ == "__main__":
    lattice = create_lattice()

    print(
        "read lattice: ",
        len(lattice.get_elements()),
        " elements, length = ",
        lattice.get_length(),
    )

    hcorr_names = ("hk1", "hk2", "hk3")
    vcorr_names = ("vk1", "vk2", "vk3")
    three_bump = Three_bump(
        lattice, "mbegin", "mend", hcorr_names, vcorr_names, "m", True
    )
    # three_bump = Three_bump(lattice, 'mbegin', 'mend', hcorr_names, None, 'm', True)

    # three_bump.information()
    # final_coords = three_bump.propagate_zero()
    # print('Final coords: ', final_coords)
    bump_settings = three_bump.set_bump((0.00075, -0.0005))

    print(
        "bump_settings: ",
        bump_settings[0],
        bump_settings[1],
        bump_settings[2],
        bump_settings[3],
        bump_settings[4],
        bump_settings[5],
    )
    three_bump.information()

    # propagate the whole lattice now
    comm = synergia.utils.Commxx()
    refpart = lattice.get_reference_particle()

    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        refpart, 8, 0.5e11
    )
    bunch = sim.get_bunch(0, 0)
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch.checkin_particles()

    # register diagnostics
    diag = synergia.bunch.Diagnostics_bulk_track("tracks.h5", 1)
    sim.reg_diag_per_step(diag)

    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.ERROR
    )
    # synergia.utils.parallel_utils.LoggerV.INFO_TURN)
    propagator.propagate(sim, simlog, 1)

    del propagator
    del simlog
    del stepper
    del diag
    del lp
    del bunch
    del sim
    # sys.exit(0)

    h5 = h5py.File("tracks.h5", "r")
    s = h5.get("track_s")[()]
    trks = h5.get("track_coords")[()]
    h5.close()

    import matplotlib.pyplot as plt

    f, ax = plt.subplots(2, 1, sharex=True)
    plt.suptitle("orbit")
    ax[0].plot(s, trks[:, 0, 0], label="X")
    ax[0].legend(loc="best")
    ax[1].plot(s, trks[:, 0, 2], label="Y")
    ax[1].legend(loc="best")

    ax[1].set_xlabel("s")

    plt.show()
