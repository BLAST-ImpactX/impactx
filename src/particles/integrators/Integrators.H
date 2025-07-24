/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_INTEGRATORS_H_
#define IMPACTX_INTEGRATORS_H_

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>       // for ParticleReal


namespace impactx::integrators
{

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is a generalization of the second-
    *  order leapfrog algorithm.  For a detailed overview:
    *
    *  E. Hairer et al, Geometric Numerical Integration:  Structure-
    *  Preserving Algorithms for Ordinary Differential Equations,
    *  2nd ed, Springer, Berlin, 2006.
    *
    * @param refpart  Reference particle data
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate (
        RefPart & refpart,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const tau1 = dz*0.5_prt;
        amrex::ParticleReal const tau2 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau1,refpart,zeval);
        }
    }

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2 + H_3.  This is the result of applying
    *  the two-terms splitting of symp2_integrate twice in succession.
    *  The form shown here appears, for example, in:
    *
    *  R. D. Ryne, Computational Methods in Accelerator Physics,
    *  USPAS lecture notes, 2009.
    *
    * @param refpart  Reference particle data
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate_split3 (
        RefPart & refpart,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const tau1 = dz*0.5_prt;
        amrex::ParticleReal const tau2 = dz*0.5_prt;
        amrex::ParticleReal const tau3 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map3(tau3,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau1,refpart,zeval);
        }
    }

   /** A fourth-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is the result of applying
    *  the two-terms splitting of symp2_integrate together with
    *  the method of Yoshida.
    *
    *  H. Yoshida, Phys. Lett. A 150, 292-268 (1990).
    *
    *  The form shown here appears, for example, in:
    *
    *  E. Forest and R. D. Ruth, Physica D 43, 105-117 (1990).
    *
    * @param refpart  Reference particle data
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp4_integrate (
        RefPart & refpart,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const alpha = 1_prt - std::pow(2_prt,1_prt/3_prt);
        amrex::ParticleReal const tau2 = dz/(1_prt + alpha);
        amrex::ParticleReal const tau1 = tau2*0.5_prt;
        amrex::ParticleReal const tau3 = alpha*tau1;
        amrex::ParticleReal const tau4 = (alpha - 1_prt)*tau2;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for (int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau3,refpart,zeval);
            element.map2(tau4,refpart,zeval);
            element.map1(tau3,refpart,zeval);
            element.map2(tau2,refpart,zeval);
            element.map1(tau1,refpart,zeval);
        }
    }

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is a generalization of the second-
    *  order leapfrog algorithm.  For a detailed overview:
    *
    *  E. Hairer et al, Geometric Numerical Integration:  Structure-
    *  Preserving Algorithms for Ordinary Differential Equations,
    *  2nd ed, Springer, Berlin, 2006.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp2_integrate_particle (
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const tau1 = dz*0.5_prt;
        amrex::ParticleReal const tau2 = dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval);
            element.map2(tau2,particle,zeval);
            element.map1(tau1,particle,zeval);
        }
    }


   /** A fourth-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is the result of applying
    *  the two-terms splitting of symp2_integrate together with
    *  the method of Yoshida.
    *
    *  H. Yoshida, Phys. Lett. A 150, 292-268 (1990).
    *
    *  The form shown here appears, for example, in:
    *
    *  E. Forest and R. D. Ruth, Physica D 43, 105-117 (1990).
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp4_integrate_particle (
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const alpha = 1_prt - std::pow(2_prt,1_prt/3_prt);
        amrex::ParticleReal const tau2 = dz/(1_prt + alpha);
        amrex::ParticleReal const tau1 = tau2*0.5_prt;
        amrex::ParticleReal const tau3 = alpha*tau1;
        amrex::ParticleReal const tau4 = (alpha - 1_prt)*tau2;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for (int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval);
            element.map2(tau2,particle,zeval);
            element.map1(tau3,particle,zeval);
            element.map2(tau4,particle,zeval);
            element.map1(tau3,particle,zeval);
            element.map2(tau2,particle,zeval);
            element.map1(tau1,particle,zeval);
        }
    }


   /** A sixth-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  The implementation below
    *  requires a reduced number of map evaluations (15 total)
    *  compared to direct use of the Yoshida construction, which
    *  requires three compositions of a fourth-order integrator step,
    *  each involving seven map evaluations (when expressed
    *  in the form obtained by Forest and Ruth, above).
    *  It is based on Section 5 of:
    *
    *  H. Yoshida, Phys. Lett. A 150, 292-268 (1990),
    *
    *  and appears also in Section V.3.2 (3.11) of:
    *
    *  E. Hairer et al, Geometric Numerical Integration, 2nd ed.,
    *  Springer-Verlag, Berlin, 2006.
    *
    * @param particle  Particle phase space data in a 6-vector
    * @param zin  Initial value of independent variable (z-location)
    * @param zout  Final value of independent variable (z-location)
    * @param nsteps  Number of integration steps
    * @param element  Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp6_integrate_particle (
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;
        amrex::ParticleReal const w0 = 1.31518632068391121888424973;
        amrex::ParticleReal const w1 = -1.17767998417887100694641568;
        amrex::ParticleReal const w2 = 0.23557321335935813368479318;
        amrex::ParticleReal const w3 = 0.78451361047755726381949763;

        amrex::ParticleReal const tau1 = w3 * dz / 2_prt;
        amrex::ParticleReal const tau2 = w3 * dz;
        amrex::ParticleReal const tau3 = (w3 + w2) * dz / 2_prt;
        amrex::ParticleReal const tau4 = w2 * dz;
        amrex::ParticleReal const tau5 = (w2 + w1) * dz / 2_prt;
        amrex::ParticleReal const tau6 = w1 * dz;
        amrex::ParticleReal const tau7 = (w1 + w0) * dz / 2_prt;
        amrex::ParticleReal const tau8 = w0 * dz;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for (int j=0; j < nsteps; ++j)
        {
            element.map1(tau1,particle,zeval);
            element.map2(tau2,particle,zeval);
            element.map1(tau3,particle,zeval);
            element.map2(tau4,particle,zeval);
            element.map1(tau5,particle,zeval);
            element.map2(tau6,particle,zeval);
            element.map1(tau7,particle,zeval);
            element.map2(tau8,particle,zeval);
            element.map1(tau7,particle,zeval);
            element.map2(tau6,particle,zeval);
            element.map1(tau5,particle,zeval);
            element.map2(tau4,particle,zeval);
            element.map1(tau3,particle,zeval);
            element.map2(tau2,particle,zeval);
            element.map1(tau1,particle,zeval);
        }
    }


   /** Applying a map on the extended phase space, for use within
    *  semiexplicit_integrator_step.  Additional documentation to be added.
    *
    *
    * @param tau Step size
    * @param particle1 Particle phase space data in a 6-vector
    * @param particle2 Particle phase space data in a 6-vector
    * @param element Element defining the two maps associated with H_1 and H_2
    * @param zeval Evaluation point in s
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void extended_phase_space_map(
        amrex::ParticleReal const tau,
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle1,
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle2,
        T_Element const & element,
        amrex::ParticleReal & zeval
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const tau1 = tau/2.0_prt;
        amrex::ParticleReal const tau2 = tau;

        // Apply the map on the extended phase space
        element.map(tau1,particle1,particle2,zeval);
        element.map(tau2,particle2,particle1,zeval);
        element.map(tau1,particle1,particle2,zeval);

    }


   /** A semiexplicit symplectic integrator for a general
    *  Hamiltonian, which does not require splitting,  based on:
    *
    *  B. Jayawardana and T. Ohsawa, "Semiexplicit symplectic
    *  integrators for non-separable Hamiltonian systems,"
    *  Math. Comput. 92, pp. 251-281 (2022).
    *
    * @param particle Particle phase space data in a 6-vector
    * @param tau Size of the step
    * @param element Element defining the two maps associated with H_1 and H_2
    * @param zeval Evaluation point in s
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void semiexplicit_integrator_step(
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
        amrex::ParticleReal const tau,
        T_Element const & element,
        amrex::ParticleReal & zeval
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // Parameters defining the implicit iteration loop
        amrex::ParticleReal const tolerance = 1.0e-12;
        int max_iterations = 100;
        int count = 0;

        // Initialize the perturbation vector
        amrex::SmallVector<amrex::ParticleReal, 6, 1> perturbation{};

        // Iteration loop (implicit step)
        for (int j = 0; j < max_iterations; j++) {

           // Create the extended phase space by creating two copies of the initial condition,
           // and applying the initial perturbation.
           amrex::SmallVector<amrex::ParticleReal, 6, 1> particle_copy1 = particle + perturbation;
           amrex::SmallVector<amrex::ParticleReal, 6, 1> particle_copy2 = particle - perturbation;

           // Apply the map for a step of size tau on the extended phase space.
           integrators::extended_phase_space_map(tau,particle_copy1,particle_copy2,element,zeval);

           // Apply the final perturbation, and evaluate the defect (difference between copies)
           particle_copy1 = particle_copy1 + perturbation;
           particle_copy2 = particle_copy2 - perturbation;
           amrex::SmallVector<amrex::ParticleReal, 6, 1> defect = particle_copy1 - particle_copy2;

           // Use the defect to update the perturbation vector.  Here, we use the simplified
           // Newton approximation (23) from Buddhika et al.
           perturbation = perturbation - (1_prt/4_prt) * defect;

           // Compute the norm of the phase space vector defect
           amrex::ParticleReal defect_norm2 = (defect.transpose()*defect).trace();
           amrex::ParticleReal defect_norm = std::sqrt(defect_norm2);

           // Update particle phase space vector once the desired tolerance is reached
           if (defect_norm < tolerance) {
              particle = particle_copy1;
              break;
           }
           count = count + 1;
       }

       // A statement for debugging:
       // amrex::Print() << "Number of iterations required: " << count << "\n";

    }

   /** A second-order symplectic integrator based on a Hamiltonian
    *  splitting H = H_1 + H_2.  This is a generalization of the second-
    *  order leapfrog algorithm.  For a detailed overview:
    *
    *  E. Hairer et al, Geometric Numerical Integration:  Structure-
    *  Preserving Algorithms for Ordinary Differential Equations,
    *  2nd ed, Springer, Berlin, 2006.
    *
    * @param particle Particle phase space data in a 6-vector
    * @param zin Initial value of independent variable (z-location)
    * @param zout Final value of independent variable (z-location)
    * @param nsteps Number of integration steps
    * @param element Element defining the two maps associated with H_1 and H_2
    */
    template <typename T_Element>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void symp_integrate_semiexplicit (
        amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const zout,
        int const nsteps,
        T_Element const & element
    )
    {
        using namespace amrex::literals; // for _rt and _prt

        // initialize numerical integration parameters
        amrex::ParticleReal const dz = (zout-zin)/nsteps;

        // initialize the value of the independent variable
        amrex::ParticleReal zeval = zin;

        // loop over integration steps
        for(int j=0; j < nsteps; ++j)
        {
            semiexplicit_integrator_step(particle,dz,element,zeval);
        }
    }



} // namespace impactx::integrators

#endif // IMPACTX_INTEGRATORS_H_
