#include "period.h"
#include "synergia/bunch/bunch.h"

namespace {
    struct Zcut_aperture {
        double half_length;

        KOKKOS_INLINE_FUNCTION
        bool
        discard(ConstParticles const& parts,
                ConstParticleMasks const&,
                int p) const
        {
            if (std::abs(parts(p, 4)) > half_length)
                return true;
            else
                return false;
        }
    };

    struct alg_z_period {
        Particles p;
        ConstParticleMasks masks;

        double length_cdt;
        double half_length;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                double tmp = p(i, 4) + half_length;
                if (tmp > 0)
                    p(i, 4) = fmod(tmp, length_cdt) - half_length;
                else
                    p(i, 4) = fmod(tmp, length_cdt) + half_length;
            }
        }
    };

    struct alg_bucket_barrier {
        Particles p;
        ConstParticleMasks masks;

        double length_cdt;
        double half_length;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                double z = p(i, 4);

                if (z > half_length || z < -half_length) {
                    double px = p(i, 1);
                    double py = p(i, 3);
                    double dpop = p(i, 5);

                    double px2 = px * px;
                    double py2 = py * py;
                    double dp2 = (dpop + 1.0) * (dpop + 1.0);

                    // flip the longitudinal momentum
                    p(i, 5) =
                        sqrt(4.0 + dp2 - 4.0 * sqrt(dp2 - px2 - py2)) - 1.0;

                    // adjust the z position
                    if (z > half_length) {
                        double off = fmod(z + half_length, length_cdt);
                        p(i, 4) = half_length - off;
                    } else {
                        double off = fmod(z - half_length, length_cdt);
                        p(i, 4) = -half_length - off;
                    }
                }
            }
        }
    };

}

void
apply_longitudinal_periodicity(Bunch& bunch, double length)
{
    // Bunch::State state=bunch.get_state();
    // bunch.convert_to_state(Bunch::fixed_z_lab);

    double beta = bunch.get_reference_particle().get_beta();
    double length_cdt = length / beta;
    double half_length = 0.5 * length_cdt;

    if (bunch.get_local_num(ParticleGroup::regular)) {
        alg_z_period alg{bunch.get_local_particles(ParticleGroup::regular),
                         bunch.get_local_particle_masks(ParticleGroup::regular),
                         length_cdt,
                         half_length};

        Kokkos::parallel_for(bunch.get_local_num(ParticleGroup::regular), alg);
    }

    if (bunch.get_local_num(ParticleGroup::spectator)) {
        alg_z_period alg{
            bunch.get_local_particles(ParticleGroup::spectator),
            bunch.get_local_particle_masks(ParticleGroup::spectator),
            length_cdt,
            half_length};

        Kokkos::parallel_for(bunch.get_local_num(ParticleGroup::spectator),
                             alg);
    }

    // bunch.convert_to_state(state);
}

void
apply_longitudinal_bucket_barrier(Bunch& bunch, double length)
{
    // Bunch::State state=bunch.get_state();
    // bunch.convert_to_state(Bunch::fixed_z_lab);

    double beta = bunch.get_reference_particle().get_beta();
    double length_cdt = length / beta;
    double half_length = 0.5 * length_cdt;

    if (bunch.get_local_num(ParticleGroup::regular)) {
        alg_bucket_barrier alg{
            bunch.get_local_particles(ParticleGroup::regular),
            bunch.get_local_particle_masks(ParticleGroup::regular),
            length_cdt,
            half_length};

        Kokkos::parallel_for(bunch.get_local_num(ParticleGroup::regular), alg);
    }

    if (bunch.get_local_num(ParticleGroup::spectator)) {
        alg_bucket_barrier alg{
            bunch.get_local_particles(ParticleGroup::spectator),
            bunch.get_local_particle_masks(ParticleGroup::spectator),
            length_cdt,
            half_length};

        Kokkos::parallel_for(bunch.get_local_num(ParticleGroup::spectator),
                             alg);
    }

    // bunch.convert_to_state(state);
}

void
apply_zcut(Bunch& bunch, double length)
{
    // Bunch::State state = bunch.get_state();
    // bunch.convert_to_state(Bunch::fixed_z_lab);

    double beta = bunch.get_reference_particle().get_beta();
    double half_length = 0.5 * length / beta;
    bunch.apply_zcut(Zcut_aperture{half_length});

    // bunch.convert_to_state(state);
}
