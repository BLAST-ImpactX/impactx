#include "fixed_t_z_converter.h"
#include "bunch.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/kokkos_views.h"
// #include "period.h"

#include <cmath>
#include <iostream>

//*********************Fixed_t_z_zeroth*****************************************

void
Fixed_t_z_zeroth::from_z_lab_to_t_bunch(Bunch& bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();

    auto particles = bunch.get_local_particles();
    auto s_particles = bunch.get_local_spectator_particles();

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    double gb = -gamma * beta;

    bool exception_flag = false;

    Kokkos::parallel_for(
        local_num,
        [particles, gamma, beta, m, p_ref, gb, &exception_flag](const int i) {
            // z in beam rest frame
            particles(i, Bunch::z) = gb * particles(i, Bunch::cdt);

            // total momentum in accelerator frame
            double p = p_ref + particles(i, Bunch::dpop) * p_ref;

            // E/c in accelerator frame
            double Eoc = std::sqrt(p * p + m * m);

            // p_{x,y,z} in accelerator frame
            double px = particles(i, Bunch::xp) * p_ref;
            double py = particles(i, Bunch::yp) * p_ref;
            double pz2 = p * p - px * px - py * py;

            if (pz2 < 0.0) exception_flag = true;
            double pz = std::sqrt(pz2);

            // zp = pz/p_{ref}^{total}
            particles(i, Bunch::zp) = gamma * (pz - beta * Eoc) / p_ref;

            // n.b. in the zeroth approximation, the transformation from
            //      t' = gamma cdt to t' = 0
            //      is a no-op.
        });

    Kokkos::parallel_for(
        local_s_num,
        [s_particles, gamma, beta, m, p_ref, gb, &exception_flag](const int i) {
            // z in beam rest frame
            s_particles(i, Bunch::z) = gb * s_particles(i, Bunch::cdt);

            // total momentum in accelerator frame
            double p = p_ref + s_particles(i, Bunch::dpop) * p_ref;

            // E/c in accelerator frame
            double Eoc = std::sqrt(p * p + m * m);

            // p_{x,y,z} in accelerator frame
            double px = s_particles(i, Bunch::xp) * p_ref;
            double py = s_particles(i, Bunch::yp) * p_ref;
            double pz2 = p * p - px * px - py * py;

            if (pz2 < 0.0) exception_flag = true;
            double pz = std::sqrt(pz2);

            // zp = pz/p_{ref}^{total}
            s_particles(i, Bunch::zp) = gamma * (pz - beta * Eoc) / p_ref;

            // n.b. in the zeroth approximation, the transformation from
            //      t' = gamma cdt to t' = 0
            //      is a no-op.
        });

#if 0
#pragma omp parallel for shared(                                               \
    local_num, gb, particles, gamma, beta, m, p_ref, exception_flag)
    for (int part = 0; part < local_num; ++part) 
    {
        // z in beam rest frame
        particles[part][Bunch::z] = gb * particles[part][Bunch::cdt];

        // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
        // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);
        // p_{x,y,z} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            exception_flag = true;
        }
        double pz = std::sqrt(pz2);
        // zp = pz/p_{ref}^{total}
        particles[part][Bunch::zp] = gamma * (pz - beta * Eoc) / p_ref;

        // n.b. in the zeroth approximation, the transformation from
        //      t' = gamma cdt to t' = 0
        //      is a no-op.
    }

#pragma omp parallel for shared(                                               \
    local_s_num, gb, s_particles, gamma, beta, m, p_ref, exception_flag)
    for (int part = 0; part < local_s_num; ++part) 
    {
        // z in beam rest frame
        s_particles[part][Bunch::z] = gb * s_particles[part][Bunch::cdt];

        // total momentum in accelerator frame
        double p = p_ref + s_particles[part][Bunch::dpop] * p_ref;
        // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);
        // p_{x,y,z} in accelerator frame
        double px = s_particles[part][Bunch::xp] * p_ref;
        double py = s_particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            exception_flag = true;
        }
        double pz = std::sqrt(pz2);
        // zp = pz/p_{ref}^{total}
        s_particles[part][Bunch::zp] = gamma * (pz - beta * Eoc) / p_ref;

        // n.b. in the zeroth approximation, the transformation from
        //      t' = gamma cdt to t' = 0
        //      is a no-op.
    }
#endif

    if (exception_flag) {
        throw std::runtime_error(
            "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
    }
}

void
Fixed_t_z_zeroth::from_t_bunch_to_z_lab(Bunch& bunch)
{
#if 0
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();

    auto particles = bunch.get_local_particles();
    auto s_particles = bunch.get_local_spectator_particles();

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    double gb = - gamma*beta;

#pragma omp parallel for shared(local_num, gb, particles, gamma, beta, m, p_ref)
    for (int part = 0; part < local_num; ++part) 
    {
        // ct in accelerator frame
        particles[part][Bunch::cdt] = particles[part][Bunch::z] / gb;

        // p'_{x,y,z} in beam frame
        double pxp = particles[part][Bunch::xp] * p_ref;
        double pyp = particles[part][Bunch::yp] * p_ref;
        double pzp = particles[part][Bunch::zp] * p_ref;
        double p_perp2 = pxp * pxp + pyp * pyp;
        // E'/c in beam frame
        double Epoc = std::sqrt(p_perp2 + pzp * pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);
        // dpop = (p - p_ref)/p_ref
        double p = std::sqrt(p_perp2 + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
    }

#pragma omp parallel for shared(                                               \
    local_s_num, gb, s_particles, gamma, beta, m, p_ref)
    for (int part = 0; part < local_s_num; ++part) 
    {
        // ct in accelerator frame
        s_particles[part][Bunch::cdt] = s_particles[part][Bunch::z] / gb;

        // p'_{x,y,z} in beam frame
        double pxp = s_particles[part][Bunch::xp] * p_ref;
        double pyp = s_particles[part][Bunch::yp] * p_ref;
        double pzp = s_particles[part][Bunch::zp] * p_ref;
        double p_perp2 = pxp * pxp + pyp * pyp;
        // E'/c in beam frame
        double Epoc = std::sqrt(p_perp2 + pzp * pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);
        // dpop = (p - p_ref)/p_ref
        double p = std::sqrt(p_perp2 + pz * pz);
        s_particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
    }
#endif
}

void
Fixed_t_z_zeroth::from_z_lab_to_t_lab(Bunch& bunch)
{
#if 0
    double beta = bunch.get_reference_particle().get_beta();
    double p_ref = bunch.get_reference_particle().get_momentum();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

#pragma omp parallel for shared(local_num, particles, beta, p_ref)
    for (int part = 0; part < local_num; ++part) 
    {
          // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
         // p_{x,y,z} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
        particles[part][Bunch::z] = - particles[part][Bunch::cdt]*beta;

        // zp = pz/p_{ref}^{total}
        particles[part][Bunch::zp] = pz/p_ref;
    }

#pragma omp parallel for shared(local_s_num, s_particles, beta, p_ref)
    for (int part = 0; part < local_s_num; ++part) 
    {
          // total momentum in accelerator frame
        double p = p_ref + s_particles[part][Bunch::dpop] * p_ref;
         // p_{x,y,z} in accelerator frame
        double px = s_particles[part][Bunch::xp] * p_ref;
        double py = s_particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
        s_particles[part][Bunch::z] = - s_particles[part][Bunch::cdt]*beta;

        // zp = pz/p_{ref}^{total}
        s_particles[part][Bunch::zp] = pz/p_ref;
    }
#endif
}

void
Fixed_t_z_zeroth::from_t_lab_to_z_lab(Bunch& bunch)
{
#if 0
    double beta = bunch.get_reference_particle().get_beta();
    double p_ref = bunch.get_reference_particle().get_momentum();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

#pragma omp parallel for shared(local_num, particles, beta, p_ref)
    for (int part = 0; part < local_num; ++part) 
    {
        // p'_{x,y,z} in beam rest frame
        double px = particles[part][Bunch::xp]* p_ref;
        double py = particles[part][Bunch::yp]* p_ref;
        double pz = particles[part][Bunch::zp]* p_ref;
        double p = std::sqrt(px * px + py * py + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
        particles[part][Bunch::cdt] = - particles[part][Bunch::z]/beta;
    }

#pragma omp parallel for shared(local_s_num, s_particles, beta, p_ref)
    for (int part = 0; part < local_s_num; ++part) 
    {
        // p'_{x,y,z} in beam rest frame
        double px = s_particles[part][Bunch::xp]* p_ref;
        double py = s_particles[part][Bunch::yp]* p_ref;
        double pz = s_particles[part][Bunch::zp]* p_ref;
        double p = std::sqrt(px * px + py * py + pz * pz);
        s_particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
        s_particles[part][Bunch::cdt] = - s_particles[part][Bunch::z]/beta;
    }
#endif
}

void
Fixed_t_z_zeroth::from_t_lab_to_t_bunch(Bunch& bunch)
{
#if 0
     double gamma=bunch.get_reference_particle().get_gamma();
     double beta = bunch.get_reference_particle().get_beta();
     double m = bunch.get_mass();
     double p_ref = bunch.get_reference_particle().get_momentum();

     MArray2d_ref particles = bunch.get_local_particles();
     MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

#pragma omp parallel for shared(local_num, particles, gamma, beta, m, p_ref)
    for (int part = 0; part < local_num; ++part) 
    {
        double px  = particles[part][Bunch::xp] * p_ref;
        double py  = particles[part][Bunch::yp] * p_ref;
        double pz  = particles[part][Bunch::zp] * p_ref;
        double p=std::sqrt(px * px + py * py + pz * pz);
        double Eoc   = std::sqrt(p * p + m * m);

        particles[part][Bunch::z]  = gamma * particles[part][Bunch::z];
        particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }

#pragma omp parallel for shared(local_s_num, s_particles, gamma, beta, m, p_ref)
    for (int part = 0; part < local_s_num; ++part) 
    {
        double px  = s_particles[part][Bunch::xp] * p_ref;
        double py  = s_particles[part][Bunch::yp] * p_ref;
        double pz  = s_particles[part][Bunch::zp] * p_ref;
        double p=std::sqrt(px * px + py * py + pz * pz);
        double Eoc   = std::sqrt(p * p + m * m);

        s_particles[part][Bunch::z]  = gamma * s_particles[part][Bunch::z];
        s_particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }
#endif
}

void
Fixed_t_z_zeroth::from_t_bunch_to_t_lab(Bunch& bunch)
{
#if 0
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

#pragma omp parallel for shared(local_num, particles, gamma, beta, m, p_ref)
    for (int part = 0; part < local_num; ++part) 
    {
          // p'_{x,y,z} in beam rest frame
        double pxp = particles[part][Bunch::xp] * p_ref;
        double pyp = particles[part][Bunch::yp] * p_ref;
        double pzp = particles[part][Bunch::zp] * p_ref;

        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);

        particles[part][Bunch::z] = particles[part][Bunch::z]/gamma;
        particles[part][Bunch::dpop] = pz / p_ref;
    }

#pragma omp parallel for shared(local_s_num, s_particles, gamma, beta, m, p_ref)
    for (int part = 0; part < local_s_num; ++part) 
    {
          // p'_{x,y,z} in beam rest frame
        double pxp = s_particles[part][Bunch::xp] * p_ref;
        double pyp = s_particles[part][Bunch::yp] * p_ref;
        double pzp = s_particles[part][Bunch::zp] * p_ref;

        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);

        s_particles[part][Bunch::z] = s_particles[part][Bunch::z]/gamma;
        s_particles[part][Bunch::dpop] = pz / p_ref;
    }
#endif
}

#if 0
template<class Archive>
    void
    Fixed_t_z_zeroth::serialize(Archive & ar, const unsigned int version)
    {
        ar &
        BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
    }

template
void
Fixed_t_z_zeroth::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Fixed_t_z_zeroth::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Fixed_t_z_zeroth::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Fixed_t_z_zeroth::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_zeroth)
#endif

//*********************Fixed_t_z_zeroth
//end*****************************************

// void
// Fixed_t_z_ballistic::from_t_bunch_to_z_lab(Bunch &bunch)
// {
//     std::cout << "stub: ballistic fixed_t_to_fixed_z\n";
// }
//
// void
// Fixed_t_z_ballistic::from_z_lab_to_t_bunch(Bunch &bunch)
// {
//     std::cout << "stub: ballistic fixed_z_to_fixed_t\n";
// }
//
// template<class Archive>
//     void
//     Fixed_t_z_ballistic::serialize(Archive & ar, const unsigned int version)
//     {
//         ar &
//         BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
//     }
//
// template
// void
// Fixed_t_z_ballistic::serialize<boost::archive::binary_oarchive >(
//         boost::archive::binary_oarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_ballistic::serialize<boost::archive::xml_oarchive >(
//         boost::archive::xml_oarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_ballistic::serialize<boost::archive::binary_iarchive >(
//         boost::archive::binary_iarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_ballistic::serialize<boost::archive::xml_iarchive >(
//         boost::archive::xml_iarchive & ar, const unsigned int version);
//
// BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_ballistic)
//
// //*********************Fixed_t_z_alex*********************************************
//
// void
// Fixed_t_z_alex::from_z_lab_to_t_bunch(Bunch &bunch)
// {
//     double gamma=bunch.get_reference_particle().get_gamma();
//     double gamma_inv = 1.0/gamma;
//     double beta = bunch.get_reference_particle().get_beta();
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//          // p_{x,y} in accelerator frame
//         double px = particles[part][Bunch::xp] * p_ref;
//         double py = particles[part][Bunch::yp] * p_ref;
//         // total momentum in accelerator frame
//         double p = p_ref + particles[part][Bunch::dpop] * p_ref;
//         double pz2 = p * p - px * px - py * py;
//         if (pz2 < 0.0) {
//             throw std::runtime_error(
//                     "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has
//                     negative pz^2");
//         }
//         double pz = std::sqrt(pz2);
//          // E/c in accelerator frame
//         double Eoc = std::sqrt(p * p + m * m);
//
//         double  betaz_part=pz/Eoc;
//         double  betax_part=px/Eoc;
//         double  betay_part=py/Eoc;
//         double  cdt=particles[part][Bunch::cdt];
//          // x,y,z in beam rest frame
//         particles[part][Bunch::z] = -
//         cdt*gamma_inv*betaz_part/(1.0-beta*betaz_part);
//         particles[part][Bunch::x] = particles[part][Bunch::x] -
//         cdt*betax_part/(1.0-beta*betaz_part); particles[part][Bunch::y] =
//         particles[part][Bunch::y] -  cdt*betay_part/(1.0-beta*betaz_part);
//
//         // zp = pz/p_{ref}^{total}, pz in the beam frame,  p_ref in the
//         accelerator frame
//        particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
//     }
// }
//
// void
// Fixed_t_z_alex::from_t_bunch_to_z_lab(Bunch &bunch)
// {
//     double gamma = bunch.get_reference_particle().get_gamma();
//     double beta = bunch.get_reference_particle().get_beta();
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//           // p'_{x,y,z} in beam rest frame
//         double pxp = particles[part][Bunch::xp]* p_ref;
//         double pyp = particles[part][Bunch::yp]* p_ref;
//         double pzp = particles[part][Bunch::zp]* p_ref;
//
//         // E'/c in beam rest frame
//         double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
//         // E/c in accelerator frame
//         double Eoc =gamma * (Epoc+ beta*pzp);
//         double pz = gamma * (pzp + beta * Epoc);
//
//
//
//         double  betaz_part=pz/Eoc;
//         double  betax_part=pxp/Eoc;
//         double  betay_part=pyp/Eoc;
//         double  cdt=-
//         particles[part][Bunch::z]*gamma*(1.0-beta*betaz_part)/betaz_part;
//         particles[part][Bunch::cdt] = cdt;
//         particles[part][Bunch::x] = particles[part][Bunch::x] +
//         cdt*betax_part/(1.0-beta*betaz_part); particles[part][Bunch::y] =
//         particles[part][Bunch::y] + cdt*betay_part/(1.0-beta*betaz_part);
//
//         double p = std::sqrt(pxp * pxp + pyp * pyp + pz * pz);
//         particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
//
//     }
//
// }
//
// void
// Fixed_t_z_alex::from_z_lab_to_t_lab(Bunch &bunch)
// {
//
//      double m = bunch.get_mass();
//      double p_ref = bunch.get_reference_particle().get_momentum();
//      MArray2d_ref particles = bunch.get_local_particles();
//      for (int part = 0; part < bunch.get_local_num(); ++part) {
//            // total momentum in accelerator frame
//          double p = p_ref + particles[part][Bunch::dpop] * p_ref;
//           // p_{x,y,z} in accelerator frame
//          double px = particles[part][Bunch::xp] * p_ref;
//          double py = particles[part][Bunch::yp] * p_ref;
//          double  cdt=particles[part][Bunch::cdt];
//          double pz2 = p * p - px * px - py * py;
//          if (pz2 < 0.0) {
//              throw std::runtime_error(
//                      "Fixed_t_z_alex::from_z_lab_to_t_lab: particle has
//                      negative pz^2");
//          }
//          double pz = std::sqrt(pz2);
//          double Eoc = std::sqrt(p * p + m * m);
//
//          double  betaz_part=pz/Eoc;
//          double  betax_part=px/Eoc;
//          double  betay_part=py/Eoc;
//            // x,y,z, zp in beam rest frame
//          particles[part][Bunch::z] = - cdt*betaz_part;
//          particles[part][Bunch::x] = particles[part][Bunch::x] -
//          cdt*betax_part; particles[part][Bunch::y] =
//          particles[part][Bunch::y] - cdt*betay_part;
//          // zp = pz/p_{ref}^{total}
//          particles[part][Bunch::zp] = pz/p_ref;
//
//      }
// }
//
//
// void
// Fixed_t_z_alex::from_t_lab_to_z_lab(Bunch &bunch)
// {
//      double m = bunch.get_mass();
//      double p_ref = bunch.get_reference_particle().get_momentum();
//      MArray2d_ref particles = bunch.get_local_particles();
//      for (int part = 0; part < bunch.get_local_num(); ++part) {
//          double px  = particles[part][Bunch::xp] * p_ref;
//          double py  = particles[part][Bunch::yp] * p_ref;
//          double pz  = particles[part][Bunch::zp]  * p_ref;
//          double p=std::sqrt(px * px + py * py + pz * pz);
//          double Eoc   = std::sqrt(p * p + m * m);
//          double  betaz_part=pz/Eoc;
//          double  betax_part=px/Eoc;
//          double  betay_part=py/Eoc;
//          double cdt = -particles[part][Bunch::z]/betaz_part;
//          particles[part][Bunch::cdt] = cdt;
//          particles[part][Bunch::x] = particles[part][Bunch::x] +
//          cdt*betax_part; particles[part][Bunch::y] =
//          particles[part][Bunch::y] +  cdt*betay_part;
//          particles[part][Bunch::dpop]=(p-p_ref)/p_ref;
//      }
//
// }
//
//
// void
// Fixed_t_z_alex::from_t_lab_to_t_bunch(Bunch &bunch)
// {
//      double gamma=bunch.get_reference_particle().get_gamma();
//      double gamma_inv = 1.0/gamma;
//      double beta = bunch.get_reference_particle().get_beta();
//      double m = bunch.get_mass();
//      double p_ref = bunch.get_reference_particle().get_momentum();
//      MArray2d_ref particles = bunch.get_local_particles();
//      for (int part = 0; part < bunch.get_local_num(); ++part) {
//          double px  = particles[part][Bunch::xp] * p_ref;
//          double py  = particles[part][Bunch::yp] * p_ref;
//          double pz  = particles[part][Bunch::zp]  * p_ref;
//          double p=std::sqrt(px * px + py * py + pz * pz);
//          double Eoc   = std::sqrt(p * p + m * m);
//          double  betaz_part=pz/Eoc;
//          double  betax_part=px/Eoc;
//          double  betay_part=py/Eoc;
//          double dz=particles[part][Bunch::z];
//          particles[part][Bunch::z] =  gamma_inv *
//          particles[part][Bunch::z]/(1.-beta* betaz_part);
//          particles[part][Bunch::x]=
//          particles[part][Bunch::x]+dz*beta*betax_part/(1.-beta* betaz_part);
//          particles[part][Bunch::y]=
//          particles[part][Bunch::y]+dz*beta*betay_part/(1.-beta* betaz_part);
//          particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
//     }
// }
//
// void
// Fixed_t_z_alex::from_t_bunch_to_t_lab(Bunch &bunch)
// {
//     double gamma = bunch.get_reference_particle().get_gamma();
//     double beta = bunch.get_reference_particle().get_beta();
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//           // p'_{x,y,z} in beam rest frame
//         double pxp = particles[part][Bunch::xp]* p_ref;
//         double pyp = particles[part][Bunch::yp]* p_ref;
//         double pzp = particles[part][Bunch::zp]* p_ref;
//
//         // E'/c in beam rest frame
//         double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
//         // E/c in accelerator frame
//         double Eoc =gamma * (Epoc+ beta*pzp);
//         double pz = gamma * (pzp + beta * Epoc);
//         double  betaz_part=pz/Eoc;
//         double  betax_part=pxp/Eoc;
//         double  betay_part=pyp/Eoc;
//         double dz = particles[part][Bunch::z]*gamma*(1.0-beta*betaz_part); //
//         in t_lab frame particles[part][Bunch::z] = dz;
//         particles[part][Bunch::x] =
//         particles[part][Bunch::x]-dz*beta*betax_part/(1.0-beta*betaz_part);
//         particles[part][Bunch::y] =
//         particles[part][Bunch::y]-dz*beta*betay_part/(1.0-beta*betaz_part);
//         particles[part][Bunch::dpop] = pz / p_ref;
//     }
// }
//
// template<class Archive>
//     void
//     Fixed_t_z_alex::serialize(Archive & ar, const unsigned int version)
//     {
//         ar &
//         BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
//     }
//
// template
// void
// Fixed_t_z_alex::serialize<boost::archive::binary_oarchive >(
//         boost::archive::binary_oarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_alex::serialize<boost::archive::xml_oarchive >(
//         boost::archive::xml_oarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_alex::serialize<boost::archive::binary_iarchive >(
//         boost::archive::binary_iarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_alex::serialize<boost::archive::xml_iarchive >(
//         boost::archive::xml_iarchive & ar, const unsigned int version);
//
// BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_alex)
//
// //*********************Fixed_t_z_alex
// end*********************************************
//
// //*********************Fixed_t_z_synergia20*****************************************
//
//
// void
// Fixed_t_z_synergia20::from_z_lab_to_t_bunch(Bunch &bunch)
// {
//     double gamma=bunch.get_reference_particle().get_gamma();
//     double beta = bunch.get_reference_particle().get_beta();
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//          // p_{x,y} in accelerator frame
//         double px = particles[part][Bunch::xp] * p_ref;
//         double py = particles[part][Bunch::yp] * p_ref;
//         // total momentum in accelerator frame
//         double p = p_ref + particles[part][Bunch::dpop] * p_ref;
//         double pz2 = p * p - px * px - py * py;
//         if (pz2 < 0.0) {
//             throw std::runtime_error(
//                     "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has
//                     negative pz^2");
//         }
//         double pz = std::sqrt(pz2);
//          // E/c in accelerator frame
//         double Eoc = std::sqrt(p * p + m * m);
//         double  betaz_part=pz/Eoc;
//
//         particles[part][Bunch::z] =
//         -particles[part][Bunch::cdt]*gamma*betaz_part;
//        // zp = pz/p_{ref}^{total}, pz in the beam frame,  p_ref in the
//        accelerator frame
//         particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
//     }
//
// }
//
//
// void Fixed_t_z_synergia20::from_t_bunch_to_z_lab(Bunch &bunch)
// {
//     double gamma = bunch.get_reference_particle().get_gamma();
//     double beta = bunch.get_reference_particle().get_beta();
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//           // p'_{x,y,z} in beam rest frame
//         double pxp = particles[part][Bunch::xp]* p_ref;
//         double pyp = particles[part][Bunch::yp]* p_ref;
//         double pzp = particles[part][Bunch::zp]* p_ref;
//
//         // E'/c in beam rest frame
//         double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
//         // E/c in accelerator frame
//         double Eoc =gamma * (Epoc+ beta*pzp);
//         double pz = gamma * (pzp + beta * Epoc);
//
//         double  betaz_part=pz/Eoc;
//         particles[part][Bunch::cdt] = -
//         particles[part][Bunch::z]/gamma/betaz_part;
//
//         double p = std::sqrt(pxp * pxp + pyp * pyp + pz * pz);
//         particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
//
//     }
// }
//
// void
// Fixed_t_z_synergia20::from_z_lab_to_t_lab(Bunch &bunch)
// {
//
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//           // total momentum in accelerator frame
//         double p = p_ref + particles[part][Bunch::dpop] * p_ref;
//          // p_{x,y,z} in accelerator frame
//         double px = particles[part][Bunch::xp] * p_ref;
//         double py = particles[part][Bunch::yp] * p_ref;
//         double pz2 = p * p - px * px - py * py;
//         if (pz2 < 0.0) {
//             throw std::runtime_error(
//                     "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has
//                     negative pz^2");
//         }
//         double pz = std::sqrt(pz2);
//          // E/c in accelerator frame
//         double Eoc = std::sqrt(p * p + m * m);
//
//         double  betaz_part=pz/Eoc;
//         particles[part][Bunch::z] = - particles[part][Bunch::cdt]*betaz_part;
//
//         // zp = pz/p_{ref}^{total}
//          particles[part][Bunch::zp] =pz/p_ref;
//
//     }
// }
//
// void Fixed_t_z_synergia20::from_t_lab_to_z_lab(Bunch &bunch)
// {
//
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//
//    for (int part = 0; part < bunch.get_local_num(); ++part) {
//           // p'_{x,y,z} in beam rest frame
//         double px = particles[part][Bunch::xp]* p_ref;
//         double py = particles[part][Bunch::yp]* p_ref;
//         double pz = particles[part][Bunch::zp]* p_ref;
//
//         double Eoc = std::sqrt(px*px + py*py + pz*pz + m * m);
//         double p = std::sqrt(px * px + py * py + pz * pz);
//         particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
//
//         double  betaz_part=pz/Eoc;
//         particles[part][Bunch::cdt] = - particles[part][Bunch::z]/betaz_part;
//     }
// }
//
//
//
// void
// Fixed_t_z_synergia20::from_t_lab_to_t_bunch(Bunch &bunch)
// {
//      double gamma=bunch.get_reference_particle().get_gamma();
//      double beta = bunch.get_reference_particle().get_beta();
//      double m = bunch.get_mass();
//      double p_ref = bunch.get_reference_particle().get_momentum();
//      MArray2d_ref particles = bunch.get_local_particles();
//      for (int part = 0; part < bunch.get_local_num(); ++part) {
//          double px  = particles[part][Bunch::xp] * p_ref;
//          double py  = particles[part][Bunch::yp] * p_ref;
//          double pz  = particles[part][Bunch::zp]  * p_ref;
//          double p=std::sqrt(px * px + py * py + pz * pz);
//          double Eoc   = std::sqrt(p * p + m * m);
//
//          particles[part][Bunch::z] =  gamma * particles[part][Bunch::z];
//
//          particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
//     }
// }
// void
// Fixed_t_z_synergia20::from_t_bunch_to_t_lab(Bunch &bunch)
// {
//     double gamma = bunch.get_reference_particle().get_gamma();
//     double beta = bunch.get_reference_particle().get_beta();
//     double m = bunch.get_mass();
//     double p_ref = bunch.get_reference_particle().get_momentum();
//     MArray2d_ref particles = bunch.get_local_particles();
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//           // p'_{x,y,z} in beam rest frame
//         double pxp = particles[part][Bunch::xp]* p_ref;
//         double pyp = particles[part][Bunch::yp]* p_ref;
//         double pzp = particles[part][Bunch::zp]* p_ref;
//         // E'/c in beam rest frame
//         double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
//         double pz = gamma * (pzp + beta * Epoc);
//         particles[part][Bunch::z] = particles[part][Bunch::z]/gamma;
//         particles[part][Bunch::dpop] = pz / p_ref;
//
//     }
// }
//
// template<class Archive>
//     void
//     Fixed_t_z_synergia20::serialize(Archive & ar, const unsigned int version)
//     {
//         ar &
//         BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
//     }
//
// template
// void
// Fixed_t_z_synergia20::serialize<boost::archive::binary_oarchive >(
//         boost::archive::binary_oarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_synergia20::serialize<boost::archive::xml_oarchive >(
//         boost::archive::xml_oarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_synergia20::serialize<boost::archive::binary_iarchive >(
//         boost::archive::binary_iarchive & ar, const unsigned int version);
//
// template
// void
// Fixed_t_z_synergia20::serialize<boost::archive::xml_iarchive >(
//         boost::archive::xml_iarchive & ar, const unsigned int version);
//
// BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_synergia20)
//
//
// //*********************Fixed_t_z_synergia20
// end*****************************************
