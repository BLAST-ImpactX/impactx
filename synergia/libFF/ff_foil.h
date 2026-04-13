#ifndef FF_FOIL_H
#define FF_FOIL_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_Random.hpp>

namespace foil_impl {
    KOKKOS_INLINE_FUNCTION
    double
    get_a()
    {
        return 12.01;
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_z()
    {
        return 6.0;
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_rho()
    {
        // return 3.52e3;
        return 2.27e3;
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_elastic_crosssection(double E)
    {
        const double energy[59] = {
            .0005, .001,  .0015, 0.002,  .0025, .003,  0.0035, .004,  .0045,
            0.005, .0055, .006,  0.0065, .007,  .0075, 0.008,  .009,  .01,
            0.011, .012,  .013,  0.014,  .015,  .0175, 0.02,   .0225, .025,
            .030,  0.035, .040,  .045,   .050,  .0510, .055,   .060,  .070,
            .080,  .090,  .100,  .110,   .120,  .140,  .160,   .180,  .200,
            .225,  .250,  .275,  .300,   .325,  .350,  .375,   .400,  .500,
            .700,  1.000, 1.500, 2.000,  2.500};

        const double z6_elastic[59] = {
            0.,    0.,    0.,    0.,    0.,    0.321, 0.335, 0.349, 0.363,
            0.377, 0.386, 0.395, 0.404, 0.413, 0.423, 0.434, 0.456, 0.479,
            0.509, 0.539, 0.569, 0.598, 0.628, 0.685, 0.743, 0.783, 0.822,
            0.878, 0.915, 0.938, 0.813, 0.688, 0.678, 0.641, 0.594, 0.515,
            0.448, 0.392, 0.345, 0.305, 0.272, 0.214, 0.167, 0.138, 0.117,
            0.098, 0.085, 0.077, .072,  0.068, 0.067, 0.066, 0.067, 0.077,
            0.09,  0.102, 0.108, 0.114, 0.113};

        double energytracker = 0.0;
        int ilow = 57;
        for (int i = 0; i < 58; i++) {
            if (E <= energy[i + 1]) {
                ilow = i;
                break;
            }
        }

        const double* ec = z6_elastic;

        double efrac = (E - energy[ilow]) / (energy[ilow + 1] - energy[ilow]);
        double cross_section = ec[ilow] + efrac * (ec[ilow + 1] - ec[ilow]);
        return cross_section;
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_inelastic_crosssection(double E)
    {
        const double energy[59] = {
            .0005, .001,  .0015, 0.002,  .0025, .003,  0.0035, .004,  .0045,
            0.005, .0055, .006,  0.0065, .007,  .0075, 0.008,  .009,  .01,
            0.011, .012,  .013,  0.014,  .015,  .0175, 0.02,   .0225, .025,
            .030,  0.035, .040,  .045,   .050,  .0510, .055,   .060,  .070,
            .080,  .090,  .100,  .110,   .120,  .140,  .160,   .180,  .200,
            .225,  .250,  .275,  .300,   .325,  .350,  .375,   .400,  .500,
            .700,  1.000, 1.500, 2.000,  2.500};

        const double z6_inelastic[59] = {
            0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
            0.013, 0.078, 0.131, 0.175, 0.212, 0.244, 0.271, 0.314, 0.346,
            0.370, 0.389, 0.403, 0.413, 0.421, 0.432, 0.434, 0.432, 0.427,
            0.412, 0.394, 0.376, 0.359, 0.344, 0.288, 0.281, 0.272, 0.257,
            0.245, 0.235, 0.228, 0.223, 0.220, 0.223, 0.213, 0.212, 0.212,
            0.212, 0.213, 0.214, 0.216, 0.217, 0.219, 0.221, 0.223, 0.229,
            0.246, 0.261, 0.261, 0.257, 0.255};

        double energytracker = 0.0;
        int ilow = 57;
        for (int i = 0; i < 58; i++) {
            if (E <= energy[i + 1]) {
                ilow = i;
                break;
            }
        }

        const double* ic = z6_inelastic;

        double efrac = (E - energy[ilow]) / (energy[ilow + 1] - energy[ilow]);
        double cross_section = ic[ilow] + efrac * (ic[ilow + 1] - ic[ilow]);
        return cross_section;
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_direction(double xp, double yp, double dpop)
    {
        double xpfac = xp / (1.0 + dpop);
        double ypfac = yp / (1.0 + dpop);

        return sqrt(1.0 + xpfac * xpfac + ypfac * ypfac);
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_beta(double dpop, double pref, double m)
    {
        double p = (dpop + 1.0) * pref;
        double E2 = p * p + m * m;
        double beta = sqrt((E2 - m * m) / E2);
        return beta;
    }

    KOKKOS_INLINE_FUNCTION
    double
    get_p(double dpop, double pref)
    {
        return (dpop + 1.0) * pref;
    }

    KOKKOS_INLINE_FUNCTION
    double
    ion_energy_loss(double beta, double z, double a)
    {
        double dE, gamma, T_max, m_e = .511, M = 938.27, I, arg;

        I = z * 10.;
        gamma = 1 / sqrt(1 - beta * beta);
        T_max = 1e6 * (2 * m_e * beta * beta * gamma * gamma) /
                (1 + 2 * gamma * m_e / M + pow(m_e / M, 2.0));
        arg = 2 * (m_e * 1e6) * beta * beta * gamma * gamma * T_max / (I * I);
        dE = 0.307075 * z / a / (beta * beta) * (log(arg) / 2.0 - beta * beta);

        // Unit conversion from MeV/(g/cm^2) to GeV/(kg/m^2)
        dE /= 10000;

        return dE;
    }

    KOKKOS_INLINE_FUNCTION
    void
    mcs_jackson(Kokkos::Random_XorShift64_Pool<>::generator_type& rand_gen,
                double& x,
                double& px,
                double& y,
                double& py,
                double dpop,
                double beta,
                double stepsize,
                double z,
                double a,
                double rho)
    {
        const double elementary_charge_CGS = 4.8032068e-10;

        // Convert to mm and mrad for this routine.	And density to g/cm3
        x *= 1000.0;
        y *= 1000.0;
        px *= 1000.0;
        py *= 1000.0;
        stepsize *= 1000.0;
        rho /= 1000.0;

        double pfac = dpop + 1.0;

        double pi = Kokkos::numbers::pi_v<double>;
        double cvel = pconstants::c * 100.0;
        double eesu = elementary_charge_CGS;
        double hbar = 1.05443e-27;
        double AMU = 1.65979e-24;
        double aBohr = 5.29172e-09;
        double AP = 1.007593;
        double gamma = 1.0 / sqrt(1.0 - beta * beta);

        double aMax = 1.4 * aBohr / pow(z, 0.333333);
        double aMin = 1.4e-13 * pow(a, 0.333333);

        double N = rho / (a * AMU);

        double T = stepsize / 10.0;

        double pmom = gamma * beta * cvel * AP * AMU;

        double thMin = hbar / (pmom * aMax);
        double thMax = hbar / (pmom * aMin);

        double thMin2i = 1.0 / (thMin * thMin);
        double thMax2i = 1.0 / (thMax * thMax);
        double th2iDiff = thMin2i - thMax2i;

        double th2s = 2.0 * log(thMax / thMin) / th2iDiff;
        double coeff = 2.0 * z * eesu * eesu / (pmom * beta * cvel);

        double sigTot = pi * coeff * coeff * th2iDiff;
        double nColl = N * T * sigTot;

        double th2Tot = nColl * th2s;

        // double probrp = Random::ran1(idum);
        // double probxy = 2.0 * pi * Random::ran1(idum);
        double probrp = rand_gen.drand();
        double probxy = 2.0 * pi * rand_gen.drand();

        double angle = sqrt(-th2Tot * log(probrp));
        double anglexMCS = angle * cos(probxy);
        double angleyMCS = angle * sin(probxy);

        double xpfac = px / (1000. * pfac);
        double ypfac = py / (1000. * pfac);

        double anglex = atan(xpfac) + anglexMCS;
        double angley = atan(ypfac) + angleyMCS;

        double tanglex = tan(anglex);
        double tangley = tan(angley);

        px = tanglex * (1000. * pfac);
        py = tangley * (1000. * pfac);

        double directionfac = sqrt(1.0 + xpfac * xpfac + ypfac * ypfac);
        double zstep = stepsize / directionfac;
        x += zstep * (xpfac + tanglex) / 2.0;
        y += zstep * (ypfac + tangley) / 2.0;

        // Convert back to m and rad.
        x /= 1000.0;
        y /= 1000.0;
        px /= 1000.0;
        py /= 1000.0;
    }

    // returns true if the particle is lost
    KOKKOS_INLINE_FUNCTION
    bool
    take_step(Kokkos::Random_XorShift64_Pool<>::generator_type& rand_gen,
              double& x,
              double& xp,
              double& y,
              double& yp,
              double& dpop,
              double pref,
              double m,
              double z,
              double a,
              double density,
              double stepsize)
    {
        double beta = get_beta(dpop, pref, m);
        double pfac = dpop + 1.0;
        double p = pfac * pref;
        double E = sqrt(p * p + m * m);

        mcs_jackson(
            rand_gen, x, xp, y, yp, dpop, beta, stepsize, z, a, density);

        double dE = ion_energy_loss(beta, z, a);

        // factors for units m->cm and MeV->GeV
        dE = -dE * density * stepsize;

        // dp/pref
        dpop = sqrt((E + dE) * (E + dE) - m * m) / pref - 1.0;

        // lost or not
        return (E + dE) < 0.02 ? true : false;
    }

    KOKKOS_INLINE_FUNCTION
    double
    ruth_scatt_jackson(
        Kokkos::Random_XorShift64_Pool<>::generator_type& rand_gen,
        double stepsize,
        double z,
        double a,
        double rho,
        double beta,
        bool trackit,
        double pfac,
        double& thetax,
        double& thetay)
    {
        const double elementary_charge_CGS = 4.8032068e-10;

        stepsize *= 1000.0; // Convert to mm
        rho /= 1000.0;      // Convert to g/cm3

        double theta[2];
        theta[0] = 0.0;
        theta[1] = 0.0;

        double pi = Kokkos::numbers::pi_v<double>;
        double cvel = pconstants::c * 100.0;
        double eesu = elementary_charge_CGS;
        double hbar = 1.05443e-27;
        double AMU = 1.65979e-24;
        double aBohr = 5.29172e-09;
        double AP = 1.007593;
        double gamma = 1.0 / sqrt(1.0 - beta * beta);

        double aMax = 1.4 * aBohr / pow(z, 0.333333);
        double aMin = 1.4e-13 * pow(a, 0.333333);

        double N = rho / (a * AMU);
        double T = stepsize / 10.0;

        double pmom = gamma * beta * cvel * AP * AMU;

        double thMin = hbar / (pmom * aMax);
        double thMax = hbar / (pmom * aMin);

        double thMin2i = 1.0 / (thMin * thMin);
        double thMax2i = 1.0 / (thMax * thMax);
        double th2iDiff = thMin2i - thMax2i;

        double th2s = 2.0 * log(thMax / thMin) / th2iDiff;
        double coeff = 2.0 * z * eesu * eesu / (pmom * beta * cvel);

        double sigTot = pi * coeff * coeff * th2iDiff;

        double nColl = N * T * sigTot;

        double th2Tot = nColl * th2s;

        if (thMin < (2.0 * sqrt(th2Tot))) thMin = 2.0 * sqrt(th2Tot);

        thMin2i = 1.0 / (thMin * thMin);
        th2iDiff = thMin2i - thMax2i;
        double rcross = 1.e+24 * pi * coeff * coeff * th2iDiff;

        if (rcross < 0.0) rcross = 0.0;

        thMin *= 1.4;
        thMin2i = 1.0 / (thMin * thMin);
        th2iDiff = thMin2i - thMax2i;

        if (trackit) {
            double thx = 0.0;
            double thy = 0.0;

            if (thMin < thMax) {
                // double probrp = Random::ran1(idum);
                // double probxy = 2.0 * pi * Random::ran1(idum);
                double probrp = rand_gen.drand();
                double probxy = 2.0 * pi * rand_gen.drand();

                double denom2 = probrp * th2iDiff + thMax2i;
                double th = sqrt(1.0 / denom2);

                thx = th * cos(probxy);
                thy = th * sin(probxy);
            }
            theta[0] = thx;
            theta[1] = thy;
        }

        return rcross;
    }

    KOKKOS_INLINE_FUNCTION
    void
    momentum_kick(Kokkos::Random_XorShift64_Pool<>::generator_type& rand_gen,
                  double t,
                  double p,
                  double& dpx,
                  double& dpy)
    {
        double va, vb, va2, vb2, r2 = 10., theta;
        // long idum = -(unsigned)time(0);
        theta = acos(1 - t / (2 * p * p));
        double dp[2] = {0.0, 0.0};

        while (r2 > 1.) {
            // va=2.*Random::ran1(idum)-1;
            // vb=Random::ran1(idum)-1;
            va = 2.0 * rand_gen.drand() - 1;
            vb = rand_gen.drand() - 1;

            va2 = va * va;
            vb2 = vb * vb;
            r2 = va2 + vb2;
        }

        dp[0] = theta * (2. * va * vb) / r2;
        dp[1] = theta * (va2 - vb2) / r2;
        dpx = dp[0];
        dpy = dp[1];
    }

    KOKKOS_INLINE_FUNCTION
    double
    bessj0(double x)
    {
        const double P1 = 1.0, P2 = -0.1098628627E-2, P3 = 0.2734510407E-4,
                     P4 = -0.2073370639E-5, P5 = 0.2093887211E-6,
                     Q1 = -0.1562499995E-1, Q2 = 0.1430488765E-3,
                     Q3 = -0.6911147651E-5, Q4 = 0.7621095161E-6,
                     Q5 = -0.9349451520E-7, R1 = 57568490574.0,
                     R2 = -13362590354.0, R3 = 651619640.7, R4 = -11214424.18,
                     R5 = 77392.33017, R6 = -184.9052456, S1 = 57568490411.0,
                     S2 = 1029532985.0, S3 = 9494680.718, S4 = 59272.64853,
                     S5 = 267.8532712, S6 = 1.0;

        double AX, FR, FS, Z, FP, FQ, XX, Y, TMP;

        if (x == 0.0) return 1.0;

        AX = fabs(x);
        if (AX < 8.0) {
            Y = x * x;
            FR = R1 + Y * (R2 + Y * (R3 + Y * (R4 + Y * (R5 + Y * R6))));
            FS = S1 + Y * (S2 + Y * (S3 + Y * (S4 + Y * (S5 + Y * S6))));
            TMP = FR / FS;
        } else {
            Z = 8. / AX;
            Y = Z * Z;
            XX = AX - 0.785398164;
            FP = P1 + Y * (P2 + Y * (P3 + Y * (P4 + Y * P5)));
            FQ = Q1 + Y * (Q2 + Y * (Q3 + Y * (Q4 + Y * Q5)));
            TMP = sqrt(0.636619772 / AX) * (FP * cos(XX) - Z * FQ * sin(XX));
        }

        return TMP;
    }

    KOKKOS_INLINE_FUNCTION
    double
    BSign(double X, double Y)
    {
        if (Y < 0.0)
            return (-fabs(X));
        else
            return (fabs(X));
    }

    KOKKOS_INLINE_FUNCTION
    double
    bessj1(double X)
    {
        const double P1 = 1.0, P2 = 0.183105E-2, P3 = -0.3516396496E-4,
                     P4 = 0.2457520174E-5, P5 = -0.240337019E-6,
                     P6 = 0.636619772, Q1 = 0.04687499995,
                     Q2 = -0.2002690873E-3, Q3 = 0.8449199096E-5,
                     Q4 = -0.88228987E-6, Q5 = 0.105787412E-6,
                     R1 = 72362614232.0, R2 = -7895059235.0, R3 = 242396853.1,
                     R4 = -2972611.439, R5 = 15704.48260, R6 = -30.16036606,
                     S1 = 144725228442.0, S2 = 2300535178.0, S3 = 18583304.74,
                     S4 = 99447.43394, S5 = 376.9991397, S6 = 1.0;

        double AX, FR, FS, Y, Z, FP, FQ, XX, TMP;

        AX = fabs(X);
        if (AX < 8.0) {
            Y = X * X;
            FR = R1 + Y * (R2 + Y * (R3 + Y * (R4 + Y * (R5 + Y * R6))));
            FS = S1 + Y * (S2 + Y * (S3 + Y * (S4 + Y * (S5 + Y * S6))));
            TMP = X * (FR / FS);
        } else {
            Z = 8.0 / AX;
            Y = Z * Z;
            XX = AX - 2.35619491;
            FP = P1 + Y * (P2 + Y * (P3 + Y * (P4 + Y * P5)));
            FQ = Q1 + Y * (Q2 + Y * (Q3 + Y * (Q4 + Y * Q5)));
            TMP = sqrt(P6 / AX) * (cos(XX) * FP - Z * sin(XX) * FQ) *
                  BSign(S6, X);
        }
        return TMP;
    }

    KOKKOS_INLINE_FUNCTION
    double
    elastic_t(Kokkos::Random_XorShift64_Pool<>::generator_type& rand_gen,
              double p,
              double a)
    {
        double c = pconstants::c;
        double PI = Kokkos::numbers::pi_v<double>;

        int found = 0;
        double R, lamda, E, u, sig, cnorm, theta;
        double random, f_x, x = 0., angle_cm = 0.;
        double M_nuc, E_cm, p_cm, m = .938, h = 4.135e-15 / 1e9;
        double t = 0.;

        M_nuc = a * .931494;
        E = sqrt(p * p + m * m);
        E_cm = sqrt(m * m + M_nuc * M_nuc + 2 * E * M_nuc); // See PDG.
        p_cm = p * M_nuc / E_cm;
        lamda = h * c / p_cm;
        R = 1.4e-15 * pow(a, 1. / 3.) + lamda;
        u = 2 * R / lamda * sin(PI / 2);
        theta = 0.0001;

        cnorm = 1. / 2. *
                (sqrt(PI) - sqrt(PI) * pow(bessj0(u), 2) -
                 sqrt(PI) * pow(bessj1(u), 2)) /
                (sqrt(PI));

        // random = Random::ran1(idum);
        random = rand_gen.drand();

        while (theta <= PI) {
            x = 2. * R * sin(theta / 2.) / lamda;
            f_x = 1. / 2. / cnorm *
                      (sqrt(PI) - sqrt(PI) * pow(bessj0(x), 2) -
                       sqrt(PI) * pow(bessj1(x), 2)) /
                      (sqrt(PI)) -
                  random;
            if ((f_x > -.01) && (f_x < .01)) {
                angle_cm = theta;
                theta = PI + 1.;
                found = 1;
            }
            theta += 1.768e-3;
        }

        // if(found==0) cout<<"Warning, never found elastic t.\n";
        t = 2. * p_cm * p_cm * (1. - cos(angle_cm));
        return t;
    }

    template <class BP>
    struct PropFoilFullScatter {
        const double nAvogadro = 6.022045e23;

        Kokkos::Random_XorShift64_Pool<> rand_pool;

        typename BP::parts_t parts;
        typename BP::masks_t masks;

        int ma_;

        double xmin;
        double xmax;
        double ymin;
        double ymax;

        double pref;
        double m;

        double z;
        double a;
        double density;
        double b_pN;

        double length;
        double dlength;

        // karray1d_dev stat;

        KOKKOS_INLINE_FUNCTION
        PropFoilFullScatter(int seed,
                            typename BP::parts_t parts,
                            typename BP::masks_t masks,
                            double xmin,
                            double xmax,
                            double ymin,
                            double ymax,
                            double thick,
                            double pref,
                            double m)
            : rand_pool(seed)
            , parts(parts)
            , masks(masks)
            , ma_(0)
            , xmin(xmin)
            , xmax(xmax)
            , ymin(ymin)
            , ymax(ymax)
            , pref(pref)
            , m(m)
            , z(get_z(/*ma_*/))
            , a(get_a(/*ma_*/))
            , density(get_rho(/*ma_*/))
            , b_pN(14.5 * pow(a, 2.0 / 3.0))
            , length(thick / (1.0e5 * density))
            , dlength(length * 1.0e-4)
        //, stat("stat", 3)
        {}

        KOKKOS_INLINE_FUNCTION
        bool
        check_foil_flag(double x, double y) const
        {
            return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i) == 0) return;

            auto x = parts(i, 0);
            auto xp = parts(i, 1);
            auto y = parts(i, 2);
            auto yp = parts(i, 3);
            auto dpop = parts(i, 5);

            int step = 0;
            double zrl = length;
            bool foil_flag = check_foil_flag(x, y);

            if (!foil_flag) {
                // TODO: drift();
                return;
            }

#if 0
            std::cout << "hit foil\n";
            std::cout << "zrl = " << zrl << "\n";
#endif

            // random number
            Kokkos::Random_XorShift64_Pool<>::generator_type rand_gen =
                rand_pool.get_state();

            while (zrl > 0) {
                // nHits++;
                double directionfac = get_direction(xp, yp, dpop);
                double rl = zrl * directionfac;

                double beta = get_beta(dpop, pref, m);

                double pfac = 1.0 + dpop;
                double p = (dpop + 1.0) * pref;
                double E = sqrt(p * p + m * m);

                double ecross = get_elastic_crosssection(E);
                double icross = get_inelastic_crosssection(E);

                double meanfreepath = 0.0;
                double stepsize = 0.0;

                if (step == 0) {
                    // If first step, do an iteration with ecross and icross to
                    // get first stepsize and first rcross
                    step++;
                    double totcross = icross + ecross;
                    meanfreepath = a / ((nAvogadro * 1e3) * density *
                                        (totcross * 1.0e-28));
                    // stepsize = -meanfreepath * log(Random::ran1(idum));
                    stepsize = -meanfreepath * log(rand_gen.drand());

#if 0
                    std::cout << "step 0\n";
                    std::cout << "  meanfreepath = " << meanfreepath 
                        << ", stepsize = " << stepsize << "\n";
#endif
                }

                double thetax = 0.0;
                double thetay = 0.0;

                double rcross = ruth_scatt_jackson(rand_gen,
                                                   stepsize,
                                                   z,
                                                   a,
                                                   density,
                                                   beta,
                                                   false,
                                                   pfac,
                                                   thetax,
                                                   thetay);

                double totcross = ecross + icross + rcross;

#if 0
                std::cout << "E = " << E 
                    << ", ecross = " << ecross
                    << ", icross = " << icross
                    << ", rcross = " << rcross
                    << ", rl = " << rl
                    << "\n";
#endif

                meanfreepath =
                    a / ((nAvogadro * 1e3) * density * (totcross * 1.0e-28));
                // stepsize = -meanfreepath * log(Random::ran1(idum));
                stepsize = -meanfreepath * log(rand_gen.drand());

#if 0
                std::cout << "  meanfreepath = " << meanfreepath 
                    << ", stepsize = " << stepsize << "\n";
#endif

                if (stepsize > rl) {
                    // Take the step but no nuclear scattering event
                    stepsize = rl + dlength;

                    if (take_step(rand_gen,
                                  x,
                                  xp,
                                  y,
                                  yp,
                                  dpop,
                                  pref,
                                  m,
                                  z,
                                  a,
                                  density,
                                  stepsize)) {
                        foil_flag = false;
                        zrl = -1.0;

                        // 0xff is marked locally when a particle is
                        // lost during the foil scattering
                        masks(i) = 0xff;
                    } else {
                        double directionfac = get_direction(xp, yp, dpop);
                        zrl -= stepsize / directionfac;
                        rl = zrl * directionfac;
                        foil_flag = check_foil_flag(x, y);
                    }
                }

                if (stepsize <= rl) {
                    // Kokkos::atomic_increment(&stat(0));

                    // Take the step and allow nuclear scatter
                    if (take_step(rand_gen,
                                  x,
                                  xp,
                                  y,
                                  yp,
                                  dpop,
                                  pref,
                                  m,
                                  z,
                                  a,
                                  density,
                                  stepsize)) {
                        foil_flag = false;
                        zrl = -1.0;
                        masks(i) = 0xff;
                    } else {
                        double directionfac = get_direction(xp, yp, dpop);
                        zrl -= stepsize / directionfac;
                        rl = zrl * directionfac;
                        foil_flag = check_foil_flag(x, y);
                    }

                    // If it still exists after MCS and energy loss, nuclear
                    // scatter
                    if (foil_flag == 1 && zrl > 0) {
                        beta = get_beta(dpop, pref, m);

                        pfac = 1.0 + dpop;

                        p = (dpop + 1.0) * pref;
                        E = sqrt(p * p + m * m);

                        ecross = get_elastic_crosssection(E);
                        icross = get_inelastic_crosssection(E);

                        double thx;
                        double thy;

                        rcross = ruth_scatt_jackson(rand_gen,
                                                    stepsize,
                                                    z,
                                                    a,
                                                    density,
                                                    beta,
                                                    false,
                                                    pfac,
                                                    thx,
                                                    thy);

                        totcross = ecross + icross + rcross;

                        double e_frac = ecross / totcross;
                        double i_frac = icross / totcross;
                        double r_frac = rcross / totcross;

#if 0
                        std::cout << "efrac = " << e_frac
                            << ", ifrac = " << e_frac
                            << ", rfrac = " << r_frac
                            << "\n";
#endif

                        // choice = Random::ran1(idum);
                        double choice = rand_gen.drand();

                        // Nuclear Elastic Scattering
                        if ((choice >= 0.) && (choice <= e_frac)) {
                            // Kokkos::atomic_increment(&stat(1));
                            // std::cout << "elastic\n";

                            double t = 0.0;

                            if (E <= 0.4) {
                                t = elastic_t(rand_gen, p, a);
                            } else {
                                // t=-log(Random::ran1(idum))/b_pN;
                                t = -log(rand_gen.drand()) / b_pN;
                            }

                            double dp_x;
                            double dp_y;

                            momentum_kick(rand_gen, t, p, dp_x, dp_y);

                            parts(i, 1) += dp_x * pfac;
                            parts(i, 3) += dp_y * pfac;
                        }

                        // Rutherford Coulomb scattering
                        if ((choice > e_frac) && (choice <= (1 - i_frac))) {
                            // Kokkos::atomic_increment(&stat(2));
                            // std::cout << "coulomb\n";

                            rcross = ruth_scatt_jackson(rand_gen,
                                                        stepsize,
                                                        z,
                                                        a,
                                                        density,
                                                        beta,
                                                        true,
                                                        pfac,
                                                        thx,
                                                        thy);

                            double xpfac = xp / pfac;
                            double ypfac = yp / pfac;

                            double anglex = atan(xpfac) + thx;
                            double angley = atan(ypfac) + thy;

                            parts(i, 1) = tan(anglex) * pfac;
                            parts(i, 3) = tan(angley) * pfac;
                        }

                        // Nuclear Inelastic absorption
                        if ((choice > (1. - i_frac)) && (choice <= 1.)) {
                            // Kokkos::atomic_increment(&stat(3));
                            // std::cout << "absorb\n";

                            foil_flag = false;
                            zrl = -1.0;
                            masks(i) = 0xff;
                        }
                    }
                }

            } // end of while(zrl>0)

            rand_pool.free_state(rand_gen);

        } // end of operator()
    };

    template <class BP>
    struct PropFoilSimpleScatter {
        Kokkos::Random_XorShift64_Pool<> rand_pool;

        typename BP::parts_t parts;
        typename BP::masks_t masks;

        const int ma_;

        double xmin;
        double xmax;
        double ymin;
        double ymax;

        double length;
        double theta_scat_min;
        double lscatter;

        // karray1d_dev stat;

        KOKKOS_INLINE_FUNCTION
        PropFoilSimpleScatter(int seed,
                              typename BP::parts_t parts,
                              typename BP::masks_t masks,
                              double xmin,
                              double xmax,
                              double ymin,
                              double ymax,
                              double thick,
                              double theta_scat_min,
                              double lscatter)
            : rand_pool(seed)
            , parts(parts)
            , masks(masks)
            , ma_(0)
            , xmin(xmin)
            , xmax(xmax)
            , ymin(ymin)
            , ymax(ymax)
            , length(thick / (1.0e3 * get_rho(/*ma_*/)))
            , theta_scat_min(theta_scat_min)
            , lscatter(lscatter)
        //, stat("stat", 101)
        {}

        KOKKOS_INLINE_FUNCTION
        bool
        check_foil_flag(double x, double y) const
        {
            return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i) == 0) return;

            auto x = parts(i, 0);
            auto y = parts(i, 2);

            bool foil_flag = check_foil_flag(x, y);

            if (!foil_flag) return;

            // If in the foil, tally the hit and start tracking
            double zrl = length; // distance remaining in foil in cm//
            double thetaX = 0.;
            double thetaY = 0.;

            // random number
            Kokkos::Random_XorShift64_Pool<>::generator_type rand_gen =
                rand_pool.get_state();

            // Generate interaction points until particle exits foil
            while (zrl >= 0.0) {
                // random1 = Random::ran1(idum);
                double random1 = rand_gen.drand();

                zrl += lscatter * log(random1);
                if (zrl < 0.0) break; // exit foil

                // Generate random angles
                // random1 = Random::ran1(idum);
                random1 = rand_gen.drand();
                double phi = 2 * Kokkos::numbers::pi_v<double> * random1;

                // random1 = Random::ran1(idum);
                random1 = rand_gen.drand();
                double theta = theta_scat_min * sqrt(random1 / (1. - random1));

                thetaX += theta * cos(phi);
                thetaY += theta * sin(phi);

#if 0
                //std::cout << theta*theta << "\n";
                Kokkos::atomic_add(&stat(100), theta*theta);

                int bin = (int)(theta/0.1e-5);
                if (bin<100) Kokkos::atomic_increment(&stat(bin));
#endif
            }

            parts(i, 1) += thetaX;
            parts(i, 3) += thetaY;

            rand_pool.free_state(rand_gen);
        }
    };

    struct Foil_aperture {
        constexpr static const char* type = "foil";

        Foil_aperture() {}

        KOKKOS_INLINE_FUNCTION
        bool
        discard(ConstParticles const& parts,
                ConstParticleMasks const& masks,
                int p) const
        {
            return masks(p) == 0xff;
        }
    };

}

namespace FF_foil {
    template <class BunchT>
    inline void
    apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        // Do notthing for trigons
        if constexpr (is_trigon<typename BunchT::part_t>::value) {
            return;
        } else {
            using namespace foil_impl;

            scoped_simple_timer timer("libFF_foil");

            // const double  length = slice.get_right() - slice.get_left();
            const double mass = bunch.get_mass();

            // element
            auto const& ele = slice.get_lattice_element();

            const double xmin = ele.get_double_attribute("xmin", 0.0);
            const double xmax = ele.get_double_attribute("xmax", 0.0);
            const double ymin = ele.get_double_attribute("ymin", 0.0);
            const double ymax = ele.get_double_attribute("ymax", 0.0);

            const double thick = ele.get_double_attribute("thick", 0.0);

            // use simple scatter or full scatter
            bool is_simple = false;
            const double simple = ele.get_double_attribute("simple", 0);
            if (simple != 0) is_simple = true;

            // random seed
            const int seed = 123;

            // ref
            Reference_particle const& ref_b = bunch.get_reference_particle();
            const double pref =
                ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

            // propagate the bunch particles
            auto apply = [&](ParticleGroup pg) {
                auto bp = bunch.get_bunch_particles(pg);
                if (!bp.num_valid()) return;

                using exec = typename BunchT::exec_space;
                auto range = Kokkos::RangePolicy<exec>(0, bp.size());

                if (is_simple) {
                    double BohrRadius =
                        0.52917706e-8;           // hydrogenic Bohr radius in cm
                    double hBar = 1.0545887e-27; // Planck's constant in erg-sec
                    double echarge = 4.803242e-10; // in esu or statcoulombs
                    double nAvogadro = 6.022045e23;
                    double deg2Rad = 1.74532925e-2;
                    double rhofoil = 2.265;
                    double muScatter = 1.35;

                    double length = thick / (1.0e3 * get_rho(/*ma_*/));

                    // Momentum in g*cm/sec
                    double pInj0 = 1.6726e-22 * mass / pconstants::proton_mass *
                                   ref_b.get_beta() * ref_b.get_gamma() *
                                   pconstants::c;

                    // Thomas-Fermi atom radius (cm):
                    double TFRadius =
                        muScatter * BohrRadius * pow(get_z(/*ma_*/), -0.33333);

                    // Minimum scattering angle:
                    double thetaScatMin = hBar / (pInj0 * TFRadius);

                    // Theta max as per Jackson (13.102)
                    double thetaScatMax =
                        274.e5 * pconstants::electron_mass * pconstants::c /
                        (pInj0 * pow(get_a(/*ma_*/), 0.33333));

                    double pv = 1.e2 * pInj0 * ref_b.get_beta() * pconstants::c;
                    double term = get_z(/*ma_*/) * echarge * echarge / pv;
                    double sigmacoul = 4 * Kokkos::numbers::pi_v<double> *
                                       term * term /
                                       (thetaScatMin * thetaScatMin);

                    // Scattering area per area
                    double nscatters = nAvogadro * (get_rho(/*ma_*/) / 1000.0) /
                                       get_a(/*ma_*/) * length * sigmacoul;

                    // Mean free path
                    double lscatter = length / nscatters;

                    PropFoilSimpleScatter<typename BunchT::bp_t> psc(
                        seed,
                        bp.parts,
                        bp.masks,
                        xmin,
                        xmax,
                        ymin,
                        ymax,
                        thick, // geometry
                        lscatter,
                        thetaScatMin);

                    Kokkos::parallel_for(range, psc);
                } else {
                    PropFoilFullScatter<typename BunchT::bp_t> pfc(
                        seed,
                        bp.parts,
                        bp.masks,
                        xmin,
                        xmax,
                        ymin,
                        ymax,
                        thick, // geometry
                        pref,
                        mass);

                    Kokkos::parallel_for(range, pfc);
                }
            };

            // apply
            apply(ParticleGroup::regular);
            apply(ParticleGroup::spectator);

            // advance the ref_part
            // bunch.get_reference_particle().increment_trajectory(length);

            if (!is_simple) {
                // aperture to fiter out all particles with mask value 0xff
                // seems only happen in the full scatter
                Foil_aperture ap;
                int ndiscarded = bunch.apply_aperture(ap);
                double charge =
                    ndiscarded * bunch.get_real_num() / bunch.get_total_num();
                slice.get_lattice_element().deposit_charge(charge);
            }

            // fencing
            Kokkos::fence();
        }
    }
}

#endif // FF_FOIL_H
