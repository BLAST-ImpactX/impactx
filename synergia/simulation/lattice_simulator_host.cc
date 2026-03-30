#include "lattice_simulator_host.h"
#include <Kokkos_MathematicalConstants.hpp>

namespace Lattice_simulator {

    using Vector = Eigen::Matrix<std::complex<double>, 6, 1>;
    using Matrix = Eigen::Matrix<std::complex<double>, 6, 6, Eigen::RowMajor>;

    Vector ev_ordering(Vector const& ev, Matrix const& B);

    std::array<double, 2>
    filter_transverse_tunes(double const* jac_arr)
    {
        constexpr static int EigenIterations = 100000;

        using MatrixD = Eigen::
            Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        std::array<double, 2> nu;
        Eigen::Matrix<double, 6, 6, Eigen::RowMajor> jac(jac_arr);

        const int ix = 0;
        const int ipx = 1;
        const int iy = 2;
        const int ipy = 3;

        if (jac(iy, ix) || jac(ix, iy) || jac(ix, ipy) || jac(iy, ipx) ||
            jac(ipy, ix) || jac(ipx, iy) || jac(ipy, ipx) || jac(ipx, ipy)) {
            Eigen::EigenSolver<MatrixD> eigensolver;

            eigensolver.setMaxIterations(EigenIterations);
            eigensolver.compute(jac);

            if (eigensolver.info() == Eigen::NoConvergence)
                throw std::runtime_error("eigensolver no convergence");

            if (eigensolver.info() != Eigen::Success)
                throw std::runtime_error("failed solving eigenvectors");

            auto eigen_val = eigensolver.eigenvalues();
            auto eigen_vec = eigensolver.eigenvectors();

            auto lambda = ev_ordering(eigen_val, eigen_vec);

            for (int i = 0; i < 6; ++i) {
                if (fabs(abs(lambda(i)) - 1.0) > 1.0e-4) {
                    std::stringstream ss;
                    ss << "filterTransverseTunes: "
                       << "the lattice is nearly unstable. "
                       << "lambda( " << i
                       << " ) has magnitude = " << abs(lambda(i));
                    throw std::runtime_error(ss.str());
                }
            }

            if ((abs(lambda(0) - std::conj(lambda(3))) > 1.0e-4) ||
                (abs(lambda(1) - std::conj(lambda(4))) > 1.0e-4)) {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "conjugacy condition has been vilated. "
                   << "The lattice may be linearly unstable. "
                   << "Eigenvalues = " << lambda;
                throw std::runtime_error(ss.str());
            }

            double csH = lambda(0).real();
            double csV = lambda(1).real();

            if (fabs(csH - csV) < 1.0e-4) {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "\"Horizontal\" and \"vertical\" tunes "
                   << "are too near each other for reasonable results. "
                   << "The calculation is meaningless.";
                throw std::runtime_error(ss.str());
            }

            double dcos, cos2phi, sin2phi, tanphi;

            MatrixD U(2, 2), S(2, 2);

            U << 1.0, 0.0, 0.0, 1.0;
            S << 0.0, 1.0, -1.0, 0.0;

            MatrixD M(2, 2), N(2, 2);
            MatrixD m(2, 2), n(2, 2);

            M << jac(ix, ix), jac(ix, ipx), jac(ipx, ix), jac(ipx, ipx);
            N << jac(iy, iy), jac(iy, ipy), jac(ipy, iy), jac(ipy, ipy);
            m << jac(iy, ix), jac(iy, ipx), jac(ipy, ix), jac(ipy, ipx);
            n << jac(ix, iy), jac(ix, ipy), jac(ipx, iy), jac(ipx, ipy);

#if 0
        std::cout << "M = " << M << "\n";
        std::cout << "N = " << N << "\n";
        std::cout << "m = " << m << "\n";
        std::cout << "n = " << n << "\n";
#endif

            dcos = csH - csV;
            cos2phi = (M - N).trace() / (2.0 * (dcos));

            if (fabs(cos2phi - 1.0) < 1.0e-4)
                cos2phi = 1.0; // ??? Rather coarse,

            if (fabs(cos2phi + 1.0) < 1.0e-4) cos2phi = -1.0; // ??? isn't it?

            if (fabs(cos2phi) > 1.0) {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "cos(2 phi) = " << std::setprecision(10) << cos2phi
                   << "; has magnitude larger than one. "
                   << "Cannot continue calculation. ";
                throw std::runtime_error(ss.str());
            }

            if (cos2phi < 0.0) {
                sin2phi = csH; // Variable used as dummy register.
                csH = csV;
                csV = sin2phi;
                dcos = -dcos;
                cos2phi = -cos2phi;
            }

            sin2phi = sqrt(1.0 - cos2phi * cos2phi);
            tanphi = sin2phi / (1.0 + cos2phi);

#if 0
        std::cout << "sin2phi = " << sin2phi << "\n";
        std::cout << "tanphi = " << tanphi << "\n";
#endif

            MatrixD D(2, 2), A(2, 2), B(2, 2);

            if (fabs(sin2phi) > 1.0e-8) {
                D = -(m + S * n.transpose() * S.transpose()) *
                    (1.0 / (dcos * sin2phi));
            } else {
                D << 1.0, 0.0, 0.0, 1.0;
            }

            if (fabs(D.determinant() - 1.0) > 1.0e-4) {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "The matrix D is non-symplectic. "
                   << "|D| = " << D.determinant();
                throw std::runtime_error(ss.str());
            }

            // ...... Edwards-Teng sign convention.
            if (D.trace() < 0.0) {
                D = -D;
                sin2phi = -sin2phi;
                tanphi = -tanphi;
            }

            A = M - D.inverse() * m * tanphi;
            B = N + D * n * tanphi;

#if 0
        std::cout << "A = " << A << "\n";
        std::cout << "B = " << B << "\n";
#endif

            // ......  First the "horizontal" ......
            MatrixD JH = A - csH * U;
            double snH = (JH(0, 1) > 0.0) ? sqrt(1.0 - csH * csH) :
                                            -sqrt(1.0 - csH * csH);

            // .......... A little test to keep everyone honest .....
            if (JH(0, 0)) {
                if (fabs((JH(0, 0) + JH(1, 1)) / (JH(0, 0) - JH(1, 1))) >
                    1.0e-4) {
                    std::cout
                        << "WARNING -- filterTransverseTunes: "
                        << "\"Horizontal\" matrix does not "
                        << "pass symplecticity test. "
                        << "JH( 0, 0 ) = " << JH(0, 0) << ", "
                        << "JH( 1, 1 ) = " << JH(1, 1) << ". "
                        << "The ratio is "
                        << fabs((JH(0, 0) + JH(1, 1)) / (JH(0, 0) - JH(1, 1)))
                        << "\n";
                }
            }

            // ......  Then  the "vertical" ......
            MatrixD JV = B - csV * U;
            double snV = (JV(0, 1) > 0.0) ? sqrt(1.0 - csV * csV) :
                                            -sqrt(1.0 - csV * csV);

            // .......... A little test to keep everyone honest .....
            if (JV(0, 0)) {
                if (fabs((JV(0, 0) + JV(1, 1)) / (JV(0, 0) - JV(1, 1))) >
                    1.0e-4) {
                    std::cout
                        << "WARNING -- filterTransverseTunes: "
                        << "\"Vertical\" matrix does not "
                        << "pass symplecticity test. "
                        << "JV( 0, 0 ) = " << JV(0, 0) << ", "
                        << "JV( 1, 1 ) = " << JV(1, 1) << ". "
                        << "The ratio is "
                        << fabs((JV(0, 0) + JV(1, 1)) / (JV(0, 0) - JV(1, 1)))
                        << "\n";
                }
            }

            constexpr double M_TWOPI = Kokkos::numbers::pi_v<double> * 2;

            double theta = atan2(snH, csH);
            if (theta < 0.0) theta += M_TWOPI;
            nu[0] = theta / M_TWOPI;

            theta = atan2(snV, csV);
            if (theta < 0.0) theta += M_TWOPI;
            nu[1] = theta / M_TWOPI;
        } else {
            double sn, cs;

            // Uncoupled calculation .....
            // (Lifted from LattFuncSage) ...
            // ... first horizontal
            cs = (jac(ix, ix) + jac(ipx, ipx)) / 2.0;

            if (fabs(cs) <= 1.0) {
                if (jac(ix, ipx) > 0.0)
                    sn = sqrt(1.0 - cs * cs);
                else
                    sn = -sqrt(1.0 - cs * cs);
            } else {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "cos( psi_H ) = " << cs << ". "
                   << "Cannot continue with calculation.";
                throw std::runtime_error(ss.str());
            }

            constexpr double M_TWOPI = Kokkos::numbers::pi_v<double> * 2;

            double theta = atan2(sn, cs);
            if (theta < 0.0) theta += M_TWOPI;
            nu[0] = theta / M_TWOPI;

            // ... then vertical.
            cs = (jac(iy, iy) + jac(ipy, ipy)) / 2.0;

            if (fabs(cs) <= 1.0) {
                if (jac(iy, ipy) > 0.0)
                    sn = sqrt(1.0 - cs * cs);
                else
                    sn = -sqrt(1.0 - cs * cs);
            } else {

                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "cos( psi_V ) = " << cs << ". "
                   << "Cannot continue with calculation.";
                throw std::runtime_error(ss.str());
            }

            theta = atan2(sn, cs);
            if (theta < 0.0) theta += M_TWOPI;
            nu[1] = theta / M_TWOPI;
        }

        return nu;
    }

    Vector
    ev_ordering(Vector const& ev, Matrix const& B)
    {
        int Lidx[6];
        int unused[6];

        for (int i = 0; i < 6; ++i) {
            unused[i] = 1;
            Lidx[i] = -1;
        }

        // first reorder to L1, L2, L3, L1*, L2*, L3*
        for (int i = 0; i < 3; ++i) {
            // find first unused
            for (int j = 0; j < 6; ++j) {
                if (unused[j]) {
                    Lidx[i] = j;
                    unused[j] = 0;
                    break;
                }
            }

            // eigenvalue of the selected index
            auto lambda1 = ev(Lidx[i]);

            // find the pairing eigenvalue
            for (int j = 0; j < 6; ++j) {
                if (unused[j]) {
                    auto lambda2 = ev(j);

                    // put the matching index to i+3
                    if (abs(lambda1 - std::conj(lambda2)) < 1e-6) {
                        Lidx[i + 3] = j;
                        unused[j] = 0;
                        break;
                    }
                }
            }

            // do we find the matching one?
            if (Lidx[i + 3] < 0)
                throw std::runtime_error("Failed to find matching eigenvalue");
        }

        // next find the largest x component and move to column 0 and 3
        double max = 0.0;
        int idx = 0;

        for (int i = 0; i < 3; ++i) {
            double x = abs(B(0, Lidx[i]));

            if (x > max) {
                max = x;
                idx = i;
            }
        }

        if (idx != 0) {
            int tmp = Lidx[0];
            Lidx[0] = Lidx[idx];
            Lidx[idx] = tmp;

            tmp = Lidx[3];
            Lidx[3] = Lidx[idx + 3];
            Lidx[idx + 3] = tmp;
        }

        // then find the largest y component and move to column 1 and 4
        max = 0.0;

        for (int i = 1; i < 3; ++i) {
            double y = abs(B(1, Lidx[i]));

            if (y > max) {
                max = y;
                idx = i;
            }
        }

        if (idx != 1) {
            int tmp = Lidx[1];
            Lidx[1] = Lidx[idx];
            Lidx[idx] = tmp;

            tmp = Lidx[4];
            Lidx[4] = Lidx[idx + 3];
            Lidx[idx + 3] = tmp;
        }

        Vector ev2;
        for (int i = 0; i < 6; ++i)
            ev2(i) = ev(Lidx[i]);

        return ev2;
    }

}
