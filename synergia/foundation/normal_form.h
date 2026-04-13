#ifndef FOUNDATION_NORMAL_FORM_H
#define FOUNDATION_NORMAL_FORM_H

#include "Eigen/Eigen"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/trigon.h"

#include <cereal/cereal.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/vector.hpp>

template <unsigned int order>
class NormalForm {
    // max iterations for the eigen solver
    constexpr static int EigenIterations = 100000;

    // matrix types
    using Vector6D = Eigen::Matrix<double, 6, 1>;
    using Vector6C = Eigen::Matrix<std::complex<double>, 6, 1>;

    using Matrix6D = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>;

    using Matrix6C = Eigen::Matrix<std::complex<double>, 6, 6, Eigen::RowMajor>;

    using Matrix3D = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

    using Matrix3C = Eigen::Matrix<std::complex<double>, 3, 3, Eigen::RowMajor>;

    using MatrixD =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using MatrixC = Eigen::Matrix<std::complex<double>,
                                  Eigen::Dynamic,
                                  Eigen::Dynamic,
                                  Eigen::RowMajor>;

  public:
    // trigon types
    using trigon_t = Trigon<double, order, 6>;
    using mapping_t = TMapping<trigon_t>;

    using trigon_c_t = Trigon<std::complex<double>, order, 6>;
    using mapping_c_t = TMapping<trigon_c_t>;

    using operators_t = std::vector<mapping_c_t>;

    // e0: reference energy
    // pc0: reference momentum
    // mass: reference particle mass
    NormalForm(mapping_t const& one_turn_map,
               double e0,
               double pc0,
               double mass);

    // default constructor for serialization only
    NormalForm() {}

    std::array<double, 3> stationaryActions(double stdx,
                                            double stdy,
                                            double stdz) const;

    std::array<std::complex<double>, 3> cnvDataToNormalForm(
        std::array<double, 6> const& hform) const;

    std::array<double, 6> cnvDataFromNormalForm(
        std::array<std::complex<double>, 3> const& nform) const;

    std::vector<mapping_c_t>&
    get_f()
    {
        return f_;
    }
    std::vector<mapping_c_t> const&
    get_f() const
    {
        return f_;
    }

    std::vector<mapping_c_t>&
    get_g()
    {
        return g_;
    }
    std::vector<mapping_c_t> const&
    get_g() const
    {
        return g_;
    }

  private:
    Matrix6C ev_ordering(Vector6C const& ev, Matrix6C const& B) const;

  private:
    constexpr static const double MLT1 = 1.0e-5;

    mapping_t CanonToSyn;
    mapping_t SynToCanon;

    Matrix6C E_;
    Matrix6C invE_;

    std::vector<mapping_c_t> f_;
    std::vector<mapping_c_t> g_;

  private:
    // serialization
    friend class cereal::access;

    template <class AR>
    void
    serialize(AR& ar)
    {
        ar(CanonToSyn);
        ar(SynToCanon);
        ar(f_);
        ar(g_);

        if (AR::is_saving::value) {
            auto size = 6 * 6 * sizeof(std::complex<double>);

            std::array<std::complex<double>, 6 * 6> ae;
            std::array<std::complex<double>, 6 * 6> ainve;

            std::memcpy(ae.data(), E_.data(), size);
            std::memcpy(ainve.data(), invE_.data(), size);

            ar(ae);
            ar(ainve);
        } else {
            auto size = 6 * 6 * sizeof(std::complex<double>);

            std::array<std::complex<double>, 6 * 6> ae;
            std::array<std::complex<double>, 6 * 6> ainve;

            ar(ae);
            ar(ainve);

            std::memcpy(E_.data(), ae.data(), size);
            std::memcpy(invE_.data(), ainve.data(), size);
        }
    }
};

template <unsigned int order>
NormalForm<order>::NormalForm(mapping_t const& one_turn_map,
                              double e0,
                              double pc0,
                              double mass)
    : E_(Matrix6C::Zero())
    , invE_(Matrix6C::Zero())
    , f_(order - 1)
    , g_(order - 1)
{

    KOKKOS_IF_ON_DEVICE((
        // empty implementation for cuda as this is only
        // supposed to run on the host
        ))
    KOKKOS_IF_ON_HOST((

        constexpr const int dim = trigon_t::dim;

        // chef index
        int c_ix = 0;
        int c_iy = 1;
        int c_it = 2;
        int c_ipx = 3;
        int c_ipy = 4;
        int c_ide = 5;

        // synergia index
        int s_ix = 0;
        int s_ipx = 1;
        int s_iy = 2;
        int s_ipy = 3;
        int s_idt = 4;
        int s_idp = 5;

        // coordinates for synergia particle
        trigon_t x(0.0, s_ix);
        trigon_t px(0.0, s_ipx);
        trigon_t y(0.0, s_iy);
        trigon_t py(0.0, s_ipy);
        trigon_t cdt(0.0, s_idt);
        trigon_t dpop(0.0, s_idp);

        // coordinates for canonical particle
        trigon_t cx(0.0, c_ix);
        trigon_t cy(0.0, c_iy);
        trigon_t ct(0.0, c_it);
        trigon_t cpx(0.0, c_ipx);
        trigon_t cpy(0.0, c_ipy);
        trigon_t cde(0.0, c_ide);

        // synergia coordinates to cannonical
        // cdt -> -dt
        // dp/p -> deltaE
        trigon_t realP = (1.0 + dpop) * pc0;
        trigon_t deltaE = sqrt(realP * realP + mass * mass) - e0;

        SynToCanon[c_ix] = x;
        SynToCanon[c_iy] = y;
        SynToCanon[c_it] = -cdt / pconstants::c;
        SynToCanon[c_ipx] = px * pc0 / pconstants::c;
        SynToCanon[c_ipy] = py * pc0 / pconstants::c;
        SynToCanon[c_ide] = deltaE;

        // canonical to synergia
        trigon_t realE = cde + e0;
        trigon_t dp = sqrt(realE * realE - mass * mass) - pc0;

        // std::cout << "realE = " << realE;
        // std::cout << "dp = " << dp;

        CanonToSyn[s_ix] = cx;
        CanonToSyn[s_ipx] = cpx * pconstants::c / pc0;
        CanonToSyn[s_iy] = cy;
        CanonToSyn[s_ipy] = cpy * pconstants::c / pc0;
        CanonToSyn[s_idt] = -ct * pconstants::c;
        CanonToSyn[s_idp] = dp / pc0;

        // set the constant part to 0
        mapping_t M = one_turn_map;
        for (int i = 0; i < mapping_t::dim; ++i) M[i].value() = 0.0;

        // The combined transformation that we will use for the normal form
        // analysis is the one turn map of a canonical particle.  To get this,
        // apply the maps that turns a canonical particle to a synergia
        // particle, one turn map of a synergia particle, syn particle to
        // canonical particle.
        mapping_t canonMap = SynToCanon(M(CanonToSyn));

        // now the normal form
        const std::complex<double> complex_0(0.0, 0.0);
        const std::complex<double> complex_1(1.0, 0.0);
        const std::complex<double> mi(0.0, -1.0);

        // establising linear normal form coordinates
        auto kjac = canonMap.jacobian();
        Matrix6D A(kjac.data());

        Eigen::EigenSolver<Matrix6D> eigensolver;

        eigensolver.setMaxIterations(EigenIterations);
        eigensolver.compute(A);

        if (eigensolver.info() == Eigen::NoConvergence) throw std::
            runtime_error("eigensolver no convergence");

        if (eigensolver.info() != Eigen::Success) throw std::runtime_error(
            "failed solving eigenvectors");

        auto ev = eigensolver.eigenvalues();
        auto B = eigensolver.eigenvectors();

        // normalizing the linear normal form coordinates
        Matrix6D J = Matrix6D::Zero();
        for (int i = dim / 2; i < dim; ++i) {
            J(i - dim / 2, i) = 1.0;
            J(i, i - dim / 2) = -1.0;
        }

        // reordering B
        Matrix6C Br = ev_ordering(ev, B);

        Matrix6C Nx = (Br.transpose() * J * Br * J) * mi;

        for (int i = 0; i < 6; ++i) Nx(i, i) = 1.0 / sqrt(abs(Nx(i, i)));

        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) if (i != j)
            Nx(i, j) = std::complex<double>(0, 0);

        // std::cout << "Nx = \n" << Nx << "\n\n";

        B = Br * Nx;

        // std::cout << "B = B * Nx = \n" << B << "\n\n";

        // try to get the phase correct
        std::complex<double> m0,
        cm0,
        m1,
        cm1,
        m2,
        cm2;
        m0 = B(0, 0) / abs(B(0, 0));
        cm0 = std::conj(m0);
        m1 = B(1, 1) / abs(B(1, 1));
        cm1 = std::conj(m1);
        m2 = B(2, 2) / abs(B(2, 2));
        cm2 = std::conj(m2);

        for (int i = 0; i < 6; ++i) {
            B(i, 0) *= cm0;
            B(i, 3) *= m0;
            B(i, 1) *= cm1;
            B(i, 4) *= m1;
            B(i, 2) *= cm2;
            B(i, 5) *= m2;
        }

        // NOTE: the variable m0 is reused here and
        // below as a dummy variable. This nullifies
        // its previous interpretation.
        if (imag(B(3, 0)) > 0.0) {
            for (int i = 0; i < 6; ++i) {
                m0 = B(i, 0);
                B(i, 0) = B(i, 3);
                B(i, 3) = m0;
            }
        }

        if (imag(B(4, 1)) > 0.0) {
            for (int i = 0; i < 6; ++i) {
                m0 = B(i, 1);
                B(i, 1) = B(i, 4);
                B(i, 4) = m0;
            }
        }

        if (imag(B(5, 2)) > 0.0) {
            for (int i = 0; i < 6; ++i) {
                m0 = B(i, 2);
                B(i, 2) = B(i, 5);
                B(i, 5) = m0;
            }
        }

        if (imag(B(5, 2)) > 0.0) {
            for (int i = 0; i < 6; ++i) {
                m0 = B(i, 2);
                B(i, 2) = B(i, 5);
                B(i, 5) = m0;
            }
        }

        // std::cout << "B after phase correct = \n" << B << "\n\n";

        E_ = B;

        // std::cout << "E = \n" << E_ << "\n\n";

        invE_ = E_.inverse();

        // std::cout << "invE = \n" << invE_ << "\n\n";

        // some useful matrices
        Matrix6C Binv = B.inverse();
        Matrix6C D = Binv * A * B;
        Matrix6C Dinv = D.inverse();

        constexpr auto dcols = Matrix6C::ColsAtCompileTime;
        std::array<std::complex<double>, dcols> lambda;
        std::array<double, dcols> nu;

        for (int i = 0; i < D.cols(); ++i) {
            lambda[i] = D(i, i);
            nu[i] = -std::arg(lambda[i]) / (Kokkos::numbers::pi_v<double> * 2);
        }

        // the following blocks are marked with "CAUTION" in CHEF
        for (int i = 0; i < 6; i++) {
            if (fabs(abs(lambda[i]) - 1.0) > MLT1) {
                std::ostringstream uic;
                uic << "NormalForm(): "
                    << "Only elliptic fixed points allowed: |lambda( " << i
                    << " )| = " << std::abs(lambda[i]) << " = 1.0 + ( "
                    << (abs(lambda[i]) - 1.0) << " )";
                throw std::runtime_error(uic.str());
            }

            if (fabs(lambda[i].imag()) < MLT1) {
                std::ostringstream uic;
                uic << "NormalForm(): Eigenvalue " << i << " = " << lambda[i]
                    << ": too close to integer or half-integer tune.";
                throw std::runtime_error(uic.str());
            }

            if (fabs(lambda[i].real()) < MLT1) {
                std::ostringstream uic;
                uic << "NormalForm(): Eigenvalue " << i << " = " << lambda[i]
                    << ": too close to quarter-integer tune.";
                throw std::runtime_error(uic.str());
            }
        }

        // A little checking and cleaning.
        for (int i = 0; i < 6; i++) {
            if (fabs(abs(D(i, i)) - 1.0) > MLT1) {
                std::ostringstream uic;
                uic << "NormalForm(): "
                    << "For now, only elliptic maps allowed: | D( " << i << ", "
                    << i << " ) | = " << std::abs(D(i, i)) << " = 1.0 + ( "
                    << (abs(D(i, i)) - 1.0) << " )";
                throw std::runtime_error(uic.str());
            }

            for (int j = 0; j < 6; ++j) {
                if (j == i) continue;

                if (abs(D(i, j)) > MLT1) {
                    std::ostringstream uic;
                    uic << "NormalForm(): "
                        << "An impossible error has occured, | D( " << i << ", "
                        << j << " ) | = " << std::abs(D(i, j)) << " > " << MLT1;
                    throw std::runtime_error(uic.str());
                }

                D(i, j) = complex_0;
            }
        }

        for (int i = 0; i < 6; i++) {
            if (fabs(abs(Dinv(i, i)) - 1.0) > MLT1) {
                std::ostringstream uic;
                uic << "NormalForm(): "
                    << "For now, only elliptic maps allowed: | Dinv( " << i
                    << ", " << i << " ) | = " << std::abs(Dinv(i, i))
                    << " = 1.0 + ( " << (abs(Dinv(i, i)) - 1.0) << " )";
                throw std::runtime_error(uic.str());
            }

            for (int j = 0; j < 6; ++j) {
                if (j == i) continue;

                if (abs(Dinv(i, j)) > MLT1) {
                    std::ostringstream uic;
                    uic << "NormalForm(): "
                        << "An impossible error has occured, | Dinv( " << i
                        << ", " << j << " ) | = " << std::abs(D(i, j)) << " > "
                        << MLT1;
                    throw std::runtime_error(uic.str());
                }

                Dinv(i, j) = complex_0;
            }
        }

        // the original near-identity transformation
        mapping_c_t id;
        for (int i = 0; i < id.dim; ++i) id[i].set(0.0, i);

        mapping_c_t CL1 = static_cast<mapping_c_t>(canonMap);
        mapping_c_t calN = Binv * CL1(B * (Dinv * id));

        for (int i = 0; i < dim; ++i) {
            calN[i].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) < 1e-10) val = 0;
            });
        }

        std::array<mapping_c_t, order> N;
        std::array<mapping_c_t, order> T;

        for (int k = 0; k < order - 1; ++k) {
            mapping_c_t reg = id;
            int ll = 0;

            while (ll < k) {
                reg = N[ll].exp_map(-complex_1, reg);
                ++ll;
            }

            reg = calN(reg);
            reg.filter(k + 2, k + 2);

            N[k] = reg;

            // N[k].filter(shear);
            N[k][0].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) == 0) return;

                arr_t<size_t, dim> exp;
                for (auto i : ind)
                    ++exp[i];

                if (exp[0] != exp[3] + 1)
                    val = 0;
                else if (exp[1] != exp[4])
                    val = 0;
                else if (exp[2] != exp[5])
                    val = 0;
            });

            N[k][1].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) == 0) return;

                arr_t<size_t, dim> exp;
                for (auto i : ind)
                    ++exp[i];

                if (exp[0] != exp[3])
                    val = 0;
                else if (exp[1] != exp[4] + 1)
                    val = 0;
                else if (exp[2] != exp[5])
                    val = 0;
            });

            N[k][2].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) == 0) return;

                arr_t<size_t, dim> exp;
                for (auto i : ind)
                    ++exp[i];

                if (exp[0] != exp[3])
                    val = 0;
                else if (exp[1] != exp[4])
                    val = 0;
                else if (exp[2] != exp[5] + 1)
                    val = 0;
            });

            N[k][3].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) == 0) return;

                arr_t<size_t, dim> exp;
                for (auto i : ind)
                    ++exp[i];

                if (exp[0] != exp[3] - 1)
                    val = 0;
                else if (exp[1] != exp[4])
                    val = 0;
                else if (exp[2] != exp[5])
                    val = 0;
            });

            N[k][4].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) == 0) return;

                arr_t<size_t, dim> exp;
                for (auto i : ind)
                    ++exp[i];

                if (exp[0] != exp[3])
                    val = 0;
                else if (exp[1] != exp[4] - 1)
                    val = 0;
                else if (exp[2] != exp[5])
                    val = 0;
            });

            N[k][5].each_term([](size_t idx, auto const& ind, auto& val) {
                if (abs(val) == 0) return;

                arr_t<size_t, dim> exp;
                for (auto i : ind)
                    ++exp[i];

                if (exp[0] != exp[3])
                    val = 0;
                else if (exp[1] != exp[4])
                    val = 0;
                else if (exp[2] != exp[5] - 1)
                    val = 0;
            });

            mapping_c_t doc = N[k] - reg;

#if 0
        std::cout << "N[" << k << "] = \n" << N[k] << "\n";
        std::cout << "reg = \n" << reg << "\n";
        std::cout << "doc[" << k << "] = \n" << doc << "\n";
#endif

            for (int d = 0; d < dim; ++d) {
                // idx is the index in the terms[] array
                // ind is the indices of the corresponding term
                // val is the coefficient of the term
                doc[d].each_term(
                    [this, k, d, &lambda, &N, &T](
                        size_t idx, auto const& ind, auto const& val) {
                        // do nothing if the term is (0,0)
                        if (!abs(val)) return;

                        const std::complex<double> complex_1(1.0, 0.0);
                        std::complex<double> factor(1.0, 0.0);

                        const int power = ind.size();

                        for (int i = 0; i < power; ++i)
                            factor *= complex_1 / lambda[ind[i]];

                        factor *= lambda[d];
                        auto denom = factor - complex_1;

                        // either absorption or resonance subtraction ...
                        if (abs(denom) < 1e-7) {
                            N[k][d].set_term(power, idx, val);

#if 0
                    // If this is not a shear term, print a warning
                    auto exponents = trigon_t::exponent(z);
                    if (!(shear[i])(exponents, std::complex<double>(0,0)))
                    {
                        // print
                    }
#endif

                        } else {
                            T[k][d].set_term(power, idx, val / denom);
                        }
                    });

            } // loop-dim

            // prepare for the next order
            mapping_c_t mapT;

            reg = Dinv * id;
            // std::cout << "reg = Dinv*id = \n" << reg << "\n";

            mapT = T[k].exp_map(complex_1, id);
            // std::cout << "mapT = \n" << mapT << "\n";

            for (int i = 0; i < 6; ++i) {
                // std::cout << "T[k]^id(" << i << ") = \n" << (T[k] ^ id[i]) <<
                // "\n";
            }

            reg = mapT(reg);
            // std::cout << "mapT(reg) = \n" << reg << "\n";

            reg = D * reg;
            // std::cout << "D*reg = \n" << reg << "\n";

            reg = calN(reg);
            // std::cout << "calN(reg) = \n" << reg << "\n";

            mapT = T[k].exp_map(-complex_1, id);
            // std::cout << "mapT = T[k].exp_map() = \n" << mapT << "\n";

            calN = mapT(reg);
            // std::cout << "calN = \n" << calN << "\n";

            // above in one line:
            // calN = T[k].exp_map( -1.0, calN( D*( T[k].exp_map( 1.0, Dinv*id )
            // ) )
            // );
        } // loop k-order

#if 0
    for(int i=0; i<order; ++i)
    {
        std::cout << "N[" << i << "] = \n" << N[i] << "\n";
        std::cout << "T[" << i << "] = \n" << T[i] << "\n";
        std::cout << "\n\n";
    }
#endif

        for (int i = 0; i < order - 1; ++i) {
            f_[i] = T[i].exp_map(-1.0, id);
            f_[i].filter(0, i + 2);

            g_[i] = T[i].exp_map(1.0, id);
            g_[i].filter(0, i + 2);

            // std::cout << "f[" << i << "] = \n" << f_[i] << "\n";
        }))
}

template <unsigned int order>
std::array<double, 3>
NormalForm<order>::stationaryActions(double stdx,
                                     double stdy,
                                     double stdz) const
{
    MatrixD bmom(3, 3);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            bmom(i, j) = 2.0 * (E_(i, j) * std::conj(E_(i, j))).real();

    // std::cout << "bmom = " << bmom << "\n";

    auto inv_bmom = bmom.inverse();

    // std::cout << "inv_bmom = " << inv_bmom << "\n";

    double stdt = stdz / pconstants::c;
    Eigen::Vector3d moments(stdx * stdx, stdy * stdy, stdt * stdt);

    // std::cout << "moments = " << moments << "\n";

    auto mact = inv_bmom * moments;

    return {mact(0), mact(1), mact(2)};
}

template <unsigned int order>
std::array<std::complex<double>, 3>
NormalForm<order>::cnvDataToNormalForm(std::array<double, 6> const& hform) const
{
    arr_t<double, 6> hf;
    for (int i = 0; i < 6; ++i)
        hf[i] = hform[i];

    auto hsymp = SynToCanon(hf);
    auto hftest = CanonToSyn(hsymp);

    for (int i = 0; i < 6; ++i) {
        // std::cout << hsymp[i] << "\t" << hftest[i] << "\n";
        if (abs(hftest[i] - hf[i]) > 1e-6)
            throw std::runtime_error("cnvDataToNormalForm: CanonToSyn error!");
    }

    Vector6C u;
    for (int i = 0; i < 6; ++i)
        u[i] = hsymp[i];
    u = invE_ * u;

    // std::cout << "u = \n" << u << "\n";

    arr_t<std::complex<double>, 6> a;

    for (int i = 0; i < order - 1; ++i) {
        for (int j = 0; j < 6; ++j)
            a[j] = u(j);
        for (int j = 0; j < 6; ++j)
            u(j) = (f_[i][j])(a);
    }

    std::array<std::complex<double>, 3> nform;
    for (int i = 0; i < 3; ++i)
        nform[i] = u[i];

    return nform;
}

template <unsigned int order>
std::array<double, 6>
NormalForm<order>::cnvDataFromNormalForm(
    std::array<std::complex<double>, 3> const& nform) const
{
    const int THROTTLE = 10;
    const int NAG = 1000;

    int nwarnings = 0;
    int nagthreshold = NAG;
    int warningsdecade = 0;
    bool throttled = false;

    Vector6C u;

    for (int i = 0; i < 3; ++i) {
        u(i) = nform[i];
        u(i + 3) = std::conj(nform[i]);
    }

    arr_t<std::complex<double>, 6> a;

    for (int i = order - 2; i >= 0; --i) {
        for (int j = 0; j < 6; ++j)
            a[j] = u(j);
        for (int j = 0; j < 6; ++j)
            u(j) = (g_[i][j])(a);
    }

    u = E_ * u;

    arr_t<double, 6> hsymp;
    const double small_thresh = 5.0e-11;

    for (int i = 0; i < 6; ++i) {
        // imaginary part of u(i) should be small.
        if (std::abs(real(u(i))) > small_thresh) {
            // if the real part is non-zero, the imaginary part should
            // be a small fraction of it.
            if (std::abs(imag(u(i)) / real(u(i))) > small_thresh) {
                ++nwarnings;
                if (!throttled) {
                    std::cout
                        << "error, imaginary part of human form coordinate "
                        << "relatively large\n"
                        << u << std::endl;

                } else if (nwarnings % nagthreshold == 0) {
                    // throttled, but maybe we'll nag
                    std::cout
                        << "error, imaginary part of human form coordinate, "
                        << nwarnings << " times." << std::endl;

                    ++warningsdecade;

                    if (warningsdecade == 9) {
                        nagthreshold *= 10;
                        warningsdecade = 0;
                    }
                }

                if (nwarnings >= THROTTLE) throttled = true;
            }
        } else {
            // the absolute value of the real part is small, the imaginary part
            // should be similarly small
            if (std::abs(imag(u(i))) > small_thresh) {
                ++nwarnings;

                if (!throttled) {
                    std::cout << "error, real and imaginary parts of human "
                                 "form coordinate "
                              << "both real and similarly small\n"
                              << u << std::endl;

                } else if (nwarnings % nagthreshold == 0) {
                    // throttled, but maybe we'll nag
                    std::cout
                        << "error, imaginary part of human form coordinate, "
                        << nwarnings << " times." << std::endl;

                    ++warningsdecade;

                    if (warningsdecade == 9) {
                        nagthreshold *= 10;
                        warningsdecade = 0;
                    }
                }

                if (nwarnings >= THROTTLE) throttled = true;
            }
        }

        hsymp[i] = u[i].real();
    }

    auto hform = CanonToSyn(hsymp);
    // hform += closed_orbit_;

    std::array<double, 6> hf;
    for (int i = 0; i < 6; ++i)
        hf[i] = hform[i];

    return hf;
}

template <unsigned int order>
typename NormalForm<order>::Matrix6C
NormalForm<order>::ev_ordering(
    typename NormalForm<order>::Vector6C const& ev,
    typename NormalForm<order>::Matrix6C const& B) const
{
    using Matrix = typename NormalForm<order>::Matrix6C;
    using Vector = typename NormalForm<order>::Vector6C;

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

    Matrix Br;

    for (int i = 0; i < 6; ++i)
        Br(i, 0) = B(i, Lidx[0]);
    for (int i = 0; i < 6; ++i)
        Br(i, 1) = B(i, Lidx[1]);
    for (int i = 0; i < 6; ++i)
        Br(i, 2) = B(i, Lidx[2]);
    for (int i = 0; i < 6; ++i)
        Br(i, 3) = B(i, Lidx[3]);
    for (int i = 0; i < 6; ++i)
        Br(i, 4) = B(i, Lidx[4]);
    for (int i = 0; i < 6; ++i)
        Br(i, 5) = B(i, Lidx[5]);

    return Br;
}

#endif
