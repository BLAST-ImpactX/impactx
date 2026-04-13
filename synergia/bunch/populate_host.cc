
#include "populate_host.h"

#include "Eigen/Cholesky"
#include "Eigen/Eigen"

using namespace Eigen;

void
adjust_moments_host(double const* means,
                    double const* covariances,
                    double const* bunch_mean,
                    double const* bunch_mom2,
                    int num_particles,
                    int num_particles_slots,
                    double* particles)
{
    Matrix<double, 6, 6, Eigen::RowMajor> C(covariances);
    Matrix<double, 6, 6, Eigen::RowMajor> G(C.llt().matrixL());
    Matrix<double, 6, 6, Eigen::RowMajor> X(bunch_mom2);
    Matrix<double, 6, 6, Eigen::RowMajor> H(X.llt().matrixL());
    Matrix<double, 6, 6, Eigen::RowMajor> A(G * H.inverse());

    // jfa: dummy exists only to work around a bad interaction betwen
    //      Eigen3 and g++ 4.1.2
    std::stringstream dummy;
    dummy << C;

    Eigen::Map<Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
        rho7(particles, num_particles_slots, 7);

    Matrix<double, 1, 6> rhobar6(bunch_mean);

    for (int part = 0; part < num_particles; ++part) {
        rho7.block<1, 6>(part, 0) -= rhobar6;
    }

    rho7.block(0, 0, num_particles, 6) *= A.transpose();

    Matrix<double, 1, 6> means6(means);

    for (int part = 0; part < num_particles; ++part) {
        rho7.block<1, 6>(part, 0) += means6;
    }
}

void
get_correlation_matrix_host(double* correlation_matrix,
                            double const* one_turn_map,
                            double arms,
                            double brms,
                            double crms,
                            double beta,
                            std::array<int, 3> const& rms_index)
{

    Matrix<double, 6, 6, Eigen::RowMajor> c_matrix;
    Matrix<double, 6, 6, Eigen::RowMajor> eigen_map(one_turn_map);

    EigenSolver<MatrixXd> es(eigen_map);
    VectorXcd evals = es.eigenvalues();
    MatrixXcd evect_matrix = es.eigenvectors();

    std::vector<MatrixXd> F;
    std::vector<int> remaining;
    for (int j = 5; j > -1; j--)
        remaining.push_back(j);

    for (int i = 0; i < 3; i++) {
        // find complex conjugate among remaining eigenvectors
        int first = remaining.back();
        remaining.pop_back();

        double best = 1.0e30;
        int conj = -1;

        for (int item = 0; item < remaining.size(); item++) {
            VectorXcd sum =
                evect_matrix.col(first) + evect_matrix.col(remaining[item]);

            if (sum.imag().cwiseAbs().maxCoeff() < best) {
                best = sum.imag().cwiseAbs().maxCoeff();
                conj = remaining[item];
            }
        }

        if (conj == -1)
            throw std::runtime_error("failed to find a conjugate pair in "
                                     "_get_correlation_matrix");

        remaining.erase(std::remove(remaining.begin(), remaining.end(), conj),
                        remaining.end());

        MatrixXd tmp = (evect_matrix.col(first) *
                            evect_matrix.col(first).conjugate().transpose() +
                        evect_matrix.col(conj) *
                            evect_matrix.col(conj).conjugate().transpose())
                           .real();

        F.push_back(tmp);
        //  F[i] is effectively 2*e[i] cross e^H[i].
    }

    // The correlation matrix is a linear combination of F[i] with
    // appropriate coefficients such that the diagonal elements C[i,i] i=(0,2,4)
    // come out to be the desired 2nd moments.
    Eigen::MatrixXd S(3, 3);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            S(i, j) = F[j](rms_index[i], rms_index[i]);

    Eigen::MatrixXd Sinv = S.inverse();

    std::array<double, 6> units = {1.0, 1.0, 1.0, 1.0, 1.0 / beta, 1.0};

    double cd1 = arms * units[rms_index[0]] * arms * units[rms_index[0]];
    double cd2 = brms * units[rms_index[1]] * brms * units[rms_index[1]];
    double cd3 = crms * units[rms_index[2]] * crms * units[rms_index[2]];

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            c_matrix(i, j) = 0.0;

            for (int k = 0; k < 3; k++)
                c_matrix(i, j) +=
                    F[k](i, j) *
                    (Sinv(k, 0) * cd1 + Sinv(k, 1) * cd2 + Sinv(k, 2) * cd3);
        }
    }

    for (int i = 0; i < 36; ++i)
        correlation_matrix[i] = c_matrix.data()[i];
}
