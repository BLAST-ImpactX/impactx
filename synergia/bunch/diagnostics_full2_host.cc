
#include "Eigen/Core"
#include "Eigen/LU"

using namespace Eigen;

namespace {
    const double tiny = 1.0e-15;

    double
    eliminate_small_negative(double x)
    {
        if ((x < 0) && (x > -tiny))
            return 0;
        else
            return x;
    }
}

void
calculate_emittances(double* mom2,
                     double& emitx,
                     double& emity,
                     double& emitz,
                     double& emitxy,
                     double& emitxyz)
{
    Matrix<double, 6, 6> mom2_matrix(mom2);

    emitx = std::sqrt(eliminate_small_negative(
        mom2_matrix.block<2, 2>(0 /*Bunch::x*/, 0 /*Bunch::x*/).determinant()));
    emity = std::sqrt(eliminate_small_negative(
        mom2_matrix.block<2, 2>(2 /*Bunch::y*/, 2 /*Bunch::y*/).determinant()));
    emitz = std::sqrt(eliminate_small_negative(
        mom2_matrix.block<2, 2>(4 /*Bunch::z*/, 4 /*Bunch::z*/).determinant()));
    emitxy = std::sqrt(eliminate_small_negative(
        mom2_matrix.block<4, 4>(0 /*Bunch::x*/, 0 /*Bunch::x*/).determinant()));

    emitxyz = std::sqrt(eliminate_small_negative(mom2_matrix.determinant()));
}
