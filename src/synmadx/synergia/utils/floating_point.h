#ifndef FLOATING_POINT_H_
#define FLOATING_POINT_H_

#include <cmath>
#include <complex>

/// check a = b within the given tolerance
inline bool
floating_point_equal(double a, double b, double tolerance)
{
    if (std::abs(a) < tolerance) {
        return (std::abs(a - b) < tolerance);
    } else {
        return (std::abs((a - b) / a) < tolerance);
    }
}

/// check a <= b within the given tolerance
inline bool
floating_point_leq(double a, double b, double tolerance)
{
    return ((b - a) > -tolerance);
}

/// check a = b within the given tolerance for complex values
inline bool
complex_floating_point_equal(std::complex<double> a,
                             std::complex<double> b,
                             double tolerance)
{
    if (std::abs(a) < tolerance) {
        return (std::abs(a - b) < tolerance);
    } else {
        return (std::abs((a - b) / a) < tolerance);
    }
}

#endif /* FLOATING_POINT_H_ */
