#ifndef APERTURE_OPERATION_H_
#define APERTURE_OPERATION_H_

#include "synergia/simulation/independent_operation.h"
#include "synergia/utils/simple_timer.h"

template <class AP>
class Aperture_operation : public Independent_operation {
  private:
    Lattice_element_slice slice;
    AP ap;

  private:
    void
    apply_impl(Bunch& bunch, Logger& logger) const override
    {
        scoped_simple_timer timer(std::string("aperture_") + ap.type);

        int ndiscarded = bunch.apply_aperture(ap);
        double charge =
            ndiscarded * bunch.get_real_num() / bunch.get_total_num();
        slice.get_lattice_element().deposit_charge(
            charge, bunch.get_bunch_index(), bunch.get_train_index());
    }

  public:
    Aperture_operation(Lattice_element_slice const& slice)
        : Independent_operation("aperture")
        , slice(slice)
        , ap(slice.get_lattice_element())
    {}

    std::string const&
    get_aperture_type() const
    {
        return ap.type;
    }
};

struct Dummy_aperture {
    constexpr static const char* type = "dummy";

    Dummy_aperture(Lattice_element const&) {}

    KOKKOS_INLINE_FUNCTION
    bool
    discard(ConstParticles const&, ConstParticleMasks const&, int p) const
    {
        return false;
    }
};

/// An aperture to remove all particles with infinite and/or NaN coordinates.
struct Finite_aperture {
    constexpr static const char* type = "finite";

    Finite_aperture(Lattice_element const&) {}

    KOKKOS_INLINE_FUNCTION
    bool
    discard(ConstParticles const& parts, ConstParticleMasks const&, int p) const
    {
#if 1
        if (!std::isfinite(parts(p, 0)) || !std::isfinite(parts(p, 1)) ||
            !std::isfinite(parts(p, 2)) || !std::isfinite(parts(p, 3)) ||
            !std::isfinite(parts(p, 4)) || !std::isfinite(parts(p, 5)))
            return true;
#endif

#if 0
        // TODO: std::isinfinite() and -ffast-math/-fno-finite-math-only issue
        if (  __isinf(parts(p, 0)) || __isnan(parts(p, 0))
           || __isinf(parts(p, 1)) || __isnan(parts(p, 1))
           || __isinf(parts(p, 2)) || __isnan(parts(p, 2))
           || __isinf(parts(p, 3)) || __isnan(parts(p, 3))
           || __isinf(parts(p, 4)) || __isnan(parts(p, 4))
           || __isinf(parts(p, 5)) || __isnan(parts(p, 5)) ) return true;
#endif

        double pt = 1.0 + parts(p, 5);
        double px = parts(p, 1);
        double py = parts(p, 3);

        return pt * pt - px * px - py * py < 0.0;
    }
};

/// A circular aperture with radius in meters determined by the
/// Lattice_element attribute "circular_aperture_radius".
/// If the radius is not defined, the default value of 1000.0 m will
/// be used.
struct Circular_aperture {
    constexpr static const char* type = "circular";
    double r2, xoff, yoff;

    Circular_aperture(Lattice_element const& ele)
        : r2(1000.0)
        , xoff(ele.get_double_attribute("hoffset", 0.0))
        , yoff(ele.get_double_attribute("voffset", 0.0))
    {
        double r = ele.get_double_attribute("circular_aperture_radius", 1000.0);
        r2 = r * r;
    }

    KOKKOS_INLINE_FUNCTION
    bool
    discard(ConstParticles const& parts, ConstParticleMasks const&, int p) const
    {
        double xrel = parts(p, 0) - xoff;
        double yrel = parts(p, 2) - yoff;

        double radius2 = xrel * xrel + yrel * yrel;
        return (radius2 > r2);
    }
};

/// An elliptical aperture with horizontal and vertical radii in meters
/// determined by the Lattice_element_attributes
/// "elliptical_aperture_horizontal_radius" and
/// "elliptical_aperture_vertical_radius", respectively.
/// Both radii must be specified. Failing to do so will cause an
/// exception.
struct Elliptical_aperture {
    constexpr static const char* type = "elliptical";
    double h2, v2, xoff, yoff;

    Elliptical_aperture(Lattice_element const& ele)
        : h2(1.0)
        , v2(1.0)
        , xoff(ele.get_double_attribute("hoffset", 0.0))
        , yoff(ele.get_double_attribute("voffset", 0.0))
    {
        double hr =
            ele.get_double_attribute("elliptical_aperture_horizontal_radius");
        double vr =
            ele.get_double_attribute("elliptical_aperture_vertical_radius");
        h2 = hr * hr;
        v2 = vr * vr;
    }

    KOKKOS_INLINE_FUNCTION
    bool
    discard(ConstParticles const& parts, ConstParticleMasks const&, int p) const
    {
        double xrel = parts(p, 0) - xoff;
        double yrel = parts(p, 2) - yoff;

        double scaled_r2 = xrel * xrel / h2 + yrel * yrel / v2;
        return (scaled_r2 > 1.0);
    }
};

/// A rectangular aperture with horizontal and vertical dimensions in meters
/// determined by the Lattice_element_attributes
/// "rectangular_aperture_width" and
/// "rectangular_aperture_height", respectively.
/// Both dimensions must be specified. Failing to do so will cause an
/// exception.
struct Rectangular_aperture {
    constexpr static const char* type = "rectangular";
    double width, height, xoff, yoff;

    Rectangular_aperture(Lattice_element const& ele)
        : width(ele.get_double_attribute("rectangular_aperture_width"))
        , height(ele.get_double_attribute("rectangular_aperture_height"))
        , xoff(ele.get_double_attribute("hoffset", 0.0))
        , yoff(ele.get_double_attribute("voffset", 0.0))
    {}

    KOKKOS_INLINE_FUNCTION
    bool
    discard(ConstParticles const& parts, ConstParticleMasks const&, int p) const
    {
        using Kokkos::fabs;

        double xrel = parts(p, 0) - xoff;
        double yrel = parts(p, 2) - yoff;

        return (fabs(xrel) > 0.5 * width) || (fabs(yrel) > 0.5 * height);
    }
};

/// A polygon aperture with vertices
/// determined by the Lattice_element_attributes
/// "pax1", "pay1", "pax2", "pay2", and so on.
/// And it also requires "the_number_of_vertices", which determines the number
/// of vertices and must be greter than and equal to 3.
/// Must have at least 3 vertcies. Failing to do so will cause an
/// exception.
struct Polygon_aperture {
    constexpr static const int max_vertices = 100;
    constexpr static const char* type = "polygon";

    Kokkos::complex<double> vertices[max_vertices];

    int num_vertices;
    double min_radius2, xoff, yoff;

    Polygon_aperture(Lattice_element const& ele)
        : num_vertices(ele.get_double_attribute("the_number_of_vertices"))
        , min_radius2(ele.get_double_attribute("min_radius2", 0.0))
        , xoff(ele.get_double_attribute("hoffset", 0.0))
        , yoff(ele.get_double_attribute("voffset", 0.0))
    {
        if (num_vertices < 3) {
            throw std::runtime_error(
                "Polygon_aperture: requires at least 3 vertices.");
        }

        if (num_vertices > max_vertices) {
            throw std::runtime_error(
                "Polygon_aperture: the_number_of_vertices exceeds the "
                "max allowed number. Please increase the max value at "
                "Polygon_aperture::max_vertices");
        }

        for (int i = 0; i < num_vertices; ++i) {
            std::stringstream sx;
            std::stringstream sy;

            sx << "pax" << (i + 1);
            sy << "pay" << (i + 1);

            vertices[i] =
                Kokkos::complex<double>(ele.get_double_attribute(sx.str()),
                                        ele.get_double_attribute(sy.str()));
        }
    }

    KOKKOS_INLINE_FUNCTION
    bool
    discard(ConstParticles const& parts, ConstParticleMasks const&, int p) const
    {
        using Kokkos::atan2;

        double xrel = parts(p, 0) - xoff;
        double yrel = parts(p, 2) - yoff;
        double r2 = xrel * xrel + yrel * yrel;

        bool keep = true;

        if (r2 >= min_radius2) {
            Kokkos::complex<double> u(xrel, yrel);
            int idx = 0;
            int size = num_vertices;
            double theta_sum = 0.0;

            while (idx < size) {
                int idx2 = idx + 1;
                if (idx2 == size) idx2 = 0;

                Kokkos::complex<double> v(vertices[idx]);
                Kokkos::complex<double> w(vertices[idx2]);

                auto r = (w - u) * conj(v - u);
                double theta = atan2(r.imag(), r.real());
                theta_sum += theta;
                ++idx;
            }

            const double tiny = 1e-12;
            if (theta_sum / (2.0 * Kokkos::numbers::pi_v<double>) < tiny)
                keep = false;
        }

        return !keep;
    }
};

#endif /* APERTURE_OPERATION_H_ */
