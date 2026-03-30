#ifndef IMPLEMENTED_COLLECTIVE_OPTIONS_H
#define IMPLEMENTED_COLLECTIVE_OPTIONS_H

#include <array>
#include <memory>

#include "synergia/utils/cereal.h"

enum class green_fn_t {
    pointlike,
    linear,
};

enum class LongitudinalDistribution {
    gaussian,
    uniform,
};

#ifdef BUILD_FD_SPACE_CHARGE_SOLVER
struct Space_charge_3d_fd_options {

    std::array<int, 3> shape;
    std::array<double, 3> size;
    std::array<double, 3> offset;
    std::array<double, 3> scale_thresholds;
    bool domain_fixed;
    double n_sigma;
    double kick_scale;
    int comm_group_size;

    Space_charge_3d_fd_options(int gridx = 32,
                               int gridy = 32,
                               int gridz = 64,
                               double scale_x_threshold = 15,
                               double scale_y_threshold = 15,
                               double scale_z_threshold = 2.5)
        : shape{gridx, gridy, gridz}
        , scale_thresholds{scale_x_threshold,
                           scale_y_threshold,
                           scale_z_threshold}
        , domain_fixed(false)
        , n_sigma(8.0)
        , kick_scale(1.0)
        , comm_group_size(1)
    {}

    void
    set_fixed_domain(std::array<double, 3> offset, std::array<double, 3> size)
    {
        this->offset = offset;
        this->size = size;

        this->domain_fixed = true;
    }

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(shape);
        ar(n_sigma);
        ar(comm_group_size);
    }
};
#endif

struct Space_charge_3d_open_hockney_options {

    std::array<int, 3> shape;
    std::array<int, 3> doubled_shape;
    std::array<double, 3> size;
    std::array<double, 3> offset;
    green_fn_t green_fn;
    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    double n_sigma;
    double kick_scale;
    bool domain_fixed;
    int comm_group_size;

    Space_charge_3d_open_hockney_options(int gridx = 32,
                                         int gridy = 32,
                                         int gridz = 64)
        : shape{gridx, gridy, gridz}
        , doubled_shape{gridx * 2, gridy * 2, gridz * 2}
        , green_fn(green_fn_t::linear)
        , periodic_z(false)
        , z_period(0.0)
        , grid_entire_period(false)
        , n_sigma(8.0)
        , kick_scale(1.0)
        , domain_fixed(false)
        , comm_group_size(4)
    {}

    void
    set_fixed_domain(std::array<double, 3> offset, std::array<double, 3> size)
    {
        if (grid_entire_period) {
            this->offset[2] = 0.0;
            this->size[2] = z_period;
        }

        this->offset = offset;
        this->size = size;

        this->domain_fixed = true;
    }

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(shape);
        ar(doubled_shape);
        ar(size);
        ar(offset);
        ar(green_fn);
        ar(periodic_z);
        ar(z_period);
        ar(grid_entire_period);
        ar(n_sigma);
        ar(comm_group_size);
    };
};

struct Space_charge_2d_open_hockney_options {

    std::array<int, 3> shape;
    std::array<int, 3> doubled_shape;
    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    double n_sigma;
    int comm_group_size;

    Space_charge_2d_open_hockney_options(int gridx = 32,
                                         int gridy = 32,
                                         int gridz = 32)
        : shape{gridx, gridy, gridz}
        , doubled_shape{gridx * 2, gridy * 2, gridz}
        , periodic_z(false)
        , z_period(0.0)
        , grid_entire_period(false)
        , n_sigma(8.0)
        , comm_group_size(4)
    {}

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(shape);
        ar(doubled_shape);
        ar(periodic_z);
        ar(z_period);
        ar(grid_entire_period);
        ar(n_sigma);
        ar(comm_group_size);
    }
};

struct Space_charge_2d_kv_options {
    using LD = LongitudinalDistribution;

    // switch to control whether the linear charge density is assumed to be
    // distributed gaussian or uniform over the bunch length.
    LD longitudinal_distribution = LD::uniform;

    bool strictly_linear = true;
    bool strictly_centered = false;

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(longitudinal_distribution);
        ar(strictly_linear);
        ar(strictly_centered);
    }
};

struct Space_charge_rectangular_options {

    std::array<int, 3> shape;
    std::array<double, 3> pipe_size;
    int comm_group_size;

    Space_charge_rectangular_options(
        std::array<int, 3> const& shape = {32, 32, 64},
        std::array<double, 3> const& pipe_size = {0.1, 0.1, 1.0})
        : shape(shape), pipe_size(pipe_size), comm_group_size(1)
    {}

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(shape);
        ar(pipe_size);
        ar(comm_group_size);
    }
};

struct Impedance_options {

    std::string wake_file;
    std::string wake_type;
    int z_grid;
    bool full_machine;
    int nstored_turns;
    int num_buckets;
    double orbit_length;
    double bunch_spacing;
    std::array<int, 3> wn;
    double mwf_xlead;
    double mwf_xtrail;
    double mwf_ylead;
    double mwf_ytrail;
    double mwf_zwake;

    Impedance_options(std::string const& wake_file = "",
                      std::string const& wake_type = "",
                      int z_grid = 1000)
        : wake_file(wake_file)
        , wake_type(wake_type)
        , z_grid(z_grid)
        , full_machine(false)
        , nstored_turns(15)
        , num_buckets(1)
        , orbit_length(1)
        , bunch_spacing(1)
        , mwf_xlead(1.0)
        , mwf_xtrail(1.0)
        , mwf_ylead(1.0)
        , mwf_ytrail(1.0)
        , mwf_zwake(1.0)
    {}

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(wake_file);
        ar(wake_type);
        ar(z_grid);
        ar(full_machine);
        ar(nstored_turns);
        ar(num_buckets);
        ar(orbit_length);
        ar(bunch_spacing);
        ar(mwf_xlead);
        ar(mwf_xtrail);
        ar(mwf_ylead);
        ar(mwf_ytrail);
        ar(mwf_zwake);
    }
};

struct Dummy_CO_options {

    template <class Archive>
    void
    serialize(Archive& ar)
    {}
};

#endif
