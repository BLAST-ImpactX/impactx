#ifndef PCG_DISTRIBUTION_H
#define PCG_DISTRIBUTION_H

#include "synergia/foundation/distribution.h"
#include "synergia/utils/pcg/pcg_random.hpp"
#include <random>

class PCG_random_distribution : public Distribution {
private:
  pcg64 rng;
  std::normal_distribution<double> normal_dist;

public:
  PCG_random_distribution(uint64_t seed, int rank)
    : rng(seed, rank), normal_dist(0.0, 1.0)
  {}

  virtual ~PCG_random_distribution() = default;

  double
  get() override
  {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
  }

  /// Fill a one-dimensional array uniformly between min and max.
  double
  get_uniform(double min, double max) override
  {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
  }

  /// Fill a one-dimensional array with Gaussian distribution of
  /// zero mean and unit standard deviation.
  double
  get_unit_gaussian() override
  {
    return normal_dist(rng);
  }

#if 0
    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    void fill_unit_disk(double* x_array, double* y_array) override
    { 
        const auto Nx = x_array.extent(0);
        const auto Ny = y_array.extent(0);

        if (Nx != Ny) {
            throw std::runtime_error(
                    "Random_distribution::fill_unit_disk: "
                    "x_array and y_array must have the same length\n");
        }

        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (auto n = 0; n < Nx; ++n) {
            double x = 1.0, y = 1.0;
            while (x * x + y * y > 1.0) {
                x = 2.0 * dist(rng) - 1.0;
                y = 2.0 * dist(rng) - 1.0;
            }
            x_array[n] = x;
            y_array[n] = y;
        }
    }
#endif

  void
  advance(uint64_t delta) override
  {
    rng.advance(delta);
  }
};

#endif
