#include <ctime>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdexcept>

#include "distribution.h"

unsigned long int
Random_distribution::get_default_seed(const char* device)
{
    unsigned long int seed;
    std::ifstream devrandom(device, std::ios::binary);
    if (devrandom) {
        devrandom.read((char*)&seed, sizeof(unsigned long int));
    } else {
        seed = std::time(0);
    }
    return seed;
}

Random_distribution::Random_distribution(unsigned long int seed,
                                         int rank,
                                         Generator generator)
    : rng(nullptr), rng_type(nullptr), rank(rank), original_seed(0)
{
    gsl_rng_env_setup();
    if (generator == ranlxd2) {
        rng_type = gsl_rng_ranlxd2;
    } else if (generator == mt19937) {
        rng_type = gsl_rng_mt19937;
    }
    rng = gsl_rng_alloc(rng_type);
    set_seed(seed); // NOTE: this resets original_seed
}

void
Random_distribution::set_seed(unsigned long int seed)
{
    if (seed == 0) {
        original_seed = get_default_seed();
    } else {
        original_seed = seed;
    }
    unsigned long int distributed_seed;
    distributed_seed =
        (1000 + 5 * (rank + original_seed)) * ((rank + original_seed) + 7) - 1;

    gsl_rng_set(rng, distributed_seed);
}

unsigned long int
Random_distribution::get_original_seed() const
{
    return original_seed;
}

double
Random_distribution::get()
{
    return gsl_rng_uniform(rng);
}

double
Random_distribution::get_uniform(double min, double max)
{
    return gsl_ran_flat(rng, min, max);
}

double
Random_distribution::get_unit_gaussian()
{
    return gsl_ran_ugaussian_ratio_method(rng);
}

#if 0
void
Random_distribution::fill_unit_disk(double* x_array, double* y_array)
{
    const auto Nx = x_array.extent(0);
    const auto Ny = y_array.extent(0);

    if (Nx != Ny) {
        throw std::runtime_error(
                "Random_distribution::fill_unit_disk: x_array and y_array must have the same length\n");
    }

    for (unsigned int n = 0; n < Nx; ++n) {
        double x = 1.0, y = 1.0;
        while (x * x + y * y > 1.0) {
            x = 2.0 * gsl_rng_uniform(rng) - 1.0;
            y = 2.0 * gsl_rng_uniform(rng) - 1.0;
        }
        x_array[n] = x;
        y_array[n] = y;
    }
}
#endif

Random_distribution::~Random_distribution()
{
    gsl_rng_free(rng);
}

void
Random_distribution::advance(uint64_t)
{
    throw std::runtime_error("Random_distribution can not be advanced.");
}
