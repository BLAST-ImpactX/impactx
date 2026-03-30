#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include <cstdint>
#include <gsl/gsl_rng.h>

/// Distribution is a virtual base class for obtaining the next number or set
/// of numbers from a sequence according to a limited set of shapes.
class Distribution {
  public:
    virtual ~Distribution() = default;

    /// Get the next number in the sequence (between 0 and 1).
    virtual double get() = 0;

    /// Fill a one-dimensional array uniformly between min and max.
    virtual double get_uniform(double min, double max) = 0;

    /// Fill a one-dimensional array with Gaussian distribution of
    /// zero mean and unit standard deviation.
    virtual double get_unit_gaussian() = 0;

    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    /// virtual void fill_unit_disk(double* x_array, double* y_array) = 0;

    /// Skip ahead the random number generator by delta.
    /// NOTE: subclasses may implement this as a no-op.
    virtual void advance(uint64_t delta) = 0;
};

/// Random_distribution provides a Distribution of random numbers. The random
/// seed is maintained across multiple processors. The implementation uses
/// random numbers from the GNU Scientific Library.
class Random_distribution : public Distribution {

  private:
    gsl_rng* rng;
    const gsl_rng_type* rng_type;
    int rank;
    unsigned long int original_seed;

  public:
    enum Generator { ranlxd2, mt19937 };

    /// Construct a Random_distribution.
    /// @param seed The random number seed. If seed == 0, the seed is
    /// obtained from Random_distribution::get_default_seed().
    /// @param rank Set the seed for this rank.
    /// @param generator The underlying random number generator to be used.
    Random_distribution(unsigned long int seed,
                        int rank,
                        Generator generator = ranlxd2);

    virtual ~Random_distribution();

    /// Generate a random seed. Attempt to read from device if present.
    /// Otherwise, use the system clock.
    /// @param device Read from pathname device.
    static unsigned long int get_default_seed(
        const char* device = "/dev/urandom");

    /// Set the random number generator seed.
    /// @param seed The seed.
    void set_seed(unsigned long int seed);

    /// Get the seed used to start the random number generator.
    unsigned long int get_original_seed() const;

    /// Get the next random number between 0 and 1.
    double get() override;

    /// Fill a one-dimensional array uniformly between min and max.
    double get_uniform(double min, double max) override;

    /// Fill a one-dimensional array with Gaussian distribution of
    /// zero mean and unit standard deviation.
    double get_unit_gaussian() override;

    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    /// void fill_unit_disk(double* x_array, double* y_array) override;

    void advance(uint64_t delta) override;

    void test();
};

#endif /* DISTRIBUTION_H_ */
