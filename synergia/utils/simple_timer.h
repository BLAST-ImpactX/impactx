#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H

#include <iomanip>
#include <map>

#include "logger.h"

struct simple_timer_counter {
  struct timing {
    double sum;
    double start;
    int count;
  };
  static std::map<std::string, timing> timings;

  static void
  start(std::string const& label, double t0)
  {
    timings.emplace(label, timing{0.0, t0, 0}).first->second.start = t0;
  }

  static void
  stop(std::string const& label, double t1)
  {
    auto iter = timings.emplace(label, timing{0.0, t1, 0}).first;
    iter->second.sum += t1 - iter->second.start;
    ++iter->second.count;
  }
};

inline void
simple_timer_start(std::string const& label)
{
#ifdef SIMPLE_TIMER
  MPI_Barrier(MPI_COMM_WORLD);
  simple_timer_counter::start(label, MPI_Wtime());
#endif
}

inline void
simple_timer_stop(std::string const& label)
{
#ifdef SIMPLE_TIMER
  MPI_Barrier(MPI_COMM_WORLD);
  simple_timer_counter::stop(label, MPI_Wtime());
#endif
}

struct scoped_simple_timer {
#ifdef SIMPLE_TIMER
  std::string const label;

  scoped_simple_timer(std::string const& label) : label(label)
  {
    simple_timer_start(label);
  }

  ~scoped_simple_timer()
  {
    simple_timer_stop(label);
  }
#else
  scoped_simple_timer(std::string const&) {}
#endif
};

void simple_timer_print(Logger& logger);

#endif
