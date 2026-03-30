
#include "synergia/utils/simple_timer.h"

std::map<std::string, simple_timer_counter::timing>
  simple_timer_counter::timings =
    std::map<std::string, simple_timer_counter::timing>();

void
simple_timer_print(Logger& logger)
{
#ifdef SIMPLE_TIMER
  using namespace std;

  logger(LoggerV::INFO)
    //<< std::scientific
    << std::setprecision(8) << left << setw(30) << "timer label" << setw(20)
    << "time(s)" << setw(12) << "cout"
    << "\n"
    << std::string(62, '-') << "\n";

  for (auto const& timing : simple_timer_counter::timings) {
    logger << left << std::setw(30) << timing.first << std::setw(20)
           << timing.second.sum << std::setw(12) << timing.second.count << "\n";
  }

  logger << "\n";
#endif
}
