#ifndef LOGGER_H_
#define LOGGER_H_

#include "synergia/utils/commxx.h"
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

enum class LoggerV {
    DEBUG,
    DINFO,

    INFO,
    INFO_OPN,  // operation
    INFO_OPR,  // operator
    INFO_STEP, // step
    INFO_TURN, // turn

    WARNING,
    ERROR,
};

class Logger {
  private:
    std::ostream* stream_ptr;
    bool have_stream;
    std::ofstream* fstream_ptr;
    bool have_fstream;

    LoggerV verbosity;
    LoggerV severity;

  public:
    /// Log to screen on a single rank
    Logger(int rank = 0, LoggerV verbosity = LoggerV::DINFO, bool log = true);

    /// Log to file and, optionally, screen on a single rank
    Logger(int rank,
           std::string const& filename,
           LoggerV verbosity = LoggerV::DINFO,
           bool screen = true,
           bool log = true);

    /// Log to file and, optionally, screen on a single rank
    /// This variation on the previous constructor is needed
    /// for technical C++ reasons...
    Logger(int rank,
           char const* filename,
           LoggerV verbosity = LoggerV::DINFO,
           bool screen = true,
           bool log = true);

    /// Log to a separate file on each rank
    Logger(std::string const& filename_base,
           LoggerV verbosity = LoggerV::DINFO,
           bool log = true);

    Logger& set_stream(std::ostream& stream);
    Logger& write(std::string const& str);

    template <typename T>
    Logger&
    operator<<(T const& t)
    {
        if (have_stream && (severity >= verbosity)) { (*stream_ptr) << t; }

        if (have_fstream && (severity >= verbosity)) { (*fstream_ptr) << t; }

        return *this;
    }

    Logger&
    operator()(LoggerV severity)
    {
        this->severity = severity;
        return *this;
    }

    Logger& operator<<(std::ostream& (*op)(std::ostream&));
    Logger& flush();

    ~Logger();
};

#endif /* LOGGER_H_ */
