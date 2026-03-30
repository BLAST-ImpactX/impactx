#ifndef DIAGNOSTICS_WORKER_H
#define DIAGNOSTICS_WORKER_H

#include <memory>
#include <string>

#include <cereal/types/memory.hpp>

#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_io.h"

class Diagnostics_worker {

  private:
    std::shared_ptr<Diagnostics> diag;
    Diagnostics_io diag_io;

  public:
    // default constructor for serialization only
    Diagnostics_worker() : diag(), diag_io() {}

    // construct a diag worker with given type of diag and filename
    // specialization is provided for s_p<Diagnostics> so the python
    // interface can register
    template <class DiagCal>
    Diagnostics_worker(DiagCal const& diag, std::shared_ptr<Commxx> const& comm)
        : diag(std::make_shared<DiagCal>(diag))
        , diag_io(diag.filename(), diag.single_file(), comm)
    {}

    // for registering from python only
    Diagnostics_worker(std::shared_ptr<Diagnostics> const& diag,
                       std::shared_ptr<Commxx> const& comm)
        : diag(diag), diag_io(diag->filename(), diag->single_file(), comm)
    {}

    std::string type() const;

    void update(Bunch const& bunch);
    void write();

    void
    update_and_write(Bunch const& bunch)
    {
        update(bunch);
        write();
    }

  private:
    friend class cereal::access;

    template <class AR>
    void
    serialize(AR& ar)
    {
        ar(CEREAL_NVP(diag));
        ar(CEREAL_NVP(diag_io));
    }
};

class Diagnostics_handler {
  private:
    Diagnostics_worker* worker;
    Bunch const* bunch;

  public:
    Diagnostics_handler() : worker(nullptr), bunch(nullptr) {}

    Diagnostics_handler(Diagnostics_worker& worker, Bunch const& bunch)
        : worker(&worker), bunch(&bunch)
    {}

    std::string
    type() const
    {
        return worker ? worker->type() : "";
    }

    void
    update()
    {
        if (worker) worker->update(*bunch);
    }

    void
    write()
    {
        if (worker) worker->write();
    }

    void
    update_and_write()
    {
        if (worker) worker->update_and_write(*bunch);
    }
};

#endif
