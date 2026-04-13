
#include "synergia/bunch/diagnostics_worker.h"
#include "synergia/bunch/bunch.h"

std::string
Diagnostics_worker::type() const
{
    return diag->type();
}

void
Diagnostics_worker::update(Bunch const& bunch)
{
    diag->update(bunch);
    // need to open file for the HDF5 backend to have
    // a meaningful root rank, previously this was done
    // implicity via get_file
    diag_io.open_file();
    diag->reduce(bunch.get_comm(), diag_io.get_root_rank());
}

void
Diagnostics_worker::write()
{
    diag->write(diag_io.get_io_device(), diag_io.get_count());
    diag_io.finish_write();
}
