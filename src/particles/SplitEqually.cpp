/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "SplitEqually.H"

namespace impactx
{
    ParticleChunk
    split_equally (amrex::Long npart, amrex::Long index, amrex::Long size)
    {
        ParticleChunk my_chunk;

        amrex::Long const navg = npart / size;  // note: integer division
        amrex::Long const nleft = npart - navg * size;
        my_chunk.size = (index < nleft) ? navg+1 : navg;  // add 1 to each chunk (proc/thread) until distributed

        if (index < nleft) { // get navg+1 items
            my_chunk.offset = index * (navg + 1);
        } else {                  // get navg items
            my_chunk.offset = index * navg + nleft;
        }

        return my_chunk;
    }

} // namespace impactx
