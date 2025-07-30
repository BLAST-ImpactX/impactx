/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SPLIT_EQUALLY_H
#define IMPACTX_SPLIT_EQUALLY_H

namespace impactx
{
    /** Offset and size of a chunk in a contigous array (e.g., of particle properties) */
    struct ParticleChunk
    {
        int offset = 0;  //! where to start in the array, in elements
        int size = 0;    //! how large is my chunk, counted from offset, in elements
    };

    /** Split an array of contiguous elements into well-balanched chunks.
     *
     * This function contains the auxiliary math to split npart contiguous particles into
     * contiguous chunks of size, returning the individual starting offset and chunk size.
     * It also deals with distributing remainders well.
     *
     * This can be used for distributed logic for MPI and OpenMP parallelism alike.
     *
     * @param npart total number of elements/particles
     * @param index my index in parallelization, e.g., MPI rank or OpenMP thread id
     * @param size the parallelization space, e.g., MPI size or OpenMP max threads
     * @return offset and size of the chunk for a specific index
     */
    ParticleChunk
    split_equally (int npart, int index, int size);

} // namespace impactx

#endif // IMPACTX_SPLIT_EQUALLY_H
