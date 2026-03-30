#include <cmath>
#include <stdexcept>

#include "parallel_utils.h"

void
decompose_1d_raw(int processors,
                 int length,
                 std::vector<int>& offsets,
                 std::vector<int>& counts)
{
    int min_counts = length / processors;
    int remainder = static_cast<int>(fmod(length, processors));
    int offset = 0;
    for (int i = 0; i < processors; ++i) {
        int count = min_counts;
        if (i >= (processors - remainder)) { count += 1; }
        offsets.at(i) = offset;
        counts.at(i) = count;
        offset += count;
    }
}

void
decompose_1d(Commxx const& comm,
             int length,
             std::vector<int>& offsets,
             std::vector<int>& counts)
{
    int size = comm.get_size();
    decompose_1d_raw(size, length, offsets, counts);
}

int
decompose_1d_local(Commxx const& comm, int length)
{
    int size = comm.get_size();
    int rank = comm.get_rank();
    std::vector<int> offsets(size), counts(size);
    decompose_1d_raw(size, length, offsets, counts);
    return counts[rank];
}

std::vector<std::vector<int>>
distribute_1d_raw(int processors, int elements)
{
    std::vector<std::vector<int>> retval(elements);
    std::vector<int> offsets(elements), counts(elements);
    if (processors < elements) {
        decompose_1d_raw(processors, elements, offsets, counts);
        int element = 0;
        for (int processor = 0; processor < processors; ++processor) {
            for (int count = 0; count < counts.at(processor); ++count) {
                retval.at(element).resize(1);
                retval.at(element).at(0) = processor;
                ++element;
            }
        }
    } else {
        decompose_1d_raw(elements, processors, offsets, counts);
        for (int element = 0; element < elements; ++element) {
            retval.at(element).resize(counts.at(element));
            for (int count = 0; count < counts.at(element); ++count) {
                retval.at(element).at(count) = offsets.at(element) + count;
            }
        }
    }
    return retval;
}

std::vector<std::vector<int>>
distribute_1d(Commxx const& comm, int elements)
{
    return distribute_1d_raw(comm.get_size(), elements);
}

void
counts_and_offsets_for_impedance_raw(unsigned int processors,
                                     int length,
                                     std::vector<int>& offsets,
                                     std::vector<int>& counts)
{

    if (offsets.size() != processors)
        throw std::runtime_error(
            "parallel utilis:counts and offsets for imped, offsets size");
    if (counts.size() != processors)
        throw std::runtime_error(
            "parallel utilis:counts and offsets for imped. counts size");
    std::vector<std::vector<int>> ranks(distribute_1d_raw(processors, length));

    int count = 0;
    int proc, proc_saved;
    proc_saved = ranks[0][0];
    for (int elm = 0; elm < length; ++elm) {
        proc = ranks[elm][0];
        counts[proc] += 1;
        if (proc != proc_saved) count += counts[proc_saved];
        offsets[proc] = count;
        proc_saved = proc;
    }
}

void
counts_and_offsets_for_impedance(Commxx const& comm,
                                 int length,
                                 std::vector<int>& offsets,
                                 std::vector<int>& counts)
{
    counts_and_offsets_for_impedance_raw(
        comm.get_size(), length, offsets, counts);
}
