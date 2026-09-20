// SPDX-License-Identifier: MIT
//
// Builds a regions x bins signal matrix from an indexed BAM, with the same
// windowing semantics as `deeptools computeMatrix reference-point`.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "profiler/bam_reader.hpp"
#include "profiler/types.hpp"

namespace profiler {

// Derives the fixed-width query window around a region's anchor. Returns false
// when the window falls entirely outside the chromosome.
bool make_window(const Region& region, const ProfileOptions& opts,
                 std::int64_t chrom_length, Region& out_window);

// Multiplier applied to every bin for the requested normalisation.
// `bin_size` is in base pairs; `total_signal` is the genome-wide denominator
// (mapped reads for CPM/RPKM, the summed per-bp rate for BPM).
double normalization_factor(Normalization norm, std::uint64_t total_mapped,
                            std::int64_t bin_size, double total_signal);

// Computes the matrix. Work is distributed one region per task across
// `opts.threads` workers; each worker leases its own BAM handle so no
// serialisation occurs on the hot path.
//
// `progress`, if set, is invoked from the calling thread with a value in
// [0, 1]. It is never called concurrently.
SignalMatrix compute_matrix(const BamReader& reader,
                            const std::vector<Region>& regions,
                            const ProfileOptions& opts,
                            const FilterOptions& filters,
                            const std::function<void(double)>& progress = {});

// Serialises `matrix` as a deeptools-compatible tab-separated table:
// one header line of bin offsets followed by one line per region.
void write_matrix_tsv(const SignalMatrix& matrix, const std::string& path);

// Writes only the meta-profile curve (bin offset, mean signal) as TSV.
void write_profile_tsv(const SignalMatrix& matrix, const std::string& path);

}  // namespace profiler
