// SPDX-License-Identifier: MIT
//
// Fragment sources: coordinate-sorted BAM, BEDPE, and 3-column fragment BED
// (the 10x / sinto style file).
//
// All three present the same interface -- a batch callback -- and none of them
// materialises the fragment set. `stream_*` reads its input once, reconstructs
// fragments, applies `FragmentFilter`, accumulates the size histogram, and
// hands out spans of at most kBatchFragments records that are invalidated as
// soon as the callback returns. A consumer that needs to keep fragments copies
// them; the three consumers in this module (duplicate marking, FRiP, size QC)
// do not, which is the entire reason the pipeline has no intermediate file.
//
// BAM support is compile-time optional (GTK_HAVE_HTSLIB). Without htslib the
// text paths still build and the module still passes its whole test suite --
// deliberate, because htslib is the one dependency that makes this module hard
// to build on a locked-down cluster, and the fragment-file path is what an
// scATAC user has anyway.
#pragma once

#include <functional>
#include <span>
#include <string>
#include <vector>

#include "toolkit/contig_dict.hpp"
#include "toolkit/types.hpp"

namespace toolkit {

// Invoked with each batch. The span is valid only for the duration of the
// call. Fragments are mutable so a consumer can set flags in place.
using BatchCallback = std::function<void(std::span<Fragment>)>;

// Streams a BEDPE file. Columns 1-6 are chrom1/start1/end1/chrom2/start2/end2;
// the fragment is the span from min(start) to max(end). Mate pairs on
// different contigs are counted in `dropped_interchrom` and skipped.
//
// Column 8, when present and numeric, is taken as MAPQ; column 9 as strand1.
// Track/comment/header lines are skipped.
StreamStats stream_bedpe(const std::string& path, const FragmentFilter& filter,
                         ContigDict& dict, const BatchCallback& on_batch);

// Streams a 3+ column fragment BED (chrom, start, end, ...). Each line is
// already a fragment, so no pairing is done, and by default each line is ONE
// fragment whatever its later columns hold.
//
// With `filter.use_count_column` (the CLI's `--with-counts`), column 5 is the
// 10x fragments.tsv read-pair count: the fragment is emitted that many times,
// and a present column that is not an integer in [1, 1,000,000] makes the
// record malformed. It is opt-in because the same column in a BED5/BED6 is a
// score, and reading a score as a count multiplies every statistic by it.
// Without the flag a 10x file counts unique fragments, not read pairs.
StreamStats stream_fragment_bed(const std::string& path, const FragmentFilter& filter,
                                ContigDict& dict, const BatchCallback& on_batch);

// Dispatches on extension: .bedpe -> stream_bedpe, .bed/.tsv/.fragments ->
// stream_fragment_bed, .bam/.cram/.sam -> stream_bam. A .gz suffix is handled
// transparently (via htslib's BGZF when available; otherwise a .gz input is a
// clear error rather than a garbled parse).
//
// `region`, when non-empty, restricts the stream to an htslib region string
// ("chr1", "chr1:1-1000000"). It requires an indexed BAM/CRAM and is an error
// on the text formats, which have no index to seek with -- see stream_bam.
StreamStats stream_fragments(const std::string& path, const FragmentFilter& filter,
                             ContigDict& dict, const BatchCallback& on_batch,
                             const std::vector<std::string>& regions = {});

#ifdef GTK_HAVE_HTSLIB
// Streams a coordinate-sorted BAM/CRAM.
//
// Only the read whose TLEN is positive emits a fragment, so a proper pair
// yields exactly one record and no mate lookup, no name hashing and no second
// pass are needed. Single-end reads emit their aligned span only when
// `filter.allow_single_end` is set, flagged kSingleEnd so a consumer can tell
// an inferred span from a measured one.
//
// `dict` is populated from the BAM header in header order, so tids are the
// BAM's own tids and stay comparable across files sharing a reference.
//
// `region` restricts the read to an htslib region string via the index. Two
// things about it are load-bearing for using this to shard a file:
//
//   * A fragment belongs to the region containing its START. htslib's iterator
//     returns every record OVERLAPPING the region, so a fragment straddling a
//     boundary comes back from the queries on BOTH sides; counting it twice is
//     the exact failure a map-reduce over regions must not have. Records whose
//     start precedes the region are therefore skipped, and they are NOT
//     counted as dropped -- they are another shard's records, not this one's.
//   * The index is required. Without one, a "region" scan would silently
//     become a full pass over the file in every shard, turning an N-way split
//     into N times the work for the same answer.
// Several regions read as one stream, their statistics summed. A shard of a
// map-reduce usually owns more than one contig, and running the tool once per
// contig would pay the process and header cost per region instead of per
// shard. Regions must not overlap: each is counted independently, so an
// overlap double-counts exactly as it would across shards.
StreamStats stream_bam(const std::string& path, const FragmentFilter& filter,
                       ContigDict& dict, const BatchCallback& on_batch,
                       const std::vector<std::string>& regions = {});

// True when this build can read BAM. Also available to tests, which skip the
// BAM cases rather than fail when htslib is absent.
inline constexpr bool kHaveBam = true;
#else
inline constexpr bool kHaveBam = false;
#endif

}  // namespace toolkit
