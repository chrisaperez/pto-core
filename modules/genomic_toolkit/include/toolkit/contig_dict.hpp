// SPDX-License-Identifier: MIT
//
// Interning of contig names to dense integer ids.
//
// Every structure downstream of the readers is keyed by an int32 tid, never by
// a chromosome name. A text pipeline compares "chr14" against "chr14" once per
// record per stage; here the string comparison happens once per *distinct
// name* in each input, and the per-record cost is an integer compare.
//
// The dictionary is also where the single most common silent failure in this
// corner of bioinformatics is caught: a peak file that says "1" and a BAM that
// says "chr1" produce a FRiP of exactly zero with no error anywhere. Binding
// two dictionaries reports the names that failed to match (see
// `unmatched_against`), and the CLI refuses to print a zero FRiP without
// saying so.
#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "toolkit/types.hpp"

namespace toolkit {

class ContigDict {
public:
    // Returns the existing id for `name`, or assigns the next one.
    std::int32_t intern(std::string_view name) {
        const auto it = ids_.find(name);
        if (it != ids_.end()) return it->second;
        const auto id = static_cast<std::int32_t>(names_.size());
        names_.emplace_back(name);
        // The map owns its own copy of the key rather than a view into
        // names_.back(). A view would dangle the moment names_ reallocates:
        // for a short name the characters live *inside* the std::string object
        // (SSO), so growing the vector moves them. Contig counts are in the
        // tens to low thousands, so the duplicate strings are irrelevant; the
        // per-record path never touches either copy.
        ids_.emplace(std::string(name), id);
        return id;
    }

    // Lookup without insertion; kNoTid when absent.
    [[nodiscard]] std::int32_t lookup(std::string_view name) const {
        const auto it = ids_.find(name);
        return it == ids_.end() ? kNoTid : it->second;
    }

    [[nodiscard]] const std::string& name(std::int32_t tid) const {
        static const std::string kUnknown = "*";
        if (tid < 0 || static_cast<std::size_t>(tid) >= names_.size()) return kUnknown;
        return names_[static_cast<std::size_t>(tid)];
    }

    [[nodiscard]] std::size_t size() const noexcept { return names_.size(); }
    [[nodiscard]] const std::vector<std::string>& names() const noexcept {
        return names_;
    }

    // Maps this dictionary's ids onto `other`'s. Result is indexed by this
    // dictionary's tid; entries are kNoTid where the name is absent from
    // `other`. Built once, then used per record.
    [[nodiscard]] std::vector<std::int32_t> map_onto(const ContigDict& other) const {
        std::vector<std::int32_t> out(names_.size(), kNoTid);
        for (std::size_t i = 0; i < names_.size(); ++i) {
            out[i] = other.lookup(names_[i]);
        }
        return out;
    }

    // Names in this dictionary with no counterpart in `other`. Empty is the
    // healthy case; a full list means the two files use different conventions.
    [[nodiscard]] std::vector<std::string> unmatched_against(
        const ContigDict& other) const {
        std::vector<std::string> out;
        for (const auto& n : names_) {
            if (other.lookup(n) == kNoTid) out.push_back(n);
        }
        return out;
    }

private:
    // Heterogeneous lookup, so a std::string_view query does not allocate a
    // std::string on every per-name probe.
    struct Hash {
        using is_transparent = void;
        std::size_t operator()(std::string_view s) const noexcept {
            return std::hash<std::string_view>{}(s);
        }
    };
    struct Eq {
        using is_transparent = void;
        bool operator()(std::string_view a, std::string_view b) const noexcept {
            return a == b;
        }
    };

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::int32_t, Hash, Eq> ids_;
};

}  // namespace toolkit
