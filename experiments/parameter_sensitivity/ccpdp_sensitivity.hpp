#ifndef JOINT_DP_CONSENSUS_V8_ADAPTIVE_HPP
#define JOINT_DP_CONSENSUS_V8_ADAPTIVE_HPP

// CCP-DP v8 adaptive release header.
// Contains only code reachable from the production CCP-DP entry point.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mcs_phase
{

static inline int base_id(char c)
{
    switch (c)
    {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
        return 3;
    default:
        return -1;
    }
}

static inline char id_base(int x)
{
    static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
    return (x >= 0 && x < 4) ? bases[x] : 'A';
}

static inline std::string normalize_seq(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s)
    {
        int b = base_id(c);
        if (b >= 0)
            out.push_back(id_base(b));
    }
    return out;
}

struct Params
{
    int target_len = 110;
    // Sensitivity-only override. The production default is 40.
    int coverage_switch_reads = 40;
    int kmer = 5;
    int max_reads = 36;
    int max_len_delta = 35;
    int band_extra = 13;
    int max_insertion = 5;
    int insertion_choices = 2;
    bool uniform_read_weight = true;

    // Channel-informed regularization, expressed in equivalent read votes.
    double no_insertion_prior = 0.85;
    double no_deletion_prior = 0.70;
    double gap_event_penalty = 0.18;
    double phase_drift_penalty = 0.0;

    // Estimate one coherent read core before building a low-coverage
    // consensus profile.
    int exact_medoid_max_reads = 20;
    double coherence_min_norm = 0.16;
    double coherence_max_norm = 0.28;
    double coherence_mad_scale = 3.0;

    // Library topology used to choose the order of one cross-cluster prior.
    int topology_k = 11;
    int topology_min_support = 2;

    // High-coverage fallback: project local k-mers to normalized strand
    // coordinates and recover one overlap-consistent phase chain.
    bool positional_phase_chain = true;
    int positional_min_reads = 60;
    int positional_k = 7;
    int positional_radius = 2;
    double positional_reference_bonus = 0.35;
    bool blind_channel = true;
    int loo_markov_k = 8;
    int loo_local_k = 5;
    double loo_local_weight = 1.25;
    double loo_reverse_weight = 1.00;
};

struct Diag
{
    int cluster_id = 0;
    int num_reads = 0;
    int used_reads = 0;
    int seed_read_index = -1;
    int seed_len = 0;
    int passes = 0;
    int aligned_reads = 0;
    int rejected_reads = 0;
    int final_len = 0;
    int seed_to_final_edits = 0;
    int weak_positions = 0;
    int coherence_filtered = 0;
    double coherence_threshold = 1.0;
    double mean_norm_dist = 1.0;
    double mean_chosen_conf = 0.0;
};

struct Result
{
    std::string sequence;
    Diag diag;
};

static inline uint32_t kmer_code(const std::string &s, int pos, int k)
{
    uint32_t x = 0;
    for (int i = 0; i < k; ++i)
        x = (x << 2) | static_cast<uint32_t>(base_id(s[pos + i]));
    return x;
}

static inline int edit_distance_linear(const std::string &a, const std::string &b)
{
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());
    std::vector<int> prev(m + 1), cur(m + 1);
    std::iota(prev.begin(), prev.end(), 0);
    for (int i = 1; i <= n; ++i)
    {
        cur[0] = i;
        for (int j = 1; j <= m; ++j)
        {
            int sub = prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
            cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, sub});
        }
        std::swap(prev, cur);
    }
    return prev[m];
}

static inline int banded_distance_only(const std::string &a,
                                       const std::string &b,
                                       int band_extra)
{
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());
    const int band = std::max(std::abs(n - m) + 2,
                              band_extra + std::abs(n - m));
    const uint16_t INF = 30000;
    std::vector<uint16_t> prev(m + 1, INF), cur(m + 1, INF);
    for (int j = 0; j <= std::min(m, band); ++j)
        prev[j] = static_cast<uint16_t>(j);
    for (int i = 1; i <= n; ++i)
    {
        std::fill(cur.begin(), cur.end(), INF);
        const int lo = std::max(0, i - band);
        const int hi = std::min(m, i + band);
        if (lo == 0)
            cur[0] = static_cast<uint16_t>(i);
        for (int j = std::max(1, lo); j <= hi; ++j)
        {
            const uint16_t sub = prev[j - 1] < INF
                                     ? static_cast<uint16_t>(
                                           prev[j - 1] + (a[i - 1] != b[j - 1]))
                                     : INF;
            const uint16_t del = prev[j] < INF
                                     ? static_cast<uint16_t>(prev[j] + 1)
                                     : INF;
            const uint16_t ins = cur[j - 1] < INF
                                     ? static_cast<uint16_t>(cur[j - 1] + 1)
                                     : INF;
            cur[j] = std::min({sub, del, ins});
        }
        std::swap(prev, cur);
    }
    return prev[m] < INF ? prev[m] : std::max(n, m);
}

static inline int choose_seed(const std::vector<std::string> &reads, const Params &p)
{
    if (reads.size() == 1)
        return 0;

    if (static_cast<int>(reads.size()) <= p.exact_medoid_max_reads)
    {
        std::vector<double> costs(reads.size(), 0.0);
        for (int i = 0; i < static_cast<int>(reads.size()); ++i)
        {
            for (int j = i + 1; j < static_cast<int>(reads.size()); ++j)
            {
                int d = banded_distance_only(reads[i], reads[j], p.band_extra);
                double norm = static_cast<double>(d) /
                              std::max(1, std::max(static_cast<int>(reads[i].size()),
                                                  static_cast<int>(reads[j].size())));
                // Capping prevents a single contaminant from dominating the
                // medoid while preserving ordering among coherent reads.
                double robust = std::min(0.32, norm);
                costs[i] += robust;
                costs[j] += robust;
            }
            costs[i] += 0.010 *
                        std::abs(static_cast<int>(reads[i].size()) - p.target_len);
        }
        return static_cast<int>(std::min_element(costs.begin(), costs.end()) -
                                costs.begin());
    }

    const int k = p.kmer;
    const int universe = 1 << (2 * k);
    std::vector<uint16_t> support(universe, 0);
    std::vector<uint32_t> seen_stamp(universe, 0);
    uint32_t stamp = 0;

    for (const auto &r : reads)
    {
        ++stamp;
        if (static_cast<int>(r.size()) < k)
            continue;
        for (int i = 0; i + k <= static_cast<int>(r.size()); ++i)
        {
            uint32_t code = kmer_code(r, i, k);
            if (seen_stamp[code] != stamp)
            {
                seen_stamp[code] = stamp;
                if (support[code] != std::numeric_limits<uint16_t>::max())
                    ++support[code];
            }
        }
    }

    int best = 0;
    double best_score = -1e100;
    for (int rid = 0; rid < static_cast<int>(reads.size()); ++rid)
    {
        const auto &r = reads[rid];
        double sum = 0.0;
        int count = 0;
        ++stamp;
        if (static_cast<int>(r.size()) >= k)
        {
            for (int i = 0; i + k <= static_cast<int>(r.size()); ++i)
            {
                uint32_t code = kmer_code(r, i, k);
                if (seen_stamp[code] == stamp)
                    continue;
                seen_stamp[code] = stamp;
                sum += std::min<int>(support[code], 12);
                ++count;
            }
        }
        const double centrality = count > 0 ? sum / count : 0.0;
        const double len_penalty = 0.055 * std::abs(static_cast<int>(r.size()) - p.target_len);
        const double score = centrality - len_penalty;
        if (score > best_score)
        {
            best_score = score;
            best = rid;
        }
    }
    return best;
}

static inline std::vector<std::string> coherent_read_core(
    const std::vector<std::string> &reads,
    int seed,
    const Params &p,
    int *filtered_out,
    double *threshold_out)
{
    if (filtered_out)
        *filtered_out = 0;
    if (threshold_out)
        *threshold_out = 1.0;
    if (reads.size() <= 4)
        return reads;

    std::vector<double> distances(reads.size(), 0.0);
    for (int i = 0; i < static_cast<int>(reads.size()); ++i)
    {
        int d = banded_distance_only(reads[i], reads[seed], p.band_extra);
        distances[i] = static_cast<double>(d) /
                       std::max(1, std::max(static_cast<int>(reads[i].size()),
                                           static_cast<int>(reads[seed].size())));
    }
    std::vector<double> sorted = distances;
    const size_t middle = sorted.size() / 2;
    std::nth_element(sorted.begin(), sorted.begin() + middle, sorted.end());
    const double center = sorted[middle];
    std::vector<double> deviations;
    deviations.reserve(distances.size());
    for (double x : distances)
        deviations.push_back(std::abs(x - center));
    std::nth_element(deviations.begin(), deviations.begin() + middle,
                     deviations.end());
    const double mad = deviations[middle];
    const double threshold = std::max(
        p.coherence_min_norm,
        std::min(p.coherence_max_norm,
                 center + p.coherence_mad_scale * std::max(0.012, mad)));

    std::vector<std::string> core;
    core.reserve(reads.size());
    for (int i = 0; i < static_cast<int>(reads.size()); ++i)
    {
        if (i == seed || distances[i] <= threshold)
            core.push_back(reads[i]);
    }
    // A core smaller than half the cluster is not trustworthy: in that case
    // retain all reads and let the robust consensus weights handle them.
    if (core.size() * 2 < reads.size() || core.size() < 3)
        return reads;
    if (filtered_out)
        *filtered_out = static_cast<int>(reads.size() - core.size());
    if (threshold_out)
        *threshold_out = threshold;
    return core;
}

struct Alignment
{
    int distance = 0;
    bool valid = false;
    std::string ref_aln;
    std::string read_aln;
};

// Global edit alignment with a diagonal band. The full parent matrix is tiny
// (normally about 110 x 140); only cells inside the band are evaluated.
static inline Alignment banded_align(const std::string &read,
                                     const std::string &ref,
                                     int band_extra)
{
    const int n = static_cast<int>(read.size());
    const int m = static_cast<int>(ref.size());
    const int band = std::max(std::abs(n - m) + 2, band_extra + std::abs(n - m));
    const int W = m + 1;
    const uint16_t INF = 30000;
    std::vector<uint16_t> dp(static_cast<size_t>(n + 1) * W, INF);
    std::vector<uint8_t> parent(static_cast<size_t>(n + 1) * W, 0);
    auto index = [W](int i, int j) -> size_t
    {
        return static_cast<size_t>(i) * W + j;
    };
    dp[0] = 0;
    for (int i = 0; i <= n; ++i)
    {
        const int lo = std::max(0, i - band);
        const int hi = std::min(m, i + band);
        for (int j = lo; j <= hi; ++j)
        {
            if (i == 0 && j == 0)
                continue;
            uint16_t best = INF;
            uint8_t op = 0;
            // Prefer a diagonal move on exact ties. This avoids turning one
            // substitution into an arbitrary insertion/deletion pair.
            if (i > 0 && j > 0)
            {
                uint16_t prev = dp[index(i - 1, j - 1)];
                if (prev < INF)
                {
                    best = static_cast<uint16_t>(prev + (read[i - 1] == ref[j - 1] ? 0 : 1));
                    op = 1;
                }
            }
            if (i > 0)
            {
                uint16_t prev = dp[index(i - 1, j)];
                uint16_t cand = prev < INF ? static_cast<uint16_t>(prev + 1) : INF;
                if (cand < best)
                {
                    best = cand;
                    op = 2;
                }
            }
            if (j > 0)
            {
                uint16_t prev = dp[index(i, j - 1)];
                uint16_t cand = prev < INF ? static_cast<uint16_t>(prev + 1) : INF;
                if (cand < best)
                {
                    best = cand;
                    op = 3;
                }
            }
            dp[index(i, j)] = best;
            parent[index(i, j)] = op;
        }
    }

    Alignment out;
    if (dp[index(n, m)] >= INF)
        return out;
    out.distance = dp[index(n, m)];
    out.valid = true;
    out.ref_aln.reserve(n + m);
    out.read_aln.reserve(n + m);
    int i = n;
    int j = m;
    while (i > 0 || j > 0)
    {
        uint8_t op = parent[index(i, j)];
        if (op == 1)
        {
            out.ref_aln.push_back(ref[j - 1]);
            out.read_aln.push_back(read[i - 1]);
            --i;
            --j;
        }
        else if (op == 2)
        {
            out.ref_aln.push_back('-');
            out.read_aln.push_back(read[i - 1]);
            --i;
        }
        else if (op == 3)
        {
            out.ref_aln.push_back(ref[j - 1]);
            out.read_aln.push_back('-');
            --j;
        }
        else
        {
            out.valid = false;
            out.ref_aln.clear();
            out.read_aln.clear();
            return out;
        }
    }
    std::reverse(out.ref_aln.begin(), out.ref_aln.end());
    std::reverse(out.read_aln.begin(), out.read_aln.end());
    return out;
}

static inline void left_normalize_gap_row(std::string &gap_row,
                                          const std::string &other_row)
{
    const int n = static_cast<int>(gap_row.size());
    for (int s = 1; s < n;)
    {
        if (gap_row[s] != '-' || gap_row[s - 1] == '-')
        {
            ++s;
            continue;
        }
        int e = s;
        while (e < n && gap_row[e] == '-')
            ++e;
        while (s > 0 && other_row[s - 1] != '-' && other_row[e - 1] != '-')
        {
            const int old_mismatch = gap_row[s - 1] != other_row[s - 1];
            const int new_mismatch = gap_row[s - 1] != other_row[e - 1];
            if (new_mismatch > old_mismatch)
                break;
            const char moved = gap_row[s - 1];
            for (int c = s - 1; c + 1 < e; ++c)
                gap_row[c] = '-';
            gap_row[e - 1] = moved;
            --s;
            --e;
        }
        s = std::max(s + 1, e + 1);
    }
}

// Different minimum-edit tracebacks can place the same indel at different
// columns inside a repeat. Canonicalizing every equal-cost gap to the left
// prevents support for one physical event from being split across slots.
static inline void left_normalize_gaps(Alignment &aln)
{
    if (!aln.valid || aln.ref_aln.size() != aln.read_aln.size())
        return;
    left_normalize_gap_row(aln.ref_aln, aln.read_aln);
    left_normalize_gap_row(aln.read_aln, aln.ref_aln);
    left_normalize_gap_row(aln.ref_aln, aln.read_aln);
}

struct ChannelCounts
{
    double matches = 0.0;
    double substitutions = 0.0;
    double insertions = 0.0;
    double deletions = 0.0;
};

struct ChannelProfile
{
    double match = 0.97;
    double substitution = 0.01;
    double insertion = 0.01;
    double deletion = 0.01;
};

static inline double median_value(std::vector<double> values)
{
    if (values.empty())
        return 0.0;
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() & 1U)
        return upper;
    std::nth_element(values.begin(), values.begin() + middle - 1,
                     values.begin() + middle);
    return 0.5 * (upper + values[middle - 1]);
}

static inline void add_alignment_channel_counts(
    const Alignment &alignment,
    ChannelCounts &counts)
{
    if (!alignment.valid ||
        alignment.ref_aln.size() != alignment.read_aln.size())
        return;
    for (size_t column = 0; column < alignment.ref_aln.size(); ++column)
    {
        const char reference_base = alignment.ref_aln[column];
        const char read_base = alignment.read_aln[column];
        if (reference_base == '-')
            counts.insertions += read_base != '-';
        else if (read_base == '-')
            counts.deletions += 1.0;
        else if (reference_base == read_base)
            counts.matches += 1.0;
        else
            counts.substitutions += 1.0;
    }
}

// Estimate an IDS channel without centers. Reads are aligned only to the
// current cluster reconstruction. A median/MAD rule rejects incoherent reads;
// no reported dataset error rate is consumed.
static inline ChannelCounts estimate_cluster_channel(
    const std::vector<std::string> &raw_reads,
    const std::string &reference,
    const Params &p)
{
    struct Observation
    {
        Alignment alignment;
        double normalized_distance = 1.0;
    };
    std::vector<Observation> observations;
    std::vector<double> distances;
    // Channel calibration needs an unbiased event sample, not every trace.
    // A small deterministic cap keeps its cost independent of exceptionally
    // deep clusters; thousands of clusters still provide a large global
    // sample.
    const int sample_count = std::min<int>(
        raw_reads.size(), std::max(1, std::min(p.max_reads, 24)));
    observations.reserve(sample_count);
    distances.reserve(sample_count);
    for (int sample = 0; sample < sample_count; ++sample)
    {
        const int read_index = static_cast<int>(
            static_cast<long long>(sample) * raw_reads.size() /
            sample_count);
        const std::string &raw_read = raw_reads[read_index];
        const std::string read = normalize_seq(raw_read);
        if (read.empty() ||
            std::abs(static_cast<int>(read.size()) - p.target_len) >
                p.max_len_delta)
            continue;
        Alignment alignment = banded_align(read, reference, p.band_extra);
        if (!alignment.valid)
            continue;
        left_normalize_gaps(alignment);
        const double normalized_distance =
            static_cast<double>(alignment.distance) /
            std::max<int>(
                1, std::max(static_cast<int>(read.size()),
                            static_cast<int>(reference.size())));
        observations.push_back(
            {std::move(alignment), normalized_distance});
        distances.push_back(normalized_distance);
    }

    ChannelCounts counts;
    if (observations.empty())
        return counts;
    const double median = median_value(distances);
    std::vector<double> deviations;
    deviations.reserve(distances.size());
    for (double distance : distances)
        deviations.push_back(std::abs(distance - median));
    const double mad = median_value(std::move(deviations));
    const double robust_limit = std::min(
        0.42, median + 3.0 * std::max(mad, 1.0 / std::max(1, p.target_len)));
    for (const Observation &observation : observations)
    {
        if (observation.normalized_distance > robust_limit)
            continue;
        add_alignment_channel_counts(observation.alignment, counts);
    }
    return counts;
}

static inline ChannelProfile channel_profile_from_alignment_counts(
    const ChannelCounts &counts)
{
    ChannelProfile profile;
    const double source_bases =
        counts.matches + counts.substitutions + counts.deletions;
    if (source_bases <= 0.0)
        return profile;
    profile.substitution = counts.substitutions / source_bases;
    profile.insertion = counts.insertions / source_bases;
    profile.deletion = counts.deletions / source_bases;
    const double total_error =
        profile.substitution + profile.insertion + profile.deletion;
    profile.match = std::max(1e-6, 1.0 - total_error);
    return profile;
}

static inline ChannelCounts add_channel_counts(
    const ChannelCounts &a,
    const ChannelCounts &b)
{
    ChannelCounts out;
    out.matches = a.matches + b.matches;
    out.substitutions = a.substitutions + b.substitutions;
    out.insertions = a.insertions + b.insertions;
    out.deletions = a.deletions + b.deletions;
    return out;
}

static inline int read_weight(int distance, int read_len, int ref_len)
{
    const double norm = static_cast<double>(distance) / std::max(1, std::max(read_len, ref_len));
    if (norm <= 0.075)
        return 256;
    if (norm <= 0.115)
        return 224;
    if (norm <= 0.165)
        return 176;
    if (norm <= 0.225)
        return 112;
    if (norm <= 0.300)
        return 56;
    return 20;
}

struct Choice
{
    std::string emit;
    double score = -1e100;
    double confidence = 0.0;
};

struct PassResult
{
    std::string sequence;
    int aligned_reads = 0;
    int rejected_reads = 0;
    int weak_positions = 0;
    double mean_norm_dist = 1.0;
    double mean_chosen_conf = 0.0;
};

static inline void accumulate_alignment_events(
    const Alignment &aln,
    int reference_len,
    int weight,
    int max_insertion,
    std::vector<std::array<int, 5>> &base_counts,
    std::vector<std::unordered_map<std::string, int>> &ins_counts)
{
    std::vector<std::string> inserted(reference_len + 1);
    std::vector<int> state(reference_len, 4);
    int rp = 0;
    for (int c = 0; c < static_cast<int>(aln.ref_aln.size()); ++c)
    {
        const char rc = aln.ref_aln[c];
        const char qc = aln.read_aln[c];
        if (rc == '-')
        {
            if (qc != '-' && rp >= 0 && rp <= reference_len &&
                static_cast<int>(inserted[rp].size()) < max_insertion)
                inserted[rp].push_back(qc);
        }
        else
        {
            if (rp < reference_len)
                state[rp] = qc == '-' ? 4 : base_id(qc);
            ++rp;
        }
    }
    for (int slot = 0; slot <= reference_len; ++slot)
        ins_counts[slot][inserted[slot]] += weight;
    for (int pos = 0; pos < reference_len; ++pos)
    {
        int b = state[pos];
        if (b < 0 || b > 4)
            b = 4;
        base_counts[pos][b] += weight;
    }
}

static inline PassResult consensus_pass(const std::vector<std::string> &reads,
                                        const std::string &ref,
                                        const Params &p)
{
    PassResult out;
    const int m = static_cast<int>(ref.size());
    if (m == 0)
        return out;

    std::vector<std::array<int, 5>> base_counts(m);
    for (auto &x : base_counts)
        x.fill(0);
    std::vector<std::unordered_map<std::string, int>> ins_counts(m + 1);
    int total_evidence_weight = 0;
    double norm_dist_sum = 0.0;

    for (const auto &read : reads)
    {
        if (std::abs(static_cast<int>(read.size()) - p.target_len) > p.max_len_delta)
        {
            ++out.rejected_reads;
            continue;
        }
        Alignment forward = banded_align(read, ref, p.band_extra);
        if (!forward.valid)
        {
            ++out.rejected_reads;
            continue;
        }
        left_normalize_gaps(forward);
        const double norm = static_cast<double>(forward.distance) /
                            std::max(1, std::max(static_cast<int>(read.size()), m));
        if (norm > 0.42)
        {
            ++out.rejected_reads;
            continue;
        }
        const int full_weight = p.uniform_read_weight
                                    ? 256
                                    : read_weight(forward.distance,
                                                  static_cast<int>(read.size()), m);
        ++out.aligned_reads;
        norm_dist_sum += norm;

        accumulate_alignment_events(forward, m, full_weight, p.max_insertion,
                                    base_counts, ins_counts);
        total_evidence_weight += full_weight;
    }

    if (out.aligned_reads == 0)
        return out;
    out.mean_norm_dist = norm_dist_sum / out.aligned_reads;

    std::vector<std::vector<Choice>> steps;
    steps.reserve(2 * m + 1);
    const double alpha = 0.35;
    const double unit_vote = 256.0;

    for (int slot = 0; slot <= m; ++slot)
    {
        auto &counts = ins_counts[slot];
        if (counts.find("") == counts.end())
            counts[""] = 0;
        counts[""] += static_cast<int>(std::llround(p.no_insertion_prior * unit_vote));
        const int slot_total = total_evidence_weight +
                               static_cast<int>(std::llround(
                                   p.no_insertion_prior * unit_vote));
        std::vector<std::pair<std::string, int>> ranked(counts.begin(), counts.end());
        std::sort(ranked.begin(), ranked.end(), [](const auto &a, const auto &b)
                  {
                      if (a.second != b.second)
                          return a.second > b.second;
                      return a.first.size() < b.first.size();
                  });

        std::vector<Choice> opts;
        opts.reserve(p.insertion_choices + 1);
        int kept_nonempty = 0;
        bool kept_empty = false;
        for (const auto &kv : ranked)
        {
            if (kv.first.empty())
            {
                if (!kept_empty)
                {
                    double prob = (kv.second + alpha) /
                                  (slot_total + alpha * std::max<size_t>(1, counts.size()));
                    opts.push_back({kv.first, std::log(prob), prob});
                    kept_empty = true;
                }
            }
            else if (kept_nonempty < p.insertion_choices)
            {
                double prob = (kv.second + alpha) /
                              (slot_total + alpha * std::max<size_t>(1, counts.size()));
                double regularized = std::log(prob) -
                                     p.gap_event_penalty * kv.first.size();
                opts.push_back({kv.first, regularized, prob});
                ++kept_nonempty;
            }
        }
        if (!kept_empty)
            opts.push_back({"", std::log(alpha / (slot_total + alpha)), 0.0});
        steps.push_back(std::move(opts));

        if (slot < m)
        {
            std::vector<Choice> bopts;
            bopts.reserve(5);
            int position_total = 0;
            for (int b = 0; b < 5; ++b)
                position_total += base_counts[slot][b];
            int best_non_gap = 0;
            for (int b = 1; b < 4; ++b)
            {
                if (base_counts[slot][b] > base_counts[slot][best_non_gap])
                    best_non_gap = b;
            }
            const int no_del_prior = static_cast<int>(std::llround(
                p.no_deletion_prior * unit_vote));
            base_counts[slot][best_non_gap] += no_del_prior;
            position_total += no_del_prior;
            int max_non_gap_count = 0;
            for (int b = 0; b < 4; ++b)
                max_non_gap_count = std::max(max_non_gap_count, base_counts[slot][b]);
            const double base_conf = static_cast<double>(max_non_gap_count) /
                                     std::max(1, position_total);
            if (base_conf < 0.72)
                ++out.weak_positions;
            for (int b = 0; b < 5; ++b)
            {
                if (base_counts[slot][b] == 0 && b != base_id(ref[slot]))
                    continue;
                const double prob = (base_counts[slot][b] + alpha) /
                                    (position_total + 5.0 * alpha);
                double score = std::log(prob);
                if (b == 4)
                    score -= p.gap_event_penalty;
                bopts.push_back({b == 4 ? "" : std::string(1, id_base(b)),
                                 score, prob});
            }
            steps.push_back(std::move(bopts));
        }
    }

    const int L = p.target_len;
    const int S = static_cast<int>(steps.size());
    const double NEG = -1e100;
    std::vector<double> prev(L + 1, NEG), cur(L + 1, NEG);
    std::vector<int16_t> parent_len(static_cast<size_t>(S + 1) * (L + 1), -1);
    std::vector<int8_t> parent_opt(static_cast<size_t>(S + 1) * (L + 1), -1);
    auto state_index = [L](int s, int len) -> size_t
    {
        return static_cast<size_t>(s) * (L + 1) + len;
    };
    prev[0] = 0.0;
    for (int s = 0; s < S; ++s)
    {
        std::fill(cur.begin(), cur.end(), NEG);
        for (int len = 0; len <= L; ++len)
        {
            if (prev[len] <= NEG / 2)
                continue;
            for (int o = 0; o < static_cast<int>(steps[s].size()); ++o)
            {
                int nl = len + static_cast<int>(steps[s][o].emit.size());
                if (nl > L)
                    continue;
                double score = prev[len] + steps[s][o].score;
                // After consuming a reference base, couple the path to the
                // expected global phase. This is a soft penalty, not a hard
                // diagonal band, so genuine separated indels remain possible.
                if ((s & 1) == 1)
                {
                    const int consumed_ref = (s + 1) / 2;
                    const double expected_phase =
                        static_cast<double>(consumed_ref) * L / std::max(1, m);
                    score -= p.phase_drift_penalty *
                             std::abs(static_cast<double>(nl) - expected_phase);
                }
                if (score > cur[nl])
                {
                    cur[nl] = score;
                    parent_len[state_index(s + 1, nl)] = static_cast<int16_t>(len);
                    parent_opt[state_index(s + 1, nl)] = static_cast<int8_t>(o);
                }
            }
        }
        std::swap(prev, cur);
    }

    if (prev[L] <= NEG / 2)
        return out;
    std::vector<std::string> reversed;
    reversed.reserve(S);
    double conf_sum = 0.0;
    int len = L;
    for (int s = S; s > 0; --s)
    {
        int o = parent_opt[state_index(s, len)];
        int pl = parent_len[state_index(s, len)];
        if (o < 0 || pl < 0)
        {
            out.sequence.clear();
            return out;
        }
        reversed.push_back(steps[s - 1][o].emit);
        conf_sum += steps[s - 1][o].confidence;
        len = pl;
    }
    std::reverse(reversed.begin(), reversed.end());
    out.sequence.reserve(L);
    for (const auto &x : reversed)
        out.sequence += x;
    out.mean_chosen_conf = conf_sum / std::max(1, S);
    return out;
}

static inline std::vector<std::string> choose_read_subset(const std::vector<std::string> &reads,
                                                          const Params &p,
                                                          int seed)
{
    if (static_cast<int>(reads.size()) <= p.max_reads)
        return reads;
    std::vector<int> order(reads.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int x, int y)
              {
                  int dx = std::abs(static_cast<int>(reads[x].size()) - p.target_len);
                  int dy = std::abs(static_cast<int>(reads[y].size()) - p.target_len);
                  if (x == seed)
                      dx -= 1000;
                  if (y == seed)
                      dy -= 1000;
                  return dx < dy;
              });
    std::vector<std::string> selected;
    selected.reserve(p.max_reads);
    for (int i = 0; i < p.max_reads; ++i)
        selected.push_back(reads[order[i]]);
    return selected;
}

static inline std::string positional_phase_chain_consensus(
    const std::vector<std::string> &reads,
    const std::string &reference,
    const Params &p)
{
    const int L = p.target_len;
    const int k = p.positional_k;
    const int layers_count = L - k + 1;
    if (!p.positional_phase_chain ||
        static_cast<int>(reads.size()) < p.positional_min_reads ||
        k < 3 || k > 12 || layers_count <= 1)
        return reference;

    std::vector<std::unordered_map<uint32_t, int>> evidence(layers_count);
    for (const auto &read : reads)
    {
        const int read_layers = static_cast<int>(read.size()) - k + 1;
        if (read_layers <= 1)
            continue;
        for (int i = 0; i < read_layers; ++i)
        {
            const int projected = static_cast<int>(std::llround(
                static_cast<double>(i) * (layers_count - 1) /
                (read_layers - 1)));
            const uint32_t code = kmer_code(read, i, k);
            for (int delta = -p.positional_radius;
                 delta <= p.positional_radius; ++delta)
            {
                const int pos = projected + delta;
                if (pos < 0 || pos >= layers_count)
                    continue;
                evidence[pos][code] +=
                    p.positional_radius + 1 - std::abs(delta);
            }
        }
    }
    if (evidence.front().empty() || evidence.back().empty())
        return reference;

    std::unordered_map<uint32_t, double> previous;
    previous.reserve(evidence.front().size());
    const uint32_t ref_first =
        static_cast<int>(reference.size()) >= k
            ? kmer_code(reference, 0, k)
            : std::numeric_limits<uint32_t>::max();
    for (const auto &entry : evidence.front())
    {
        double score = std::log1p(static_cast<double>(entry.second));
        if (entry.first == ref_first)
            score += p.positional_reference_bonus;
        previous[entry.first] = score;
    }

    std::vector<std::unordered_map<uint32_t, uint32_t>> parents(layers_count);
    for (int pos = 1; pos < layers_count; ++pos)
    {
        std::unordered_map<uint32_t, double> current;
        current.reserve(evidence[pos].size());
        parents[pos].reserve(evidence[pos].size());
        const int ref_base =
            pos + k - 1 < static_cast<int>(reference.size())
                ? base_id(reference[pos + k - 1])
                : -1;
        for (const auto &entry : evidence[pos])
        {
            const uint32_t code = entry.first;
            const uint32_t prefix = code >> 2;
            double best = -1e100;
            uint32_t best_parent = 0;
            bool found = false;
            for (uint32_t first = 0; first < 4; ++first)
            {
                const uint32_t predecessor =
                    (first << (2 * (k - 1))) | prefix;
                auto it = previous.find(predecessor);
                if (it != previous.end() && it->second > best)
                {
                    best = it->second;
                    best_parent = predecessor;
                    found = true;
                }
            }
            if (!found)
                continue;
            double score = best + std::log1p(static_cast<double>(entry.second));
            if (static_cast<int>(code & 3U) == ref_base)
                score += p.positional_reference_bonus;
            current[code] = score;
            parents[pos][code] = best_parent;
        }
        if (current.empty())
            return reference;
        previous = std::move(current);
    }

    auto best_it = std::max_element(
        previous.begin(), previous.end(),
        [](const auto &a, const auto &b)
        { return a.second < b.second; });
    if (best_it == previous.end())
        return reference;
    std::vector<uint32_t> path(layers_count);
    path.back() = best_it->first;
    for (int pos = layers_count - 1; pos > 0; --pos)
    {
        auto it = parents[pos].find(path[pos]);
        if (it == parents[pos].end())
            return reference;
        path[pos - 1] = it->second;
    }
    std::string sequence(k, 'A');
    uint32_t first = path.front();
    for (int i = k - 1; i >= 0; --i)
    {
        sequence[i] = id_base(static_cast<int>(first & 3U));
        first >>= 2;
    }
    sequence.reserve(L);
    for (int pos = 1; pos < layers_count; ++pos)
        sequence.push_back(id_base(static_cast<int>(path[pos] & 3U)));
    return sequence;
}

// Stage 1: construct one fixed-length local coordinate guide. High-coverage
// clusters use the positional k-mer chain; the others use one banded-evidence
// DP. The output is never fed back for re-alignment.
static inline Result build_local_coordinate_once(
    const std::vector<std::string> &raw_reads,
    int cluster_id,
    const Params &p)
{
    Result result;
    result.diag.cluster_id = cluster_id;
    std::vector<std::string> reads;
    reads.reserve(raw_reads.size());
    for (const auto &raw : raw_reads)
    {
        std::string read = normalize_seq(raw);
        if (!read.empty())
            reads.push_back(std::move(read));
    }
    result.diag.num_reads = static_cast<int>(reads.size());
    if (reads.empty())
        return result;

    const int seed = choose_seed(reads, p);
    result.diag.seed_read_index = seed;
    result.diag.seed_len = static_cast<int>(reads[seed].size());
    int coherence_filtered = 0;
    double coherence_threshold = 1.0;
    std::vector<std::string> coherent = coherent_read_core(
        reads, seed, p, &coherence_filtered, &coherence_threshold);
    int coherent_seed = 0;
    for (int i = 0; i < static_cast<int>(coherent.size()); ++i)
        if (coherent[i] == reads[seed])
        {
            coherent_seed = i;
            break;
        }
    const std::vector<std::string> selected =
        choose_read_subset(coherent, p, coherent_seed);
    PassResult pass;
    if (static_cast<int>(reads.size()) >= p.coverage_switch_reads)
    {
        Params positional_params = p;
        positional_params.positional_phase_chain = true;
        positional_params.positional_min_reads = 6;
        positional_params.positional_reference_bonus = 0.0;
        result.sequence = positional_phase_chain_consensus(
            reads, reads[seed], positional_params);
        pass.aligned_reads = static_cast<int>(selected.size());
        pass.mean_norm_dist = 0.0;
        pass.mean_chosen_conf = 1.0;
        pass.weak_positions = 1;
    }
    else
    {
        pass = consensus_pass(selected, reads[seed], p);
        result.sequence = pass.sequence;
    }
    if (result.sequence.empty())
    {
        result.sequence = reads[seed];
        if (static_cast<int>(result.sequence.size()) > p.target_len)
            result.sequence.resize(p.target_len);
        while (static_cast<int>(result.sequence.size()) < p.target_len)
            result.sequence.push_back('A');
    }
    result.diag.used_reads = static_cast<int>(selected.size());
    result.diag.coherence_filtered = coherence_filtered;
    result.diag.coherence_threshold = coherence_threshold;
    result.diag.passes = 1;
    result.diag.aligned_reads = pass.aligned_reads;
    result.diag.rejected_reads = pass.rejected_reads;
    result.diag.mean_norm_dist = pass.mean_norm_dist;
    result.diag.mean_chosen_conf = pass.mean_chosen_conf;
    result.diag.weak_positions = pass.weak_positions;
    result.diag.final_len = static_cast<int>(result.sequence.size());
    result.diag.seed_to_final_edits =
        edit_distance_linear(reads[seed], result.sequence);
    return result;
}

struct CrossClusterReuseModel
{
    int k = 11;
    int edge_len = 12;
    int min_support = 2;
    std::vector<uint16_t> counts;
};

static inline uint32_t encode_reuse_edge(const std::string &s,
                                         int pos,
                                         int edge_len)
{
    uint32_t code = 0;
    for (int i = 0; i < edge_len; ++i)
        code = (code << 2) | static_cast<uint32_t>(base_id(s[pos + i]));
    return code;
}

static inline CrossClusterReuseModel build_cross_cluster_reuse_model(
    const std::vector<std::string> &sequences,
    const Params &p)
{
    CrossClusterReuseModel model;
    model.k = p.topology_k;
    model.edge_len = p.topology_k + 1;
    model.min_support = p.topology_min_support;
    if (model.edge_len <= 0 || model.edge_len > 13)
        return model;
    const size_t universe = static_cast<size_t>(1) << (2 * model.edge_len);
    const uint32_t edge_mask =
        (static_cast<uint32_t>(1) << (2 * model.edge_len)) - 1U;
    model.counts.assign(universe, 0);
    for (const auto &s : sequences)
    {
        if (static_cast<int>(s.size()) < model.edge_len)
            continue;
        // One guide contributes at most one vote to each edge. Repeated
        // occurrences inside the same guide are not cross-cluster reuse.
        const int windows =
            static_cast<int>(s.size()) - model.edge_len + 1;
        if (windows <= 192)
        {
            // Allocation-free per-guide set for the short strands used by
            // this decoder. Load is at most 0.5, so linear probing is cheap.
            std::array<uint32_t, 256> seen;
            seen.fill(std::numeric_limits<uint32_t>::max());
            uint32_t code = encode_reuse_edge(
                s, 0, model.edge_len);
            for (int i = 0; i < windows; ++i)
            {
                if (i > 0)
                    code = ((code << 2) & edge_mask) |
                           static_cast<uint32_t>(base_id(
                               s[i + model.edge_len - 1]));
                size_t slot =
                    (static_cast<size_t>(code) * 2654435761U) & 255U;
                while (seen[slot] !=
                           std::numeric_limits<uint32_t>::max() &&
                       seen[slot] != code)
                    slot = (slot + 1) & 255U;
                if (seen[slot] == code)
                    continue;
                seen[slot] = code;
                uint16_t &count = model.counts[code];
                if (count != std::numeric_limits<uint16_t>::max())
                    ++count;
            }
        }
        else
        {
            // General fallback for unusually long targets.
            std::vector<uint32_t> edges_in_guide;
            edges_in_guide.reserve(windows);
            uint32_t code = encode_reuse_edge(
                s, 0, model.edge_len);
            for (int i = 0; i < windows; ++i)
            {
                if (i > 0)
                    code = ((code << 2) & edge_mask) |
                           static_cast<uint32_t>(base_id(
                               s[i + model.edge_len - 1]));
                edges_in_guide.push_back(code);
            }
            std::sort(edges_in_guide.begin(), edges_in_guide.end());
            edges_in_guide.erase(
                std::unique(edges_in_guide.begin(),
                            edges_in_guide.end()),
                edges_in_guide.end());
            for (uint32_t code : edges_in_guide)
            {
                uint16_t &count = model.counts[code];
                if (count != std::numeric_limits<uint16_t>::max())
                    ++count;
            }
        }
    }
    return model;
}

// Fraction of sequence windows whose edge is observed in at least
// min_support distinct guides. For min_support=2, a low value means that
// inter-cluster transitions are mostly unique and should receive little or no
// cross-cluster weight.
static inline double cross_cluster_reuse_fraction(
    const std::vector<std::string> &sequences,
    const CrossClusterReuseModel &model)
{
    uint64_t supported = 0;
    uint64_t total = 0;
    if (model.counts.empty())
        return 0.0;
    const uint32_t edge_mask =
        (static_cast<uint32_t>(1) << (2 * model.edge_len)) - 1U;
    for (const auto &s : sequences)
    {
        if (static_cast<int>(s.size()) < model.edge_len)
            continue;
        uint32_t code = encode_reuse_edge(s, 0, model.edge_len);
        for (int pos = 0;
             pos + model.edge_len <= static_cast<int>(s.size()); ++pos)
        {
            if (pos > 0)
                code = ((code << 2) & edge_mask) |
                       static_cast<uint32_t>(base_id(
                           s[pos + model.edge_len - 1]));
            ++total;
            if (model.counts[code] >= model.min_support)
                ++supported;
        }
    }
    return total == 0
               ? 0.0
               : static_cast<double>(supported) /
                     static_cast<double>(total);
}

struct LooMarkovModel
{
    int k = 8;
    std::vector<uint16_t> edge_counts;
    std::vector<uint32_t> node_counts;
    std::vector<uint32_t> suffix_node_counts;
};

static inline LooMarkovModel build_loo_markov_model(
    const std::vector<std::string> &sequences,
    int k)
{
    LooMarkovModel model;
    model.k = k;
    if (k < 2 || k > 11)
        return model;
    model.edge_counts.assign(
        static_cast<size_t>(1) << (2 * (k + 1)), 0);
    model.node_counts.assign(
        static_cast<size_t>(1) << (2 * k), 0);
    model.suffix_node_counts.assign(
        static_cast<size_t>(1) << (2 * k), 0);
    const uint32_t node_mask =
        (static_cast<uint32_t>(1) << (2 * k)) - 1U;
    for (const auto &sequence : sequences)
    {
        for (int pos = 0;
             pos + k < static_cast<int>(sequence.size()); ++pos)
        {
            const uint32_t edge =
                kmer_code(sequence, pos, k + 1);
            const uint32_t node = edge >> 2;
            uint16_t &edge_count = model.edge_counts[edge];
            if (edge_count != std::numeric_limits<uint16_t>::max())
                ++edge_count;
            if (model.node_counts[node] !=
                std::numeric_limits<uint32_t>::max())
                ++model.node_counts[node];
            const uint32_t suffix_node = edge & node_mask;
            if (model.suffix_node_counts[suffix_node] !=
                std::numeric_limits<uint32_t>::max())
                ++model.suffix_node_counts[suffix_node];
        }
    }
    return model;
}

// Dataset-level reliability of the cross-cluster Markov field.  Calibration
// uses only frozen provisional guides.  It never sees the true centers or any
// sequence produced by the final trellis.
struct CrossClusterCalibration
{
    double trust = 0.0;
    double train_gain = 0.0;
    double validation_gain = 0.0;
    double validation_se = 0.0;
    int sampled_guides = 0;
    int sampled_windows = 0;
    bool accepted = false;
};

// Select a shrinkage factor in [0,1] by deterministic train/validation
// splitting.  For trust=t, each leave-one-out transition probability is
// mixed with the uniform DNA probability: (1-t)/4 + t*p_LOO.  The selected
// factor is accepted only if its validation gain is positive by at least two
// standard errors.  Otherwise trust is exactly zero, which is the safe local
// evidence fallback for libraries without reusable cross-strand structure.
static inline CrossClusterCalibration calibrate_cross_cluster_trust(
    const std::vector<std::string> &sequences,
    const LooMarkovModel &model,
    const Params &p,
    int maximum_sample_guides = 512)
{
    CrossClusterCalibration result;
    const int n = static_cast<int>(sequences.size());
    const int k = model.k;
    if (n < 16 || k < 2 || model.edge_counts.empty())
        return result;

    const int sample_count = std::min(n, std::max(16, maximum_sample_guides));
    result.sampled_guides = sample_count;
    const uint32_t node_mask =
        (static_cast<uint32_t>(1) << (2 * k)) - 1U;
    constexpr double alpha = 0.5;
    constexpr std::array<double, 9> candidates{
        0.0, 0.125, 0.25, 0.375, 0.50,
        0.625, 0.75, 0.875, 1.0};

    std::array<double, candidates.size()> train_sum{};
    std::array<int, candidates.size()> train_guides{};
    std::vector<std::array<double, candidates.size()>> validation_gains;
    validation_gains.reserve(sample_count / 2);

    auto count_of = [](const auto &table, uint32_t code) -> int
    {
        const auto found = table.find(code);
        return found == table.end() ? 0 : found->second;
    };

    for (int sample = 0; sample < sample_count; ++sample)
    {
        const int sequence_index = static_cast<int>(
            static_cast<long long>(sample) * n / sample_count);
        const std::string &sequence = sequences[sequence_index];
        const int windows = static_cast<int>(sequence.size()) - k;
        if (windows <= 0)
            continue;

        std::unordered_map<uint32_t, int> own_edges;
        std::unordered_map<uint32_t, int> own_nodes;
        std::unordered_map<uint32_t, int> own_suffix_nodes;
        own_edges.reserve(windows * 2);
        own_nodes.reserve(windows * 2);
        own_suffix_nodes.reserve(windows * 2);
        for (int pos = 0; pos < windows; ++pos)
        {
            const uint32_t edge = kmer_code(sequence, pos, k + 1);
            ++own_edges[edge];
            ++own_nodes[edge >> 2];
            ++own_suffix_nodes[edge & node_mask];
        }

        std::array<double, candidates.size()> guide_gain{};
        int valid_windows = 0;
        for (int pos = 0; pos < windows; ++pos)
        {
            const uint32_t edge = kmer_code(sequence, pos, k + 1);
            const uint32_t prefix = edge >> 2;
            const uint32_t suffix = edge & node_mask;
            const int edge_count = std::max(
                0, static_cast<int>(model.edge_counts[edge]) -
                       count_of(own_edges, edge));
            const int prefix_count = std::max(
                0, static_cast<int>(model.node_counts[prefix]) -
                       count_of(own_nodes, prefix));
            const int suffix_count = std::max(
                0, static_cast<int>(model.suffix_node_counts[suffix]) -
                       count_of(own_suffix_nodes, suffix));
            const double forward_probability =
                (edge_count + alpha) /
                (prefix_count + 4.0 * alpha);
            const double reverse_probability =
                (edge_count + alpha) /
                (suffix_count + 4.0 * alpha);
            for (size_t candidate = 0;
                 candidate < candidates.size(); ++candidate)
            {
                const double trust = candidates[candidate];
                const double mixed_forward =
                    (1.0 - trust) * 0.25 +
                    trust * forward_probability;
                double gain = std::log(mixed_forward / 0.25);
                if (p.loo_reverse_weight > 0.0)
                {
                    const double mixed_reverse =
                        (1.0 - trust) * 0.25 +
                        trust * reverse_probability;
                    gain += p.loo_reverse_weight *
                            std::log(mixed_reverse / 0.25);
                }
                guide_gain[candidate] += gain;
            }
            ++valid_windows;
        }
        if (valid_windows == 0)
            continue;
        result.sampled_windows += valid_windows;
        for (double &gain : guide_gain)
            gain /= valid_windows;
        if ((sample & 1) == 0)
            for (size_t candidate = 0;
                 candidate < candidates.size(); ++candidate)
            {
                train_sum[candidate] += guide_gain[candidate];
                ++train_guides[candidate];
            }
        else
            validation_gains.push_back(guide_gain);
    }

    size_t best_candidate = 0;
    double best_train_gain = 0.0;
    for (size_t candidate = 1; candidate < candidates.size(); ++candidate)
    {
        if (train_guides[candidate] == 0)
            continue;
        const double mean = train_sum[candidate] /
                            train_guides[candidate];
        // A strict improvement keeps ties at the simpler, smaller trust.
        if (mean > best_train_gain + 1e-12)
        {
            best_train_gain = mean;
            best_candidate = candidate;
        }
    }
    result.train_gain = best_train_gain;
    if (best_candidate == 0 || validation_gains.size() < 8)
        return result;

    std::vector<double> selected_validation;
    selected_validation.reserve(validation_gains.size());
    for (const auto &gain : validation_gains)
        selected_validation.push_back(gain[best_candidate]);
    result.validation_gain = std::accumulate(
        selected_validation.begin(), selected_validation.end(), 0.0) /
        selected_validation.size();
    double squared_deviation = 0.0;
    for (double gain : selected_validation)
    {
        const double difference = gain - result.validation_gain;
        squared_deviation += difference * difference;
    }
    if (selected_validation.size() > 1)
        result.validation_se = std::sqrt(
            squared_deviation /
            (selected_validation.size() - 1) /
            selected_validation.size());

    const double acceptance_margin = std::max(
        0.005, 2.0 * result.validation_se);
    result.accepted =
        result.validation_gain > acceptance_margin;
    result.trust = result.accepted
                       ? candidates[best_candidate]
                       : 0.0;
    return result;
}

struct JointPhaseTrellisEvidence
{
    int target_len = 0;
    int reference_len = 0;
    int cross_k = 0;
    int local_k = 0;
    int positional_k = 0;
    double cross_weight = 1.0;
    std::string anchor;
    ChannelProfile channel;
    const LooMarkovModel *cross_model = nullptr;
    std::vector<std::array<double, 5>> reference_action_score;
    std::vector<std::array<double, 4>> insertion_action_score;
    std::vector<std::array<double, 4>> insertion_continuation_score;
    std::vector<std::array<double, 2>> insertion_stop_score;
    std::unordered_map<uint32_t, int> own_cross_edges;
    std::unordered_map<uint32_t, int> own_cross_nodes;
    std::unordered_map<uint32_t, int> own_cross_suffix_nodes;
    std::vector<double> local_edge_log_score;
    std::vector<std::unordered_map<uint32_t, double>> positional_log_score;
};

// Construct the complete, immutable scoring field for one phase trellis.
// Reads are aligned only to the frozen anchor. No decoded sequence is fed
// back into this function.
static inline JointPhaseTrellisEvidence build_joint_phase_trellis_evidence(
    const std::string &anchor,
    const std::vector<std::string> &raw_reads,
    const LooMarkovModel &model,
    const ChannelProfile &channel,
    const Params &p,
    double cross_weight = 1.0)
{
    JointPhaseTrellisEvidence field;
    const int R = static_cast<int>(anchor.size());
    const int L = p.target_len;
    field.target_len = L;
    field.reference_len = R;
    field.cross_k = model.k;
    field.local_k = p.loo_local_k;
    field.positional_k = p.positional_k;
    field.cross_weight = std::max(0.0, cross_weight);
    field.anchor = anchor;
    field.channel = channel;
    field.cross_model = &model;
    if (R <= 0 || L <= 0 || model.edge_counts.empty())
        return field;

    std::vector<std::string> reads;
    reads.reserve(raw_reads.size());
    for (const auto &raw : raw_reads)
    {
        std::string read = normalize_seq(raw);
        if (!read.empty())
            reads.push_back(std::move(read));
    }
    if (reads.empty())
        return field;
    const int seed = choose_seed(reads, p);
    int coherence_filtered = 0;
    double coherence_threshold = 1.0;
    std::vector<std::string> coherent = coherent_read_core(
        reads, seed, p, &coherence_filtered, &coherence_threshold);
    int coherent_seed = 0;
    for (int i = 0; i < static_cast<int>(coherent.size()); ++i)
        if (coherent[i] == reads[seed])
        {
            coherent_seed = i;
            break;
        }
    const std::vector<std::string> selected =
        choose_read_subset(coherent, p, coherent_seed);

    std::vector<std::array<int, 5>> base_counts(R);
    for (auto &row : base_counts)
        row.fill(0);
    std::vector<std::unordered_map<std::string, int>> insertion_counts(R + 1);
    int total_weight = 0;
    for (const auto &read : selected)
    {
        Alignment alignment = banded_align(read, anchor, p.band_extra);
        if (!alignment.valid)
            continue;
        left_normalize_gaps(alignment);
        const double normalized_distance =
            static_cast<double>(alignment.distance) /
            std::max(1, std::max<int>(read.size(), anchor.size()));
        if (normalized_distance > 0.42)
            continue;
        const int weight = p.uniform_read_weight
                               ? 256
                               : read_weight(
            alignment.distance, read.size(), R);
        accumulate_alignment_events(
            alignment, R, weight, p.max_insertion,
            base_counts, insertion_counts);
        total_weight += weight;
    }
    if (total_weight <= 0)
        return field;

    constexpr double alpha = 0.5;
    field.reference_action_score.resize(R);
    for (int pos = 0; pos < R; ++pos)
    {
        double total = 5.0 * alpha;
        for (int count : base_counts[pos])
            total += count;
        for (int action = 0; action < 5; ++action)
        {
            const bool anchor_action =
                action == base_id(anchor[pos]);
            field.reference_action_score[pos][action] =
                base_counts[pos][action] == 0 && !anchor_action
                    ? -1e100
                    : std::log(
                          (base_counts[pos][action] + alpha) / total);
        }
    }
    field.insertion_action_score.resize(R + 1);
    field.insertion_continuation_score.resize(R + 1);
    field.insertion_stop_score.resize(R + 1);
    for (int slot = 0; slot <= R; ++slot)
    {
        std::array<double, 4> counts{alpha, alpha, alpha, alpha};
        double total = 4.0 * alpha;
        double insertion_events = 0.0;
        double continuation_events = 0.0;
        for (const auto &entry : insertion_counts[slot])
        {
            if (entry.first.empty())
                continue;
            insertion_events += entry.second;
            continuation_events +=
                entry.second * std::max<int>(
                                   0, entry.first.size() - 1);
            for (char base_char : entry.first)
            {
                const int base = base_id(base_char);
                if (base >= 0)
                {
                    counts[base] += entry.second;
                    total += entry.second;
                }
            }
        }
        const double no_insertion_events = std::max(
            0.0, static_cast<double>(total_weight) - insertion_events);
        const double event_probability =
            (insertion_events + alpha) /
            (insertion_events + no_insertion_events + 2.0 * alpha);
        const double continuation_probability =
            (continuation_events + alpha) /
            (insertion_events + continuation_events + 2.0 * alpha);
        field.insertion_stop_score[slot][0] =
            std::log(1.0 - event_probability);
        field.insertion_stop_score[slot][1] =
            std::log(1.0 - continuation_probability);
        for (int base = 0; base < 4; ++base)
        {
            const double base_score = std::log(counts[base] / total);
            field.insertion_action_score[slot][base] =
                std::log(event_probability) + base_score;
            field.insertion_continuation_score[slot][base] =
                std::log(continuation_probability) + base_score;
        }
    }

    const int cross_k = field.cross_k;
    const uint32_t cross_node_mask =
        (static_cast<uint32_t>(1) << (2 * cross_k)) - 1U;
    for (int pos = 0; pos + cross_k < R; ++pos)
    {
        const uint32_t edge = kmer_code(anchor, pos, cross_k + 1);
        ++field.own_cross_edges[edge];
        ++field.own_cross_nodes[edge >> 2];
        ++field.own_cross_suffix_nodes[edge & cross_node_mask];
    }
    const int local_k = field.local_k;
    const size_t local_edge_count =
        static_cast<size_t>(1) << (2 * (local_k + 1));
    const size_t local_node_count =
        static_cast<size_t>(1) << (2 * local_k);
    std::vector<uint32_t> local_edges(local_edge_count, 0);
    std::vector<uint32_t> local_nodes(local_node_count, 0);
    std::vector<uint32_t> local_suffix_nodes(local_node_count, 0);
    const uint32_t local_node_mask =
        static_cast<uint32_t>(local_node_count - 1);
    for (const auto &read : reads)
        for (int pos = 0;
             pos + local_k < static_cast<int>(read.size()); ++pos)
        {
            const uint32_t edge = kmer_code(read, pos, local_k + 1);
            ++local_edges[edge];
            ++local_nodes[edge >> 2];
            ++local_suffix_nodes[edge & local_node_mask];
        }
    field.local_edge_log_score.resize(local_edge_count);
    for (size_t edge_index = 0; edge_index < local_edge_count; ++edge_index)
    {
        const uint32_t edge = static_cast<uint32_t>(edge_index);
        const uint32_t prefix = edge >> 2;
        field.local_edge_log_score[edge] = std::log(
            (local_edges[edge] + alpha) /
            (local_nodes[prefix] + 4.0 * alpha));
        field.local_edge_log_score[edge] +=
            p.loo_reverse_weight * std::log(
                (local_edges[edge] + alpha) /
                (local_suffix_nodes[edge & local_node_mask] + 4.0 * alpha));
    }

    const int positional_k = field.positional_k;
    const int positional_layers = L - positional_k + 1;
    field.positional_log_score.resize(std::max(0, positional_layers));
    std::vector<std::unordered_map<uint32_t, double>>
        positional_counts(std::max(0, positional_layers));
    if (positional_layers > 1)
        for (const auto &read : reads)
        {
            const int read_layers =
                static_cast<int>(read.size()) - positional_k + 1;
            if (read_layers <= 1)
                continue;
            for (int read_pos = 0; read_pos < read_layers; ++read_pos)
            {
                const int projected = static_cast<int>(std::llround(
                    static_cast<double>(read_pos) *
                    (positional_layers - 1) / (read_layers - 1)));
                const uint32_t code =
                    kmer_code(read, read_pos, positional_k);
                for (int delta = -p.positional_radius;
                     delta <= p.positional_radius; ++delta)
                {
                    const int layer = projected + delta;
                    if (layer < 0 || layer >= positional_layers)
                        continue;
                    positional_counts[layer][code] +=
                        positional_k - std::abs(delta);
                }
            }
        }
    for (int layer = 0; layer < positional_layers; ++layer)
    {
        double maximum = 0.0;
        for (const auto &entry : positional_counts[layer])
            maximum = std::max(maximum, std::log1p(entry.second));
        auto &scores = field.positional_log_score[layer];
        scores.reserve(positional_counts[layer].size());
        for (const auto &entry : positional_counts[layer])
            scores[entry.first] = std::log1p(entry.second) - maximum;
    }
    return field;
}

// One banded Viterbi decode over (reference position, output length, k-mer
// context). Match/substitution, insertion, and deletion transitions all share
// the same path and simultaneously consume the four evidence terms.
static inline std::string decode_joint_phase_trellis(
    const JointPhaseTrellisEvidence &field,
    const Params &p,
    int phase_band = 6,
    int context_states_per_cell = 8)
{
    const int L = field.target_len;
    const int R = field.reference_len;
    const int k = field.cross_k;
    if (L <= 0 || R <= 0 ||
        static_cast<int>(field.anchor.size()) != R ||
        field.reference_action_score.size() != static_cast<size_t>(R) ||
        field.cross_model == nullptr ||
        field.cross_model->edge_counts.empty() || k < field.local_k)
        return field.anchor;

    struct Node
    {
        double score = -1e100;
        uint32_t suffix = 0;
        int parent = -1;
        char emit = '\0';
        bool inserted = false;
    };
    std::vector<Node> nodes;
    nodes.reserve(static_cast<size_t>(R) * (2 * phase_band + 1) *
                  context_states_per_cell);
    nodes.push_back({0.0, 0U, -1, '\0', false});
    std::vector<std::vector<int>> cells(
        static_cast<size_t>(R + 1) * (L + 1));
    auto cell = [L](int reference_pos, int output_len) -> size_t
    {
        return static_cast<size_t>(reference_pos) * (L + 1) + output_len;
    };
    cells[cell(0, 0)].push_back(0);
    const uint32_t suffix_mask =
        (static_cast<uint32_t>(1) << (2 * k)) - 1U;
    const uint32_t local_edge_mask =
        (static_cast<uint32_t>(1) << (2 * (field.local_k + 1))) - 1U;
    const uint32_t positional_edge_mask =
        (static_cast<uint32_t>(1) << (2 * field.positional_k)) - 1U;
    const double epsilon = 1e-9;
    const uint32_t cross_node_mask = suffix_mask;
    std::unordered_map<uint32_t, double> cross_score_cache;
    auto cross_edge_score = [&](uint32_t edge) -> double
    {
        const auto cached = cross_score_cache.find(edge);
        if (cached != cross_score_cache.end())
            return cached->second;
        auto own_count = [](const auto &table, uint32_t code) -> int
        {
            const auto found = table.find(code);
            return found == table.end() ? 0 : found->second;
        };
        const auto &model = *field.cross_model;
        const uint32_t prefix = edge >> 2;
        const uint32_t suffix = edge & cross_node_mask;
        const int edge_count = std::max(
            0, static_cast<int>(model.edge_counts[edge]) -
                   own_count(field.own_cross_edges, edge));
        const int prefix_count = std::max(
            0, static_cast<int>(model.node_counts[prefix]) -
                   own_count(field.own_cross_nodes, prefix));
        const int suffix_count = std::max(
            0, static_cast<int>(model.suffix_node_counts[suffix]) -
                   own_count(field.own_cross_suffix_nodes, suffix));
        constexpr double alpha = 0.5;
        double score = std::log(
            (edge_count + alpha) /
            (prefix_count + 4.0 * alpha));
        score += p.loo_reverse_weight * std::log(
            (edge_count + alpha) /
            (suffix_count + 4.0 * alpha));
        cross_score_cache.emplace(edge, score);
        return score;
    };
    auto push_state = [&](int reference_pos, int output_len,
                          double score, uint32_t suffix,
                          int parent, char emit, bool inserted)
    {
        const int expected_output = static_cast<int>(std::llround(
            static_cast<double>(reference_pos) * L / R));
        if (reference_pos < 0 || reference_pos > R ||
            output_len < 0 || output_len > L ||
            std::abs(output_len - expected_output) > phase_band)
            return;
        auto &target = cells[cell(reference_pos, output_len)];
        for (int node_index : target)
            if (nodes[node_index].suffix == suffix &&
                nodes[node_index].inserted == inserted)
            {
                if (score > nodes[node_index].score)
                    nodes[node_index] =
                        {score, suffix, parent, emit, inserted};
                return;
            }
        if (static_cast<int>(target.size()) < context_states_per_cell)
        {
            target.push_back(static_cast<int>(nodes.size()));
            nodes.push_back({score, suffix, parent, emit, inserted});
            return;
        }
        int worst_slot = 0;
        for (int slot = 1; slot < static_cast<int>(target.size()); ++slot)
            if (nodes[target[slot]].score <
                nodes[target[worst_slot]].score)
                worst_slot = slot;
        if (score > nodes[target[worst_slot]].score)
            nodes[target[worst_slot]] =
                {score, suffix, parent, emit, inserted};
    };
    auto emission_score = [&](uint32_t suffix, int output_len,
                              int base) -> double
    {
        double score = 0.0;
        const uint32_t extended =
            ((suffix << 2) | static_cast<uint32_t>(base));
        // When calibration rejects the cross-cluster field, its weight is
        // exactly zero.  Skip the hash lookups entirely instead of computing
        // a score that would then be multiplied by zero.
        if (output_len >= k && field.cross_weight > 0.0)
            score += field.cross_weight * cross_edge_score(
                extended & ((static_cast<uint32_t>(1) << (2 * (k + 1))) - 1U));
        if (output_len >= field.local_k)
            score += p.loo_local_weight *
                     field.local_edge_log_score[
                         extended & local_edge_mask];
        if (output_len + 1 >= field.positional_k)
        {
            const int layer = output_len + 1 - field.positional_k;
            const uint32_t code = extended & positional_edge_mask;
            const auto found = field.positional_log_score[layer].find(code);
            score += 0.10 *
                     (found == field.positional_log_score[layer].end()
                          ? -4.0
                          : found->second);
        }
        return score;
    };

    // Score the explicit all-KEEP path under the same objective. This is the
    // conservative path inside the trellis, not a second post-hoc candidate
    // test. It also protects KEEP from finite context-state pruning.
    double keep_path_score = -1e100;
    if (R == L)
    {
        keep_path_score = 0.0;
        uint32_t keep_suffix = 0U;
        for (int pos = 0; pos < R; ++pos)
        {
            const int base = base_id(field.anchor[pos]);
            keep_path_score +=
                field.insertion_stop_score[pos][0] +
                field.reference_action_score[pos][base] +
                emission_score(keep_suffix, pos, base) +
                0.10 * std::log(std::max(
                           epsilon, field.channel.match));
            keep_suffix =
                ((keep_suffix << 2) |
                 static_cast<uint32_t>(base)) &
                suffix_mask;
        }
    }

    for (int reference_pos = 0; reference_pos <= R; ++reference_pos)
    {
        const int expected_output = static_cast<int>(std::llround(
            static_cast<double>(reference_pos) * L / R));
        for (int output_len = std::max(0, expected_output - phase_band);
             output_len <= std::min(L, expected_output + phase_band);
             ++output_len)
        {
            const std::vector<int> current = cells[
                cell(reference_pos, output_len)];
            for (int node_index : current)
            {
                const Node node = nodes[node_index];
                if (reference_pos < R)
                {
                    const double deletion_score =
                        node.score +
                        field.insertion_stop_score[reference_pos]
                                                  [node.inserted ? 1 : 0] +
                        field.reference_action_score[reference_pos][4] +
                        0.10 * std::log(std::max(
                                   epsilon, field.channel.deletion));
                    push_state(reference_pos + 1, output_len,
                               deletion_score, node.suffix,
                               node_index, '\0', false);
                    if (output_len < L)
                        for (int base = 0; base < 4; ++base)
                        {
                            const bool matches =
                                base == base_id(field.anchor[reference_pos]);
                            const double channel_score = std::log(std::max(
                                epsilon,
                                matches ? field.channel.match
                                        : field.channel.substitution / 3.0));
                            const double score =
                                node.score +
                                field.insertion_stop_score[reference_pos]
                                                          [node.inserted ? 1 : 0] +
                                field.reference_action_score[reference_pos][base] +
                                emission_score(node.suffix, output_len, base) +
                                0.10 * channel_score;
                            const uint32_t suffix =
                                ((node.suffix << 2) |
                                 static_cast<uint32_t>(base)) &
                                suffix_mask;
                            push_state(reference_pos + 1, output_len + 1,
                                       score, suffix, node_index,
                                       id_base(base), false);
                        }
                }
                if (output_len < L)
                    for (int base = 0; base < 4; ++base)
                    {
                        const double score =
                            node.score +
                            (node.inserted
                                 ? field.insertion_continuation_score
                                       [reference_pos][base]
                                 : field.insertion_action_score
                                       [reference_pos][base]) +
                            emission_score(node.suffix, output_len, base) +
                            0.10 * std::log(std::max(
                                       epsilon,
                                       field.channel.insertion / 4.0));
                        const uint32_t suffix =
                            ((node.suffix << 2) |
                             static_cast<uint32_t>(base)) & suffix_mask;
                        push_state(reference_pos, output_len + 1,
                                   score, suffix, node_index,
                                   id_base(base), true);
                    }
            }
        }
    }
    const auto &terminal = cells[cell(R, L)];
    if (terminal.empty())
        return field.anchor;
    const int best_node = *std::max_element(
        terminal.begin(), terminal.end(),
        [&](int left, int right)
        { return nodes[left].score < nodes[right].score; });
    if (nodes[best_node].score <= keep_path_score)
        return field.anchor;
    std::string reversed;
    for (int node = best_node; node > 0; node = nodes[node].parent)
        if (nodes[node].emit != '\0')
            reversed.push_back(nodes[node].emit);
    if (static_cast<int>(reversed.size()) != L)
        return field.anchor;
    std::reverse(reversed.begin(), reversed.end());
    return reversed;
}

struct ConservedAnchorModel
{
    std::string left;
    std::string right;
};

// A routing predicate, not a decoder: identify the unique boundary action
// family from frozen anchors before entering an action subgraph.
static inline int shifted_boundary_delete_position(
    const std::string &sequence,
    const ConservedAnchorModel &anchors)
{
    const int L = static_cast<int>(sequence.size());
    const int left_length = static_cast<int>(anchors.left.size());
    const int right_length = static_cast<int>(anchors.right.size());
    if (left_length < 8 || right_length < 8 ||
        left_length + right_length + 4 >= L)
        return -1;

    const bool left_normal =
        sequence.compare(0, left_length, anchors.left) == 0;
    const bool right_normal =
        sequence.compare(L - right_length, right_length,
                         anchors.right) == 0;
    const bool left_shifted_right =
        sequence.compare(1, left_length, anchors.left) == 0;
    const bool right_shifted_left =
        sequence.compare(L - right_length - 1, right_length,
                         anchors.right) == 0;
    if (left_shifted_right && !left_normal && right_normal)
        return 0;
    if (right_shifted_left && !right_normal && left_normal)
        return L - 1;
    return -1;
}

// Learn boundary anchors from high-coverage provisional guides only.
// No primer sequence, payload boundary, center answer, or dataset name is
// supplied by the caller.
static inline ConservedAnchorModel build_conserved_anchor_model(
    const std::vector<std::string> &high_coverage_guides,
    double minimum_fraction = 0.985)
{
    ConservedAnchorModel model;
    if (high_coverage_guides.size() < 64)
        return model;
    const int L = static_cast<int>(high_coverage_guides.front().size());
    if (L < 24)
        return model;
    std::vector<char> consensus(L, 'N');
    std::vector<double> confidence(L, 0.0);
    for (int pos = 0; pos < L; ++pos)
    {
        std::array<int, 4> counts{0, 0, 0, 0};
        int observations = 0;
        for (const auto &sequence : high_coverage_guides)
        {
            if (static_cast<int>(sequence.size()) != L)
                continue;
            ++counts[base_id(sequence[pos])];
            ++observations;
        }
        const int best =
            static_cast<int>(std::max_element(
                                 counts.begin(), counts.end()) -
                             counts.begin());
        if (observations > 0)
        {
            static constexpr char alphabet[] = "ACGT";
            consensus[pos] = alphabet[best];
            confidence[pos] =
                static_cast<double>(counts[best]) / observations;
        }
    }
    const int maximum_anchor = L / 3;
    int left_length = 0;
    while (left_length < maximum_anchor &&
           confidence[left_length] >= minimum_fraction)
        ++left_length;
    int right_length = 0;
    while (right_length < maximum_anchor &&
           confidence[L - 1 - right_length] >= minimum_fraction)
        ++right_length;
    if (left_length >= 8)
        model.left.assign(consensus.begin(),
                          consensus.begin() + left_length);
    if (right_length >= 8)
        model.right.assign(consensus.end() - right_length,
                           consensus.end());
    return model;
}

// Repair a one-base phase excursion whose unmatched event was absorbed by a
// sequence boundary. The conserved anchor makes the boundary event
// identifiable; only the missing payload position/base remains latent and is
// selected by trimmed read energy.
static inline std::string decode_boundary_action_graph(
    const std::string &sequence,
    const std::vector<std::string> &reads,
    const ConservedAnchorModel &anchors,
    const Params &p,
    int *read_gain_out = nullptr)
{
    if (read_gain_out)
        *read_gain_out = 0;
    const int L = static_cast<int>(sequence.size());
    const int left_length = static_cast<int>(anchors.left.size());
    const int right_length = static_cast<int>(anchors.right.size());
    if (left_length < 8 || right_length < 8 ||
        left_length + right_length + 4 >= L)
        return sequence;

    // Exactly one shifted boundary is required. This keeps the decoder
    // dormant on datasets without a phase-boundary inconsistency.
    const int delete_position =
        shifted_boundary_delete_position(sequence, anchors);
    if (delete_position < 0)
        return sequence;

    std::string shortened = sequence;
    shortened.erase(shortened.begin() + delete_position);
    const int insertion_begin = left_length;
    const int insertion_end = L - right_length;
    const int shortened_length = static_cast<int>(shortened.size());
    std::vector<int> action_energy(
        static_cast<size_t>(shortened_length + 1) * 4, 0);
    int base_energy = 0;
    for (const auto &raw_read : reads)
    {
        const std::string read = normalize_seq(raw_read);
        if (read.empty())
            continue;
        base_energy += banded_distance_only(
            sequence, read, p.band_extra);
        const int read_length = static_cast<int>(read.size());
        const int stride = read_length + 1;
        std::vector<int> forward(
            static_cast<size_t>(shortened_length + 1) * stride);
        std::vector<int> backward(forward.size());
        auto cell = [stride](int i, int j)
        {
            return static_cast<size_t>(i) * stride + j;
        };
        for (int i = 0; i <= shortened_length; ++i)
            forward[cell(i, 0)] = i;
        for (int j = 0; j <= read_length; ++j)
            forward[cell(0, j)] = j;
        for (int i = 1; i <= shortened_length; ++i)
            for (int j = 1; j <= read_length; ++j)
                forward[cell(i, j)] = std::min(
                    {forward[cell(i - 1, j)] + 1,
                     forward[cell(i, j - 1)] + 1,
                     forward[cell(i - 1, j - 1)] +
                         (shortened[i - 1] != read[j - 1])});
        for (int i = 0; i <= shortened_length; ++i)
            backward[cell(i, read_length)] =
                shortened_length - i;
        for (int j = 0; j <= read_length; ++j)
            backward[cell(shortened_length, j)] =
                read_length - j;
        for (int i = shortened_length - 1; i >= 0; --i)
            for (int j = read_length - 1; j >= 0; --j)
                backward[cell(i, j)] = std::min(
                    {backward[cell(i + 1, j)] + 1,
                     backward[cell(i, j + 1)] + 1,
                     backward[cell(i + 1, j + 1)] +
                         (shortened[i] != read[j])});

        for (int position = insertion_begin;
             position <= insertion_end; ++position)
        {
            int deleted_insert_cost =
                std::numeric_limits<int>::max();
            for (int j = 0; j <= read_length; ++j)
                deleted_insert_cost = std::min(
                    deleted_insert_cost,
                    forward[cell(position, j)] + 1 +
                        backward[cell(position, j)]);
            for (int base_code = 0; base_code < 4; ++base_code)
            {
                static constexpr char alphabet[] = "ACGT";
                int best = deleted_insert_cost;
                for (int j = 0; j < read_length; ++j)
                    best = std::min(
                        best,
                        forward[cell(position, j)] +
                            (alphabet[base_code] != read[j]) +
                            backward[cell(position, j + 1)]);
                action_energy[
                    static_cast<size_t>(position) * 4 +
                    base_code] += best;
            }
        }
    }

    int best_energy = std::numeric_limits<int>::max();
    std::string best = sequence;
    for (int position = insertion_begin;
         position <= insertion_end; ++position)
    {
        static constexpr char alphabet[] = "ACGT";
        for (int base_code = 0; base_code < 4; ++base_code)
        {
            const int energy = action_energy[
                static_cast<size_t>(position) * 4 + base_code];
            if (energy < best_energy)
            {
                best_energy = energy;
                best = shortened;
                best.insert(best.begin() + position,
                            alphabet[base_code]);
            }
        }
    }
    // The original guide is the all-KEEP path of this boundary action graph.
    // A shifted guide violates one high-confidence conserved anchor, so KEEP
    // carries one edit-equivalent penalty per violated anchor base. The move
    // is committed only when read energy plus anchor evidence beats KEEP.
    const int read_gain = base_energy - best_energy;
    if (read_gain_out)
        *read_gain_out = read_gain;
    const int restored_anchor_support =
        delete_position == 0 ? left_length : right_length;
    if (best == sequence ||
        read_gain + restored_anchor_support <= 0)
        return sequence;
    return best;
}

// Rank clusters without using center answers. A sequence is suspicious only
// when its weakest phase windows are surprising both to the leave-one-out
// library field and to its own raw-read field.
static inline double loo_markov_anomaly_score(
    const std::string &sequence,
    const std::vector<std::string> &reads,
    const LooMarkovModel &model,
    const Params &p)
{
    const int k = model.k;
    const int local_k = p.loo_local_k;
    const int L = static_cast<int>(sequence.size());
    if (model.edge_counts.empty() || L <= k ||
        local_k < 2 || local_k > 8)
        return 0.0;

    std::unordered_map<uint32_t, int> own_edges;
    std::unordered_map<uint32_t, int> own_nodes;
    std::unordered_map<uint32_t, int> own_suffix_nodes;
    own_edges.reserve(sequence.size());
    own_nodes.reserve(sequence.size());
    own_suffix_nodes.reserve(sequence.size());
    const uint32_t node_mask =
        (static_cast<uint32_t>(1) << (2 * k)) - 1U;
    for (int pos = 0; pos + k < L; ++pos)
    {
        const uint32_t edge = kmer_code(sequence, pos, k + 1);
        ++own_edges[edge];
        ++own_nodes[edge >> 2];
        ++own_suffix_nodes[edge & node_mask];
    }
    auto own_count = [](const auto &counts, uint32_t code) -> int
    {
        const auto found = counts.find(code);
        return found == counts.end() ? 0 : found->second;
    };
    constexpr double alpha = 0.5;
    std::vector<double> anomaly(L, 0.0);
    for (int pos = 0; pos + k < L; ++pos)
    {
        const uint32_t edge = kmer_code(sequence, pos, k + 1);
        const uint32_t node = edge >> 2;
        const int edge_count = std::max(
            0, static_cast<int>(model.edge_counts[edge]) -
                   own_count(own_edges, edge));
        const int node_count = std::max(
            0, static_cast<int>(model.node_counts[node]) -
                   own_count(own_nodes, node));
        double log_probability = std::log(
            (edge_count + alpha) /
            (node_count + 4.0 * alpha));
        if (p.loo_reverse_weight > 0.0)
        {
            const uint32_t suffix_node = edge & node_mask;
            const int suffix_count = std::max(
                0, static_cast<int>(
                       model.suffix_node_counts[suffix_node]) -
                       own_count(own_suffix_nodes, suffix_node));
            log_probability += p.loo_reverse_weight * std::log(
                (edge_count + alpha) /
                (suffix_count + 4.0 * alpha));
        }
        const double surprise = -log_probability;
        for (int base_pos = pos; base_pos <= pos + k; ++base_pos)
            anomaly[base_pos] =
                std::max(anomaly[base_pos], surprise);
    }

    const size_t local_edge_universe =
        static_cast<size_t>(1) << (2 * (local_k + 1));
    const size_t local_node_universe =
        static_cast<size_t>(1) << (2 * local_k);
    std::vector<uint32_t> local_edges(local_edge_universe, 0);
    std::vector<uint32_t> local_nodes(local_node_universe, 0);
    std::vector<uint32_t> local_suffix_nodes(
        local_node_universe, 0);
    const uint32_t local_node_mask =
        static_cast<uint32_t>(local_node_universe - 1);
    for (const auto &read : reads)
    {
        for (int pos = 0;
             pos + local_k < static_cast<int>(read.size()); ++pos)
        {
            const uint32_t edge =
                kmer_code(read, pos, local_k + 1);
            ++local_edges[edge];
            ++local_nodes[edge >> 2];
            ++local_suffix_nodes[edge & local_node_mask];
        }
    }
    std::vector<double> local_anomaly(L, 0.0);
    for (int pos = 0; pos + local_k < L; ++pos)
    {
        const uint32_t edge =
            kmer_code(sequence, pos, local_k + 1);
        const uint32_t node = edge >> 2;
        double log_probability = std::log(
            (local_edges[edge] + alpha) /
            (local_nodes[node] + 4.0 * alpha));
        if (p.loo_reverse_weight > 0.0)
            log_probability += p.loo_reverse_weight * std::log(
                (local_edges[edge] + alpha) /
                (local_suffix_nodes[edge & local_node_mask] +
                 4.0 * alpha));
        const double surprise = -log_probability;
        for (int base_pos = pos;
             base_pos <= pos + local_k; ++base_pos)
            local_anomaly[base_pos] =
                std::max(local_anomaly[base_pos], surprise);
    }
    for (int pos = 0; pos < L; ++pos)
        anomaly[pos] +=
            p.loo_local_weight * local_anomaly[pos];

    constexpr int top_count = 24;
    const int keep = std::min<int>(top_count, anomaly.size());
    if (keep < static_cast<int>(anomaly.size()))
        std::nth_element(
            anomaly.begin(), anomaly.begin() + keep, anomaly.end(),
            std::greater<double>());
    else
        std::sort(anomaly.begin(), anomaly.end(),
                  std::greater<double>());
    return std::accumulate(
               anomaly.begin(), anomaly.begin() + keep, 0.0) /
           std::max(1, keep);
}

} // namespace mcs_phase

#endif
