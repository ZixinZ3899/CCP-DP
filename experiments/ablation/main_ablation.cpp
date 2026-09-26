#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ccpdp_ablation.hpp"

struct Args
{
    std::string input_path;
    std::string output_path;
    std::string diag_path;
    std::string guide_output_path;
    std::string separator = "====";
    std::string mode = "auto";
    int target_len = 110;
    int limit = -1;
    int jobs = 1;
    bool enable_unified_graph = true;
    bool decode_all = false;
    bool disable_shallow_dense_fallback = false;
    bool disable_cross_loo_score = false;
    bool disable_blind_ids_likelihood = false;
    double route_fraction = -1.0;
    double graph_cross_weight = -1.0;
    double graph_local_weight = -1.0;
    int graph_phase_band = -1;
    int graph_context_states = -1;
    mcs_phase::Params params;
};

static inline std::string trim(const std::string &s)
{
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a])))
        ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1])))
        --b;
    return s.substr(a, b - a);
}

static void apply_auto_mode(Args &args)
{
    args.params = mcs_phase::Params{};
    args.mode = "auto";
    // One dataset-independent pipeline: one local guide solve followed by at
    // most one evidence-gated structured decode.
    args.params.band_extra = 13;
    args.params.max_reads = 36;
    args.params.uniform_read_weight = true;
    args.params.no_insertion_prior = 0.0;
    args.params.no_deletion_prior = 0.0;
    args.params.gap_event_penalty = 0.0;
    args.params.topology_k = 12;
    args.params.topology_min_support = 2;
    args.params.loo_local_weight = 0.4;
    args.params.blind_channel = true;
}

static void usage(const char *prog)
{
    std::cerr
        << "CCP-DP v8.2 (adaptive calibration and data-driven routing)\n\n"
        << "Usage:\n"
        << "  " << prog << " -i Clusters.txt -l 110 -s \"====\" -o result.txt"
        << " --diag diag.csv --mode auto --jobs 32\n\n"
        << "Options:\n"
        << "  --mode auto            blind channel calibration and evidence gate\n"
        << "  --jobs INT             OpenMP threads\n"
        << "  --limit INT            process first INT clusters (debug only)\n"
        << "  --disable-unified-graph         output coordinate guides\n"
        << "  --decode-all                    ablation: decode every cluster\n"
        << "  --disable-cross-loo-score       ablation: omit Cross-LOO edge score only\n"
        << "  --disable-blind-ids-likelihood  ablation: omit IDS likelihood from DP only\n"
        << "  --route-fraction FLOAT          diagnostic: decode top fraction by anomaly\n"
        << "  --disable-shallow-dense-fallback  debug: retain sparse 8-state route\n"
        << "  --graph-cross-weight FLOAT      debug/ablation override\n"
        << "  --graph-local-weight FLOAT      debug/ablation override\n"
        << "  --graph-phase-band INT          debug/ablation override\n"
        << "  --graph-context-states INT      paths retained per trellis cell\n"
        << "  --guide-output PATH             write frozen coordinate guides\n";
}

static Args parse_args(int argc, char **argv)
{
    Args args;
    apply_auto_mode(args);
    auto value = [&](int &i) -> std::string
    {
        if (i + 1 >= argc)
            throw std::runtime_error(std::string("Missing value for ") + argv[i]);
        return argv[++i];
    };
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (a == "-h" || a == "--help")
        {
            usage(argv[0]);
            std::exit(0);
        }
        else if (a == "-i" || a == "--input")
            args.input_path = value(i);
        else if (a == "-o" || a == "--output")
            args.output_path = value(i);
        else if (a == "-l" || a == "--length")
            args.target_len = std::stoi(value(i));
        else if (a == "-s" || a == "--separator")
            args.separator = value(i);
        else if (a == "--diag")
            args.diag_path = value(i);
        else if (a == "--guide-output")
            args.guide_output_path = value(i);
        else if (a == "--mode")
        {
            const std::string mode = value(i);
            if (mode != "auto")
                throw std::runtime_error("Only --mode auto is supported");
            apply_auto_mode(args);
        }
        else if (a == "--jobs")
            args.jobs = std::max(1, std::stoi(value(i)));
        else if (a == "--limit")
            args.limit = std::stoi(value(i));
        else if (a == "--disable-unified-graph")
            args.enable_unified_graph = false;
        else if (a == "--decode-all")
            args.decode_all = true;
        else if (a == "--disable-cross-loo-score")
            args.disable_cross_loo_score = true;
        else if (a == "--disable-blind-ids-likelihood")
            args.disable_blind_ids_likelihood = true;
        else if (a == "--route-fraction")
        {
            args.route_fraction = std::stod(value(i));
            if (args.route_fraction < 0.0 || args.route_fraction > 1.0)
                throw std::runtime_error(
                    "--route-fraction must be between 0 and 1");
        }
        else if (a == "--disable-shallow-dense-fallback")
            args.disable_shallow_dense_fallback = true;
        else if (a == "--graph-cross-weight")
            args.graph_cross_weight = std::stod(value(i));
        else if (a == "--graph-local-weight")
            args.graph_local_weight = std::stod(value(i));
        else if (a == "--graph-phase-band")
            args.graph_phase_band = std::max(1, std::stoi(value(i)));
        else if (a == "--graph-context-states")
            args.graph_context_states = std::max(1, std::stoi(value(i)));
        else
            throw std::runtime_error("Unknown argument: " + a);
    }
    if (args.input_path.empty() || args.output_path.empty())
        throw std::runtime_error("Both -i and -o are required");
    if (args.target_len <= 0)
        throw std::runtime_error("Invalid target length");
    args.params.target_len = args.target_len;
    args.params.use_blind_ids_likelihood =
        !args.disable_blind_ids_likelihood;
    if (args.diag_path.empty())
        args.diag_path = args.output_path + ".diag.csv";
    return args;
}

static std::vector<std::vector<std::string>> read_clusters(const std::string &path,
                                                           const std::string &separator,
                                                           int limit)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open input: " + path);
    std::vector<std::vector<std::string>> clusters;
    std::vector<std::string> cur;
    auto push = [&]()
    {
        if (!cur.empty())
        {
            clusters.push_back(std::move(cur));
            cur.clear();
        }
    };
    std::string line;
    while (std::getline(in, line))
    {
        std::string t = trim(line);
        if (t.empty())
            continue;
        const bool equals_separator =
            t.size() >= 3 &&
            std::all_of(t.begin(), t.end(),
                        [](char c) { return c == '='; });
        const bool named_cluster_separator =
            t.rfind("CLUSTER ", 0) == 0;
        if ((!separator.empty() &&
             t.find(separator) != std::string::npos) ||
            equals_separator || named_cluster_separator)
        {
            push();
            if (limit > 0 && static_cast<int>(clusters.size()) >= limit)
                break;
            continue;
        }
        std::stringstream ss(t);
        std::string token;
        while (ss >> token)
        {
            if (!separator.empty() && token.find(separator) != std::string::npos)
                push();
            else
                cur.push_back(token);
        }
        if (limit > 0 && static_cast<int>(clusters.size()) >= limit)
            break;
    }
    if (limit <= 0 || static_cast<int>(clusters.size()) < limit)
        push();
    if (limit > 0 && static_cast<int>(clusters.size()) > limit)
        clusters.resize(limit);
    return clusters;
}

struct AnomalyMixture
{
    double low_mean = 0.0;
    double high_mean = 0.0;
    double separation = 0.0;
    size_t high_count = 0;
    bool separated = false;
};

// A deterministic one-dimensional two-means split.  It is used only to
// estimate how much of a high-reuse library belongs to the high-risk anomaly
// component; center answers are never consulted.
static AnomalyMixture estimate_anomaly_mixture(
    const std::vector<double> &scores)
{
    AnomalyMixture result;
    if (scores.size() < 2)
        return result;
    double low = *std::min_element(scores.begin(), scores.end());
    double high = *std::max_element(scores.begin(), scores.end());
    if (!(high > low))
        return result;
    for (int iteration = 0; iteration < 64; ++iteration)
    {
        const double boundary = 0.5 * (low + high);
        double low_sum = 0.0;
        double high_sum = 0.0;
        size_t low_count = 0;
        size_t high_count = 0;
        for (double score : scores)
            if (score <= boundary)
            {
                low_sum += score;
                ++low_count;
            }
            else
            {
                high_sum += score;
                ++high_count;
            }
        if (low_count == 0 || high_count == 0)
            return result;
        const double next_low = low_sum / low_count;
        const double next_high = high_sum / high_count;
        if (std::max(std::abs(next_low - low),
                     std::abs(next_high - high)) < 1e-10)
        {
            low = next_low;
            high = next_high;
            break;
        }
        low = next_low;
        high = next_high;
    }
    const double boundary = 0.5 * (low + high);
    double low_sq = 0.0;
    double high_sq = 0.0;
    size_t low_count = 0;
    size_t high_count = 0;
    for (double score : scores)
        if (score <= boundary)
        {
            low_sq += (score - low) * (score - low);
            ++low_count;
        }
        else
        {
            high_sq += (score - high) * (score - high);
            ++high_count;
        }
    if (low_count == 0 || high_count == 0)
        return result;
    const double low_sd = std::sqrt(low_sq / low_count);
    const double high_sd = std::sqrt(high_sq / high_count);
    result.low_mean = low;
    result.high_mean = high;
    result.high_count = high_count;
    result.separation =
        (high - low) / std::max(1e-12, low_sd + high_sd);
    // A single approximately Gaussian population is split by two-means with
    // separation near 1.3.  Requiring 2.0 avoids expanding the route on a
    // merely broad but unimodal easy-library distribution.
    result.separated = result.separation >= 2.0;
    return result;
}

static void activate_top_anomalies(const std::vector<double> &scores,
                                   size_t active_count,
                                   std::vector<uint8_t> &active)
{
    active_count = std::min(active_count, scores.size());
    if (active_count == 0)
        return;
    if (active_count == scores.size())
    {
        std::fill(active.begin(), active.end(), 1);
        return;
    }
    std::vector<double> distribution = scores;
    const size_t gate_index = distribution.size() - active_count;
    std::nth_element(distribution.begin(),
                     distribution.begin() + gate_index,
                     distribution.end());
    const double gate = distribution[gate_index];
    size_t selected = 0;
    for (size_t i = 0; i < scores.size(); ++i)
        if (scores[i] > gate)
        {
            active[i] = 1;
            ++selected;
        }
    // Deterministically resolve the usually empty threshold tie while keeping
    // the requested compute budget exact.
    for (size_t i = 0; i < scores.size() && selected < active_count; ++i)
        if (scores[i] == gate)
        {
            active[i] = 1;
            ++selected;
        }
}

int main(int argc, char **argv)
{
    try
    {
        Args args = parse_args(argc, argv);
#ifdef _OPENMP
        omp_set_num_threads(args.jobs);
#endif
        auto t0 = std::chrono::steady_clock::now();
        auto clusters = read_clusters(args.input_path, args.separator, args.limit);
        const int n = static_cast<int>(clusters.size());
        if (n == 0)
            throw std::runtime_error(
                "No non-empty clusters were found in the input");
        // Stage 1: exactly one local fixed-length DP builds a coordinate
        // estimate.  It is not iteratively refined.
        args.params.positional_phase_chain = false;
        std::vector<std::string> results(n);
        std::vector<std::string> coordinate_guides(n);
        std::vector<mcs_phase::Diag> diags(n);
        double cross_cluster_reuse = -1.0;
        int graph_decoded_clusters = 0;
        mcs_phase::ChannelProfile global_channel;
        int unified_graph_changes = 0;
        double preliminary_phase_load = 0.0;
        double preliminary_total_error = 0.0;
        mcs_phase::CrossClusterCalibration cross_calibration;
        std::string decode_mode =
            "one local guide solve + at most one gated structured decode";

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int i = 0; i < n; ++i)
        {
            mcs_phase::Result local =
                mcs_phase::build_local_coordinate_once(
                    clusters[i], i + 1, args.params);
            coordinate_guides[i] = std::move(local.sequence);
            diags[i] = std::move(local.diag);
        }
        if (!args.guide_output_path.empty())
        {
            std::ofstream guide_out(args.guide_output_path);
            if (!guide_out)
                throw std::runtime_error(
                    "Cannot write guide output: " + args.guide_output_path);
            for (const auto &guide : coordinate_guides)
                guide_out << guide << '\n';
        }

        // Build the blind channel evidence once against the unchanged
        // observed guides.  It is a probability field for the final DP, not
        // a post-hoc sequence correction.
        mcs_phase::ChannelCounts preliminary_counts;
        const int preliminary_sample_count = std::min(
            n, std::max(128, static_cast<int>(
                                  std::llround(4.0 * std::sqrt(
                                      static_cast<double>(n))))));
        std::vector<mcs_phase::ChannelCounts>
            preliminary_cluster_counts(preliminary_sample_count);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int sample = 0; sample < preliminary_sample_count; ++sample)
        {
            const int i = static_cast<int>(
                static_cast<long long>(sample) * n /
                preliminary_sample_count);
            preliminary_cluster_counts[sample] =
                mcs_phase::estimate_cluster_channel(
                    clusters[i], coordinate_guides[i], args.params);
        }
        for (const auto &counts : preliminary_cluster_counts)
            preliminary_counts = mcs_phase::add_channel_counts(
                preliminary_counts, counts);
        const mcs_phase::ChannelProfile preliminary_channel =
            mcs_phase::channel_profile_from_alignment_counts(
                preliminary_counts);
        global_channel = preliminary_channel;
        preliminary_phase_load =
            args.target_len *
            (preliminary_channel.insertion + preliminary_channel.deletion);
        preliminary_total_error =
            preliminary_channel.substitution +
            preliminary_channel.insertion +
            preliminary_channel.deletion;

        // Build the cross-cluster Context prior once from observed guides.
        // A second, center-free leave-one-out calibration below determines
        // whether this prior predicts held-out guide transitions.  A prior
        // that fails validation receives exactly zero weight.
        mcs_phase::Params reuse_params = args.params;
        reuse_params.topology_k = 11;
        const mcs_phase::CrossClusterReuseModel reuse_model =
            mcs_phase::build_cross_cluster_reuse_model(
                coordinate_guides, reuse_params);
        cross_cluster_reuse =
            mcs_phase::cross_cluster_reuse_fraction(
                coordinate_guides, reuse_model);
        args.params.positional_min_reads = 6;
        args.params.positional_k = 7;
        args.params.positional_radius = 2;
        args.params.loo_markov_k =
            cross_cluster_reuse >= 0.5 ? 11 : 8;
        const mcs_phase::LooMarkovModel cross_model =
            mcs_phase::build_loo_markov_model(
                coordinate_guides, args.params.loo_markov_k);
        if (args.graph_cross_weight < 0.0)
            cross_calibration =
                mcs_phase::calibrate_cross_cluster_trust(
                    coordinate_guides, cross_model, args.params);

        // If the cross-cluster model is rejected on a moderate/high-indel
        // random library, a sparse 8-state trellis wastes work on context
        // hypotheses that have zero weight and can miss correctable strands.
        // Use a dense but shallow 4-state trellis instead.  This condition is
        // derived entirely from frozen guides and blind channel statistics.
        const bool use_shallow_dense_fallback =
            args.enable_unified_graph &&
            !args.disable_shallow_dense_fallback &&
            args.graph_cross_weight < 0.0 &&
            !cross_calibration.accepted &&
            cross_cluster_reuse < 0.10 &&
            preliminary_phase_load >= 3.0 &&
            preliminary_phase_load < 6.0;
        std::vector<double> joint_anomaly(n, 0.0);
        std::vector<uint8_t> joint_active(n, 0);
        std::string routing_policy = "legacy";
        size_t planned_active_count = 0;
        AnomalyMixture anomaly_mixture;
        if (args.decode_all || use_shallow_dense_fallback)
        {
            // Anomaly scores are routing-only quantities.  Do not calculate
            // them when the route has already been determined as dense.
            std::fill(joint_active.begin(), joint_active.end(), 1);
            planned_active_count = n;
            routing_policy = args.decode_all
                                 ? "explicit-decode-all"
                                 : "shallow-dense-fallback";
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
            for (int i = 0; i < n; ++i)
                joint_anomaly[i] = mcs_phase::loo_markov_anomaly_score(
                    coordinate_guides[i], clusters[i], cross_model,
                    args.params);
            if (args.route_fraction >= 0.0)
            {
                const size_t active_count = std::min<size_t>(
                    n, static_cast<size_t>(std::ceil(
                           args.route_fraction * static_cast<double>(n))));
                activate_top_anomalies(
                    joint_anomaly, active_count, joint_active);
                planned_active_count = active_count;
                routing_policy = "manual-anomaly-fraction";
            }
            else if (cross_cluster_reuse >= 0.5)
            {
                anomaly_mixture = estimate_anomaly_mixture(joint_anomaly);
                const size_t minimum_count = static_cast<size_t>(
                    std::ceil(0.10 * static_cast<double>(n)));
                const size_t maximum_count = static_cast<size_t>(
                    std::ceil(0.50 * static_cast<double>(n)));
                size_t active_count = minimum_count;
                if (anomaly_mixture.separated)
                    active_count = std::clamp(
                        anomaly_mixture.high_count,
                        minimum_count, maximum_count);
                activate_top_anomalies(
                    joint_anomaly, active_count, joint_active);
                planned_active_count = active_count;
                routing_policy = anomaly_mixture.separated
                                     ? "adaptive-bimodal-anomaly"
                                     : "sparse-10pct-fallback";
            }
            else
            {
                std::vector<double> anomaly_distribution = joint_anomaly;
                const size_t anomaly_gate_index =
                    anomaly_distribution.size() / 2;
                std::nth_element(
                    anomaly_distribution.begin(),
                    anomaly_distribution.begin() + anomaly_gate_index,
                    anomaly_distribution.end());
                const double joint_anomaly_gate =
                    anomaly_distribution[anomaly_gate_index];
                for (int i = 0; i < n; ++i)
                    if (joint_anomaly[i] >= joint_anomaly_gate)
                        joint_active[i] = 1;
                std::vector<int> coverages(n);
                for (int i = 0; i < n; ++i)
                    coverages[i] = static_cast<int>(clusters[i].size());
                const size_t coverage_gate_index =
                    coverages.size() * 3 / 4;
                std::nth_element(
                    coverages.begin(),
                    coverages.begin() + coverage_gate_index,
                    coverages.end());
                const int coverage_gate = coverages[coverage_gate_index];
                for (int i = 0; i < n; ++i)
                    if (static_cast<int>(clusters[i].size()) >=
                        coverage_gate)
                        joint_active[i] = 1;
                if (preliminary_phase_load >= 6.0)
                {
                    std::fill(joint_active.begin(),
                              joint_active.end(), 1);
                    routing_policy = "low-reuse-high-phase-dense";
                }
                if (preliminary_phase_load < 3.0)
                {
                    const int minimum_weak_positions =
                        cross_cluster_reuse >= 0.25 ? 2 : 1;
                    for (int i = 0; i < n; ++i)
                        joint_active[i] =
                            diags[i].weak_positions >=
                                    minimum_weak_positions
                                ? 1
                                : 0;
                    routing_policy = "low-phase-weak-position";
                }
                if (routing_policy == "legacy")
                    routing_policy = "low-reuse-moderate-phase";
                planned_active_count = static_cast<size_t>(std::count(
                    joint_active.begin(), joint_active.end(), uint8_t{1}));
            }
        }

        // Freeze the boundary evidence before any strand is decoded.  This
        // removes the v7.8 dependency on primary decoded outputs and makes a
        // single generalized-graph decode possible.
        const bool use_unified_boundary_edges =
            args.enable_unified_graph &&
            cross_cluster_reuse < 0.5 && preliminary_phase_load < 3.0;
        mcs_phase::ConservedAnchorModel unified_boundary_model;
        if (use_unified_boundary_edges)
        {
            std::vector<int> coverages(n);
            for (int i = 0; i < n; ++i)
                coverages[i] = static_cast<int>(clusters[i].size());
            const size_t q3_index = coverages.size() * 3 / 4;
            std::nth_element(coverages.begin(),
                             coverages.begin() + q3_index,
                             coverages.end());
            const int coverage_q3 = coverages[q3_index];
            std::vector<std::string> high_confidence_guides;
            high_confidence_guides.reserve(n / 3);
            for (int i = 0; i < n; ++i)
                if (static_cast<int>(clusters[i].size()) >= coverage_q3)
                    high_confidence_guides.push_back(
                        coordinate_guides[i]);
            unified_boundary_model =
                mcs_phase::build_conserved_anchor_model(
                    high_confidence_guides);
        }

        // Dataset-blind trellis configuration selected once from measured
        // reuse and phase load. These are the actual values used below and
        // reported in the run log.
        const bool low_phase_load = preliminary_phase_load < 3.0;
        const double effective_graph_local_weight =
            args.graph_local_weight >= 0.0
                ? args.graph_local_weight
                : (cross_cluster_reuse >= 0.5
                       ? 0.0
                       : (low_phase_load ? 0.80 : 0.40));
        const double base_graph_cross_weight =
            cross_cluster_reuse < 0.5
                ? (low_phase_load ? 0.0 : 0.50)
                : 0.68;
        const double calibrated_graph_cross_weight =
            args.graph_cross_weight >= 0.0
                ? args.graph_cross_weight
                : base_graph_cross_weight * cross_calibration.trust;
        const double effective_graph_cross_weight =
            args.disable_cross_loo_score
                ? 0.0
                : calibrated_graph_cross_weight;
        const int effective_graph_phase_band =
            args.graph_phase_band >= 1
                ? args.graph_phase_band
                : (cross_cluster_reuse >= 0.5
                       ? 6
                       : (low_phase_load ? 3 : 6));
        const int effective_graph_context_states =
            args.graph_context_states >= 1
                ? args.graph_context_states
                : (use_shallow_dense_fallback ? 4 : 8);

        // The only sequence-producing operation after coordinate setup: one
        // graph path consumes local alignment, cross-cluster Context,
        // blind-channel evidence, and optional boundary-shift actions.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : graph_decoded_clusters, unified_graph_changes)
#endif
        for (int i = 0; i < n; ++i)
        {
            if (coordinate_guides[i].empty())
            {
                results[i].assign(args.target_len, 'A');
                continue;
            }
            if (!args.enable_unified_graph)
            {
                results[i] = coordinate_guides[i];
                continue;
            }
            if (use_unified_boundary_edges)
            {
                // This signature test is only graph routing; it does not
                // decode a sequence. A detected strand enters the boundary
                // boundary-action graph, then exits after that single solve.
                const bool boundary_active =
                    mcs_phase::shifted_boundary_delete_position(
                        coordinate_guides[i],
                        unified_boundary_model) >= 0;
                if (boundary_active)
                {
                    results[i] =
                        mcs_phase::decode_boundary_action_graph(
                        coordinate_guides[i], clusters[i],
                        unified_boundary_model, args.params);
                    ++graph_decoded_clusters;
                    if (results[i] != coordinate_guides[i])
                        ++unified_graph_changes;
                    diags[i].passes = 1;
                    diags[i].final_len =
                        static_cast<int>(results[i].size());
                    diags[i].seed_to_final_edits =
                        mcs_phase::edit_distance_linear(
                            coordinate_guides[i], results[i]);
                    continue;
                }
            }
            if (!joint_active[i])
            {
                results[i] = coordinate_guides[i];
                continue;
            }
            ++graph_decoded_clusters;
            mcs_phase::Params joint_params = args.params;
            if (cross_cluster_reuse < 0.5 && !low_phase_load)
                joint_params.max_reads = std::max(
                    joint_params.max_reads, 60);
            joint_params.loo_local_weight =
                effective_graph_local_weight;
            const mcs_phase::JointPhaseTrellisEvidence evidence =
                mcs_phase::build_joint_phase_trellis_evidence(
                    coordinate_guides[i], clusters[i], cross_model,
                    preliminary_channel, joint_params,
                    effective_graph_cross_weight);
            results[i] = mcs_phase::decode_joint_phase_trellis(
                evidence, joint_params,
                effective_graph_phase_band,
                effective_graph_context_states);
            if (results[i] != coordinate_guides[i])
                ++unified_graph_changes;
            diags[i].passes = 1;
            diags[i].final_len = static_cast<int>(results[i].size());
            diags[i].seed_to_final_edits =
                mcs_phase::edit_distance_linear(
                    coordinate_guides[i], results[i]);
        }

        std::ofstream out(args.output_path);
        if (!out)
            throw std::runtime_error("Cannot write output: " + args.output_path);
        for (const auto &seq : results)
            out << seq << '\n';

        std::ofstream diag(args.diag_path);
        if (!diag)
            throw std::runtime_error("Cannot write diag: " + args.diag_path);
        diag << "cluster_id,num_reads,used_reads,seed_read_index,seed_len,passes,"
             << "aligned_reads,rejected_reads,final_len,seed_to_final_edits,weak_positions,"
             << "coherence_filtered,coherence_threshold,"
             << "mean_norm_dist,mean_chosen_conf,joint_anomaly,joint_active\n";
        for (int i = 0; i < n; ++i)
        {
            const auto &d = diags[i];
            diag << d.cluster_id << ',' << d.num_reads << ',' << d.used_reads << ','
                 << d.seed_read_index << ',' << d.seed_len << ',' << d.passes << ','
                 << d.aligned_reads << ',' << d.rejected_reads << ',' << d.final_len << ','
                 << d.seed_to_final_edits << ',' << d.weak_positions << ','
                 << d.coherence_filtered << ',' << d.coherence_threshold << ','
                 << d.mean_norm_dist << ','
                 << d.mean_chosen_conf << ','
                 << joint_anomaly[i] << ','
                 << static_cast<int>(joint_active[i]) << '\n';
        }
        const double sec = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - t0)
                               .count();
        std::cerr << "build: CCP-DP-v8.2-ablation-20260916\n"
                  << "clusters: " << n << '\n'
                  << "mode: " << args.mode << '\n'
                  << "decode mode: " << decode_mode << '\n'
                  << "unified graph: "
                  << (args.enable_unified_graph ? "enabled" : "disabled")
                  << '\n'
                  << "preliminary phase load: "
                  << preliminary_phase_load << '\n'
                  << "preliminary total error: "
                  << preliminary_total_error << '\n'
                  << (graph_decoded_clusters > 0
                          ? "graph-decoded clusters: " +
                                std::to_string(graph_decoded_clusters) +
                                "/" + std::to_string(n) + "\n"
                          : "")
                  << "unified-graph changes: "
                  << unified_graph_changes << "/" << n << '\n'
                  << (use_unified_boundary_edges
                          ? "boundary anchor lengths (left,right): " +
                                std::to_string(unified_boundary_model.left.size()) + "," +
                                std::to_string(unified_boundary_model.right.size()) + "\n"
                          : "")
                  << "cross-cluster reuse (distinct guides): "
                  << cross_cluster_reuse << '\n'
                  << "cross calibration (trust,accepted): "
                  << cross_calibration.trust << ','
                  << (cross_calibration.accepted ? "yes" : "no") << '\n'
                  << "cross calibration gain (train,validation,se): "
                  << cross_calibration.train_gain << ','
                  << cross_calibration.validation_gain << ','
                  << cross_calibration.validation_se << '\n'
                  << "cross calibration sample (guides,windows): "
                  << cross_calibration.sampled_guides << ','
                  << cross_calibration.sampled_windows << '\n'
                  << "trellis weights (cross,local): "
                  << effective_graph_cross_weight << ','
                  << effective_graph_local_weight << '\n'
                  << "ablation cross-LOO score: "
                  << (args.disable_cross_loo_score ? "disabled" : "enabled")
                  << '\n'
                  << "ablation blind-IDS likelihood: "
                  << (args.disable_blind_ids_likelihood ? "disabled" : "enabled")
                  << '\n'
                  << "trellis phase band: "
                  << effective_graph_phase_band << '\n'
                  << "trellis context states: "
                  << effective_graph_context_states << '\n'
                  << "shallow dense fallback: "
                  << (use_shallow_dense_fallback ? "yes" : "no") << '\n'
                  << "routing policy: " << routing_policy << '\n'
                  << "planned structured decoding: "
                  << planned_active_count << "/" << n << '\n'
                  << "anomaly mixture (low,high,separation,high-count): "
                  << anomaly_mixture.low_mean << ','
                  << anomaly_mixture.high_mean << ','
                  << anomaly_mixture.separation << ','
                  << anomaly_mixture.high_count << '\n'
                  << "route fraction override: "
                  << args.route_fraction << '\n'
                  << (args.params.blind_channel
                          ? "calibration IDS (sub,ins,del): " +
                                std::to_string(global_channel.substitution) + "," +
                                std::to_string(global_channel.insertion) + "," +
                                std::to_string(global_channel.deletion) + "\n"
                          : "")
                  << "written output: " << args.output_path << '\n'
                  << "written diag: " << args.diag_path << '\n'
                  << "elapsed seconds: " << sec << std::endl;
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << '\n';
        usage(argv[0]);
        return 1;
    }
}
