#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "Cluster2.hpp"

using namespace std;

struct Config {
    string input_file;
    string output_file;
    string separator = "====";

    int target_len = 0;
    int threads = 32;
    int max_reads = 25;
    int limit = -1;

    uint32_t seed = 20260802u;
};

static string trim(const string& input) {
    size_t begin = 0;
    size_t end = input.size();

    while (begin < end &&
           isspace(static_cast<unsigned char>(input[begin]))) {
        ++begin;
    }

    while (end > begin &&
           isspace(static_cast<unsigned char>(input[end - 1]))) {
        --end;
    }

    return input.substr(begin, end - begin);
}

static vector<vector<string>> load_clusters(
    const string& path,
    const string& separator
) {
    ifstream input(path);

    if (!input.is_open()) {
        throw runtime_error("cannot open input file: " + path);
    }

    vector<vector<string>> clusters;
    vector<string> current_cluster;
    string line;

    while (getline(input, line)) {
        line = trim(line);

        if (line.empty()) {
            continue;
        }

        if (line.find(separator) != string::npos) {
            if (!current_cluster.empty()) {
                clusters.push_back(current_cluster);
                current_cluster.clear();
            }
        } else {
            current_cluster.push_back(line);
        }
    }

    if (!current_cluster.empty()) {
        clusters.push_back(current_cluster);
    }

    return clusters;
}

static void print_help(const char* program) {
    cerr
        << "Usage:\n"
        << program
        << " -i Clusters.txt"
        << " -l target_len"
        << " -o result.txt"
        << " [-s separator]"
        << " [--threads 32]"
        << " [--max_reads 25]"
        << " [--seed 20260802]"
        << " [--limit N]\n";
}

static Config parse_arguments(int argc, char* argv[]) {
    Config config;

    for (int i = 1; i < argc; ++i) {
        string argument = argv[i];

        auto require_value = [&](const string& option) -> string {
            if (i + 1 >= argc) {
                throw runtime_error("missing value for " + option);
            }

            return argv[++i];
        };

        if (argument == "-i" || argument == "--input") {
            config.input_file = require_value(argument);
        } else if (argument == "-o" || argument == "--output") {
            config.output_file = require_value(argument);
        } else if (argument == "-l" || argument == "--target_len") {
            config.target_len = stoi(require_value(argument));
        } else if (argument == "-s" || argument == "--separator") {
            config.separator = require_value(argument);
        } else if (argument == "--threads" || argument == "--jobs") {
            config.threads = stoi(require_value(argument));
        } else if (argument == "--max_reads") {
            config.max_reads = stoi(require_value(argument));
        } else if (argument == "--seed") {
            config.seed =
                static_cast<uint32_t>(stoul(require_value(argument)));
        } else if (argument == "--limit") {
            config.limit = stoi(require_value(argument));
        } else if (argument == "-h" || argument == "--help") {
            print_help(argv[0]);
            exit(0);
        } else {
            throw runtime_error("unknown argument: " + argument);
        }
    }

    if (config.input_file.empty()) {
        throw runtime_error("input file is required");
    }

    if (config.output_file.empty()) {
        throw runtime_error("output file is required");
    }

    if (config.target_len <= 0) {
        throw runtime_error("target length must be greater than zero");
    }

    if (config.threads <= 0) {
        throw runtime_error("threads must be greater than zero");
    }

    if (config.max_reads <= 0) {
        throw runtime_error("max_reads must be greater than zero");
    }

    return config;
}

int main(int argc, char* argv[]) {
    try {
        Config config = parse_arguments(argc, argv);

        vector<vector<string>> clusters =
            load_clusters(config.input_file, config.separator);

        if (clusters.empty()) {
            throw runtime_error("no clusters were loaded");
        }

        size_t run_count = clusters.size();

        if (config.limit > 0) {
            run_count = min(
                run_count,
                static_cast<size_t>(config.limit)
            );
        }

        vector<string> predictions(run_count);
        vector<string> error_messages(run_count);

        int failed_clusters = 0;

        /*
         * 作者的主程序参数：
         *
         * delPatternLen = 3
         * subPriority   = 0
         * delPriority   = 0
         * insPriority   = 0
         * maxReps       = 2
         * maxCopies     = 25
         */
        const int del_pattern_len = 3;
        const int sub_priority = 0;
        const int del_priority = 0;
        const int ins_priority = 0;
        const int max_reps = 2;

        omp_set_dynamic(0);
        omp_set_num_threads(config.threads);

        auto start_time = chrono::steady_clock::now();

        #pragma omp parallel for schedule(dynamic, 1) \
            reduction(+:failed_clusters)
        for (long long cluster_index = 0;
             cluster_index < static_cast<long long>(run_count);
             ++cluster_index) {
            try {
                const vector<string>& all_reads =
                    clusters[cluster_index];

                size_t reads_to_use = min(
                    all_reads.size(),
                    static_cast<size_t>(config.max_reads)
                );

                if (reads_to_use == 0) {
                    throw runtime_error("cluster has no reads");
                }

                vector<string> reads(
                    all_reads.begin(),
                    all_reads.begin() + reads_to_use
                );

                /*
                 * Cluster2 类要求传入 original。
                 * 这里传入长度为 target_len 的虚拟序列。
                 *
                 * ITR 的 TestBest 主重建路径只使用 original.size()
                 * 约束候选长度，不使用真实中心碱基。
                 */
                string dummy_original(config.target_len, 'A');

                Cluster2 cluster(dummy_original, reads);

                uint32_t local_seed =
                    config.seed +
                    static_cast<uint32_t>(cluster_index) * 1000003u;

                mt19937 generator(local_seed);

                int final_guess_edit_distance = 0;
                string final_guess;

                if (reads.size() == 1) {
                    final_guess = reads[0];
                } else {
                    final_guess = cluster.TestBest(
                        del_pattern_len,
                        final_guess_edit_distance,
                        sub_priority,
                        del_priority,
                        ins_priority,
                        generator,
                        max_reps
                    );
                }

                predictions[cluster_index] = final_guess;
            } catch (const exception& error) {
                error_messages[cluster_index] = error.what();
                predictions[cluster_index].clear();
                failed_clusters += 1;
            }
        }

        auto end_time = chrono::steady_clock::now();

        double elapsed_seconds =
            chrono::duration<double>(
                end_time - start_time
            ).count();

        ofstream output(config.output_file);

        if (!output.is_open()) {
            throw runtime_error(
                "cannot create output file: " +
                config.output_file
            );
        }

        for (const string& prediction : predictions) {
            output << prediction << '\n';
        }

        output.close();

        cout << "algorithm: ITR-OpenMP" << '\n';
        cout << "clusters loaded: " << clusters.size() << '\n';
        cout << "clusters processed: " << run_count << '\n';
        cout << "target length: " << config.target_len << '\n';
        cout << "threads: " << config.threads << '\n';
        cout << "max reads per cluster: "
             << config.max_reads << '\n';
        cout << "seed: " << config.seed << '\n';
        cout << "failed clusters: " << failed_clusters << '\n';
        cout << "written output: " << config.output_file << '\n';
        cout << "elapsed seconds: " << elapsed_seconds << '\n';

        if (failed_clusters > 0) {
            int shown = 0;

            for (size_t i = 0;
                 i < error_messages.size() && shown < 10;
                 ++i) {
                if (!error_messages[i].empty()) {
                    cerr << "cluster " << (i + 1)
                         << " failed: "
                         << error_messages[i] << '\n';
                    ++shown;
                }
            }

            return 2;
        }

        return 0;
    } catch (const exception& error) {
        cerr << "error: " << error.what() << '\n';
        print_help(argv[0]);
        return 1;
    }
}
