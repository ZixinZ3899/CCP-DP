#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

#include <omp.h>

// Reuse the authors' CPL implementation, while replacing only its serial I/O
// driver. Renaming its main avoids modifying the original DNAV8.cpp.
#define main cpl_original_main_disabled
#include "DNAV8.cpp"
#undef main

namespace
{

    struct Options
    {
        std::string input;
        std::string output;
        std::string separator = "====";
        int target_len = 0;
        int threads = 32;
        int max_copies = 32; // Same limit used by the authors' main program.
        int batch_size = 256;
        std::uint32_t seed = 20260803u;
    };

    void print_usage(const char *program)
    {
        std::cerr
            << "Usage: " << program
            << " -i CLUSTERS -l TARGET_LEN -o OUTPUT"
            << " [-s SEPARATOR] [-t THREADS]"
            << " [--max-copies N] [--batch-size N] [--seed N]\n";
    }

    Options parse_args(int argc, char **argv)
    {
        Options opt;
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            auto value = [&](const std::string &name) -> std::string
            {
                if (i + 1 >= argc)
                {
                    throw std::runtime_error("Missing value after " + name);
                }
                return argv[++i];
            };

            if (arg == "-i" || arg == "--input")
            {
                opt.input = value(arg);
            }
            else if (arg == "-o" || arg == "--output")
            {
                opt.output = value(arg);
            }
            else if (arg == "-l" || arg == "--target-len")
            {
                opt.target_len = std::stoi(value(arg));
            }
            else if (arg == "-s" || arg == "--separator")
            {
                opt.separator = value(arg);
            }
            else if (arg == "-t" || arg == "--threads")
            {
                opt.threads = std::stoi(value(arg));
            }
            else if (arg == "--max-copies")
            {
                opt.max_copies = std::stoi(value(arg));
            }
            else if (arg == "--batch-size")
            {
                opt.batch_size = std::stoi(value(arg));
            }
            else if (arg == "--seed")
            {
                opt.seed = static_cast<std::uint32_t>(std::stoul(value(arg)));
            }
            else if (arg == "-h" || arg == "--help")
            {
                print_usage(argv[0]);
                std::exit(0);
            }
            else
            {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

        if (opt.input.empty() || opt.output.empty() || opt.target_len <= 0)
        {
            throw std::runtime_error("-i, -o and a positive -l are required");
        }
        if (opt.separator.empty() || opt.threads <= 0 || opt.max_copies <= 0 ||
            opt.batch_size <= 0)
        {
            throw std::runtime_error(
                "separator must be nonempty and numeric options must be positive");
        }
        return opt;
    }

    std::string trim(std::string text)
    {
        auto not_space = [](unsigned char c)
        { return !std::isspace(c); };
        text.erase(text.begin(),
                   std::find_if(text.begin(), text.end(), not_space));
        text.erase(std::find_if(text.rbegin(), text.rend(), not_space).base(),
                   text.end());
        return text;
    }

    bool is_separator(const std::string &line, const std::string &separator)
    {
        // Accept both "====" and a longer all-equals separator when -s "====" is
        // supplied, matching the user's existing datasets.
        return line == separator || line.rfind(separator, 0) == 0;
    }

    std::uint32_t cluster_seed(std::uint32_t base, std::uint64_t cluster_id)
    {
        std::uint64_t x = static_cast<std::uint64_t>(base) +
                          0x9E3779B97F4A7C15ULL * (cluster_id + 1);
        x ^= x >> 30;
        x *= 0xBF58476D1CE4E5B9ULL;
        x ^= x >> 27;
        x *= 0x94D049BB133111EBULL;
        x ^= x >> 31;
        return static_cast<std::uint32_t>(x);
    }

    std::vector<std::string> limit_copies(const std::vector<std::string> &reads,
                                          int max_copies,
                                          std::uint32_t seed)
    {
        if (static_cast<int>(reads.size()) <= max_copies)
        {
            return reads;
        }
        std::vector<std::string> selected;
        selected.reserve(static_cast<std::size_t>(max_copies));
        std::mt19937 generator(seed);
        std::sample(reads.begin(), reads.end(), std::back_inserter(selected),
                    static_cast<std::size_t>(max_copies), generator);
        return selected;
    }

    std::string reconstruct_cluster(const std::vector<std::string> &input_reads,
                                    int target_len,
                                    int max_copies,
                                    std::uint32_t seed)
    {
        std::vector<std::string> reads = limit_copies(input_reads, max_copies, seed);
        if (reads.empty())
        {
            throw std::runtime_error("empty cluster");
        }
        if (reads.size() == 1)
        {
            return reads.front();
        }

        std::string unused_original;
        Cluster cluster(unused_original, reads);
        std::mt19937 generator(seed ^ 0xA5A5A5A5u);
        constexpr int rounds = 1; // Same R value used by the authors.
        return MinSumEDSimpleCorrectedClusterTwiceJoinR(
            cluster, target_len, generator, FixedLenMarkov, rounds);
    }

    bool process_batch(const std::vector<std::vector<std::string>> &clusters,
                       std::uint64_t first_cluster_id,
                       const Options &opt,
                       std::ofstream &output,
                       std::string &error_message)
    {
        std::vector<std::string> predictions(clusters.size());
        std::vector<std::string> errors(clusters.size());
        std::atomic<bool> failed{false};

#pragma omp parallel for schedule(dynamic, 1)
        for (long long i = 0; i < static_cast<long long>(clusters.size()); ++i)
        {
            try
            {
                const std::uint64_t cluster_id =
                    first_cluster_id + static_cast<std::uint64_t>(i);
                predictions[static_cast<std::size_t>(i)] = reconstruct_cluster(
                    clusters[static_cast<std::size_t>(i)], opt.target_len,
                    opt.max_copies, cluster_seed(opt.seed, cluster_id));
            }
            catch (const std::exception &ex)
            {
                failed.store(true, std::memory_order_relaxed);
                errors[static_cast<std::size_t>(i)] = ex.what();
            }
            catch (...)
            {
                failed.store(true, std::memory_order_relaxed);
                errors[static_cast<std::size_t>(i)] = "unknown exception";
            }
        }

        if (failed.load(std::memory_order_relaxed))
        {
            for (std::size_t i = 0; i < errors.size(); ++i)
            {
                if (!errors[i].empty())
                {
                    error_message = "cluster " +
                                    std::to_string(first_cluster_id + i + 1) +
                                    ": " + errors[i];
                    break;
                }
            }
            return false;
        }

        for (const std::string &prediction : predictions)
        {
            output << prediction << '\n';
        }
        return static_cast<bool>(output);
    }

} // namespace

int main(int argc, char **argv)
{
    try
    {
        const Options opt = parse_args(argc, argv);
        const auto start = std::chrono::steady_clock::now();

        std::ifstream input(opt.input);
        if (!input)
        {
            throw std::runtime_error("Cannot open input file: " + opt.input);
        }
        std::ofstream output(opt.output);
        if (!output)
        {
            throw std::runtime_error("Cannot open output file: " + opt.output);
        }

        omp_set_dynamic(0);
        omp_set_num_threads(opt.threads);
        int actual_threads = 1;
#pragma omp parallel
        {
#pragma omp single
            actual_threads = omp_get_num_threads();
        }

        std::vector<std::vector<std::string>> batch;
        batch.reserve(static_cast<std::size_t>(opt.batch_size));
        std::vector<std::string> current;
        std::uint64_t total_clusters = 0;

        auto flush_batch = [&]()
        {
            if (batch.empty())
            {
                return;
            }
            std::string error;
            if (!process_batch(batch, total_clusters, opt, output, error))
            {
                if (error.empty())
                {
                    error = "failed while writing predictions";
                }
                throw std::runtime_error(error);
            }
            total_clusters += batch.size();
            batch.clear();
        };

        auto finish_cluster = [&]()
        {
            if (current.empty())
            {
                return;
            }
            batch.push_back(std::move(current));
            current.clear();
            if (static_cast<int>(batch.size()) >= opt.batch_size)
            {
                flush_batch();
            }
        };

        std::string line;
        while (std::getline(input, line))
        {
            line = trim(std::move(line));
            if (is_separator(line, opt.separator))
            {
                finish_cluster();
            }
            else if (!line.empty())
            {
                current.push_back(std::move(line));
            }
        }
        finish_cluster();
        flush_batch();
        output.close();

        const double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          start)
                .count();
        rusage usage{};
        getrusage(RUSAGE_SELF, &usage);

        std::cout << "clusters: " << total_clusters << '\n';
        std::cout << "threads: " << actual_threads << '\n';
        std::cout << "max copies per cluster: " << opt.max_copies << '\n';
        std::cout << "written output: " << opt.output << '\n';
        std::cout << "elapsed seconds: " << elapsed << '\n';
        std::cout << "Maximum resident set size (kbytes): "
                  << usage.ru_maxrss << '\n';
        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "error: " << ex.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
