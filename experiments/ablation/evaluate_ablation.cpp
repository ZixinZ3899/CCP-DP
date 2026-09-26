#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

static std::string clean(const std::string &line)
{
    std::string result;
    for (char c : line)
    {
        const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' || upper == 'N')
            result.push_back(upper);
    }
    return result;
}

static std::vector<std::string> read_sequences(const std::string &path, int limit)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open " + path);
    std::vector<std::string> sequences;
    std::string line, current;
    bool fasta = false;
    while (std::getline(in, line))
    {
        if (!line.empty() && line[0] == '>')
        {
            fasta = true;
            if (!current.empty())
            {
                sequences.push_back(current);
                current.clear();
                if (limit > 0 && static_cast<int>(sequences.size()) >= limit)
                    break;
            }
            continue;
        }
        const std::string sequence = clean(line);
        if (sequence.empty())
            continue;
        if (fasta)
            current += sequence;
        else
        {
            sequences.push_back(sequence);
            if (limit > 0 && static_cast<int>(sequences.size()) >= limit)
                break;
        }
    }
    if (fasta && !current.empty() &&
        (limit <= 0 || static_cast<int>(sequences.size()) < limit))
        sequences.push_back(current);
    return sequences;
}

static int edit_distance(const std::string &left, const std::string &right)
{
    const std::string *a = &left;
    const std::string *b = &right;
    if (a->size() < b->size())
        std::swap(a, b);
    std::vector<int> previous(b->size() + 1), current(b->size() + 1);
    for (size_t j = 0; j <= b->size(); ++j)
        previous[j] = static_cast<int>(j);
    for (size_t i = 1; i <= a->size(); ++i)
    {
        current[0] = static_cast<int>(i);
        for (size_t j = 1; j <= b->size(); ++j)
            current[j] = std::min({
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + ((*a)[i - 1] != (*b)[j - 1])});
        previous.swap(current);
    }
    return previous[b->size()];
}

int main(int argc, char **argv)
{
    try
    {
        if (argc < 3 || argc > 4)
            throw std::runtime_error("Usage: evaluate_ablation CENTERS PREDICTIONS [LIMIT]");
        const int limit = argc == 4 ? std::stoi(argv[3]) : -1;
        const auto centers = read_sequences(argv[1], limit);
        const auto predictions = read_sequences(argv[2], limit);
        if (centers.size() != predictions.size() || centers.empty())
            throw std::runtime_error("Sequence-count mismatch or empty input");
        long long exact = 0;
        long long total_ed = 0;
        long long total_bases = 0;
        for (size_t i = 0; i < centers.size(); ++i)
        {
            exact += centers[i] == predictions[i];
            total_ed += edit_distance(centers[i], predictions[i]);
            total_bases += static_cast<long long>(centers[i].size());
        }
        const double n = static_cast<double>(centers.size());
        std::cout << std::setprecision(15)
                  << "sequences=" << centers.size() << '\n'
                  << "exact_count=" << exact << '\n'
                  << "success_percent=" << 100.0 * exact / n << '\n'
                  << "mean_ed=" << total_ed / n << '\n'
                  << "reconstruction_percent="
                  << 100.0 * (1.0 - static_cast<double>(total_ed) / total_bases)
                  << '\n';
    }
    catch (const std::exception &error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
