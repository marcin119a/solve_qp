#ifndef FIT_HPP
#define FIT_HPP

#include <string>
#include <vector>
#include <utility>

std::pair<int, std::vector<int>> process_sample(
    int index,
    const std::vector<double>& col,
    const std::vector<std::vector<double>>& sigs,
    double threshold,
    int mutation_count,
    int R,
    double significance_level,
    std::pair<std::vector<std::vector<double>>, std::vector<double>>& estimation_exposures
);

void fit(
    const std::string& samples_file,
    const std::string& output_folder,
    double threshold = 0.01,
    int mutation_count = -1,
    int R = 100,
    double significance_level = 0.01,
    const std::string& signatures_file = "data/COSMIC_v3.4_SBS_GRCh37.txt",
    bool drop_zeros_columns = false
);

#endif // FIT_HPP