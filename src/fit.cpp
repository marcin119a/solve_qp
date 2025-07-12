#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <iomanip>
#include "bootstrap.hpp"
#include "csv_utils.hpp"

using namespace std;

std::pair<int, std::vector<int>> process_sample(
    int index,
    const std::vector<double>& col,
    const std::vector<std::vector<double>>& sigs,
    double threshold,
    int mutation_count,
    int R,
    double significance_level,
    std::pair<std::vector<std::vector<double>>, std::vector<double>>& estimation_exposures
) {
    try {
        auto [best_columns, result] = backward_elimination(col, sigs, R, threshold, mutation_count, significance_level);
        estimation_exposures = result;
        return {index, best_columns};
    } catch (const std::exception& e) {
        std::cerr << "Error processing sample " << index << ": " << e.what() << std::endl;
        return {index, {}};
    }
}

void fit(
    const std::string& samples_file,
    const std::string& output_folder,
    double threshold = 0.01,
    int mutation_count = -1,
    int R = 100,
    double significance_level = 0.01,
    const std::string& signatures_file = "data/COSMIC_v2_SBS_GRCh37.txt",
    bool drop_zeros_columns = false
) {
    auto [samples, patient_names] = load_samples_csv(samples_file);
    auto [sigs, sig_names] = load_signatures_csv(signatures_file);

    size_t G = samples[0].size();    // liczba pacjentów
    size_t N = sigs[0].size();       // liczba sygnatur

    std::vector<std::vector<std::string>> output(G + 1, std::vector<std::string>(N + 1, "0"));

    // Nagłówki
    output[0][0] = "Patient";
    for (size_t j = 0; j < N; ++j) output[0][j + 1] = sig_names[j];

    for (size_t i = 0; i < G; ++i) {
        std::vector<double> column(samples.size());
        for (size_t k = 0; k < samples.size(); ++k) column[k] = samples[k][i];

        std::pair<std::vector<std::vector<double>>, std::vector<double>> exposures_result;
        auto [_, best_columns] = process_sample(i, column, sigs, threshold, mutation_count, R, significance_level, exposures_result);

        output[i + 1][0] = patient_names[i];
        if (!best_columns.empty()) {
            for (size_t j = 0; j < best_columns.size(); ++j) {
                int col_index = best_columns[j];
                double val = exposures_result.first[j][0];
                output[i + 1][col_index + 1] = std::to_string(val);
            }
        }

        // Pasek postępu
        float percent = static_cast<float>(i + 1) / G;
        printf("\r[%-20s] %d%%", std::string(int(20 * percent), '=').c_str(), int(100 * percent));
        fflush(stdout);
    }
    std::cout << std::endl;

    // Zapis CSV
    std::filesystem::create_directories(output_folder);
    std::ofstream file(output_folder + "/Assignment_Solution_Activities.csv");

    if (!file) {
        std::cerr << "Failed to open output file!" << std::endl;
        return;
    }

    for (const auto& row : output) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j < row.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();
    std::cout << "✅ Saved: " << output_folder + "/Assignment_Solution_Activities.csv" << std::endl;
}