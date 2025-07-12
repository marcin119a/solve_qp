#include <iostream>
#include <chrono>
#include "fit.hpp"

int main() {
    std::string samples_path = "/Users/mw/CLionProjects/solver_test/data/tumorBRCA.txt";
    std::string signatures_path = "/Users/mw/CLionProjects/solver_test/data/COSMIC_v2_SBS_GRCh37.txt";
    std::string output_folder = "output";

    double threshold = 0.01;
    int mutation_count = 1000; // -1: automatycznie wykrywany
    int R = 100;
    double significance_level = 0.01;

    auto start = std::chrono::high_resolution_clock::now();

    fit(
        samples_path,
        output_folder,
        threshold,
        mutation_count,
        R,
        significance_level,
        signatures_path,
        false  // drop_zeros_columns
    );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "⏱️ Execution time: " << elapsed_seconds.count() << "s\n";

    return 0;
}