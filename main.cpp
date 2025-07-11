#include <iostream>
#include "include/bootstrap.hpp"

int main() {
    std::vector<double> m = {50, 90, 70};  // surowe mutacje
    std::vector<std::vector<double>> P = {
        {0.2, 0.3, 0.5},
        {0.1, 0.4, 0.5},
        {0.3, 0.1, 0.6}
    };

    int R = 100;
    double threshold = 0.01;
    int mutation_count = 210;
    double significance_level = 0.1;

    auto [selected_columns, result] = backward_elimination(
        m, P, R, threshold, mutation_count, significance_level
    );

    std::cout << "Selected signatures: ";
    for (int idx : selected_columns) std::cout << idx << " ";
    std::cout << std::endl;

    std::cout << "Final exposures: ";
    for (const auto& val : result.first) {
        std::cout << val[0] << " ";
    }
    std::cout << std::endl;

    std::cout << "Final error: " << result.second[0] << std::endl;

    return 0;
}