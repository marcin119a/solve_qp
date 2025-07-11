#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "findSigExposures.hpp"
#include "decomposeQP.hpp"

bool almost_equal(double a, double b, double tol = 1e-7) {
    return std::fabs(a - b) <= tol;
}

int main() {
    // Dane testowe (odpowiednik numpy array)
    std::vector<std::vector<double>> M = {
        {0.5, 0.3, 0.2},
        {0.9, 0.05, 0.05},
        {0.7, 0.1, 0.2}
    };

    std::vector<std::vector<double>> P = {
        {0.2, 0.3, 0.5},
        {0.1, 0.4, 0.5},
        {0.3, 0.1, 0.6}
    };

    // Uruchomienie funkcji
    auto [exposures, errors] = findSigExposures(M, P);

    // Oczekiwane wartości
    std::vector<std::vector<double>> expected_exposures = {
        {0.2007233, 0.4317862, 0.6437908},
        {0.4755877, 0.2883263, 0.0000000},
        {0.3236890, 0.2798875, 0.3562092}
    };

    std::vector<double> expected_errors = {
        0.1245907, 0.4137049, 0.1939068
    };

    // Sprawdzenie exposures
    bool success = true;
    for (size_t i = 0; i < exposures.size(); ++i) {
        for (size_t j = 0; j < exposures[0].size(); ++j) {
            if (!almost_equal(exposures[i][j], expected_exposures[i][j])) {
                std::cerr << "Exposure mismatch at (" << i << "," << j << "): "
                          << exposures[i][j] << " != " << expected_exposures[i][j] << std::endl;
                success = false;
            }
        }
    }

    // Sprawdzenie errors
    for (size_t i = 0; i < errors.size(); ++i) {
        if (!almost_equal(errors[i], expected_errors[i])) {
            std::cerr << "Error mismatch at [" << i << "]: "
                      << errors[i] << " != " << expected_errors[i] << std::endl;
            success = false;
        }
    }

    if (success) {
        std::cout << "✅ test_findSigExposures PASSED" << std::endl;
    } else {
        std::cerr << "❌ test_findSigExposures FAILED" << std::endl;
    }

    return success ? 0 : 1;
}