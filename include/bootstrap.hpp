#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <vector>
#include <random>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include "findSigExposures.hpp"

inline bool is_wholenumber(double val, double tol = 1e-8) {
    return std::fabs(val - std::round(val)) < tol;
}

std::vector<std::vector<double>> bootstraped_patient(
    const std::vector<double>& m, int mutation_count, int R
);

std::vector<double> compute_p_value(
    const std::vector<std::vector<double>>& exposures, double threshold
);

std::pair<std::vector<int>, std::pair<std::vector<std::vector<double>>, std::vector<double>>>
backward_elimination(
    const std::vector<double>& m,
    const std::vector<std::vector<double>>& P,
    int R,
    double threshold,
    int mutation_count,
    double significance_level,
    DecompositionMethod decomposition_method = decomposeQP
);

#endif