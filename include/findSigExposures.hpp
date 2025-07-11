#ifndef FIND_SIG_EXPOSURES_HPP
#define FIND_SIG_EXPOSURES_HPP

#include <vector>
#include <functional>
#include "solve_qp.hpp"
#include "decomposeQP.hpp"

using DecompositionMethod = std::function<std::vector<double>(
    const std::vector<double>&,
    const std::vector<std::vector<double>>&
)>;

double FrobeniusNorm(const std::vector<double>& m,
                     const std::vector<std::vector<double>>& P,
                     const std::vector<double>& x);

std::pair<std::vector<std::vector<double>>, std::vector<double>> findSigExposures(
    const std::vector<std::vector<double>>& M,
    const std::vector<std::vector<double>>& P,
    DecompositionMethod decomposition_method = decomposeQP
);

#endif