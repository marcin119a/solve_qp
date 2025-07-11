#ifndef SOLVE_QP_HPP
#define SOLVE_QP_HPP

#include <vector>

struct QPResult {
    std::vector<double> x;
    double f;
    std::vector<double> xu;
    std::vector<int> iterations;
    std::vector<double> lagrangian;
    std::vector<int> iact;
};

QPResult solve_qp(
    const std::vector<std::vector<double>>& G,
    const std::vector<double>& a,
    const std::vector<std::vector<double>>& C = {},
    const std::vector<double>& b = {},
    int meq = 0,
    bool factorized = false
);

#endif