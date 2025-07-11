#ifndef DECOMPOSE_QP_HPP
#define DECOMPOSE_QP_HPP

#include <vector>
#include "solve_qp.hpp"

std::vector<double> decomposeQP(const std::vector<double>& m, const std::vector<std::vector<double>>& P);

#endif