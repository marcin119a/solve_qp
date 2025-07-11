#include <vector>
#include <stdexcept>
#include <functional>
#include <numeric>
#include <cmath>
#include <iostream>
#include "solve_qp.hpp"
#include "decomposeQP.hpp"

using DecompositionMethod = std::function<std::vector<double>(const std::vector<double>&, const std::vector<std::vector<double>>&)>;

// Oblicz normę Frobeniusa (błąd rekonstrukcji między m a P * x)
double FrobeniusNorm(const std::vector<double>& m, const std::vector<std::vector<double>>& P, const std::vector<double>& x) {
    int rows = m.size();
    int cols = x.size();
    double err = 0.0;

    for (int i = 0; i < rows; ++i) {
        double approx = 0.0;
        for (int j = 0; j < cols; ++j) {
            approx += P[i][j] * x[j];
        }
        double diff = m[i] - approx;
        err += diff * diff;
    }

    return std::sqrt(err);
}

// Zwraca: exposures (N x G), errors (rozmiar G)
std::pair<std::vector<std::vector<double>>, std::vector<double>> findSigExposures(
    const std::vector<std::vector<double>>& M,  // 96 x G
    const std::vector<std::vector<double>>& P,  // 96 x N
    DecompositionMethod decomposition_method = decomposeQP
) {
    int rows = M.size();
    int G = M[0].size();     // liczba pacjentów
    int N = P[0].size();     // liczba sygnatur

    if (P.size() != rows)
        throw std::invalid_argument("Matrices 'M' and 'P' must have the same number of rows (mutation types).");

    if (N < 2)
        throw std::invalid_argument("Matrix 'P' must have at least 2 columns (signatures).");

    // Normalizacja kolumn M
    std::vector<std::vector<double>> Mnorm = M;
    for (int g = 0; g < G; ++g) {
        double col_sum = 0.0;
        for (int i = 0; i < rows; ++i)
            col_sum += M[i][g];
        if (col_sum > 0.0) {
            for (int i = 0; i < rows; ++i)
                Mnorm[i][g] /= col_sum;
        }
    }

    // Oblicz ekspozycje i błędy
    std::vector<std::vector<double>> exposures(N, std::vector<double>(G, 0.0));
    std::vector<double> errors(G, 0.0);

    for (int g = 0; g < G; ++g) {
        std::vector<double> m_g(rows);
        for (int i = 0; i < rows; ++i)
            m_g[i] = Mnorm[i][g];

        std::vector<double> x = decomposition_method(m_g, P);
        for (int j = 0; j < N; ++j)
            exposures[j][g] = x[j];

        errors[g] = FrobeniusNorm(m_g, P, x);
    }

    return {exposures, errors};
}