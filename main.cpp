#include <iostream>
#include <vector>
#include <cmath>
#include "solve_qp.hpp"

std::vector<double> decomposeQP(const std::vector<double>& m, const std::vector<std::vector<double>>& P) {
    int N = P[0].size();  // liczba sygnatur
    int M = P.size();     // długość wektora m

    // Oblicz G = P^T * P
    std::vector<std::vector<double>> G(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < M; ++k)
                G[i][j] += P[k][i] * P[k][j];

    // Oblicz d = P^T * m
    std::vector<double> d(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int k = 0; k < M; ++k)
            d[i] += m[k] * P[k][i];

    // Ograniczenia: C = [1, I]
    std::vector<std::vector<double>> C(N, std::vector<double>(N + 1, 0.0));
    for (int i = 0; i < N; ++i) {
        C[i][0] = 1.0;        // pierwsza kolumna: 1
        C[i][i + 1] = 1.0;    // kolejne: jednostkowa
    }

    // b = [1, 0, 0, ..., 0]
    std::vector<double> b(N + 1, 0.0);
    b[0] = 1.0;

    QPResult res = solve_qp(G, d, C, b, 1);

    // Poprawka na bardzo małe ujemne liczby
    for (double& x : res.x) {
        if (x < 0) x = 0.0;
    }

    // Normalizacja: sum(x) = 1
    double total = 0.0;
    for (double x : res.x) total += x;
    for (double& x : res.x) x /= total;

    return res.x;
}

int main() {
    std::vector<double> m = {0.3, 0.4, 0.3};
    std::vector<std::vector<double>> P = {
        {0.2, 0.4},
        {0.3, 0.3},
        {0.5, 0.3}
    };
    std::vector<double> exposures = decomposeQP(m, P);

    std::cout << "Decomposed exposures:" << std::endl;
    for (double x : exposures) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    return 0;
}