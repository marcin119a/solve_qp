#include "decomposeQP.hpp"

std::vector<double> decomposeQP(const std::vector<double>& m, const std::vector<std::vector<double>>& P) {
    int N = P[0].size();  // liczba sygnatur (kolumn)
    int M = P.size();     // liczba typów mutacji (rzędów)

    // G = P^T * P
    std::vector<std::vector<double>> G(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < M; ++k)
                G[i][j] += P[k][i] * P[k][j];

    // d = P^T * m
    std::vector<double> d(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int k = 0; k < M; ++k)
            d[i] += m[k] * P[k][i];

    // Macierz ograniczeń C = [1; I]
    std::vector<std::vector<double>> C(N, std::vector<double>(N + 1, 0.0));
    for (int i = 0; i < N; ++i) {
        C[i][0] = 1.0;
        C[i][i + 1] = 1.0;
    }

    // b = [1, 0, ..., 0]
    std::vector<double> b(N + 1, 0.0);
    b[0] = 1.0;

    // Rozwiązanie QP
    QPResult res = solve_qp(G, d, C, b, 1);

    // Korekta ujemnych wartości i normalizacja
    for (double& x : res.x) {
        if (x < 0.0) x = 0.0;
    }

    double sum = 0.0;
    for (double x : res.x) sum += x;
    for (double& x : res.x) x /= sum;

    return res.x;
}