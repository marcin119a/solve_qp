#include "solve_qp.hpp"
#include <stdexcept>

extern "C" int qpgen2_(double *G, double *av, int n,
                       double *xv, double *lagr, double *obj,
                       double *C, double *bv, int q, int meq,
                       int* iact, int* nact, int* iter,
                       double* work, int factorized);

QPResult solve_qp(
    const std::vector<std::vector<double>>& G,
    const std::vector<double>& a,
    const std::vector<std::vector<double>>& C,
    const std::vector<double>& b,
    int meq,
    bool factorized
) {
    int n = G.size();
    if (G[0].size() != n) throw std::invalid_argument("G must be square");

    std::vector<std::vector<double>> C_matrix = C;
    std::vector<double> b_vector = b;
    int m = C.empty() ? 0 : C[0].size();

    if (C.empty()) {
        C_matrix = std::vector<std::vector<double>>(n, std::vector<double>(1, 0.0));
        b_vector = std::vector<double>(1, -1.0);
        m = 1;
    }

    if ((int)a.size() != n)
        throw std::invalid_argument("a must have same dimension as G");
    if ((int)C_matrix.size() != n)
        throw std::invalid_argument("C must have same row count as G");
    if ((int)b_vector.size() != m)
        throw std::invalid_argument("b must match columns of C");

    std::vector<double> G_flat(n * n);
    std::vector<double> C_flat(n * m);
    std::vector<double> a_copy = a;
    std::vector<double> x(n, 0.0);
    std::vector<double> lagr(m, 0.0);
    std::vector<double> xu(n);
    std::vector<int> iact(m, 0);
    int nact = 0;
    int iter[2] = {0, 0};
    double obj = 0.0;
    int r = std::min(n, m);
    std::vector<double> work(2 * n + 2 * m + r * (r + 5) / 2, 0.0);

    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            G_flat[j * n + i] = G[i][j];

    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            C_flat[j * n + i] = C_matrix[i][j];

    int factorized_ = factorized ? 1 : 0;

    int result = qpgen2_(
        G_flat.data(), a_copy.data(), n,
        x.data(), lagr.data(), &obj,
        C_flat.data(), const_cast<double*>(b_vector.data()), m, meq,
        iact.data(), &nact, iter,
        work.data(), factorized_
    );

    if (result == 1)
        throw std::runtime_error("Infeasible constraints");
    if (result == 2)
        throw std::runtime_error("Matrix G is not positive definite");

    QPResult res;
    res.x = x;
    res.f = obj;
    res.xu = a_copy;
    res.iterations = {iter[0], iter[1]};
    res.lagrangian = lagr;
    res.iact = std::vector<int>(iact.begin(), iact.begin() + nact);
    return res;
}