#include "bootstrap.hpp"

std::vector<std::vector<double>> bootstraped_patient(const std::vector<double>& m, int mutation_count, int R) {
    int K = m.size();
    double sum = std::accumulate(m.begin(), m.end(), 0.0);
    auto is_whole_with_tol = [](double val) {
        return std::fabs(val - std::round(val)) < 1e-8;
    };

    if (mutation_count <= 0) {
        if (std::all_of(m.begin(), m.end(), is_whole_with_tol)) {
            mutation_count = static_cast<int>(std::round(sum));
        } else {
            throw std::invalid_argument("Provide mutation_count or ensure m contains counts.");
        }
    }

    std::vector<double> prob(m);
    for (double& p : prob) p /= sum;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> dist(prob.begin(), prob.end());

    std::vector<std::vector<double>> M(K, std::vector<double>(R, 0.0));

    for (int r = 0; r < R; ++r) {
        std::vector<int> counts(K, 0);
        for (int i = 0; i < mutation_count; ++i) {
            counts[dist(gen)]++;
        }
        for (int k = 0; k < K; ++k) {
            M[k][r] = static_cast<double>(counts[k]) / mutation_count;
        }
    }

    return M;
}

std::vector<double> compute_p_value(const std::vector<std::vector<double>>& exposures, double threshold) {
    int N = exposures.size();     // sygnatury
    int R = exposures[0].size();  // próbki bootstrapowe

    std::vector<double> p_vals(N, 0.0);
    for (int i = 0; i < N; ++i) {
        int count = 0;
        for (int j = 0; j < R; ++j) {
            if (exposures[i][j] > threshold) count++;
        }
        p_vals[i] = 1.0 - static_cast<double>(count) / R;
    }

    return p_vals;
}

std::pair<std::vector<int>, std::pair<std::vector<std::vector<double>>, std::vector<double>>>
backward_elimination(const std::vector<double>& m,
                     const std::vector<std::vector<double>>& P,
                     int R,
                     double threshold,
                     int mutation_count,
                     double significance_level,
                     DecompositionMethod decomposition_method) {
    int N = P[0].size();  // liczba sygnatur

    std::vector<int> best_columns(N);
    std::iota(best_columns.begin(), best_columns.end(), 0);

    std::vector<std::vector<double>> P_temp = P;
    std::vector<std::vector<double>> M = bootstraped_patient(m, mutation_count, R);

    while (true) {
        auto [exposures, errors] = findSigExposures(M, P_temp, decomposition_method);
        std::vector<double> p_values = compute_p_value(exposures, threshold);

        auto max_it = std::max_element(p_values.begin(), p_values.end());
        double max_p_value = *max_it;

        if (max_p_value > significance_level) {
            int idx = std::distance(p_values.begin(), max_it);
            best_columns.erase(best_columns.begin() + idx);

            std::vector<std::vector<double>> P_new(P.size(), std::vector<double>(best_columns.size()));
            for (size_t i = 0; i < P.size(); ++i)
                for (size_t j = 0; j < best_columns.size(); ++j)
                    P_new[i][j] = P[i][best_columns[j]];
            P_temp = P_new;
        } else {
            break;
        }
    }

    // Final exposure estimation
    std::vector<std::vector<double>> m_matrix(m.size(), std::vector<double>(1));
    for (size_t i = 0; i < m.size(); ++i)
        m_matrix[i][0] = m[i];

    auto final_result = findSigExposures(m_matrix, P_temp, decomposition_method);
    return {best_columns, final_result};
}