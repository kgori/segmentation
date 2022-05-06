#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
std::vector<int> exact_pcf_(const std::vector<double> &y, unsigned int kmin, double gamma) {
    std::size_t N = y.size();
    std::vector<double> A(N, 0);
    std::vector<double> D(N, 0);
    std::vector<double> S(N, 0);// Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
            A[j] += y[k];
            if (((j > 0 && j < kmin) || j > k + 1 - kmin)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                D[j] = -A[j] * A[j] / (k - j + 1);
                S[j] = D[j] + E[j] + gamma;
            }
        }
        
        auto min_element = std::min_element(S.begin(), S.begin() + k);
        auto min_position = static_cast<int>(std::distance(S.begin(), min_element));
        auto min_value = *min_element;
        E[k + 1] = min_value;
        T[k] = min_position;
    }
    
    // Find start positions
    std::vector<int> starts;
    int pos = T.back();
    while (pos > 0) {
        starts.push_back(pos);
        pos = T[pos - 1];
    }
    starts.push_back(0);
    std::reverse(starts.begin(), starts.end());
    return starts;
}

struct Aggregates {
    std::vector<double> aggregates;
    std::vector<int> pos;
};

Aggregates make_aggregates(const std::vector<double> &y, std::vector<int> r) {
    auto min_r = std::min_element(r.begin(), r.end());
    auto max_r = std::max_element(r.begin(), r.end());
    
    if (*min_r < 0) {
        throw std::runtime_error("Negative indices are not allowed");
    }
    
    if (*max_r > y.size()) {
        throw std::runtime_error("Max index in r is out of bounds of y");
    }
    
    if (*min_r > 0) {
        r.push_back(0);
    }
    
    if (*max_r < y.size() - 1) {
        r.push_back(y.size() - 1);
    }
    
    if (*max_r < y.size()) {
        r.push_back(y.size());
    }
    
    if (!std::is_sorted(r.begin(), r.end())) {
        std::sort(r.begin(), r.end());
    }
    
    std::vector<double> agg;
    agg.reserve(r.size());
    
    for (auto it = r.begin(), it2 = std::next(it);
         it != r.end() && it2 != r.end();
         ++it, ++it2) {
        double a{0};
        int start = *it;
        int end = *it2;
        for (int i = start; i < end; ++i) {
            a += y[i];
        }
        agg.push_back(a);
    }
    
    return Aggregates{agg, r};
}

// [[Rcpp::export]]
std::vector<int> fast_pcf_(const std::vector<double> &y, const std::vector<int> &available_breakpoints, int kmin, double gamma) {
    Aggregates agg = make_aggregates(y, available_breakpoints);
    const std::vector<double> &u = agg.aggregates;
    const std::vector<int> &r = agg.pos;
    
    std::size_t N = u.size();
    std::vector<double> A(N, 0);
    std::vector<int> C(N, 0);
    std::vector<double> D(N, 0);
    std::vector<double> S(N, 0);// Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
            A[j] += u[k];
            C[j] += r[k+1] - r[k];
            if (r[j] > 0 && (r[j] < kmin || r[k] - r[j] + 1 < kmin)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                D[j] = -A[j] * A[j] / C[j];
                S[j] = D[j] + E[j] + gamma;
            }
        }
        
        auto min_element = std::min_element(S.begin(), S.begin() + k);
        auto min_position = static_cast<int>(std::distance(S.begin(), min_element));
        auto min_value = *min_element;
        E[k + 1] = min_value;
        T[k] = min_position;
    }
    
    // Find start positions
    std::vector<int> starts;
    int pos = T.back();
    while (pos > 0) {
        starts.push_back(r[pos]);
        pos = T[pos - 1];
    }
    starts.push_back(0);
    std::reverse(starts.begin(), starts.end());
    return starts;
}