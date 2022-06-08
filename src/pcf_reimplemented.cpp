#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
struct Aggregates {
    std::vector<double> aggregates;
    std::vector<int> pos;
    std::vector<int> pos_diff;
};

struct MultiAggregates {
    NumericMatrix aggregates;
    std::vector<int> pos;
    std::vector<int> pos_diff;
};

/* 
 * Checks that the index list 'r' is valid for the data 'y'.
 * Returns an index list that is in bounds, and fully spans the data.
 */
void sanitise_indices(const std::vector<double> &y, std::vector<int> &r) {
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
}

Aggregates make_aggregates(const std::vector<double> &y, std::vector<int> &r) {
    
    sanitise_indices(y, r);
    
    std::vector<double> agg;
    agg.reserve(r.size() - 1);
    
    std::vector<int> pos_diff;
    pos_diff.reserve(r.size() - 1);
    
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
        pos_diff.push_back(end - start);
    }
    
    return Aggregates{agg, r, pos_diff};
}

MultiAggregates make_multi_aggregates(const NumericMatrix &y, std::vector<int> &r) {
    
    NumericVector y0 = y(_, 0);
    sanitise_indices(Rcpp::as<std::vector<double>>(y0), r);
    
    int samples = y.ncol();
    NumericMatrix agg(r.size() - 1, samples);
    
    std::vector<int> pos_diff;
    pos_diff.reserve(r.size() - 1);
    
    int j = 0; // row index in agg, for insertion
    for (auto it = r.begin(), it2 = std::next(it);
         it != r.end() && it2 != r.end();
         ++it, ++it2) {
        NumericVector a(samples);
        int start = *it;
        int end = *it2;
        for (int i = start; i < end; ++i) {
            a += y(i, _);
        }
        agg(j++, _) = a;
        pos_diff.push_back(end - start);
    }
    
    return MultiAggregates{agg, r, pos_diff};
}


// [[Rcpp::export]]
std::vector<int> exact_multipcf_(const NumericMatrix &y, unsigned int kmin, double gamma) {
    std::size_t N = y.nrow();
    std::size_t samples = y.ncol();
    NumericMatrix A(N, samples);
    double D = 0;
    std::vector<double> S(N, 0); // Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
            A(j, _) = A(j, _) + y(k, _);
            if (j > 0 && (j < kmin || k + 1 - j < kmin)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                D = sum(-1 * A(j, _) * A(j, _) / (k - j + 1));
                S[j] = D + E[j] + gamma;
            }
        }
        
        auto min_element = std::min_element(S.begin(), S.begin() + k + 1);
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


void print_matrix(const NumericMatrix &m) {
    for (int r = 0; r < m.nrow(); ++r) {
        for (int c = 0; c < m.ncol(); ++c) {
            double elem = m(r, c);
            std::cout << elem << ',';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

// [[Rcpp::export]]
std::vector<int> fast_multipcf_(const NumericMatrix &y, std::vector<int> &available_breakpoints, int kmin, double gamma) {
    
    int samples = y.ncol();
    MultiAggregates agg = make_multi_aggregates(y, available_breakpoints);
    const NumericMatrix &u = agg.aggregates;
    const std::vector<int> &r = agg.pos;
    
    std::size_t N = u.nrow();
    NumericMatrix A(N, samples);
    std::vector<int> C(N, 0);
    double D;
    std::vector<double> S(N, 0); // Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
            A(j, _) = A(j, _) + u(k, _);
            C[j] += agg.pos_diff[k];
            if (r[j] > 0 && (r[j] < kmin || r[k + 1] - r[j] < kmin)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                D = sum(-1 * A(j, _) * A(j, _) / C[j]);
                S[j] = D + E[j] + gamma;
            }
        }
        
        auto min_element = std::min_element(S.begin(), S.begin() + k + 1);
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


// [[Rcpp::export]]
std::vector<int> exact_pcf_(const std::vector<double> &y, unsigned int kmin, double gamma) {
    std::size_t N = y.size();
    std::vector<double> A(N, 0);
    double D = 0;
    std::vector<double> S(N, 0);// Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
            A[j] += y[k];
            if (j > 0 && (j < kmin || k + 1 - j < kmin)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                D = -A[j] * A[j] / (k - j + 1);
                S[j] = D + E[j] + gamma;
            }
        }
        
        auto min_element = std::min_element(S.begin(), S.begin() + k + 1);
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

template<typename T>
void print_vec(std::string label, const std::vector<T> &v) {
    std::cout << label << " [";
    for (size_t i = 0; i < v.size() - 1; ++i) {
        std::cout << v[i] << ", ";
    }
    std::cout << v.back() << "]" << std::endl;
}

// [[Rcpp::export]]
std::vector<int> fast_pcf_(const std::vector<double> &y, std::vector<int> &available_breakpoints, int kmin, double gamma) {
    Aggregates agg = make_aggregates(y, available_breakpoints);
    const std::vector<double> &u = agg.aggregates;
    const std::vector<int> &r = agg.pos;
    
    std::size_t N = u.size();
    std::vector<double> A(N, 0);
    std::vector<int> C(N, 0);
    double D;
    std::vector<double> S(N, 0); // Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
            A[j] += u[k];
            C[j] += agg.pos_diff[k];
            if (r[j] > 0 && (r[j] < kmin || r[k + 1] - r[j] < kmin)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                D = -A[j] * A[j] / C[j];
                S[j] = D + E[j] + gamma;
            }
        }
        
        auto min_element = std::min_element(S.begin(), S.begin() + k + 1);
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

// [[Rcpp::export]]
std::vector<double> convolve_(const std::vector<double>& x, const std::vector<double>& k) {
    int nx = x.size();
    int nk = k.size();
    std::vector<double> out(nx + nk - 1, 0);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < nk; ++j) {
            out[i+j] += x[i] * k[j];
        }
    }
    return std::vector<double>(out.begin() + nk / 2 - 1, out.end() - nk / 2);
}

// [[Rcpp::export]]
double median_(NumericVector x) {
    NumericVector y = clone(x);
    int n, half;
    double y1, y2;
    n = y.size();
    half = n / 2;
    if(n % 2 == 1) {
        // median for odd length vector
        std::nth_element(y.begin(), y.begin()+half, y.end());
        return y[half];
    } else {
        // median for even length vector
        std::nth_element(y.begin(), y.begin()+half, y.end());
        y1 = y[half];
        std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
        y2 = y[half-1];
        return (y1 + y2) / 2.0;
    }
}

// [[Rcpp::export]]
double mad_(NumericVector x, double scale_factor = 1.4826) {
    // scale_factor = 1.4826; default for normal distribution consistent with R
    return median_(abs(x - median_(x))) * scale_factor;
}

// [[Rcpp::export]]
std::vector<int> mark_(const std::vector<double>& x, double nmad = 1.0, int filter_size = 4) {
    // Make the smoothing sawtooth filter for edge detection
    std::vector<double> k;
    for (int i = 0; i < filter_size; ++i) {
        k.push_back(-1);
    }
    for (int i = 0; i < 2 * filter_size; ++i) {
        k.push_back(-2);
    }
    for (int i = 0; i < 2 * filter_size; ++i) {
        k.push_back(2);
    }
    for (int i = 0; i < filter_size; ++i) {
        k.push_back(1);
    }
    
    NumericVector hpf = Rcpp::wrap(convolve_(x, k));
    NumericVector abshpf = abs(hpf);
    double threshold = median_(abshpf) + nmad * mad_(hpf);
    std::vector<int> out;
    for (int i = 0; i < abshpf.size(); ++i) {
        if (abshpf[i] > threshold) {
            out.push_back(i);
        }
    }
    return out;
}