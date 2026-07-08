#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
struct Aggregates {
    std::vector<double> aggregates;
    std::vector<std::size_t> pos;
    std::vector<std::size_t> pos_diff;
};

struct MultiAggregates {
    std::vector<double> aggregates;
    std::size_t samples;
    std::vector<std::size_t> pos;
    std::vector<std::size_t> pos_diff;
};

/*
 * Checks that the index list 'r' is valid for data of size n
 * Returns a validated an index list that is in bounds,
 * and fully spans the data.
 */
std::vector<std::size_t> sanitise_indices(std::size_t n, const std::vector<int> &r) {
    if (r.empty()) {
        if (n == 0) {
            return {0};
        }
        return {0, n};
    }

    const auto [min_it, max_it] = std::minmax_element(r.begin(), r.end());
    const int min_r = *min_it;
    const int max_r = *max_it;

    if (min_r < 0) {
        Rcpp::stop("Negative indices are not allowed");
    }

    if (static_cast<std::size_t>(max_r) > n) {
        Rcpp::stop("Max index in r is out of bounds of y");
    }

    std::vector<std::size_t> r_sanitised(r.begin(), r.end());

    if (min_r > 0) {
        r_sanitised.push_back(0);
    }

    if (static_cast<std::size_t>(max_r) < n) {
        r_sanitised.push_back(n);
    }

    std::sort(r_sanitised.begin(), r_sanitised.end());
    auto last = std::unique(r_sanitised.begin(), r_sanitised.end());
    r_sanitised.erase(last, r_sanitised.end());
    return r_sanitised;
}

Aggregates make_aggregates(const std::vector<double> &y, const std::vector<int> &breakpoints) {

    auto r = sanitise_indices(y.size(), breakpoints);

    std::vector<double> agg;
    agg.reserve(r.size() - 1);

    std::vector<std::size_t> pos_diff;
    pos_diff.reserve(r.size() - 1);

    for (auto it = r.begin(), it2 = std::next(it); it != r.end() && it2 != r.end(); ++it, ++it2) {
        double a{0};
        std::size_t start = *it;
        std::size_t end = *it2;
        for (std::size_t i = start; i < end; ++i) {
            a += y[i];
        }
        agg.push_back(a);
        pos_diff.push_back(end - start);
    }

    return Aggregates{agg, r, pos_diff};
}

/*
 * Aggregates the data in y over segments defined by breakpoint vector r
 */
MultiAggregates make_multi_aggregates(const NumericMatrix &y, const std::vector<int> &breakpoints) {
    const std::size_t nrow = y.nrow();
    const std::size_t samples = y.ncol();

    // Make sure r fully spans the data and is in bounds
    auto r = sanitise_indices(nrow, breakpoints);

    const std::size_t n_segs = r.size() - 1;
    std::vector<double> agg(n_segs * samples, 0.0); // row-major storage for the aggregated matrix
    std::vector<std::size_t> pos_diff;
    pos_diff.reserve(n_segs);

    // Iterate over columns
    for (std::size_t j = 0; j < samples; ++j) {
        auto col = y.column(j);
        for (std::size_t seg = 0; seg < n_segs; ++seg) {
            double segsum = 0.0;
            for (std::size_t i = r[seg]; i < r[seg + 1]; ++i) {
                segsum += col[i];
            }
            agg[seg * samples + j] = segsum;
        }
    }

    // An extra loop over nsegs doesn't seem like too much overhead
    for (std::size_t seg = 0; seg < n_segs; ++seg) {
        pos_diff.push_back(r[seg + 1] - r[seg]);
    }

    return MultiAggregates{agg, samples, r, pos_diff};
}

// [[Rcpp::export]]
std::vector<int> exact_multipcf_(const NumericMatrix &y, int kmin, double gamma) {
    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    const std::size_t N = y.nrow();
    const std::size_t samples = y.ncol();

    if (N == 0 || samples == 0) {
        return {};
    }

    // Transpose y into row-major format
    std::vector<double> Y(N * samples);
    for (std::size_t j = 0; j < samples; ++j) {
        auto col = y.column(j);
        for (std::size_t i = 0; i < N; ++i) {
            Y[i * samples + j] = col[i];
        }
    }

    // Matrix A is an accumulator for the sums of the rows of y, used to compute the score function.
    std::vector<double> A(N * samples, 0.0);
    std::vector<double> S(N, 0.0); // Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t j = 0; j <= k; ++j) {
            if (j > 0 && (j < kmin_size || k + 1 - j < kmin_size)) {
                // Do the sum in a manual loop instead of relying on Rcpp temporaries
                for (std::size_t s = 0; s < samples; ++s) {
                    A[j * samples + s] += Y[k * samples + s];
                }
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                const double inv = -1.0 / (k - j + 1);
                double D = 0.0;
                for (std::size_t s = 0; s < samples; ++s) {
                    A[j * samples + s] += Y[k * samples + s];
                    D += A[j * samples + s] * A[j * samples + s];
                }
                S[j] = inv * D + E[j] + gamma;
            }
        }

        const auto min_element = std::min_element(S.begin(), S.begin() + k + 1);
        const auto min_position = static_cast<int>(std::distance(S.begin(), min_element));
        const auto min_value = *min_element;
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

// [[Rcpp::export]]
std::vector<int> fast_multipcf_(const NumericMatrix &y,
                                const std::vector<int> &available_breakpoints,
                                int kmin, double gamma) {

    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    MultiAggregates agg = make_multi_aggregates(y, available_breakpoints);
    const std::vector<double> &u = agg.aggregates;
    const std::vector<std::size_t> &r = agg.pos;
    const std::size_t samples = agg.samples;
    const std::size_t N = r.size() - 1;

    if (N == 0 || samples == 0) {
        return {};
    }

    std::vector<double> A(N * samples, 0.0);
    std::vector<std::size_t> C(N, 0);
    std::vector<double> S(N, 0.0); // Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t j = 0; j <= k; ++j) {
            C[j] += agg.pos_diff[k];

            if (r[j] > 0 && (r[j] < kmin_size || r[k + 1] - r[j] < kmin_size)) {
                // Do the sum
                for (std::size_t s = 0; s < samples; ++s) {
                    A[j * samples + s] += u[k * samples + s];
                }
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                const double inv = -1.0 / C[j];
                double D = 0.0;
                for (std::size_t s = 0; s < samples; ++s) {
                    A[j * samples + s] += u[k * samples + s];
                    D += A[j * samples + s] * A[j * samples + s];
                }

                S[j] = inv * D + E[j] + gamma;
            }
        }

        const auto min_element = std::min_element(S.begin(), S.begin() + k + 1);
        const auto min_position = static_cast<int>(std::distance(S.begin(), min_element));
        const auto min_value = *min_element;
        E[k + 1] = min_value;
        T[k] = min_position;
    }

    // Find start positions
    std::vector<int> starts;
    int pos = T.back();
    while (pos > 0) {
        starts.push_back(static_cast<int>(r[pos])); // To be returned to R, so cast to int
        pos = T[pos - 1];
    }
    starts.push_back(0);
    std::reverse(starts.begin(), starts.end());
    return starts;
}

// [[Rcpp::export]]
std::vector<int> exact_pcf_(const std::vector<double> &y, int kmin, double gamma) {
    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    std::size_t N = y.size();

    if (N == 0) {
        return {};
    }

    std::vector<double> A(N, 0);
    std::vector<double> S(N, 0);// Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t j = 0; j <= k; ++j) {
            A[j] += y[k];
            if (j > 0 && (j < kmin_size || k + 1 - j < kmin_size)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                double D = -A[j] * A[j] / (k - j + 1);
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

// [[Rcpp::export]]
std::vector<int> fast_pcf_(const std::vector<double> &y, const std::vector<int> &available_breakpoints, int kmin, double gamma) {
    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    Aggregates agg = make_aggregates(y, available_breakpoints);
    const std::vector<double> &u = agg.aggregates;
    const std::vector<std::size_t> &r = agg.pos;

    std::size_t N = u.size();
    std::vector<double> A(N, 0);
    std::vector<std::size_t> C(N, 0);
    std::vector<double> S(N, 0); // Score
    std::vector<double> E(N + 1, 0);
    std::vector<int> T(N, -1);
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    if (N == 0) {
        return {};
    }

    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t j = 0; j <= k; ++j) {
            A[j] += u[k];
            C[j] += agg.pos_diff[k];
            if (r[j] > 0 && (r[j] < kmin_size || r[k + 1] - r[j] < kmin_size)) {
                S[j] = std::numeric_limits<double>::infinity();
            } else {
                double D = -A[j] * A[j] / C[j];
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
        starts.push_back(static_cast<int>(r[pos])); // To be returned to R, so cast to int
        pos = T[pos - 1];
    }
    starts.push_back(0);
    std::reverse(starts.begin(), starts.end());
    return starts;
}

// [[Rcpp::export]]
std::vector<double> convolve_(const std::vector<double>& x, const std::vector<double>& k) {
    const std::size_t nx = x.size();
    const std::size_t nk = k.size();
    if (nx == 0 || nk == 0) {
        return x;
    }
    std::vector<double> out(nx + nk - 1, 0);
    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < nk; ++j) {
            out[i+j] += x[i] * k[j];
        }
    }
    return out;
}

// [[Rcpp::export]]
double median_(NumericVector x) {
    if (x.size() == 0) {
        return NA_REAL;
    }
    NumericVector y = clone(x);
    R_xlen_t n, half;
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
    if (filter_size < 1) {
        Rcpp::stop("filter_size must be at least 1");
    }
    if (x.size() < static_cast<std::size_t>(filter_size * 6)) {
        Rcpp::stop("Input vector x is too short for the specified filter_size");
    }

    // Make the smoothing sawtooth filter for edge detection
    std::vector<double> k;
    k.reserve(6 * filter_size);
    k.insert(k.end(), filter_size, -1);
    k.insert(k.end(), 2 * filter_size, -2);
    k.insert(k.end(), 2 * filter_size, 2);
    k.insert(k.end(), filter_size, 1);

    std::vector<double> convolved = convolve_(x, k);

    // The edge detection signal is strongest at the sign-change point of the filter; this offset aligns the convolved output with this signal.
    auto offset = filter_size * 3 - 1;
    NumericVector hpf = Rcpp::wrap(std::vector<double>(convolved.begin() + offset, convolved.begin() + offset + x.size()));
    NumericVector abshpf = abs(hpf);
    double threshold = median_(abshpf) + nmad * mad_(hpf);
    std::vector<int> out;
    for (std::size_t i = 0; i < abshpf.size(); ++i) {
        if (abshpf[i] > threshold) {
            out.push_back(static_cast<int>(i));
        }
    }
    return out;
}

#ifdef DEBUG_BUILD
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

template<typename T>
void print_vec(std::string label, const std::vector<T> &v) {
    if (v.empty()) {
        std::cout << label << " []" << std::endl;
        return;
    }
    std::cout << label << " [";
    for (std::size_t i = 0; i < v.size() - 1; ++i) {
        std::cout << v[i] << ", ";
    }
    std::cout << v.back() << "]" << std::endl;
}
#endif
