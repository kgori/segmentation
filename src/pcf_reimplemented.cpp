#include <Rcpp.h>
using namespace Rcpp;

// Utility functions
std::vector<double> make_prefix_sums(const std::vector<double>& x) {
    const std::size_t N = x.size();
    std::vector<double> cumulative_sum(N + 1, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        cumulative_sum[i + 1] = cumulative_sum[i] + x[i];
    }
    return cumulative_sum;
}

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
std::vector<int> pelt_multipcf_(const NumericMatrix &y, int kmin, double gamma) {
    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    const std::size_t N = y.nrow();
    const std::size_t samples = y.ncol();
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    if (N == 0 || samples == 0) {
        return {};
    }

    // A is a matrix of per-sample prefix sums (stored contiguously as a vector)
    std::vector<double> A((N + 1) * samples, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t s = 0; s < samples; ++s) {
            A[(i + 1) * samples + s] = A[i * samples + s] + y(i, s);
        }
    }
    std::vector<double> E(N + 1, std::numeric_limits<double>::infinity());
    std::vector<int> T(N, -1);

    E[0] = 0.0;

    std::vector<std::size_t> R{0};
    std::vector<std::size_t> R_new;

    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t end = k + 1;

        double min_value = std::numeric_limits<double>::infinity();
        int min_position = -1;

        for (const std::size_t j : R) {
            if (j > 0 && end - j < kmin_size) {
                continue;
            }

            double sum_of_squares = 0.0;

            for (std::size_t s = 0; s < samples; ++s) {
                const double sum = A[end * samples + s] - A[j * samples + s];
                sum_of_squares += sum * sum;
            }

            const double D = -sum_of_squares / static_cast<double>(end - j);
            const double score = D + E[j] + gamma;

            if (score < min_value) {
                min_value = score;
                min_position = static_cast<int>(j);
            }
        }

        E[end] = min_value;
        T[k] = min_position;
        const double EPSILON = 1e-10;  // Small tolerance for floating-point comparisons

        // Prune the candidate set with kmin adjustment "delayed pruning"
        R_new.clear();
        R_new.reserve(R.size() + 1);

        for (const std::size_t j : R) {
            if (end - j < kmin_size) {
                R_new.push_back(j); // Keep it for future iterations
                continue;
            }

            // Fix the pruning reference point to be kmin size away from the end
            const std::size_t ref = end - kmin_size + 1;

            double sum_of_squares = 0.0;
            for (std::size_t s = 0; s < samples; ++s) {
                const double sum = A[ref * samples + s] - A[j * samples + s];
                sum_of_squares += sum * sum;
            }

            const double D = -sum_of_squares / static_cast<double>(ref - j);

            if (E[j] + D <= E[ref] + EPSILON) {
                R_new.push_back(j);
            }
        }

        if (end >= kmin_size && end < N) {
            R_new.push_back(end);
        }

        R.swap(R_new);
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
std::vector<int> pelt_pcf_(const std::vector<double> &y, int kmin, double gamma) {
    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    std::size_t N = y.size();
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    if (N == 0) {
        return {};
    }

    // Difference from exact: A is now a prefix sum vector so we can compute the sum for any size segment in O(1) time.
    std::vector<double> A = make_prefix_sums(y);
    // Difference from exact: no need for score vector S
    // Difference from exact: E is now initialized to infinity so that we can prune the search space.
    std::vector<double> E(N + 1, std::numeric_limits<double>::infinity());
    std::vector<int> T(N, -1);

    // Initialise E[0] to 0
    E[0] = 0.0;

    // Difference from exact: R is a vector of candidate change points, initialized with 0 (the start of the series).
    std::vector<std::size_t> R{0};

    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t end = k + 1;
        double min_value = std::numeric_limits<double>::infinity();
        int min_position = -1;

        // Difference from exact: Iterate only over candidates, not all j
        for (const std::size_t j : R) {
            if (j > 0 && end - j < kmin_size) {
                continue; // Skip if the segment is too short
            }

            const double sum = A[end] - A[j];
            const double D = -sum * sum / static_cast<double>(end - j);
            const double score = D + E[j] + gamma;

            if (score < min_value) {
                min_value = score;
                min_position = static_cast<int>(j);
            }
        }

        E[end] = min_value;
        T[k] = min_position;
        const double EPSILON = 1e-10;  // Small tolerance for floating-point comparisons

        // Difference from exact: Prune the candidate set R based on the pruning
        // condition: only keep candidates j for which E[j] + D <= E[end]
        std::vector<std::size_t> R_new;
        R_new.reserve(R.size() + 1);

        for (const std::size_t j : R) {
            if (end - j < kmin_size) {
                // Segment is too short, but we still need to keep it in R_new for future iterations
                R_new.push_back(j);
                continue;
            }

            // Fix the pruning reference point to be kmin size away from the end "delayed pruning"
            const std::size_t ref = end - kmin_size + 1;

            const double sum = A[ref] - A[j];
            const double D = -sum * sum / static_cast<double>(ref - j);

            if (E[j] + D <= E[ref] + EPSILON) {
                R_new.push_back(j);
            }
        }

        if (end >= kmin_size && end < N) {
            R_new.push_back(end); // Add the current end as a new candidate
        }

        R.swap(R_new); // Update R to the new candidate set
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


/*
 * fast_* functions need a set of candidate breakpoints. The following mark_ and mark_multi_
 * functions provide a method to generate a breakpoint set using a sawtooth kernel and local
 * thresholding.
 */


// Takes v by value deliberately (nth_element modifies in place).
double quantile_p(std::vector<double> v, double p) {
    if (v.empty()) {
        return 0.0;
    }
    std::size_t idx = static_cast<std::size_t>(p * static_cast<double>(v.size() - 1));
    if (idx >= v.size()) {
        idx = v.size() - 1;
    }
    std::nth_element(v.begin(), v.begin() + idx, v.end());
    return v[idx];
};

// [[Rcpp::export]]
std::vector<double> sliding_max_7(const std::vector<double>& v) {
    const std::size_t n = v.size();
    std::vector<double> out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t lo = (i >= 3) ? i - 3 : 0;
        const std::size_t hi = std::min(n, i + 4);
        double mx = 0.0;
        for (std::size_t w = lo; w < hi; ++w)
            if (v[w] > mx) {
                mx = v[w];
            }
        out[i] = mx;
    }
    return out;
};

// [[Rcpp::export]]
std::vector<double> make_cost_vector_(const std::vector<double>& data, std::size_t kernel_size) {
    std::size_t N = data.size();
    std::size_t size_check = static_cast<std::size_t>(6 * kernel_size);
    if (kernel_size < 1) {
        Rcpp::stop("kernel_size must be >= 1");
    }
    if (N < size_check) {
        Rcpp::stop("Input too short for filter size (need >= 6*L = %d points)", 6 * kernel_size);
    }

    // The cost function is faster to calculate if we collect prefix sums up front
    std::vector<double> prefix_sums = make_prefix_sums(data);

    std::vector<double> cost(N, 0.0);
    const std::size_t n_valid1 = N - size_check + 1;
    for (std::size_t i = 0; i < n_valid1; ++i) {
        cost[i + 3 * kernel_size - 1] = std::abs(
            4.0 * prefix_sums[i + 3 * kernel_size]
            - prefix_sums[i]
            - prefix_sums[i + kernel_size]
            - prefix_sums[i + 5 * kernel_size]
            - prefix_sums[i + 6 * kernel_size]);
    }
    return cost;
}

// [[Rcpp::export]]
std::vector<double> make_cost_vector_padded_(const std::vector<double>& data, std::size_t kernel_size) {
    std::size_t N = data.size();
    std::size_t size_check = static_cast<std::size_t>(6 * kernel_size);
    if (kernel_size < 1) {
        Rcpp::stop("kernel_size must be >= 1");
    }
    if (N < size_check) {
        Rcpp::stop("Input too short for filter size (need >= 6*L = %d points)", 6 * kernel_size);
    }

    // The cost function is faster to calculate if we collect prefix sums up front
    std::vector<double> prefix_sums = make_prefix_sums(data);

    std::vector<double> cost(N, 0.0);

    // lambda to get the prefix sum with zero padding for out-of-bounds indices
    auto get_prefix_sum = [&](int idx) -> double {
        if (idx < 0) {
            return 0.0; // padding with zeros
        } else if (static_cast<std::size_t>(idx) >= prefix_sums.size()) {
            return prefix_sums.back(); // return the last valid prefix sum
        } else {
            return prefix_sums[idx];
        }
    };

    // Simulate left and right padding - when the lookup index is out of bounds, it will return 0.0
    for (std::size_t j = 0; j < N; ++j) {
        // Transformation to index into the prefix sum array, accounting for the kernel size and padding
        int i = static_cast<int>(j) - (3 * kernel_size - 1);
        cost[j] = std::abs(
            4.0 * get_prefix_sum(i + 3 * kernel_size)
            - get_prefix_sum(i)
            - get_prefix_sum(i + kernel_size)
            - get_prefix_sum(i + 5 * kernel_size)
            - get_prefix_sum(i + 6 * kernel_size));
    }

    return cost;
}


/*
    * Marks likely breakpoints in a signal using a sawtooth kernel and local thresholding.
    *
    * @param x The input signal as a vector of doubles.
    * @param frac The fraction of the signal to mark (default 0.12).
    * @param kernel_size The size of the sawtooth kernel (default 8). The kernel structure is [-1, -2, -2, +2, +2, +1]
    * with each element repeated 'kernel_size' times, so the actual kernel length is 6 * kernel_size.
    * @param thres The threshold for marking peaks relative to local maxima (default 0.9).
    * @return A vector of *ZERO-BASED* indices in x that are marked as likely breakpoints.
    */
// [[Rcpp::export]]
std::vector<int> mark_(const std::vector<double>& x,
                        double frac = 0.12,
                        int kernel_size = 8,
                        double thres = 0.9) {
    const std::size_t N = x.size();

    // Extra careful size validation
    const std::size_t size_check = static_cast<std::size_t>(6 * kernel_size);
    if (kernel_size < 1) Rcpp::stop("kernel_size must be >= 1");
    if (N < size_check) Rcpp::stop("Input too short for filter size (need >= 6*kernel_size = %d points)", 6 * kernel_size);

    // Convolution of data with a sawtooth kernel structured as [-1l, -2l, -2l, 2l, 2l, 1l],
    // where l is kernel size, and '-1l' means repeat -1 l times.
    std::vector<double> cost = make_cost_vector_(x, kernel_size);

    // Make a local upper threshold vector using sliding window local max
    const auto local_max = sliding_max_7(cost);

    // Select elements of the cost vector when they are within 'thres' of the local max
    std::vector<double> peaks;
    for (std::size_t i = 0; i < N; ++i) {
        if (cost[i] > 0.0 && cost[i] >= thres * local_max[i]) {
            peaks.push_back(cost[i]);
        }
    }

    // If in the unexpected case no peaks are found, return an empty vector
    if (peaks.empty()) {
        return {};
    }

    // Find the quantile of the selected peaks that will result in marking approximately 'frac' of the total signal
    const double adjusted_frac = std::min(1 - frac, frac * static_cast<double>(N) / static_cast<double>(peaks.size()));
    const double limit = quantile_p(peaks, 1.0 - adjusted_frac);

    // Mark indices where the cost exceeds both the limit and the thresholded local max
    std::vector<int> marked_indices;
    for (std::size_t i = 0; i < N; ++i) {
        if (cost[i] > limit && cost[i] > thres * local_max[i]) {
            marked_indices.push_back(static_cast<int>(i));
        }
    }
    return marked_indices;
}

/*
    * Marks likely breakpoints jointly in multiple samples using a sawtooth kernel and local thresholding.
    *
    * @param x The input signal as a matrix of doubles, samples in columns.
    * @param frac The fraction of the signal to mark (default 0.12).
    * @param kernel_size The size of the sawtooth kernel (default 8). The kernel structure is [-1, -2, -2, +2, +2, +1]
    * with each element repeated 'kernel_size' times, so the actual kernel length is 6 * kernel_size.
    * @param thres The threshold for marking peaks relative to local maxima (default 0.9).
    * @return A vector of *ZERO-BASED* indices in x that are marked as likely breakpoints.
    */
// [[Rcpp::export]]
std::vector<int> mark_multi_(const NumericMatrix& x,
                        double frac = 0.12,
                        int kernel_size = 8,
                        double thres = 0.9) {
    const std::size_t N = x.nrow();
    const std::size_t S = x.ncol();
    if (N == 0 || S == 0) {
        return {};
    }

    // Extra careful size validation
    const std::size_t size_check = static_cast<std::size_t>(6 * kernel_size);
    if (kernel_size < 1) Rcpp::stop("kernel_size must be >= 1");
    if (N < size_check) Rcpp::stop("Input too short for filter size (need >= 6*kernel_size = %d points)", 6 * kernel_size);

    // Convolution of each sample with a sawtooth kernel structured as [-1l, -2l, -2l, 2l, 2l, 1l],
    // where l is kernel size, and '-1l' means repeat -1 l times. Take the max across all samples.
    std::vector<double> cost_joint(N, 0.0);
    std::vector<double> samplebuf(N, 0.0);
    for (std::size_t s = 0; s < S; ++s) {
        const auto sample = x.column(s);
        samplebuf.assign(sample.begin(), sample.end());
        std::vector<double> cost_sample = make_cost_vector_(samplebuf, kernel_size);
        for (std::size_t i = 0; i < N; ++i) {
            if (cost_joint[i] < cost_sample[i]) {
                cost_joint[i] = cost_sample[i];
            }
        }
    }

    // Make a local upper threshold vector using sliding window local max
    const auto local_max = sliding_max_7(cost_joint);

    // Select elements of the cost vector when they are within 'thres' of the local max
    std::vector<double> peaks;
    for (std::size_t i = 0; i < N; ++i) {
        if (cost_joint[i] > 0.0 && cost_joint[i] >= thres * local_max[i]) {
            peaks.push_back(cost_joint[i]);
        }
    }

    // If in the unexpected case no peaks are found, return an empty vector
    if (peaks.empty()) {
        return {};
    }

    // Find the quantile of the selected peaks that will result in marking approximately 'frac' of the total signal
    const double adjusted_frac = std::min(1 - frac, frac * static_cast<double>(N) / static_cast<double>(peaks.size()));
    const double limit = quantile_p(peaks, 1.0 - adjusted_frac);

    // Mark indices where the cost exceeds both the limit and the thresholded local max
    std::vector<int> marked_indices;
    for (std::size_t i = 0; i < N; ++i) {
        if (cost_joint[i] > limit && cost_joint[i] > thres * local_max[i]) {
            marked_indices.push_back(static_cast<int>(i));
        }
    }
    return marked_indices;
}

/*
 * The median absolute deviation (MAD) of a numeric vector.
 *
 * @param x A numeric vector.
 * @param scale_factor A scaling factor to adjust the MAD (default is 1.4826).
 * @return The scaled MAD of the input vector.
 */
// [[Rcpp::export]]
double mad_(NumericVector x, double scale_factor = 1.4826) {
    // scale_factor = 1.4826; default for normal distribution consistent with R
    return median_(abs(x - median_(x))) * scale_factor;
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
