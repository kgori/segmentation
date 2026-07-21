#include <Rcpp.h>
using namespace Rcpp;

constexpr double EPSILON = 1e-10;  // Small tolerance for floating-point comparisons

// Utility functions
std::vector<double> make_prefix_sums(const std::vector<double>& x) {
    const std::size_t N = x.size();
    std::vector<double> cumulative_sum(N + 1, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        cumulative_sum[i + 1] = cumulative_sum[i] + x[i];
    }
    return cumulative_sum;
}

std::vector<double> make_prefix_sums_of_squares(const std::vector<double> &x) {
    const std::size_t N = x.size();
    std::vector<double> cumulative_sum(N + 1, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        cumulative_sum[i + 1] = cumulative_sum[i] + x[i] * x[i];
    }
    return cumulative_sum;
}

// [[Rcpp::export]]
std::vector<int> pelt_pcf_(const std::vector<double> &y, int kmin, double gamma,
                           const std::vector<int>& allowed_breakpoints, bool constrain_integer = false) {
    if (kmin < 1) {
        Rcpp::stop("kmin must be at least 1");
    }

    std::size_t N = y.size();
    std::size_t kmin_size = static_cast<std::size_t>(kmin);

    if (N == 0) {
        return {};
    }

    for (const auto& bp : allowed_breakpoints) {
        if (bp < 0 || static_cast<std::size_t>(bp) >= y.size()) {
            Rcpp::stop("Allowed breakpoints must be within the range of the data");
        }
    }

    // Difference from exact: A is now a prefix sum vector so we can compute the sum for any size segment in O(1) time.
    std::vector<double> A = make_prefix_sums(y);
    std::vector<double> A2 = make_prefix_sums_of_squares(y);
    // Difference from exact: no need for score vector S
    // Difference from exact: E is now initialized to infinity so that we can prune the search space.
    std::vector<double> E(N + 1, std::numeric_limits<double>::infinity());
    std::vector<int> T(N, -1);

    // Construct the allowed set of breakpoints
    bool restrict_breakpoints = !allowed_breakpoints.empty();
    std::set<std::size_t> allowed_set(allowed_breakpoints.begin(), allowed_breakpoints.end());
    // Always allow the start and end of the data as a breakpoint
    allowed_set.insert(0);
    allowed_set.insert(N);

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
            if (restrict_breakpoints && !allowed_set.count(j)) {
                continue; // Skip if j is not in the allowed breakpoints
            }

            if (j > end) {
                Rcpp::stop("Candidate j is greater than end k+1, which should not happen");
            }

            if (end - j < kmin_size) {
                continue; // Skip if the segment is too short
            }

            const double sum = A[end] - A[j];
            std::size_t len = end - j;
            double D;

            if (constrain_integer) {
                // If we want to constrain the mean to be an integer, we can round the mean and compute the cost accordingly
                const double sum2 = A2[end] - A2[j];
                double mean = sum / static_cast<double>(len);
                double rounded_mean = std::round(mean);
                D = sum2 - 2 * rounded_mean * sum + len * rounded_mean * rounded_mean;
            } else {
                D = -sum * sum / static_cast<double>(end - j);
            }
            
            const double score = D + E[j] + gamma;

            if (score < min_value) {
                min_value = score;
                min_position = static_cast<int>(j);
            }
        }

        E[end] = min_value;
        T[k] = min_position;

        // Difference from exact: Prune the candidate set R based on the pruning
        // condition: only keep candidates j for which E[j] + D <= E[end]
        std::vector<std::size_t> R_new;
        R_new.reserve(R.size() + 1);

        for (const std::size_t j : R) {
            if (restrict_breakpoints && !allowed_set.count(j)) {
                continue; // Skip if j is not in the allowed breakpoints
            }

            if (end - j < kmin_size) {
                // Segment is too short, but we still need to keep it in R_new for future iterations
                R_new.push_back(j);
                continue;
            }

            // Fix the pruning reference point to be kmin size away from the end ("delayed pruning")
            const std::size_t ref = end - kmin_size + 1;
            const std::size_t len = ref - j;

            const double sum = A[ref] - A[j];
            double D;

            if (constrain_integer) {
                const double sum2 = A2[ref] - A2[j];
                double mean = sum / static_cast<double>(len);
                double rounded_mean = std::round(mean);
                D = sum2 - 2 * rounded_mean * sum + len * rounded_mean * rounded_mean;
            } else {
                D = -sum * sum / static_cast<double>(len);
            }

            if (E[j] + D <= E[ref] + EPSILON) {
                R_new.push_back(j);
            }
        }

        if (end >= kmin_size && end < N) {
            if (!restrict_breakpoints || allowed_set.count(end)) {
                R_new.push_back(end); // Add the current end as a new candidate
            }
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
std::vector<int> pelt_multipcf_(const NumericMatrix &y, int kmin, double gamma,
                                const std::vector<int>& allowed_breakpoints, bool constrain_integer = false) {
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
    // A2 is a matrix of per-sample prefix sums of squares (stored contiguously as a vector)
    std::vector<double> A((N + 1) * samples, 0.0);
    std::vector<double> A2((N + 1) * samples, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        auto row = y.row(i);
        for (std::size_t s = 0; s < samples; ++s) {
            const double value = row[s];
            A[(i + 1) * samples + s] = A[i * samples + s] + value;
            A2[(i + 1) * samples + s] = A2[i * samples + s] + value * value;
        }
    }
    std::vector<double> E(N + 1, std::numeric_limits<double>::infinity());
    std::vector<int> T(N, -1);

    // Construct the allowed set of breakpoints
    bool restrict_breakpoints = !allowed_breakpoints.empty();
    std::set<std::size_t> allowed_set(allowed_breakpoints.begin(), allowed_breakpoints.end());
    // Always allow the start and end of the data as a breakpoint
    allowed_set.insert(0);
    allowed_set.insert(N);

    E[0] = 0.0;

    std::vector<std::size_t> R{0};

    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t end = k + 1;

        double min_value = std::numeric_limits<double>::infinity();
        int min_position = -1;

        for (const std::size_t j : R) {
            if (restrict_breakpoints && !allowed_set.count(j)) {
                continue; // Skip if j is not in the allowed breakpoints
            }

            if (j > end) {
                Rcpp::stop("Candidate j is greater than end k+1, which should not happen");
            }
          
            if (end - j < kmin_size) {
                continue; // Skip if the segment is too short
            }

            double cost = 0.0;
            double score;

            const std::size_t len = end - j;

            if (constrain_integer) {
                for (std::size_t s = 0; s < samples; ++s) {
                    const double sum = A[end * samples + s] - A[j * samples + s];
                    const double sum2 = A2[end * samples + s] - A2[j * samples + s];
                    double mean = sum / static_cast<double>(len);
                    double rounded_mean = std::round(mean);
                    cost += sum2 - 2 * rounded_mean * sum + len * rounded_mean * rounded_mean;
                }
                score = cost + E[j] + gamma;
            } else {
                for (std::size_t s = 0; s < samples; ++s) {
                    const double sum = A[end * samples + s] - A[j * samples + s];
                    cost += sum * sum;
                }
                score = -cost / static_cast<double>(len) + E[j] + gamma;
            }
            
            if (score < min_value) {
                min_value = score;
                min_position = static_cast<int>(j);
            }
        }

        E[end] = min_value;
        T[k] = min_position;

        std::vector<std::size_t> R_new;
        R_new.reserve(R.size() + 1);

        for (const std::size_t j : R) {
            if (restrict_breakpoints && !allowed_set.count(j)) {
                continue; // Skip if j is not in the allowed breakpoints
            }
          
            if (end - j < kmin_size) {
                R_new.push_back(j); // Keep it for future iterations
                continue;
            }

            // Fix the pruning reference point to be kmin size away from the end ("delayed pruning")
            const std::size_t ref = end - kmin_size + 1;
            const std::size_t len = ref - j;

            double cost = 0.0;
            double D;
            if (constrain_integer) {
                for (std::size_t s = 0; s < samples; ++s) {
                    const double sum = A[ref * samples + s] - A[j * samples + s];
                    const double sum2 = A2[ref * samples + s] - A2[j * samples + s];
                    double mean = sum / static_cast<double>(len);
                    double rounded_mean = std::round(mean);
                    cost += sum2 - 2 * rounded_mean * sum + len * rounded_mean * rounded_mean;
                }
                D = cost;
            } else {
                for (std::size_t s = 0; s < samples; ++s) {
                    const double sum = A[ref * samples + s] - A[j * samples + s];
                    cost += sum * sum;
                }
                D = -cost / static_cast<double>(len);
            }

            if (E[j] + D <= E[ref] + EPSILON) {
                R_new.push_back(j);
            }
        }

        if (end >= kmin_size && end < N) {
            if (!restrict_breakpoints || allowed_set.count(end)) {
                R_new.push_back(end); // Add the current end as a new candidate
            }
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


