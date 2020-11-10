#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/*
 * Gives the variance of the data in vector v
 */
template <typename T>
double var(const std::vector<T>& v) {
    double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
    double m =  sum / v.size();
    
    double accum = 0.0;
    std::for_each (std::begin(v), std::end(v), [&](const double d) {
        accum += (d - m) * (d - m);
    });
    
    return accum / (v.size() - 1);
}

/*
 * Cross product of x and y (NumericVector)
 */
// [[Rcpp::export]]
NumericVector::elem_type Crossprod(const NumericVector& x, const NumericVector& y) {
    R_xlen_t xlen = x.size();
    if (xlen != y.size()) {
        stop ("Vectors are not the same size (x=%d, y=%d)", xlen, y.size());
    }
    NumericVector::elem_type out = 0;
    
    for (R_xlen_t i = 0; i < xlen; ++i) {
        out += x[i]*y[i];
    }
    return out;
}

/*
 * Cross product of x and y (std::vector<double>)
 */
double Crossprod(const std::vector<double>& x, const std::vector<double>& y) {
    size_t xlen = x.size();
    if (xlen != y.size()) {
        stop ("Vectors are not the same size");
    }
    double out = 0;
    
    for (size_t i = 0; i < xlen; ++i) {
        out += x[i]*y[i];
    }
    
    return out;
}

/*
 * Cumulative sum of vector v
 */
std::vector<double> cumsum(const std::vector<double>& v) {
    std::vector<double> out;
    double acc = 0;
    std::transform(v.begin(), v.end(), std::back_inserter(out),
                   [&acc](double x) {
                       acc += x;
                       return acc;
                   });
    return out;
}

/*
 * Cumulative sum of vector v, result stored in vector out
 */
void cumsum(const NumericVector& v, std::vector<double>& out) {
    double acc = 0;
    std::transform(v.begin(), v.end(), std::back_inserter(out),
                   [&acc](double x) {
                       acc += x;
                       return acc;
                   });
}

/*
 * Returns `p`-percentage quantile of data in x
 */
// [[Rcpp::export]]
double Quantile(const std::vector<double>& x, double p) {
    if (p < 0) {
        return x.front();
    }
    if (p > 1) {
        return x.back();
    }
    
    std::vector<double> y(x);
    std::sort(y.begin(), y.end());
    
    double index = (y.size()-1.0) * p;
    size_t lo = std::floor(index);
    size_t hi = std::ceil(index);
    
    double qs = y[lo];
    double x_hi = y[hi];
    
    if ((index > lo) && (x_hi != qs)) {
        double h = index - lo;
        qs = (1.0 - h) * qs + h * x_hi;
    }
    return qs;
}

/*
 Finds maximum element within a sliding window in input, for windows of
 size 'windowSize'. Result goes in 'output' vector (modified in-place).
 */
void slidingWindowMax(size_t windowSize, const std::vector<double>& input, std::vector<double>& output) {
    for (int l = 0; l <= input.size() - windowSize; ++l) {
        int r = l + windowSize;
        double windowMax = *std::max_element(input.begin() + l, input.begin() + r);
        output.push_back(windowMax);
    }
}

template <typename T>
std::vector<T> Diff(const std::vector<T>& input) {
    std::vector<T> output;
    if (input.size()==0) {
        return output;
    }
    output.push_back(input[0]);
    for (int i=1; i < input.size(); ++i) {
        output.push_back(input[i] - input[i-1]);
    }
    return output;
}

/*
 * Compliance with Artistic License 2.0 - modifications
 * All functions below are derived from R code in the R package "copynumber",
 * but are now written in C++.
 */

/*
 * Sawtooth filter for fast approximate identification of breakpoints
 */
// [[Rcpp::export]]
std::vector<bool> sawMarkMCpp(const NumericMatrix& x, int L, double frac1, double frac2) {
    auto nrProbes = x.nrow();
    auto nrSample = x.ncol();
    
    std::vector<bool> mark(nrProbes, false);
    std::vector<double> sawValue(nrProbes, 0);
    NumericVector filter(2 * L);
    std::vector<double> sawValue2(nrProbes, 0);
    NumericVector filter2(6);
    
    for (int k = 0; k < L; ++k) {
        double d = k;
        filter[k] = (d+1) / L;
        filter[2 * L - k - 1] = -(d+1) / L;
    }
    
    for (int k = 0; k < 3; ++k) {
        double d = k;
        filter2[k] = (d + 1) / 3;
        filter2[6 - k - 1] = -(d + 1) / 3;
    }
    
    for (int l = 0; l < (nrProbes - 2 * L + 1); ++l) {
        for (int m = 0; m < nrSample; ++m) {
            double diff = 0.0;
            for (int i = 0; i < 2*L; ++i) {
                diff += filter[i] * x(i+l, m);
            }
            sawValue[l + L - 1] = sawValue[l + L - 1] + abs(diff);
        }
    }
    
    double limit = Quantile(sawValue, (1 - frac1));
    for (int l = 0; l < (nrProbes - 2 * L); ++l) {
        if (sawValue[l + L - 1] > limit) {
            mark[l + L - 1] = 1;
        }
    }
    
    for (int l = (L - 2); l < (nrProbes - L - 2); ++l) {
        for (int m = 0; m < nrSample; ++m) {
            double diff2 = 0.0;
            for (int i = 0; i < 6; ++ i) {
                diff2 += filter2[i] * x(i+l, m);
            }
            sawValue2[l + 2] = sawValue2[l + 2] + abs(diff2);
        }
    }
    
    double limit2 = Quantile(sawValue2, (1 - frac2));
    for (int l = (L - 2); l < (nrProbes - L - 2); ++l) {
        if (sawValue2[l + 2] > limit2) {
            mark[l + 2] = 1;
        }
    }
    
    for (int l = 0; l < L; ++l) {
        mark[l] = 1;
        mark[nrProbes - 1 - l] = 1;
    }
    return(mark);
}

// [[Rcpp::export]]
List compactMultiCpp(const NumericMatrix& y, const std::vector<bool>& mark) {
    size_t antGen = y.ncol();
    size_t antSample = y.nrow();
    int antMark = 0;
    for (bool val : mark) {
        antMark += val ? 1 : 0;
    }
    
    std::vector<int> ant(antMark, 0);
    NumericMatrix sum(antSample, antMark);
    
    int pos = 1;
    int oldPos = 0;
    int count = 1;
    
    while (pos <= antGen) {
        std::vector<double> delSum(antSample, 0);
        while (!mark[pos-1]) {
            for (int i = 0; i < antSample; ++i) {
                delSum[i] += y(i, pos-1);
            }
            pos++;
        }
        ant[count-1] = pos - oldPos;
        
        for (int i = 0; i < antSample; ++i) {
            sum(i, count-1) = delSum[i] + y(i, pos-1);
        }
        oldPos = pos;
        pos++;
        count++;
    }
    
    List L = List::create(Named("Nr") = ant,
                          Named("Sum") = sum);
    return L;
}

std::vector<bool> findMarksMultiCpp(const std::vector<bool>& markSub, const std::vector<int>& Nr, int subSize) {
    std::vector<bool> mark(subSize, false);
    bool any = false;
    for (bool val : markSub) {
        if (val) {
            any = true;
            break;
        }
    }
    if (!any) {
        return mark;
    }
    size_t N = markSub.size();
    std::vector<int> help;
    for (int i = 0; i < N; ++i) {
        if (markSub[i]) {
            help.push_back(i+1);
        }
    }
    
    size_t lengdeHelp = help.size();
    auto lengde = Diff(help);
    
    // _1b = reminder that these are 1-based indices ported from R
    int start0_1b = 1;
    int oldStart_1b = 1;
    int startOrig_1b = 1;
    
    for (int i=0; i < lengdeHelp; ++i) {
        start0_1b += lengde[i];
        double lengdeOrig = std::accumulate(Nr.begin() + oldStart_1b - 1,
                                            Nr.begin() + start0_1b - 1,
                                            0);
        startOrig_1b += lengdeOrig;
        mark[startOrig_1b - 2] = true;
        oldStart_1b = start0_1b;
    }
    
    return mark;
}

// [[Rcpp::export]]
std::vector<bool> markMultiPottsCpp(const std::vector<int>& nr, const NumericMatrix& sum, double gamma, int subSize) {
    size_t N = nr.size();
    std::vector<bool> markSub(N, false);
    size_t nSamples = sum.nrow();
    std::vector<double> bestCost(N, 0);
    std::vector<int> bestSplit(N+1, 0);
    NumericMatrix bestAver(nSamples, N);
    NumericMatrix Sum(nSamples, N);
    NumericMatrix Nevner(nSamples, N);
    NumericMatrix eachCost(nSamples, N);
    std::vector<double> Cost(N, 0);
    
    Sum.column(0) = sum.column(0);
    std::fill(Nevner.column(0).begin(), Nevner.column(0).end(), nr[0]);
    bestAver.column(0) = sum.column(0) / nr[0];
    bestCost[0] = Crossprod(-Sum.column(0), bestAver.column(0));
    
    std::vector<size_t> lengde(N, 0);
    
    for (int n = 1; n < N; ++n) {
        // This loop replaces several vectorised operations done in the original R code
        for (int col = 0; col <= n; ++col) {
            Sum.column(col) = Sum.column(col) + sum.column(n);
            Nevner.column(col) = Nevner.column(col) + nr[n];
            eachCost.column(col) = -(Sum.column(col)*Sum.column(col))/Nevner.column(col);
            Cost[col] = std::accumulate(eachCost.column(col).begin(),
                                        eachCost.column(col).end(), 0.0);
            if (col > 0) {
                Cost[col] += bestCost[col-1] + gamma;
            }
        }
        
        int Pos = std::distance(Cost.begin(),
                                std::min_element(Cost.begin(),
                                                 Cost.begin() + n + 1));
        bestCost[n] = Cost[Pos];
        bestAver.column(n) = Sum.column(Pos) / Nevner.column(Pos);
        bestSplit[n] = Pos; // Check 0/1-based indexing
        if (Pos > 0) {
            markSub[Pos-1] = true;
        }
    }
    std::vector<bool> help = findMarksMultiCpp(markSub, nr, subSize);
    return help;
}

// [[Rcpp::export]]
List multiPCFcompactCpp(const std::vector<int>& nr, const NumericMatrix& sum, double gamma) {
    size_t N = nr.size();
    size_t nSamples = sum.nrow();
    NumericMatrix yhat(nSamples, N);
    std::vector<double> bestCost(N, 0);
    std::vector<int> bestSplit(N+1, 0);
    NumericMatrix bestAver(nSamples, N);
    NumericMatrix Sum(nSamples, N);
    NumericMatrix Nevner(nSamples, N);
    NumericMatrix eachCost(nSamples, N);
    std::vector<double> Cost(N, 0);
    
    Sum.column(0) = sum.column(0);
    std::fill(Nevner.column(0).begin(), Nevner.column(0).end(), nr[0]);
    bestAver.column(0) = sum.column(0) / nr[0];
    bestCost[0] = Crossprod(-Sum.column(0), bestAver.column(0));
    std::vector<int> lengde(N, 0);
    
    for (int n = 1; n < N; ++n) {
        // This loop replaces several vectorised operations done in the original R code
        for (int col = 0; col <= n; ++col) {
            Sum.column(col) = Sum.column(col) + sum.column(n);
            Nevner.column(col) = Nevner.column(col) + nr[n];
            eachCost.column(col) = -(Sum.column(col)*Sum.column(col))/Nevner.column(col);
            Cost[col] = std::accumulate(eachCost.column(col).begin(),
                                        eachCost.column(col).end(), 0.0);
            if (col > 0) {
                Cost[col] += bestCost[col-1] + gamma;
            }
        }
        int Pos = std::distance(Cost.begin(),
                                std::min_element(Cost.begin(),
                                                 Cost.begin() + n + 1));
        
        bestCost[n] = Cost[Pos];
        bestAver.column(n) = Sum.column(Pos) / Nevner.column(Pos);
        bestSplit[n] = Pos; // Check 0/1-based indexing
    }
    
    // n is a 1-based index, used to look up splits in bestSplits,
    // which is a collection of 1-based indices
    int n = N;
    int antInt = 0;
    
    while(n > 0) {
        int min = n < bestSplit[n-1] ? n : bestSplit[n-1];
        int max = n < bestSplit[n-1] ? bestSplit[n-1] : n;
        int nrsum = 0;
        for (int col = min; col < max; ++col) {
            yhat.column(col) = bestAver.column(n - 1);
            nrsum += nr[col];
        }
        
        antInt++;
        lengde[antInt-1] = nrsum;
        n = bestSplit[n - 1];
    }
    
    std::vector<int> lengdeRev(lengde.rbegin() + lengde.size() - antInt, lengde.rend());
    std::vector<int> init(antInt, 0);
    init[0] = 1;
    if (antInt >= 1) {
        for (int k = 1; k < antInt; ++k) {
            init[k] = init[k-1] + lengdeRev[k-1];
        }
    }
    
    n = N;
    NumericMatrix verdi(nSamples, antInt);
    bestSplit[n] = n;
    int antall = antInt;
    
    while (n > 0) {
        verdi.column(antall-1) = bestAver.column(n-1);
        n = bestSplit[n-1];
        antall--;
    }
    return List::create(Named("Lengde") = lengdeRev, Named("sta") = init,
                        Named("mean") = verdi, Named("nIntervals") = antInt);
}

// [[Rcpp::export]]
NumericMatrix expandMultiCpp(size_t nProbes, size_t nSamples, const std::vector<int>& lengthInt, const NumericMatrix& mean) {
    NumericMatrix Potts(nSamples, nProbes);
    size_t lengthCompArr = lengthInt.size();
    int k = 0;
    
    for (int i = 0; i < lengthCompArr; ++i) {
        for (int j = 0; j < lengthInt[i]; ++j) {
            Potts.column(k) = mean.column(i);
            k++;
        }
    }
    return Potts;
}

/* marks potential breakpoints, partially by a two 6*L and 6*L2 highpass
 filters (L>L2), then by a filter seaching for potential kmin long segments */
// [[Rcpp::export]]
std::vector<bool> filterMarkCpp(const NumericVector& x, int kmin, int L, int L2, double frac1, double frac2, double frac3, double thres) {
    bool fail = any(is_nan(x));
    if (fail) stop("Input must not contain NaN values");
    fail = any(is_na(x));
    if (fail) stop("Input must not contain NA values");
    fail = any(is_infinite(x));
    if (fail) stop("Input must not contain Inf values");
    
    int lengdeArr = x.size();
    std::vector<double> xc{0};
    cumsum(x, xc);
    
    size_t padding = 3*L;
    std::vector<double> cost1(padding-1, 0);
    std::vector<double> zeros(padding, 0);
    
    int ind11, ind12, ind13, ind14, ind15;
    for (int i = 0; i <= lengdeArr - 6 * L; ++i) {
        ind11 = i;
        ind12 = i+L;
        ind13 = i+3*L;
        ind14 = i+5*L;
        ind15 = i+6*L;
        cost1.push_back(fabs(4 * xc[ind13] - xc[ind11] - xc[ind12] - xc[ind14] - xc[ind15]));
    }
    
    cost1.insert(cost1.end(), zeros.begin(), zeros.end());
    
    // Find sliding-window maximum element in cost1, for windows of size 7. Result goes in 'test' vector.
    // test vector has padding at front and back.
    // Add front padding here, back padding after the loop
    std::vector<double> test{0, 0, 0};
    slidingWindowMax(7, cost1, test);
    
    // Add padding to 'test' vector
    for (int i = 0; i < 3; ++i) {
        test.push_back(0);
    }
    
    // Create cost1B, which is cost1[cost1 >= thres * test]
    assert (cost1.size() == test.size());
    std::vector<double> cost1B;
    for (int i = 0; i < cost1.size(); ++i) {
        if (cost1[i] >= thres * test[i]) {
            cost1B.push_back(cost1[i]);
        }
    }
    
    double frac1B = std::min(0.8, frac1 * (double)cost1.size()/(double)cost1B.size());
    auto limit = Quantile(cost1B, (1 - frac1B));
    
    std::vector<bool> mark;
    for (int i=0; i < cost1.size(); ++i) {
        mark.push_back((cost1[i] > limit) & (cost1[i] > 0.9 * test[i]));
    }
    
    std::vector<double> cost2;
    int ind21, ind22, ind23, ind24, ind25;
    for (int i = 0; i <= lengdeArr - 6 * L2; ++i) {
        ind21 = i;
        ind22 = i+L2;
        ind23 = i+3*L2;
        ind24 = i+5*L2;
        ind25 = i+6*L2;
        cost2.push_back(fabs(4 * xc[ind23] - xc[ind21] - xc[ind22] - xc[ind24] - xc[ind25]));
    }
    
    double limit2 = Quantile(cost2, (1 - frac2));
    
    // 'mark2' is padded
    padding = 3*L2;
    std::vector<bool> mark2(padding-1, 0);
    
    for (int i=0; i < cost2.size(); ++i) {
        mark2.push_back((cost2[i] > limit2));
    }
    zeros = std::vector<double>(padding, 0);
    mark2.insert(mark2.end(), zeros.begin(), zeros.end());
    
    if (3 * L > kmin) {
        for (int i = kmin-1; i < 3*L-1; ++i) {
            mark[i] = true;
        }
        for (int j = lengdeArr - 3 * L; j < (lengdeArr - kmin); ++j ) {
            mark[j] = true;
        }
    } else {
        mark[kmin-1] = true;
        mark[lengdeArr - kmin - 1] = true;
    }
    
    if (kmin > 1) {
        int ind1, ind2, ind3, ind4;
        std::vector<double> shortAb;
        // TODO: potential index error here, when lengdeArr < (3*kmin+1).
        //       Crashes R if not detected.
        for (int i = 0; i < (lengdeArr - 3 * kmin + 1); ++i) {
            ind1 = i;
            ind2 = ind1 + 3 * kmin;
            ind3 = ind1 + kmin;
            ind4 = ind1 + 2*kmin;
            shortAb.push_back(fabs(3 * (xc[ind4] - xc[ind3]) - (xc[ind2] - xc[ind1])));
        }
        
        // Repeat of sliding window routine, on 'shortAb' input
        // reset test vector (with padding)
        test = std::vector<double>{0, 0, 0};
        slidingWindowMax(7, shortAb, test);
        // Add padding to 'test' vector
        for (int i = 0; i < 3; ++i) {
            test.push_back(0);
        }
        
        // Create cost1C, which is shortAb[shortAb >= thres * test]
        assert (shortAb.size() == test.size());
        std::vector<double> cost1C;
        for (int i = 0; i < shortAb.size(); ++i) {
            if (shortAb[i] >= thres * test[i]) {
                cost1C.push_back(shortAb[i]);
            }
        }
        
        double frac1C = std::min(0.8, frac3 * (double)shortAb.size()/(double)cost1C.size());
        double limit3 = Quantile(cost1C, (1 - frac1C));
        
        std::vector<bool> markH1, markH2, markH3;
        for (int i=0; i < shortAb.size(); ++i) {
            bool q = (shortAb[i] > limit3) && (shortAb[i] > thres * test[i]);
            markH1.push_back(q);
        }
        
        for (int i=0; i < kmin-1; ++i) markH2.push_back(false);
        std::copy(markH1.begin(), markH1.end(), std::back_inserter(markH2));
        for (int i=0; i < 2*kmin; ++i) markH2.push_back(false);
        
        for (int i=0; i < 2*kmin-1; ++i) markH3.push_back(false);
        std::copy(markH1.begin(), markH1.end(), std::back_inserter(markH3));
        for (int i=0; i < kmin; ++i) markH3.push_back(false);
        
        assert (mark.size() == mark2.size() == markH2.size() == markH3.size());
        for (int i=0; i < mark.size(); ++i) {
            mark[i] = (mark[i] || mark2[i] || markH2[i] || markH3[i]);
        }
    } else {
        assert (mark.size() == mark2.size());
        for (int i=0; i < mark.size(); ++i) {
            mark[i] = (mark[i] || mark2[i]);
        }
    }
    
    for (int i=0; i < kmin-1; ++i) {
        mark[i] = false;
    }
    
    if (3*L > kmin) {
        for (int i=kmin-1; i < 3*L-1; ++i) {
            mark[i] = true;
        }
        for (int i=lengdeArr - 3*L; i < (lengdeArr-kmin); ++i){
            mark[i] = true;
        }
        for (int i=lengdeArr - kmin; i < (lengdeArr-1); ++i){
            mark[i] = false;
        }
        mark[lengdeArr-1] = true;
    } else {
        for (int i=lengdeArr - kmin; i < (lengdeArr-1); ++i){
            mark[i] = false;
        }
        mark[lengdeArr-1] = true;
        mark[kmin-1] = true;
        mark[lengdeArr - kmin - 1] = true;
    }
    return mark;
}

// [[Rcpp::export]]
List compactCpp(const NumericVector& y, const std::vector<bool>& mark) {
    std::vector<int> cCTell;
    for (int i=0; i < y.size(); ++i) {
        if (mark[i]) cCTell.push_back(i+1);
    }
    auto ant = Diff(cCTell);
    
    std::vector<double> cy;
    cumsum(y, cy);
    std::vector<double> cCcy;
    for (int i=0; i < mark.size(); ++i) {
        if (mark[i]) cCcy.push_back(cy[i]);
    }
    
    auto sum = Diff(cCcy);
    List L = List::create(Named("Nr") = ant, Named("Sum") = sum);
    return L;
}

// [[Rcpp::export]]
List findEstCpp(const std::vector<int>& bestSplit, size_t N, const std::vector<int>& Nr, const std::vector<double>& Sum, bool yest) {
    size_t n = N;
    std::vector<int> lengde(N, 0);
    size_t antInt = 0;
    while(n > 0) {
        if (n > bestSplit.size()) {
            break;
        }
        if (antInt < 0 || antInt >= lengde.size()) {
            break;
        }
        lengde[antInt] = n - bestSplit[n-1];
        n = bestSplit[n-1];
        antInt++;
    }
    std::reverse(lengde.begin(), lengde.begin()+antInt);
    lengde.erase(lengde.begin() + antInt, lengde.end());
    
    std::vector<int> lengdeOrig(antInt, 0);
    std::vector<int> startOrig(antInt+1, 1);
    std::vector<double> verdi;
    std::copy(lengdeOrig.begin(), lengdeOrig.end(), std::back_inserter(verdi));
    std::vector<int> start(startOrig);
    
    for (int i = 0; i < antInt; ++i) {
        start[i+1] = start[i] + lengde[i];
        lengdeOrig[i] = std::accumulate(Nr.begin() + start[i] - 1,
                                        Nr.begin() + start[i+1] - 1,
                                        0);
        startOrig[i + 1] = startOrig[i] + lengdeOrig[i];
        verdi[i] = std::accumulate(Sum.begin() + start[i] - 1,
                                   Sum.begin() + start[i+1] - 1,
                                   0.0) / lengdeOrig[i];
    }
    
    if (yest) {
        std::vector<double> yhat(startOrig[antInt] - 1, 0.0);
        for (int i = 0; i < antInt; ++i) {
            for (int j = startOrig[i] - 1; j < startOrig[i+1] - 1; ++j) {
                yhat[j] = verdi[i];
            }
        }
        startOrig.erase(startOrig.begin() + antInt, startOrig.end());
        return List::create(Named("Lengde") = lengdeOrig,
                            Named("sta") = startOrig,
                            Named("mean") = verdi,
                            Named("nIntervals") = antInt,
                            Named("yhat") = yhat);
    } else {
        startOrig.erase(startOrig.begin() + antInt, startOrig.end());
        return List::create(Named("Lengde") = lengdeOrig,
                            Named("sta") = startOrig,
                            Named("mean") = verdi,
                            Named("nIntervals") = antInt);
    }
}

// [[Rcpp::export]]
List PottsCompactCpp(int kmin, double gamma, const std::vector<int>& nr, const std::vector<double>& res, bool yest) {
    size_t N = nr.size();
    std::vector<double> Ant(N, 0);
    std::vector<double> Sum(N, 0);
    std::vector<double> Cost(N, 0);
    std::vector<double> bestCost(N, 0);
    std::vector<int> bestSplit(N, 0);
    
    double nrsum = std::accumulate(nr.begin(), nr.end(), 0.0);
    if (nrsum < 2 * kmin) {
        double ressum = std::accumulate(res.begin(), res.end(), 0.0);
        double estim = ressum / nrsum;
        return List::create(Named("estim") = estim);
    }
    
    int initAnt = nr.front();
    double initSum = res.front();
    double initAve = initSum / initAnt;
    
    bestCost[0] = -initSum * initAve;
    
    int k = 1; // adjusted to 0-based indexing
    double sum_nr_to_k = nr[0] + nr[1];
    
    while (sum_nr_to_k < 2 * kmin) {
        for (int i = 1; i <= k; ++i) {
            Ant[i] += nr[k];
            Sum[i] += res[k];
        }
        bestCost[k] = (-pow(initSum + Sum[1], 2)/(initAnt + Ant[1]));
        k += 1;
        sum_nr_to_k += nr[k];
    }
    
    for (int n = k; n < N; ++n) {
        for (int i = 1; i <= n; ++i) {
            Ant[i] += nr[n];
            Sum[i] += res[n];
        }
        
        // Ported from R, limit is a 1-based index. I'm leaving it as 1-based because it defines
        // an inclusive upper limit in the R code, and I can use it as an exclusive upper limit
        // here and get the same effect.
        int limit = n+1;
        while (limit > 2 && Ant[limit-1] < kmin) {
            limit -= 1;
        }
        
        for (int i = 1; i < limit; ++i) {
            Cost[i] = bestCost[i-1] - pow(Sum[i], 2) / Ant[i];
        }
        int Pos = std::distance(Cost.begin(), std::min_element(Cost.begin() + 1, Cost.begin() + limit));
        double cost = Cost[Pos] + gamma;
        double totCost = (-pow(Sum[1] + initSum,2)/(Ant[1] + initAnt));
        if (totCost < cost) {
            Pos = 0;
            cost = totCost;
        }
        bestCost[n] = cost;
        bestSplit[n] = Pos;
    }
    
    if (yest) {
        std::vector<double> yhat(N, 0);
        return(findEstCpp(bestSplit, N, nr, res, true));
    } else {
        return(findEstCpp(bestSplit, N, nr, res, false));
    }
}

// [[Rcpp::export]]
std::vector<bool> findMarksCpp(const std::vector<bool> markSub, const std::vector<int>& Nr, int subsize) {
    std::vector<bool> mark(subsize, false);
    
    bool any = false;
    for (bool item : markSub) {
        if (item) {
            any = true;
            break;
        }
    }
    if (!any) return mark;
    
    size_t N = markSub.size();
    std::vector<int> help;
    
    for (int i = 0; i < N; ++i) {
        if (markSub[i]) {
            help.push_back(i+1);
        }
    }
    
    size_t lengdehelp = help.size();
    auto lengde = Diff(help);
    
    int start = 0;
    int oldStart = 0;
    int startOrig = 0;
    
    for (int i = 0; i < lengdehelp; ++i) {
        start = start + lengde[i];
        int lengdeOrig = std::accumulate(Nr.begin() + oldStart, Nr.begin() + start, 0);
        startOrig = startOrig + lengdeOrig;
        if (startOrig > 0) {
            mark[startOrig-1] = true;
        }
        oldStart = start;
    }
    return mark;
}

// [[Rcpp::export]]
std::vector<bool> markWithPottsCpp(int kmin, double gamma, const std::vector<int>& nr, const std::vector<double>& res, int subsize) {
    size_t N = nr.size();
    std::vector<double> Ant(N, 0);
    std::vector<double> Sum(N, 0);
    std::vector<double> Cost(N, 0);
    std::vector<double> bestCost(N, 0);
    std::vector<int> bestSplit(N, 0);
    std::vector<bool> markSub(N, false);
    
    int initAnt = nr.front();
    double initSum = res.front();
    double initAve = initSum / (double)initAnt;
    
    
    bestCost[0] = -initSum * initAve;
    
    int k = 1; // adjusted to 0-based indexing
    double sum_nr_to_k = nr[0] + nr[1];
    
    while (sum_nr_to_k < 2 * kmin) {
        for (int i = 1; i <= k; ++i) {
            Ant[i] += nr[k];
            Sum[i] += res[k];
        }
        bestCost[k] = (-pow(initSum + Sum[1], 2)/(initAnt + Ant[1]));
        k += 1;
        sum_nr_to_k += nr[k];
    }
    
    for (int n = k; n < N; ++n) {
        for (int i = 1; i <= n; ++i) {
            Ant[i] += nr[n];
            Sum[i] += res[n];
        }
        
        // Ported from R, limit is a 1-based index, so need to be careful here on the C++ side,
        // where indices are 0-based. (See note in PottsCompact for details)
        int limit = n+1;
        while (limit > 2 && Ant[limit-1] < kmin) {
            limit -= 1;
        }
        
        for (int i = 1; i < limit; ++i) {
            Cost[i] = bestCost[i-1] - pow(Sum[i], 2) / Ant[i];
        }
        int Pos = std::distance(Cost.begin(), std::min_element(Cost.begin() + 1, Cost.begin() + limit));
        double cost = Cost[Pos] + gamma;
        double totCost = (-pow(Sum[1] + initSum,2)/(Ant[1] + initAnt));
        if (totCost < cost) {
            Pos = 0;
            cost = totCost;
        }
        bestCost[n] = cost;
        bestSplit[n] = Pos;
        if (Pos > 0) {
            markSub[Pos-1] = true;
        }
    }
    std::vector<bool> help = findMarksCpp(markSub, nr, subsize);
    return help;
}
