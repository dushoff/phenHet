

#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cmath>

using namespace Rcpp;

// ---------- Fenwick (Binary Indexed) Tree ----------
struct Fenwick {
    int n;
    std::vector<double> bit; // 1-based internal indexing

    Fenwick(int n_=0) : n(n_), bit(n_+1, 0.0) {}

    void reset(int n_) {
        n = n_;
        bit.assign(n+1, 0.0);
    }

    // point update: add delta to index i (0-based external)
    inline void add(int i, double delta) {
        if (delta == 0.0) return;
        for (int idx = i + 1; idx <= n; idx += (idx & -idx)) {
            bit[idx] += delta;
        }
    }

    // sum of [0..i], i is 0-based
    inline double sumPrefix(int i) const {
        double res = 0.0;
        for (int idx = i + 1; idx > 0; idx -= (idx & -idx)) {
            res += bit[idx];
        }
        return res;
    }

    // total sum
    inline double total() const {
        return sumPrefix(n - 1);
    }

    // return smallest i (0-based) such that sumPrefix(i) >= target
    // assumes target in (0, total()]
    inline int findByCumulative(double target) const {
        int idx = 0;
        double acc = 0.0;
        int bitMask = 1;
        while ((bitMask << 1) <= n) bitMask <<= 1;

        for (int step = bitMask; step > 0; step >>= 1) {
            int next = idx + step;
            if (next <= n && acc + bit[next] < target) {
                idx = next;
                acc += bit[next];
            }
        }
        // idx is the largest 1-based index with sum < target.
        // Answer in 1-based is idx+1; zero-based is (idx).
        if (idx >= n) idx = n - 1; // guard
        return idx;
    }
};

// [[Rcpp::export]]
List GilAlgoCpp(    List adjList,
                    int size,
                    double beta,
                    double gamma,
                    double MaxTime,
                    int InitInfSize = 1,
                    bool TrackDyn = true,
                    bool debug = false,
                    int debug_freq = 1,
                    int debug_low = 500,
                    int debug_up = 600) {

    long int debug_ctr = 0;
    long int event_ctr = 0;
    int N = size;
    double t = 0.0;

    // States and per-node rates
    IntegerVector Status(N, 0); // 0=S, 1=I, 2=R
    NumericVector Rate(N, 0.0);
    IntegerVector Deg_vec(N);

    // Pre-compute degrees and probabilities for selecting initial infected
    NumericVector prob(N);
    double deg_sum = 0.0;
    for (int i = 0; i < N; ++i) {
        IntegerVector neighbors = adjList[i];
        int deg = neighbors.size();
        Deg_vec[i] = deg;
        prob[i] = deg;
        deg_sum += deg;
    }
    prob = prob / sum(prob);

    // Initial situation
    IntegerVector noseq = seq(0, N - 1);
    IntegerVector nodes = seq(1, N);
    IntegerVector InitIndex = Rcpp::sample(noseq, InitInfSize, false, prob);

    Rprintf("Initial index %d \n", InitIndex[0]);

    // Outputs
    NumericVector Infect_time(N, NA_REAL);
    NumericVector Recovery_time(N, NA_REAL);
    IntegerVector Infect_num_rnd(N, 0);
    IntegerVector S_NbrDeg(N, 0);
    IntegerVector Infector_rnd(N, NA_INTEGER);

    NumericVector t_vec, S_vec, I_vec, R_vec;

    // Initialize counts
    int S_cnt = N - InitInfSize;
    int I_cnt = InitInfSize;
    int R_cnt = 0;

    // Initialize statuses and gamma recovery rates for initial infected
    for (int i = 0; i < InitInfSize; ++i) {
        int idx = InitIndex[i];
        Status[idx] = 1;
        Rate[idx] = gamma;
        if (TrackDyn) Infect_time[idx] = 0.0;
    }

    // Build Fenwick tree and SumRate after initial neighbor infection rates are applied
    Fenwick fw(N);
    double SumRate = 0.0;

    // Apply neighbor infection rate contributions from initially infected
    for (int i = 0; i < InitInfSize; ++i) {
        int x = InitIndex[i];
        IntegerVector neighbors = adjList[x];
        for (int j : neighbors) {
            int nbr = j - 1; // adjList is 1-based
            if (Status[nbr] == 0) {
                Rate[nbr] += beta;
                ++S_NbrDeg[x];
            }
        }
    }

    // Load Fenwick from Rate and compute SumRate
    for (int i = 0; i < N; ++i) {
        if (Rate[i] != 0.0) {
            fw.add(i, Rate[i]);
            SumRate += Rate[i];
        }
    }

    if (TrackDyn) {
        t_vec.push_back(t);
        S_vec.push_back(static_cast<double>(S_cnt) / N);
        I_vec.push_back(static_cast<double>(I_cnt) / N);
        R_vec.push_back(static_cast<double>(R_cnt) / N);
    }

    int Istep = I_cnt;

    while (t < MaxTime && Istep > 0) {
        if (SumRate <= 0.0) break;

        // Draw two uniforms (avoid temporary NumericVector)
        double r1 = R::runif(0.0, 1.0);
        double r2 = R::runif(0.0, 1.0);

        // Guard against target==0 choosing a zero-rate index
        double target = r1 * SumRate;
        if (target <= 0.0) {
            target = std::numeric_limits<double>::min();
        }

        // Sample event index via Fenwick inverse CDF
        int Event = fw.findByCumulative(target);

        // Update status (0->1 infection, 1->2 recovery)
        Status[Event] += 1;
        event_ctr++;

        if (debug && (debug_ctr % debug_freq == 0) && (debug_ctr > debug_low) && (debug_ctr < debug_up)) {
            Rprintf("%ld %f %f %f %d", event_ctr, t, r1, r2, Event);
        }

        // Neighbors of Event
        IntegerVector neighbors = adjList[Event];

        // Partition neighbors into susceptible contacts and infected infector candidates
        std::vector<int> Contact;
        std::vector<int> Infector; // std::vector to avoid Rcpp overhead
        Contact.reserve(neighbors.size());
        Infector.reserve(neighbors.size());

        for (int j : neighbors) {
            int nbr = j - 1;
            if (Status[nbr] == 0) {
                Contact.push_back(nbr);
            } else if (Status[nbr] == 1) {
                Infector.push_back(nbr);
            }
        }

        // Advance time
        double Tstep = -std::log(r2) / SumRate;
        t += Tstep;

        if (Status[Event] == 2) {
            // Recovery event: I -> R

            // Remove gamma recovery rate for Event
            if (Rate[Event] != 0.0) {
                fw.add(Event, -Rate[Event]);
                SumRate -= Rate[Event];
                Rate[Event] = 0.0;
            }

            // For each susceptible neighbor, reduce infection pressure by beta
            for (int nbr : Contact) {
                // nbr is susceptible, so its infection hazard must drop by beta
                fw.add(nbr, -beta);
                SumRate -= beta;
                Rate[nbr] -= beta;
            }

            // Bookkeeping
            if (TrackDyn) {
                Recovery_time[Event] = t;
                if (debug && (debug_ctr % debug_freq == 0) && (debug_ctr > debug_low) && (debug_ctr < debug_up)) {
                    Rprintf(", %d, \n", 0);
                }
            }

            // Counts
            I_cnt--;
            R_cnt++;

        } else if (Status[Event] == 1) {
            // Infection event: S -> I

            // Event's rate changes from (sum of beta's from infected neighbors) to gamma
            {
                double oldRate = Rate[Event]; // could be zero or >= beta
                double delta = gamma - oldRate;
                fw.add(Event, delta);
                SumRate += delta;
                Rate[Event] = gamma;
            }

            // For each susceptible neighbor, add beta
            for (int nbr : Contact) {
                fw.add(nbr, beta);
                SumRate += beta;
                Rate[nbr] += beta;
            }

            // Bookkeeping
            if (TrackDyn) {
                Infect_time[Event] = t;
                S_NbrDeg[Event] = static_cast<int>(Contact.size());

                int infsize = static_cast<int>(Infector.size());
                int samp_inf = (infsize > 0 ? Infector[0] : Event); // safe fallback

                if (debug && (debug_ctr % debug_freq == 0) && (debug_ctr > debug_low) && (debug_ctr < debug_up)) {
                    Rprintf(", %d, \n", infsize);
                }

                if (infsize > 1) {
                    // Simple uniform draw without Rcpp::sample
                    double r3 = R::runif(0.0, 1.0);
                    int pos = static_cast<int>(std::floor(r3 * infsize));
                    if (pos >= infsize) pos = infsize - 1;
                    samp_inf = Infector[pos];
                    if (debug && (debug_ctr % debug_freq == 0) && (debug_ctr > debug_low) && (debug_ctr < debug_up)) {
                        Rprintf("call samp, %d \n", samp_inf);
                    }
                }

                Infect_num_rnd[samp_inf] += 1;
                Infector_rnd[Event] = samp_inf + 1; // store 1-based id
            }

            // Counts
            S_cnt--;
            I_cnt++;
        }

        Istep = I_cnt;

        if (TrackDyn) {
            t_vec.push_back(t);
            S_vec.push_back(static_cast<double>(S_cnt) / N);
            I_vec.push_back(static_cast<double>(I_cnt) / N);
            R_vec.push_back(static_cast<double>(R_cnt) / N);
        }

        debug_ctr++;
    }

    // Final stats
    DataFrame FinalStat = DataFrame::create(
        Named("FinishTime") = t,
        Named("Ssize") = static_cast<double>(S_cnt) / N,
        Named("Isize") = static_cast<double>(I_cnt) / N,
        Named("Rsize") = static_cast<double>(R_cnt) / N
    );

    if (TrackDyn) {
        DataFrame Details = DataFrame::create(
            Named("t_vec") = t_vec,
            Named("S_vec") = S_vec,
            Named("I_vec") = I_vec,
            Named("R_vec") = R_vec
        );

        DataFrame Reff = DataFrame::create(
            Named("Node") = nodes,
            Named("Degree") = Deg_vec,
            Named("Infect_time") = Infect_time,
            Named("Recovery_time") = Recovery_time,
            Named("S_NbrDeg") = S_NbrDeg,
            Named("Infect_num_rnd") = Infect_num_rnd,
            Named("Infector_rnd") = Infector_rnd
        );

        return List::create(
            Named("FinalStat") = FinalStat,
            Named("Details") = Details,
            Named("Reff") = Reff,
            Named("Init") = InitIndex
        );
    } else {
        return List::create(Named("FinalStat") = FinalStat);
    }
}

