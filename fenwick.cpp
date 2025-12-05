#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cmath>

using namespace Rcpp;

// ---------- Fenwick (Binary Indexed) Tree over doubles ----------
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

    // State: 0=S, 1=I, 2=R
    IntegerVector Status(N, 0);
    IntegerVector Deg_vec(N);

    // Outputs
    IntegerVector nodes = seq(1, N);
    NumericVector Infect_time(N, NA_REAL);
    NumericVector Recovery_time(N, NA_REAL);
    IntegerVector Infect_num_rnd(N, 0);
    IntegerVector S_NbrDeg(N, 0);
    IntegerVector Infector_rnd(N, NA_INTEGER);

    NumericVector t_vec, S_vec, I_vec, R_vec;

    // Pre-compute degrees and probabilities for initial infected
    NumericVector prob(N);
    for (int i = 0; i < N; ++i) {
        IntegerVector neighbors = adjList[i];
        int deg = neighbors.size();
        Deg_vec[i] = deg;
        prob[i] = deg;
    }
    prob = prob / sum(prob);

    // Pick initial infected nodes (0-based indices)
    IntegerVector noseq = seq(0, N - 1);
    IntegerVector InitIndex = Rcpp::sample(noseq, InitInfSize, false, prob);
    Rprintf("Initial index %d \n", InitIndex[0]);

    // Initialize counts
    int S_cnt = N - InitInfSize;
    int I_cnt = InitInfSize;
    int R_cnt = 0;

    // Mark initial infected
    for (int k = 0; k < InitInfSize; ++k) {
        int idx = InitIndex[k];
        Status[idx] = 1;
        Infect_time[idx] = 0.0;
    }

    // Infected set with O(1) remove via swap
    std::vector<int> infected;
    infected.reserve(N);
    std::vector<int> pos_infected(N, -1);
    for (int k = 0; k < InitInfSize; ++k) {
        int idx = InitIndex[k];
        pos_infected[idx] = (int)infected.size();
        infected.push_back(idx);
    }

    // kS[i] = number of susceptible neighbors of i (defined only for infected nodes)
    std::vector<int> kS(N, 0);
    Fenwick fw_kS(N);
    double SI_edges = 0.0;

    // Compute initial kS for infected nodes and load Fenwick
    for (int k = 0; k < (int)infected.size(); ++k) {
        int i = infected[k];
        int countS = 0;
        IntegerVector neighbors = adjList[i];
        for (int j1 : neighbors) {
            int j = j1 - 1; // 1-based -> 0-based
            if (Status[j] == 0) countS++;
        }
        kS[i] = countS;
        S_NbrDeg[i] = countS; // optional: record baseline susceptible-degree for initial infected
        if (countS > 0) {
            fw_kS.add(i, (double)countS);
            SI_edges += countS;
        }
    }

    // Dynamics tracking
    if (TrackDyn) {
        t_vec.push_back(t);
        S_vec.push_back((double)S_cnt / N);
        I_vec.push_back((double)I_cnt / N);
        R_vec.push_back((double)R_cnt / N);
    }

    // Main loop
    while (t < MaxTime && I_cnt > 0) {
        // Total rate
        double lambda = beta * SI_edges + gamma * I_cnt;
        if (lambda <= 0.0) break;

        // Draw waiting time and event type
        double r1 = R::runif(0.0, 1.0);
        double r2 = R::runif(0.0, 1.0);
        double Tstep = -std::log(r2) / lambda;
        t += Tstep;

        bool infection_event = (r1 * lambda < beta * SI_edges);

        if (infection_event) {
            // If no SI edges, fallback to recovery
            if (SI_edges <= 0.0) {
                // treat as recovery
                infection_event = false;
            }
        }

        if (!infection_event) {
            // ---------- Recovery: pick infected uniformly ----------
            int idx_infected = (int)std::floor(R::runif(0.0, 1.0) * I_cnt);
            if (idx_infected >= I_cnt) idx_infected = I_cnt - 1;
            int i = infected[idx_infected];

            // Update outputs/counts
            Status[i] = 2;
            Recovery_time[i] = t;
            I_cnt--;
            R_cnt++;

            // Remove i from Fenwick and SI count
            if (kS[i] > 0) {
                fw_kS.add(i, -(double)kS[i]);
                SI_edges -= kS[i];
            }
            kS[i] = 0;

            // Remove i from infected vector (swap-delete)
            int last = infected.back();
            infected[idx_infected] = last;
            pos_infected[last] = idx_infected;
            infected.pop_back();
            pos_infected[i] = -1;

            if (debug && (debug_ctr % debug_freq == 0) && (debug_ctr > debug_low) && (debug_ctr < debug_up)) {
                Rprintf("REC %ld t=%.6f node=%d I=%d SI=%.0f\n",
                        ++event_ctr, t, i, I_cnt, SI_edges);
            }

        } else {
            // ---------- Infection: pick infected by kS weight, then a susceptible neighbor ----------
            double target = R::runif(0.0, 1.0) * SI_edges;
            if (target <= 0.0) target = std::numeric_limits<double>::min();
            int i = fw_kS.findByCumulative(target);

            // Build list of susceptible neighbors of i
            IntegerVector neighbors_i = adjList[i];
            std::vector<int> susNbrs;
            susNbrs.reserve(neighbors_i.size());
            for (int j1 : neighbors_i) {
                int j = j1 - 1;
                if (Status[j] == 0) susNbrs.push_back(j);
            }

            // Guard (shouldn't happen if kS[i] is maintained correctly)
            if (susNbrs.empty()) {
                continue; // skip silently; in well-maintained state this won't trigger
            }

            // Choose recipient j uniformly among susceptible neighbors of i
            int pick = (int)std::floor(R::runif(0.0, 1.0) * (int)susNbrs.size());
            if (pick >= (int)susNbrs.size()) pick = (int)susNbrs.size() - 1;
            int j = susNbrs[pick];

            // Infect j
            Status[j] = 1;
            Infect_time[j] = t;
            Infector_rnd[j] = i + 1;      // store 1-based id
            Infect_num_rnd[i] += 1;
            S_cnt--;
            I_cnt++;

            // Add j to infected set
            pos_infected[j] = (int)infected.size();
            infected.push_back(j);

            // Compute j's susceptible-degree (kS[j]) at infection time
            int kSj = 0;
            IntegerVector neighbors_j = adjList[j];
            for (int u1 : neighbors_j) {
                int u = u1 - 1;
                if (Status[u] == 0) kSj++;
            }
            kS[j] = kSj;
            S_NbrDeg[j] = kSj;

            // Update Fenwick & SI_edges: add j's kS
            if (kSj > 0) {
                fw_kS.add(j, (double)kSj);
                SI_edges += kSj;
            }

            // For each infected neighbor v of j, its kS[v] decreases by 1 (edge v-j goes SI -> II)
            for (int v1 : neighbors_j) {
                int v = v1 - 1;
                if (Status[v] == 1 && v != j) {
                    if (kS[v] > 0) {
                        kS[v] -= 1;
                        fw_kS.add(v, -1.0);
                        SI_edges -= 1.0;
                    }
                }
            }

            if (debug && (debug_ctr % debug_freq == 0) && (debug_ctr > debug_low) && (debug_ctr < debug_up)) {
                Rprintf("INF %ld t=%.6f inf=%d sus=%d I=%d SI=%.0f\n",
                        ++event_ctr, t, i, j, I_cnt, SI_edges);
            }
        }

        if (TrackDyn) {
            t_vec.push_back(t);
            S_vec.push_back((double)S_cnt / N);
            I_vec.push_back((double)I_cnt / N);
            R_vec.push_back((double)R_cnt / N);
        }

        debug_ctr++;
    }

    // Final stats
    DataFrame FinalStat = DataFrame::create(
        Named("FinishTime") = t,
        Named("Ssize") = (double)S_cnt / N,
        Named("Isize") = (double)I_cnt / N,
        Named("Rsize") = (double)R_cnt / N
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
            Named("Infector_rnd") = Infector_rnd,
            Named("Init") = InitIndex
        );

        return List::create(
            Named("FinalStat") = FinalStat,
            Named("Details") = Details,
            Named("Reff") = Reff
        );
    } else {
        return List::create(Named("FinalStat") = FinalStat,
                            Named("Init") = InitIndex);
    }
}

