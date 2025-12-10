#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <cstdint>

using namespace Rcpp;

struct Edge { int src; int dst; };

inline int rand_index(int n) {
	int idx = (int)std::floor(R::runif(0.0, 1.0) * (double)n);
	if (idx >= n) idx = n - 1;
	return idx;
}

inline void vuln_add(int *vuln, int &vsize, int vcap, std::vector<int> &pos_map, int eidx) {
	if (pos_map[eidx] >= 0) return;
	if (vsize >= vcap) return; // capacity guard; should not trigger for undirected inputs
	pos_map[eidx] = vsize;
	vuln[vsize] = eidx;
	vsize += 1;
}
inline void vuln_remove(int *vuln, int &vsize, std::vector<int> &pos_map, int eidx) {
	int pos = pos_map[eidx];
	if (pos < 0) return;
	int last_eidx = vuln[vsize - 1];
	vuln[pos] = last_eidx;
	pos_map[last_eidx] = pos;
	vsize -= 1;
	pos_map[eidx] = -1;
}

// [[Rcpp::export]]
List GilAlgoCpp(
	List adjList,
	int size,
	double beta,
	double gamma,
	double MaxTime,
	int InitInfSize = 1
) {
	int N = size;
	double t = 0.0;

	// --- Node state & per-node outputs ---
	IntegerVector Status(N, 0); // 0=S, 1=I, 2=R
	NumericVector Infect_time(N, NA_REAL);
	NumericVector Recovery_time(N, NA_REAL);
	IntegerVector Infect_num(N, 0);

	// --- Initial sampling weights: degree-based ---
	NumericVector prob(N);
	for (int i = 0; i < N; ++i) {
		IntegerVector neighbors = adjList[i];
		prob[i] = neighbors.size();
	}
	double prob_sum = sum(prob);
	if (prob_sum > 0.0) prob = prob / prob_sum;

	// --- Pick initial infected (0-based inside C++) ---
	IntegerVector noseq = seq(0, N - 1);
	IntegerVector InitIndex = Rcpp::sample(noseq, InitInfSize, false, prob);

	// --- Build unique undirected edge set ---
	struct PairHash {
		size_t operator()(const std::pair<int,int>& p) const noexcept {
			return (static_cast<uint64_t>(p.first) << 32) ^ static_cast<uint64_t>(p.second);
		}
	};
	std::unordered_set<std::pair<int,int>, PairHash> seen_pairs; seen_pairs.reserve(N * 4);
	std::vector< std::pair<int,int> > undirected; undirected.reserve(N * 2);
	for (int i = 0; i < N; ++i) {
		IntegerVector neighbors = adjList[i];
		for (int j1 : neighbors) {
			int j = j1 - 1; // 1-based -> 0-based
			if (j == i) continue; // skip self-loops
			int u = (i < j ? i : j);
			int v = (i < j ? j : i);
			if (seen_pairs.insert({u, v}).second) undirected.push_back({u, v});
		}
	}
	int M_undirected = (int)undirected.size();

	// --- Emit directed edges both ways ---
	std::vector<Edge> edges; edges.reserve(std::max(2 * M_undirected, 1));
	std::vector< std::vector<int> > out_edges(N), in_edges(N);
	int ecount = 0;
	for (auto &uv : undirected) {
		int u = uv.first, v = uv.second;
		// u -> v
		edges.push_back({u, v});
		out_edges[u].push_back(ecount);
		in_edges[v].push_back(ecount);
		ecount++;
		// v -> u
		edges.push_back({v, u});
		out_edges[v].push_back(ecount);
		in_edges[u].push_back(ecount);
		ecount++;
	}

	// --- Vulnerable edges preallocated array ---
	std::vector<int> pos_in_vuln(ecount, -1);
	std::vector<int> vuln_storage(std::max(M_undirected, 1), -1);
	int *vuln = vuln_storage.data();
	int vuln_size = 0; // equals |SI| (number of vulnerable edges)

	// --- Infected set (for uniform recovery) ---
	std::vector<int> infected; infected.reserve(N);
	std::vector<int> pos_infected(N, -1);

	// --- Initialize counts ---
	int S_cnt = N - InitInfSize;
	int I_cnt = InitInfSize;
	int R_cnt = 0;
	for (int k = 0; k < InitInfSize; ++k) {
		int idx = InitIndex[k];
		Status[idx] = 1;
		Infect_time[idx] = 0.0;
		pos_infected[idx] = (int)infected.size();
		infected.push_back(idx);
	}

	// Populate initial vulnerable edges from infected sources
	for (int k = 0; k < (int)infected.size(); ++k) {
		int i = infected[k];
		for (int eidx : out_edges[i]) {
			int j = edges[eidx].dst;
			if (Status[j] == 0) vuln_add(vuln, vuln_size, M_undirected, pos_in_vuln, eidx);
		}
	}

	// --- Integer-time logging setup (no TMAX; allocate to ceil(MaxTime)+1) ---
	int Kalloc = (int)std::ceil(MaxTime);
	std::vector<double> t_series(Kalloc + 1);
	std::vector<double> S_series(Kalloc + 1);
	std::vector<double> I_series(Kalloc + 1);
	std::vector<double> R_series(Kalloc + 1);
	std::vector<double> VE_series(Kalloc + 1);
	for (int k = 0; k <= Kalloc; ++k) t_series[k] = (double)k;
	int last_logged = 0;
	S_series[0] = (double)S_cnt / N;
	I_series[0] = (double)I_cnt / N;
	R_series[0] = (double)R_cnt / N;
	VE_series[0] = (double)vuln_size;

	// --- Main loop ---
	while (t < MaxTime && I_cnt > 0) {
		// Total rate
		double lambda = beta * (double)vuln_size + gamma * (double)I_cnt;
		if (lambda <= 0.0) break;

		double r1 = R::runif(0.0, 1.0);
		double r2 = R::runif(0.0, 1.0);
		double Tstep = -std::log(r2) / lambda;
		double t_old = t;
		t += Tstep;
		if (t > MaxTime) t = MaxTime; // enforce ceiling

		// Log integer marks in (t_old, t]
		int start_k = (int)std::floor(t_old) + 1;
		int end_k = (int)std::floor(t);
		if (end_k > Kalloc) end_k = Kalloc;
		for (int k = start_k; k <= end_k; ++k) {
			S_series[k] = (double)S_cnt / N;
			I_series[k] = (double)I_cnt / N;
			R_series[k] = (double)R_cnt / N;
			VE_series[k] = (double)vuln_size;
			last_logged = k;
		}

		bool infection_event = (r1 * lambda < beta * (double)vuln_size);
		if (!infection_event) {
			// Recovery: pick infected uniformly
			int idx_pos = rand_index(I_cnt);
			int i = infected[idx_pos];
			Status[i] = 2;
			Recovery_time[i] = t;
			I_cnt--; R_cnt++;

			// Remove vulnerable edges originating at i
			for (int eidx : out_edges[i]) {
				int j = edges[eidx].dst;
				if (Status[j] == 0) vuln_remove(vuln, vuln_size, pos_in_vuln, eidx);
			}

			// Remove i from infected set (swap-delete)
			int last = infected.back();
			infected[idx_pos] = last;
			pos_infected[last] = idx_pos;
			infected.pop_back();
			pos_infected[i] = -1;

		} else {
			// Infection: pick vulnerable edge uniformly
			if (vuln_size <= 0) continue; // safety
			int vpos = rand_index(vuln_size);
			int eidx = vuln[vpos];
			int i = edges[eidx].src;
			int j = edges[eidx].dst;

			// Infect j
			Status[j] = 1;
			Infect_time[j] = t;
			S_cnt--; I_cnt++;
			Infect_num[i] += 1;

			// Remove vulnerable edges incoming to j (from infected neighbors)
			for (int ein : in_edges[j]) {
				int src = edges[ein].src;
				if (Status[src] == 1) vuln_remove(vuln, vuln_size, pos_in_vuln, ein);
			}

			// Add j to infected set
			pos_infected[j] = (int)infected.size();
			infected.push_back(j);

			// Add vulnerable edges from j -> susceptible neighbors
			for (int eout : out_edges[j]) {
				int dst = edges[eout].dst;
				if (Status[dst] == 0) vuln_add(vuln, vuln_size, M_undirected, pos_in_vuln, eout);
			}
		}
	}

	// Fill remaining marks up to ceil(FinishTime) with final state
	int Kret = (int)std::ceil(t);
	if (Kret > Kalloc) Kret = Kalloc; // guard if rounding pushes over
	for (int k = last_logged + 1; k <= Kret; ++k) {
		S_series[k] = (double)S_cnt / N;
		I_series[k] = (double)I_cnt / N;
		R_series[k] = (double)R_cnt / N;
		VE_series[k] = (double)vuln_size;
	}

	// Wrap outputs: truncate to 0..Kret
	NumericVector t_vec(Kret + 1), S_vec(Kret + 1), I_vec(Kret + 1), R_vec(Kret + 1), VE_vec(Kret + 1);
	for (int k = 0; k <= Kret; ++k) {
		t_vec[k] = t_series[k];
		S_vec[k] = S_series[k];
		I_vec[k] = I_series[k];
		R_vec[k] = R_series[k];
		VE_vec[k] = VE_series[k];
	}

	DataFrame State = DataFrame::create(
		Named("t") = t_vec,
		Named("S") = S_vec,
		Named("I") = I_vec,
		Named("R") = R_vec,
		Named("VE") = VE_vec
	);

	IntegerVector nodes = seq(1, N);
	DataFrame Infector = DataFrame::create(
		Named("Node") = nodes,
		Named("InfectTime") = Infect_time,
		Named("RecoveryTime") = Recovery_time,
		Named("NumInfected") = Infect_num
	);

	DataFrame FinalStat = DataFrame::create(
		Named("FinishTime") = t,
		Named("Ssize") = (double)S_cnt / N,
		Named("Isize") = (double)I_cnt / N,
		Named("Rsize") = (double)R_cnt / N
	);

	return List::create(
		Named("FinalStat") = FinalStat,
		Named("State") = State,
		Named("Infector") = Infector,
		Named("Init") = InitIndex
	);
}
