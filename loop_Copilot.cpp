#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void SSA_loop_optimized(  NumericVector Rate
                        , IntegerVector Status
                        , List G
                        , double MaxTime
                        , double b
                        , double g
                        , int Istep
                        , bool TrackDyn
                        , NumericVector Recovery_time
                        , NumericVector Infect_time
                        , IntegerVector S_NbrDeg
                        , IntegerVector Infect_num_rnd
                        , IntegerVector Infector_rnd) {
  
  double t = 0.0;
  int N = Rate.size();
  
  NumericVector Cum(N);
  
  while (t < MaxTime && Istep != 0) {
    double Sum = 0.0;
    
    // Manual sum and cumulative sum for better performance
    Cum[0] = Rate[0];
    Sum += Rate[0];
    for (int i = 1; i < N; ++i) {
      Sum += Rate[i];
      Cum[i] = Cum[i - 1] + Rate[i];
    }
    
    if (Sum <= 0.0) break;
    
    NumericVector r = runif(2);
    double threshold = r[0] * Sum;
    
    // Binary search instead of linear search
    int Event = std::lower_bound(Cum.begin(), Cum.end(), threshold) - Cum.begin();
    Status[Event] += 1;
    
    IntegerVector Neighbor = G[Event];
    std::vector<int> Contact;
    std::vector<int> Infector;
    
    for (int i = 0; i < Neighbor.size(); ++i) {
      int nbr = Neighbor[i];
      int stat = Status[nbr];
      if (stat == 0) {
        Contact.push_back(nbr);
      } else if (stat == 1) {
        Infector.push_back(nbr);
      }
    }
    
    double Tstep = -std::log(r[1]) / Sum;
    t += Tstep;
    
    if (Status[Event] == 2) { // Recovery
      Rate[Event] = 0;
      for (int i : Contact) {
        Rate[i] -= b;
      }
      if (TrackDyn) {
        Recovery_time[Event] = t;
      }
    } else if (Status[Event] == 1) { // Infection
      Rate[Event] = g;
      for (int i : Contact) {
        Rate[i] += b;
      }
      if (TrackDyn) {
        Infect_time[Event] = t;
        S_NbrDeg[Event] = Contact.size();
        
        if (!Infector.empty()) {
          int samp_inf = Infector.size() == 1 ? Infector[0] :
          Infector[floor(R::runif(0, Infector.size()))];
          Infect_num_rnd[samp_inf] += 1;
          Infector_rnd[Event] = samp_inf;
        }
      }
    }
  }
}