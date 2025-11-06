#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List GilAlgoCpp(List adjList, int size, double beta, double gamma, double MaxTime,
                int InitInfSize = 1, bool TrackDyn = true) {
  int N = size;
  double t = 0.0;
  NumericVector Status(N, 0); // 0=S, 1=I, 2=R
  NumericVector Rate(N, 0.0);
  NumericVector Deg_vec(N);
  
  // Compute degrees
  for (int i = 0; i < N; ++i) {
    IntegerVector neighbors = adjList[i];
    Deg_vec[i] = neighbors.size();
  }
  
  // Sample initial infected nodes with probability proportional to degree
  NumericVector prob = Deg_vec / sum(Deg_vec);
  IntegerVector nodes = seq(0, N - 1);
  IntegerVector InitIndex = Rcpp::sample(nodes, InitInfSize, false, prob);
  
  for (int i = 0; i < InitInfSize; ++i) {
    Status[InitIndex[i]] = 1;
    Rate[InitIndex[i]] = gamma;
  }
  
  NumericVector Infect_time(N, NA_REAL);
  NumericVector Recovery_time(N, NA_REAL);
  IntegerVector Infect_num_rnd(N, 0);
  IntegerVector S_NbrDeg(N, 0);
  IntegerVector Infector_rnd(N, NA_INTEGER);
  
  NumericVector t_vec, S_vec, I_vec, R_vec;
  if (TrackDyn) {
    t_vec.push_back(t);
    S_vec.push_back(sum(Status == 0) / N);
    I_vec.push_back(sum(Status == 1) / N);
    R_vec.push_back(0.0);
    for (int i = 0; i < InitInfSize; ++i) {
      Infect_time[InitIndex[i]] = 0.0;
    }
  }
  
  // Update rates for susceptible neighbors of initially infected
  for (int i = 0; i < InitInfSize; ++i) {
    int x = InitIndex[i];
    IntegerVector neighbors = adjList[x];
    IntegerVector Contact;
    for (int j = 0; j < neighbors.size(); ++j) {
      int nbr = neighbors[j];
      if (Status[nbr] == 0) {
        Contact.push_back(nbr);
      }
    }
    S_NbrDeg[x] = Contact.size();
    for (int j = 0; j < Contact.size(); ++j) {
      Rate[Contact[j]] += beta;
    }
  }
  
  int Istep = sum(Status == 1);
  while (t < MaxTime && Istep != 0) {
    double SumRate = sum(Rate);
    if (SumRate <= 0.0) break;
    
    NumericVector CumRate = cumsum(Rate);
    double r1 = R::runif(0, 1);
    double r2 = R::runif(0, 1);
    int Event = 0;
    while (Event < N && CumRate[Event] < r1 * SumRate) ++Event;
    
    Status[Event] += 1;
    IntegerVector neighbors = adjList[Event];
    IntegerVector Contact;
    IntegerVector Infector;
    for (int j = 0; j < neighbors.size(); ++j) {
      int nbr = neighbors[j];
      if (Status[nbr] == 0) Contact.push_back(nbr);
      if (Status[nbr] == 1) Infector.push_back(nbr);
    }
    
    double Tstep = -log(r2) / SumRate;
    t += Tstep;
    
    if (Status[Event] == 2) { // Recovery
      Rate[Event] = 0.0;
      for (int j = 0; j < Contact.size(); ++j) {
        Rate[Contact[j]] -= beta;
      }
      if (TrackDyn) Recovery_time[Event] = t;
    } else if (Status[Event] == 1) { // Infection
      Rate[Event] = gamma;
      for (int j = 0; j < Contact.size(); ++j) {
        Rate[Contact[j]] += beta;
      }
      if (TrackDyn) {
        Infect_time[Event] = t;
        S_NbrDeg[Event] = Contact.size();
        if (Infector.size() > 0) {
          int samp_inf = Infector.size() == 1 ? Infector[0] : Infector[floor(R::runif(0, Infector.size()))];
          Infect_num_rnd[samp_inf] += 1;
          Infector_rnd[Event] = samp_inf;
        }
      }
    }
    
    Istep = sum(Status == 1);
    if (TrackDyn) {
      t_vec.push_back(t);
      S_vec.push_back(sum(Status == 0) / N);
      I_vec.push_back(sum(Status == 1) / N);
      R_vec.push_back(sum(Status == 2) / N);
    }
  }
  
  DataFrame FinalStat = DataFrame::create(
    Named("FinishTime") = t,
    Named("Ssize") = sum(Status == 0) / N,
    Named("Isize") = sum(Status == 1) / N,
    Named("Rsize") = sum(Status == 2) / N
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
    
    return List::create(Named("FinalStat") = FinalStat,
                        Named("Details") = Details,
                        Named("Reff") = Reff);
  } else {
    return List::create(Named("FinalStat") = FinalStat);
  }
}
