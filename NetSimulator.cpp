#include <Rcpp.h>
#include <RcppClock.h>

using namespace Rcpp;
//[[Rcpp::depends(RcppClock)]]

// [[Rcpp::export]]
List GilAlgoCpp(  List adjList
                , int size
                , double beta
                , double gamma
                , double MaxTime
                , int InitInfSize = 1
		, bool TrackDyn = true
		, bool debug = false
		, int debug_freq = 1
    , int debug_low = 500
    , int debug_up = 600
                    ) {

  //Profiling
  //Rcpp::Clock clock;
  
  //clock.tick("Init");
  
  long int debug_ctr = 0;
  long int event_ctr = 0;
  int N = size;
  double t = 0.0;
  IntegerVector Status(N, 0); 
  // 0=S, 1=I, 2=R
  
  NumericVector Rate(N, 0.0);
  IntegerVector Deg_vec(N);
  
  // Pre-compute degrees and probabilities
  NumericVector prob(N);
  for (int i = 0; i < N; ++i) {
    IntegerVector neighbors = adjList[i];
    int deg = neighbors.size();
    Deg_vec[i] = deg;
    prob[i] = deg;
  }
  
  prob = prob/sum(prob);
  
  // Initial situation
  IntegerVector noseq = seq(0, N - 1);
  IntegerVector nodes = seq(1, N);
  IntegerVector InitIndex = Rcpp::sample(noseq, InitInfSize, false, prob);
  
  Rprintf("Initial index %d \n", InitIndex[0]);
  // Initialize output vectors
  NumericVector Infect_time(N, NA_REAL);
  NumericVector Recovery_time(N, NA_REAL);
  IntegerVector Infect_num_rnd(N, 0);
  IntegerVector S_NbrDeg(N, 0);
  IntegerVector Infector_rnd(N, NA_INTEGER);
  
  NumericVector t_vec, S_vec, I_vec, R_vec;

  for (int i = 0; i < InitInfSize; ++i) {
    int idx = InitIndex[i];
    Status[idx] = 1;
    Rate[idx] = gamma;
    if (TrackDyn) {
      Infect_time[idx] = 0.0;
    }
  }
  
  if (TrackDyn) {
    t_vec.push_back(t);
    S_vec.push_back((double)(N - InitInfSize) / N);
    I_vec.push_back((double)InitInfSize / N);
    R_vec.push_back(0.0);
  }
  
  // Update rates for susceptible neighbors of initially infected
  for (int i = 0; i < InitInfSize; ++i) {
    int x = InitIndex[i];
    IntegerVector neighbors = adjList[x];
    for (int j : neighbors) {
      int nbr = j - 1;
      if (Status[nbr] == 0) {
        Rate[nbr] += beta;
        ++S_NbrDeg[x];
      }
    }
  }
  //clock.tock("Init");
  
  int Istep = InitInfSize;
  //clock.tick("Loop");
  while (t < MaxTime && Istep > 0) {
    double SumRate = sum(Rate);
    if (SumRate <= 0.0) break;
    
    //clock.tick("Draw_r");
    NumericVector CumRate = cumsum(Rate);
    //double r1 = R::runif(0, 1);
    //double r2 = R::runif(0, 1);
    NumericVector r = Rcpp::runif(2,0,1);
    double r1 = r[0];
    double r2 = r[1];
    //Rprintf("%ld, %f %f \n", event_ctr, r1, r2);
    int Event = std::lower_bound(CumRate.begin(), CumRate.end(), r1*SumRate) - CumRate.begin();
    //clock.tock("Draw_r");
    
    //clock.tick("Assign+Calc");
    Status[Event] += 1;

    event_ctr++;
    if (debug 
          & (debug_ctr % debug_freq == 0) 
          & (debug_ctr > debug_low)
          & (debug_ctr < debug_up)) {
      Rprintf("%ld %f %f %f %d", event_ctr, t, r1, r2, Event);
    }
    //clock.tock("Assign+Calc");
      
    //clock.tick("Call_Net");
    IntegerVector neighbors = adjList[Event];
    //clock.tock("Call_Net");
    
    //clock.tick("Assign+Calc");
    std::vector<int> Contact; 
    IntegerVector Infector;
    
    for (int j : neighbors) {
      int nbr = j - 1;
      if (Status[nbr] == 0) Contact.push_back(nbr);
      else if (Status[nbr] == 1) Infector.push_back(nbr);
    }
    
    double Tstep = -std::log(r2) / SumRate;
    t += Tstep;
    //clock.tock("Assign+Calc");
    
    if (Status[Event] == 2) { // Recovery
      //clock.tick("Assign+Calc");
      Rate[Event] = 0.0;
      for (int nbr : Contact) {
        Rate[nbr] -= beta;
      }
      if (TrackDyn) {
        Recovery_time[Event] = t;
        if (debug 
              & (debug_ctr % debug_freq == 0) 
              & (debug_ctr > debug_low)
              & (debug_ctr < debug_up)) {
          int infsize = 0;
          Rprintf(", %d, \n", infsize);
        }
      }
      //clock.tock("Assign+Calc");
    } else if (Status[Event] == 1) { // Infection
      
      //clock.tick("Assign+Calc");
      Rate[Event] = gamma;
      
      for (int nbr : Contact) {
        Rate[nbr] += beta; 
      }
      //clock.tock("Assign+Calc");
      
      if (TrackDyn) {
        //clock.tick("Assign+Calc");
        Infect_time[Event] = t;
        S_NbrDeg[Event] = Contact.size();
        
        int infsize = Infector.size();
	      int samp_inf = Infector[0];
	      //clock.tock("Assign+Calc");
	      
        if (debug 
              & (debug_ctr % debug_freq == 0) 
              & (debug_ctr > debug_low)
              & (debug_ctr < debug_up)) {
          Rprintf(", %d, \n", infsize);
        }
        
        if (infsize>1) {
          //clock.tick("Draw_Infector");
          samp_inf = Rcpp::sample(Infector, 1, false)[0];
          //clock.tock("Draw_Infector");
          if (debug 
                & (debug_ctr % debug_freq == 0) 
                & (debug_ctr > debug_low)
                & (debug_ctr < debug_up)){
            Rprintf("call samp, %d \n", samp_inf);
            }
          }
        
        //clock.tick("Assign+Calc");
	      Infect_num_rnd[samp_inf] += 1;
	      Infector_rnd[Event] = samp_inf + 1;
	      //clock.tock("Assign+Calc");
	      } // TrackDyn
      } // infection event
    
    //clock.tick("Assign+Calc");
    Istep = std::count(Status.begin(), Status.end(), 1);
    if (TrackDyn) {
      t_vec.push_back(t);
      S_vec.push_back(std::count(Status.begin(), Status.end(), 0) / (double)N);
      I_vec.push_back(Istep / (double)N);
      R_vec.push_back(std::count(Status.begin(), Status.end(), 2) / (double)N);
    }
    debug_ctr++;
    //clock.tock("Assign+Calc");
  }  // time loop
  //clock.tock("Loop");
  
  //clock.tick("Output");
  DataFrame FinalStat = DataFrame::create(
    Named("FinishTime") = t,
    Named("Ssize") = std::count(Status.begin(), Status.end(), 0) / (double)N,
    Named("Isize") = std::count(Status.begin(), Status.end(), 1) / (double)N,
    Named("Rsize") = std::count(Status.begin(), Status.end(), 2) / (double)N
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
    //clock.tock("Output");
    //clock.stop("profile_cpp");
    return List::create(Named("FinalStat") = FinalStat,
                        Named("Details") = Details,
                        Named("Reff") = Reff,
                        Named("Init") = InitIndex);
  } else {
    //clock.tock("Output");
    //clock.stop("profile_cpp");
    return List::create(Named("FinalStat") = FinalStat);
  }
}