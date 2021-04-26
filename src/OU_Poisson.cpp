#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector stl_sortb(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}
// [[Rcpp::export]]
NumericVector OU_Poisson_Rcpp_(double t0, int Nsim,double x0, double alpha, double mu)
{
  NumericVector NN=rpois(Nsim,round(mu*t0));
  NumericVector x_sim(Nsim);
  double x, loc;
  int N;
  NumericVector time_events;
  for(int jsim=0;jsim<Nsim;jsim++){
    x=x0;
    N=NN[jsim];
    if(N>0){
    time_events=stl_sortb(runif(N,0,t0));
    loc=exp(-alpha*time_events[0]);
    x=x*loc+1;
    for(int jj=1;jj<N;jj++) {
      loc=exp(-alpha*(time_events[jj]-time_events[jj-1]));
      x=x*loc+1;
      }
    x_sim[jsim]=x*exp(-alpha*(t0-time_events[N-1]));
    }
    else x_sim[jsim]=x*exp(-alpha*t0);
    }
  return x_sim;};

