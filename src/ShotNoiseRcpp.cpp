#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
NumericVector cumsum_sug(NumericVector x){
  // initialize the result vector
  NumericVector res = cumsum(x);
  return res;
}

// [[Rcpp::export]]
NumericVector shotnoise_JnP_Rcpp_(double t0, int Nsim,double x0, double alpha, double mu,double lambda)
{
  NumericVector NN=rpois(Nsim,round(lambda*t0));
  NumericVector x_sim(Nsim);
  double x, loc;
  int N;
  NumericVector time_events;
  for(int jsim=0;jsim<Nsim;jsim++){
    x=x0;
    N=NN[jsim];
    if(N>0){
      NumericVector rpo(N);
      time_events=stl_sort(runif(N,0,t0));
      loc=exp(-alpha*time_events[0]);
      rpo=rpois(N,mu/lambda);
      x=x*loc+rpo[0];
      for(int jj=1;jj<N;jj++) {
        loc=exp(-alpha*(time_events[jj]-time_events[jj-1]));
        x=x*loc+rpo[jj];
      }
      x_sim[jsim]=x*exp(-alpha*(t0-time_events[N-1]));
    }
    else x_sim[jsim]=x*exp(-alpha*t0);}
  return x_sim;};

// [[Rcpp::export]]
NumericVector shotnoise_JnBer_Rcpp_(double t0, int Nsim,double x0, double alpha, double mu,double lambda)
{
  NumericVector NN=rpois(Nsim,round(lambda*t0));
  NumericVector x_sim(Nsim);
  double x, loc;
  int N;
  NumericVector time_events;
  for(int jsim=0;jsim<Nsim;jsim++){
    x=x0;
    N=NN[jsim];
    if(N>0){
      NumericVector rpo(N);
      time_events=stl_sort(runif(N,0,t0));
      loc=exp(-alpha*time_events[0]);
      rpo=rbinom(N,1,mu/lambda);
      x=x*loc+rpo[0];
      for(int jj=1;jj<N;jj++) {
        loc=exp(-alpha*(time_events[jj]-time_events[jj-1]));
        x=x*loc+rpo[jj];
      }
      x_sim[jsim]=x*exp(-alpha*(t0-time_events[N-1]));
    }
    else x_sim[jsim]=x*exp(-alpha*t0);}
  return x_sim;};


// [[Rcpp::export]]
NumericVector shotnoise_JnG_Rcpp_(double t0, int Nsim,double x0, double alpha, double mu,double sigma2,double lambda)
{
  NumericVector NN=rpois(Nsim,round(lambda*t0));
  NumericVector x_sim(Nsim);
  double x, loc;
  int N;
  NumericVector time_events;
  for(int jsim=0;jsim<Nsim;jsim++){
    x=x0;
    N=NN[jsim];
    if(N>0){
      NumericVector rpo(N);
      time_events=stl_sort(runif(N,0,t0));
      loc=exp(-alpha*time_events[0]);
      rpo=rgamma(N,(mu/lambda)/(sigma2/mu-mu/lambda),(sigma2/mu-mu/lambda));
      x=x*loc+rpo[0];
      for(int jj=1;jj<N;jj++) {
        loc=exp(-alpha*(time_events[jj]-time_events[jj-1]));
        x=x*loc+rpo[jj];
      }
      x_sim[jsim]=x*exp(-alpha*(t0-time_events[N-1]));
    }
    else x_sim[jsim]=x*exp(-alpha*t0);}
  return x_sim;};
//

// [[Rcpp::export]]
NumericVector shotnoise_JnG_trajectory_Rcpp_(double t0,int N_time, double x0, double alpha, double mu,double sigma2,double lambda)
{
  int N=rpois(1,round(lambda*t0))[0];
  NumericVector x_sim(N_time+1);
  //  NumericVector time_events(1);
  double tau=t0/N_time, cumtau;
  double loc=exp(-alpha*tau);
  int index=0;
  x_sim[0]=x0;
  if(N>0){
    NumericVector time_events(N);
    NumericVector rpo(N);
    time_events=stl_sort(runif(N,0,t0));
    rpo=rgamma(N,mu*mu/(lambda*sigma2),sigma2/mu);
    for(int i=1;i<=N_time;i++){
      cumtau=i*tau;
      x_sim[i]=x_sim[i-1]*loc;
      if(cumtau>=time_events[index]&index<N) {x_sim[i]=x_sim[i]+rpo[index];
        index=index+1;}
    }}
  return x_sim;};

// [[Rcpp::export]]
NumericVector shotnoise_JnChiSq_Rcpp_(double t0, int Nsim,double x0, double alpha, double mu,double lambda)
{
  NumericVector NN=rpois(Nsim,round(lambda*t0));
  NumericVector x_sim(Nsim);
  double x, loc;
  int N;
  NumericVector time_events;
  for(int jsim=0;jsim<Nsim;jsim++){
    x=x0;
    N=NN[jsim];
    if(N>0){
      NumericVector rpo(N);
      time_events=stl_sort(runif(N,0,t0));
      loc=exp(-alpha*time_events[0]);
      rpo=rchisq(N,mu/lambda);
      x=x*loc+rpo[0];
      for(int jj=1;jj<N;jj++) {
        loc=exp(-alpha*(time_events[jj]-time_events[jj-1]));
        x=x*loc+rpo[jj];
      }
      x_sim[jsim]=x*exp(-alpha*(t0-time_events[N-1]));
    }
    else x_sim[jsim]=x*exp(-alpha*t0);}
  return x_sim;};


// [[Rcpp::export]]
NumericVector shotnoise_JnIG_Rcpp_(double t0, int Nsim,double x0, double alpha, double mu,double sigma2,double lambda)
{
  Environment statmod("package:statmod");
  Function f = statmod["rinvgauss"];
  NumericVector NN=rpois(Nsim,round(lambda*t0));
  NumericVector x_sim(Nsim);
  double x, loc;
  int N;
  NumericVector exp_j(N),time_events;
  for(int jsim=0;jsim<Nsim;jsim++){
    x=x0;
    N=NN[jsim];
    if(N>0){
      NumericVector rpo(N);
      time_events=stl_sort(runif(N,0,t0));
      loc=exp(-alpha*time_events[0]);
      rpo=f(N,mu/lambda,pow(mu,3)/(lambda*lambda*sigma2));
      x=x*loc+rpo[0];
      for(int jj=1;jj<N;jj++) {
        loc=exp(-alpha*(time_events[jj]-time_events[jj-1]));
        x=x*loc+rpo[jj];
      }
      x_sim[jsim]=x*exp(-alpha*(t0-time_events[N-1]));
    }
    else x_sim[jsim]=x*exp(-alpha*t0);}
  return x_sim;};


// [[Rcpp::export]]
NumericVector shotnoise_JnIG_trajectory_Rcpp_(double t0,double N_time,double x0, double alpha, double mu,double sigma2,double lambda)
{
  Environment statmod("package:statmod");
  Function f = statmod["rinvgauss"];
  int N=rpois(1,round(lambda*t0))[0],index=0;
  NumericVector x_sim(N_time+1);
  double tau=t0/N_time,cumtau;
  double loc=exp(-alpha*tau);
  NumericVector exp_j(N);
  x_sim[0]=x0;
  if(N>0){
    NumericVector rpo(N),time_events(N);
    time_events=stl_sort(runif(N,0,t0));
    rpo=f(N,mu/lambda,pow(mu,3)/(lambda*lambda*sigma2));
   for(int i=1;i<=N_time;i++){
     cumtau=i*tau;
     x_sim[i]=x_sim[i-1]*loc;
     if(cumtau>=time_events[index]&index<N) {x_sim[i]=x_sim[i]+rpo[index];
       index=index+1;}
   }}
  return x_sim;};


