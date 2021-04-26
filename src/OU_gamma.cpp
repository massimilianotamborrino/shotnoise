#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector OU_Gamma_Rcpp_(double T_terminal,int N_time, int N_path, double Y_0, double delta, double mu, double sigma2,int scenario){
  double rho=1;
  double alpha,beta,tau,w,Poiss_rate,gamma_1,gamma_2;
  if(scenario<2) {alpha=mu*0.5;
    beta=0.5;}
  if(scenario>1) {alpha=pow(mu,2)/sigma2;
    beta=mu/sigma2;}
  tau=T_terminal/N_time;// length of each small time interval
  w=exp(-delta*tau);
  Poiss_rate=0.5*alpha*rho*delta*pow(tau,2.0);
  gamma_1=alpha*rho*tau; //% Gamma shape parameter
  gamma_2=beta*exp(delta*tau); //% Gamma rate parameter
  //Simulations
  NumericVector   Y_T_terminal(N_path);  //% record intensity Y at time T_terminal
  //       %% Simulation Loop for N_path Paths
  double Y_process,V;  //%Y process
  NumericVector gammarn(N_time);
  for(int j_path=0;j_path<N_path;++j_path){
    //          %% Simulation Loop for One Path
    Y_process=Y_0;
    double Poisson_rv;
    //return(Poisson_rv);
    gammarn=rgamma(N_time,gamma_1,1/gamma_2); //% (n,shape,scale)
    double a;
    for(int i_time=0;i_time<N_time;++i_time){
      Poisson_rv=rpois(1,Poiss_rate)[0];
      {if(Poisson_rv>0){
        a=sqrt(runif(1))[0];
        V=1/pow(w,a);
      NumericVector log2(Poisson_rv), S_vector(Poisson_rv);
      log2=log(runif(Poisson_rv));
      S_vector=-1/(beta*V)*log2; //% exponential r.v. of rate beta*V
      Y_process = w*Y_process +  gammarn[i_time] + sum(S_vector);}
      else Y_process = w*Y_process +  gammarn[i_time];}
    }
    Y_T_terminal[j_path]=Y_process;
  }
  return(Y_T_terminal);
}

// [[Rcpp::export]]
NumericVector OU_Gamma_trajectory_Rcpp_(double T_terminal,int N_time, double Y_0, double delta, double mu, double sigma2,int scenario){
  double rho=1;
  double alpha,beta,tau,w,Poiss_rate,gamma_1,gamma_2;
  if(scenario<2) {alpha=mu*0.5;
    beta=0.5;}
  if(scenario>1) {alpha=pow(mu,2)/sigma2;
    beta=mu/sigma2;}
  tau=T_terminal/N_time;// length of each small time interval
  w=exp(-delta*tau);
  Poiss_rate=0.5*alpha*rho*delta*pow(tau,2.0);
  gamma_1=alpha*rho*tau; //% Gamma shape parameter
  gamma_2=beta*exp(delta*tau); //% Gamma rate parameter
  //Simulations
  NumericVector   Y_process(N_time+1);  //% record intensity Y at time T_terminal
  //       %% Simulation Loop for N_path Paths
  double V;  //%Y process
  NumericVector gammarn(N_time),Poisson_rv(N_time);
    Y_process[0]=Y_0;
    gammarn=rgamma(N_time,gamma_1,1/gamma_2); //% (n,shape,scale)
    double a;
    Poisson_rv=rpois(N_time,Poiss_rate);
    for(int i_time=0;i_time<N_time;++i_time){
      {if(Poisson_rv[i_time]>0){
        a=sqrt(runif(1))[0];
        V=1/pow(w,a);
        NumericVector log2(Poisson_rv[i_time]), S_vector(Poisson_rv[i_time]);
        log2=log(runif(Poisson_rv[i_time]));
        S_vector=-1/(beta*V)*log2; //% exponential r.v. of rate beta*V
        Y_process[i_time+1] = w*Y_process[i_time] +  gammarn[i_time] + sum(S_vector);}
      else Y_process[i_time+1] = w*Y_process[i_time] +  gammarn[i_time];}
    }
  return(Y_process);
}


