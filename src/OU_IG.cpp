#include <Rcpp.h>
using namespace Rcpp;

//%%%%%%%%%%% Generating one random variable V based on A/R Non-uniform envelop 1 %%%%%%%%%%%%%%%%%
//  function V=OU_IG_V_generator2(w) %
// [[Rcpp::export]]

double OU_IG_V_generator2(double w){
int accept=0;
double E_e,p,V=0;
while(accept<1){
E_e= pow( sqrt(pow(pow(w,-0.25) -1,2)*runif(1)[0])  +1,4);
p= (pow(E_e,.5)-1 )/(pow(E_e,.5)-pow(E_e,.25))/2;
if (runif(1)[0] <= p){ accept=1;
V=E_e;}
};
return(V);}


// //%%%%%%%%%%% Generating n random variables from an inverse-Gaussian (IG) r.v. by Michael 1976, see Wiki  %%%%%%%%%%%%%%%%%
// [[Rcpp::export]]
NumericVector  function_IG_generator(int n, double mu,double  lambda){
NumericVector   rinv_gaussian(n),y=pow(rnorm(n),2),x(n);
 x=mu + (pow(mu,2)*y)/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*y + pow(mu*y,2));
 NumericVector test=runif(n);
 for (int i=0;i<n;++i){
 if (test[i] <= (mu)/(mu + x[i]))   rinv_gaussian[i] = x[i];
 else  rinv_gaussian[i] = pow(mu,2)/x[i];}
return(rinv_gaussian);}

// [[Rcpp::export]]
NumericVector OU_IG_Rcpp_(double T_terminal,int N_time,int N_path,double X_0,double alpha,double mu,double sigma2){
////%% Parameter Setting
double rho=mu*sqrt(mu/sigma2);
double IG_c, tau,w,Poiss_rate,IG_mu,IG_lambda;
IG_c=sqrt(mu/sigma2); //  % IG parameter
tau=T_terminal/N_time; //% length of each small time interval
w=exp(-alpha*tau);
Poiss_rate=rho*IG_c/alpha*( 2/sqrt(w) -2 -alpha*tau );
IG_mu=2*rho/(alpha*IG_c)*(sqrt(w)-w);// % IG parameter mu
IG_lambda=pow(2*rho/alpha*(1-sqrt(w)),2);// % IG parameter lambda
//  %% Simulations
NumericVector X_T_terminal(N_path); // % record intensity X at time T_terminal
// %% Simulation Loop for N_path Paths
double X_process,V,gamma_2;
NumericVector Poisson_rv(N_time),rinv(N_time);
for(int j_path=0;j_path<N_path;++j_path){
//  %% Simulation Loop for One Path
X_process=X_0; // %lambda process
Poisson_rv=rpois(N_time,Poiss_rate);
rinv=function_IG_generator(N_time,IG_mu,IG_lambda);
for (int i_time=0;i_time<N_time;++i_time){
V=OU_IG_V_generator2(w); //% simulate V by A/R
gamma_2=0.5*pow(IG_c,2)*V; //% Gamma rate parameter
NumericVector S_vector(Poisson_rv[i_time]);
S_vector=rgamma(Poisson_rv[i_time],0.5,1/gamma_2); //% (n,shape,scale)
//S_vector=gamrnd(0.5,1/gamma_2,[Poisson_rv[i_time] 1]); //% gamrnd(shape,scale)
X_process = w*X_process + rinv[i_time] +sum(S_vector);
 }
X_T_terminal[j_path]=X_process;
}
return(X_T_terminal);}

// [[Rcpp::export]]
NumericVector OU_IG_trajectory_Rcpp_(double T_terminal,int N_time,double X_0,double alpha,double mu,double sigma2){
  ////%% Parameter Setting
  double rho=mu*sqrt(mu/sigma2);
  double IG_c, tau,w,Poiss_rate,IG_mu,IG_lambda;
  IG_c=sqrt(mu/sigma2); //  % IG parameter
  tau=T_terminal/N_time; //% length of each small time interval
  w=exp(-alpha*tau);
  Poiss_rate=rho*IG_c/alpha*( 2/sqrt(w) -2 -alpha*tau );
  IG_mu=2*rho/(alpha*IG_c)*(sqrt(w)-w);// % IG parameter mu
  IG_lambda=pow(2*rho/alpha*(1-sqrt(w)),2);// % IG parameter lambda
  //  %% Simulations
  NumericVector X_process(N_time+1); // % record intensity X at time T_terminal
  // %% Simulation Loop for N_path Paths
  double V,gamma_2;
  NumericVector Poisson_rv(N_time),rinv(N_time);
    //  %% Simulation Loop for One Path
    X_process[0]=X_0; // %lambda process
    Poisson_rv=rpois(N_time,Poiss_rate);
    rinv=function_IG_generator(N_time,IG_mu,IG_lambda);
    for (int i_time=0;i_time<N_time;++i_time){
      V=OU_IG_V_generator2(w); //% simulate V by A/R
      gamma_2=0.5*pow(IG_c,2)*V; //% Gamma rate parameter
      NumericVector S_vector(Poisson_rv[i_time]);
      S_vector=rgamma(Poisson_rv[i_time],0.5,1/gamma_2); //% (n,shape,scale)
      //S_vector=gamrnd(0.5,1/gamma_2,[Poisson_rv[i_time] 1]); //% gamrnd(shape,scale)
      X_process[i_time+1] = w*X_process[i_time] + rinv[i_time] +sum(S_vector);
    }
  return(X_process);}

