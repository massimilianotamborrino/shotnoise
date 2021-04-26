#' @useDynLib shotnoise
#' @importFrom Rcpp sourceCpp
NULL


#'@rdname errorOUIG
#'@title errorOUIG
#'@description Integrate Absolute error (IAE) of the shot noise vs the OU and OU-Inverse Gaussian (OU-IG) process
#' with mean x0*exp(-alpha*t) + mu*(1-exp(-alpha*t))/alpha and standard deviation sqrt(tilde.sigma2*(1-exp(-2*alpha*t))/(2*alpha)) at time t.
#'@param t time at which to evaluate the IAE;
#'@param lambda positive constant rate of the underlying Poisson process;
#'@param alpha positive exponential decay rate of the shot noise process
#'@param mu positive quantity entering into the mean and shape parameters of the IG jump distribution Jn (Table 1, Tamborrino and Lansky, Physica D, 418 (2021) 132845, TL2021), corresponding to the limiting first infinitesimal moment
#' eq. (23);
#'@param tilde.sigma2 positive quantity entering into the shape parameter of the IG Jn, corresponding also to the
#' limiting second moment of eq. (24);
#'@param x0 non-negative starting condition of the processes;
#'@param deltax simulation time step;
#'@param NpointsOUIG  number of points to simulate in the OU-IG process.
#'@param Nsim number of simulation of the shot noise process;
#'@return IAEs of the shot noise process at time t with respect to the OU process and the limiting non-Gaussian OU-IG process
#'@export
errorOUIG<-function(t,lambda,alpha,mu,tilde.sigma2,x0,deltax,NpointsOUIG,Nsim){
  simshot<-shotnoise_JnIG(t,Nsim,x0,alpha,mu,tilde.sigma2,lambda)
  simOUIG<-OU_IG(t,NpointsOUIG,Nsim,x0,alpha,mu,tilde.sigma2)
  a<-min(simshot,simOUIG)
  b<-max(simshot,simOUIG)
  npoints<-round((b-a)/deltax)
  D_IG<-density(simshot,from=a,to=b,n=npoints)
  X<-D_IG$x
  D_IGy<-D_IG$y
  N<-length(D_IGy)
  OUIG<-density(simOUIG,from=a,to=b,n=npoints)$y
  asn<-dnorm(X,mean=x0*exp(-alpha*t) + mu*(1-exp(-alpha*t))/alpha,sd=sqrt(tilde.sigma2*(1-exp(-2*alpha*t))/(2*alpha)))
  deltax<-diff(X)[1]
    c(deltax*c(sum(abs(D_IGy-asn))-abs(D_IGy[1]-asn[1])/2-abs(D_IGy[N]-asn[N])/2),
    deltax*c(sum(abs(D_IGy-OUIG))-abs(D_IGy[1]-OUIG[1])/2-abs(D_IGy[N]-OUIG[N])/2))
}


#'@rdname errorOUGamma
#'@title errorOUGamma
#'@description Integrate Absolute error (IAE) of the shot noise vs the OU and OU-Gamma  process
#' with mean x0*exp(-alpha*t) + mu*(1-exp(-alpha*t))/alpha and standard deviation sqrt(tilde.sigma2*(1-exp(-2*alpha*t))/(2*alpha)) at time t.
#' Scenario 1 for Chi-square distributed jumps Jn, scenario 2 for Gamma distributed jumps Jn
#'@param t time at which to evaluate the IAE;
#'@param lambda positive constant rate of the underlying Poisson process;
#'@param alpha positive exponential decay rate of the shot noise process
#'@param mu  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#'@param x0 non-negative starting condition of the processes;
#'@param deltax simulation time step;
#'@param NpointsOUG number of points to simulate in the OU-Gamma process.
#'@param Nsim number of simulation of the shot noise process;
#'@param scenario 1 if Jn are chi-square distributed, 2 if Gamma
#'@return IAEs of the shot noise process at time t with respect to the OU process and the limiting non-Gaussian OU-Gamma process
#'@export
errorOUGamma<-function(t,lambda,alpha,mu,tilde.sigma2,x0,deltax,NpointsOUG,Nsim,scenario){
  if(scenario==1) {simshot<-shotnoise_JnChiSq(t,Nsim,x0,alpha,mu,lambda)
  sigma2<-2*mu}
  if(scenario==2) simshot<-shotnoise_JnG(t,Nsim,x0,alpha,mu,tilde.sigma2,lambda)
  Xsim<-OU_Gamma(t,NpointsOUG,Nsim,x0,alpha,mu,tilde.sigma2,scenario)
  a<-min(simshot,Xsim)
  b<-max(simshot,Xsim)
  npoints<-round((b-a)/deltax)
  D_G<-density(simshot,from=a,to=b,n=npoints)
  X<-D_G$x
  shotsim<-D_G$y
  N<-length(shotsim)
  OUG<-density(Xsim,from=a,to=b,n=npoints)$y
  asn<-dnorm(X,mean=x0*exp(-alpha*t) + mu*(1-exp(-alpha*t))/alpha,sd=sqrt(tilde.sigma2*(1-exp(-2*alpha*t))/(2*alpha)))
  deltax<-diff(X)[1]
  c(deltax*c(sum(abs(shotsim-asn))-abs(shotsim[1]-asn[1])/2-abs(shotsim[N]-asn[N])/2),
    deltax*c(sum(abs(shotsim-OUG))-abs(shotsim[1]-OUG[1])/2-abs(shotsim[N]-OUG[N])/2))
}

#'@rdname errorOUPoisson
#'@title errorOUPoisson
#'@description Integrate Absolute error (IAE) of the shot noise vs the OU and the OU-Poisson process
#' with mean x0*exp(-alpha*t) + mu*(1-exp(-alpha*t))/alpha and standard deviation sqrt(mu*(1-exp(-2*alpha*t))/(2*alpha)) at time t.
#' Scenario 1 for Bernoulli distributed jumps Jn, scenario 2 for Poisson distributed jumps Jn
#'@param t time at which to evaluate the IAE;
#'@param lambda positive constant rate of the underlying Poisson process;
#'@param alpha positive exponential decay rate of the shot noise process
#'@param mu positive quantity entering into the jump distribution Jn;
#'@param x0 non-negative starting condition of the processes;
#'@param deltax simulation time step;
#'@param Nsim number of simulation of the shot noise;
#'@param scenario 1 if Jn are Bernoulli distributed, 2 if Poisson
#'@return IAEs of the shot noise process at time t with respect to the OU process and the limiting non-Gaussian OU-Poisson
#'@export
errorOUPoisson<-function(t,lambda,alpha,mu,x0,deltax,Nsim,scenario){
  if(scenario==1)   simshot<-shotnoise_JnBer(t,Nsim,x0,alpha,mu,lambda)
  if(scenario==2)   simshot<-shotnoise_JnP(t,Nsim,x0,alpha,mu,lambda)
  simOUP<-OU_Poisson(t,Nsim,x0,alpha,mu)
  a<-min(simshot,simOUP)
  b<-max(simshot,simOUP)
  npoints<-round((b-a)/deltax)
  Sim<-density(simshot,from=a,to=b,n=npoints)
  X<-Sim$x
  SN<-Sim$y
  asn<-dnorm(X,mean=x0*exp(-alpha*t) + mu*(1-exp(-alpha*t))/alpha,sd=sqrt(mu*(1-exp(-2*alpha*t))/(2*alpha)))
  N<-length(SN)
  OUP<-density(simOUP,from=a,to=b,n=npoints)$y
  deltax<-diff(X)[1]
  c(deltax*c(sum(abs(SN-asn))-abs(SN[1]-asn[1])/2-abs(SN[N]-asn[N])/2),
    deltax*c(sum(abs(SN-OUP))-abs(SN[1]-OUP[1])/2-abs(SN[N]-OUP[N])/2))
}
