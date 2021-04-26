#'@rdname shotnoise_JnP
#'@title shotnoise_JnP
#'@description Simulation of Nsim values from the shot noise process Xn(t) with Poissonian distributed jumps Jn at time Tmax
#'@param Tmax time at which we want to sample from;
#'@param Nsim number of simulated values from Xn(Tmax)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the rate of the jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021)
#'@param lambda positive intensity of the underlying Poisson process
#'@return Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being Poisson distributed
#'@export
shotnoise_JnP <- function(Tmax,Nsim,x0,alpha,mu,lambda)   return(shotnoise_JnP_Rcpp_(Tmax,Nsim,x0,alpha,mu,lambda))

#'@rdname shotnoise_JnBer
#'@title shotnoise_JnBer
#'@description Simulation of Nsim values from the shot noise process Xn(t) with Bernoulli distributed jumps Jn at time Tmax
#'@param Tmax time at which we want to sample from;
#'@param Nsim number of simulated values from Xn(Tmax)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the probability of success of the Bernoulli jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021)
#'@param lambda positive intensity of the underlying Poisson process
#'@return Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being Bernoulli distributed
#'@export
shotnoise_JnBer <- function(Tmax,Nsim,x0,alpha,mu,lambda)   return(shotnoise_JnBer_Rcpp_(Tmax,Nsim,x0,alpha,mu,lambda))


#'@rdname shotnoise_JnG
#'@title shotnoise_JnG
#'@description Simulation of Nsim values from the shot noise process Xn(t) with Gamma distributed jumps Jn at time Tmax
#'@param Tmax time at which we want to sample from;
#'@param Nsim number of simulated values from Xn(Tmax)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#'@param lambda positive intensity of the underlying Poisson process
#'@return Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being Gamma distributed
#'@export
shotnoise_JnG <- function(Tmax,Nsim,x0,alpha,mu,tilde.sigma2,lambda)   return(shotnoise_JnG_Rcpp_(Tmax,Nsim,x0,alpha,mu,tilde.sigma2,lambda))


#'@rdname shotnoise_JnChiSq
#'@title shotnoise_JnChiSq
#'@description Simulation of Nsim values from the shot noise process Xn(t) with Chi-square distributed jumps Jn at time Tmax
#'@param Tmax time at which we want to sample from;
#'@param Nsim number of simulated values from Xn(Tmax)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the shape parameter of the Gamma jump amplitude Jn, see point 2 of Section 3.1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#'@param lambda positive intensity of the underlying Poisson process
#'@return Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being chi-square distributed
#'@export
shotnoise_JnChiSq  <- function(Tmax,Nsim,x0,alpha,mu,lambda)   return(shotnoise_JnChiSq_Rcpp_(Tmax,Nsim,x0,alpha,mu,lambda))


#'@rdname shotnoise_JnIG
#'@title shotnoise_JnIG
#'@description Simulation of Nsim values from the shot noise process Xn(t) with Inverse Gaussian distributed jumps Jn at time Tmax
#'@param Tmax time at which we want to sample from;
#'@param Nsim number of simulated values from Xn(Tmax)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the mean and shape parameters of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the shape parameter  of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#'@param lambda positive intensity of the underlying Poisson process
#'@return Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being IG distributed
#'@export
shotnoise_JnIG  <- function(Tmax,Nsim,x0,alpha,mu,tilde.sigma2,lambda)   return(shotnoise_JnIG_Rcpp_(Tmax,Nsim,x0,alpha,mu,tilde.sigma2,lambda))


#'@rdname shotnoise_JnG_trajectory
#'@title shotnoise_JnG_trajectory
#'@description Simulation of a trajectory of the shot noise process Xn(t) with Gamma distributed jumps in [0,Tmax]
#'@param Tmax time at which we want to sample from;
#'@param N_time number of equidistant points to simulate in [0,Tmax]
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#'@param lambda positive intensity of the underlying Poisson process
#'@return N_time equidistant points of the shot noise process Xn(t) with jumps Jn being Gamma distributed, for t in [0,Tmax]
#'@export
shotnoise_JnG_trajectory <- function(Tmax,N_time,x0,alpha,mu,tilde.sigma2,lambda)   return(shotnoise_JnG_trajectory_Rcpp_(Tmax,N_time,x0,alpha,mu,tilde.sigma2,lambda))


#'@rdname shotnoise_JnIG_trajectory
#'@title shotnoise_JnIG_trajectory
#'@description Simulation of a trajectory of the shot noise process Xn(t) with IG distributed jumps in [0,Tmax]
#'@param Tmax time at which we want to sample from;
#'@param N_time number of equidistant points to simulate in [0,Tmax]
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the mean and shape parameters of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the shape parameter  of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#'@param lambda positive intensity of the underlying Poisson process
#'@return Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being IG distributed
#'@export
shotnoise_JnIG_trajectory <- function(Tmax,N_time,x0,alpha,mu,tilde.sigma2,lambda)   return(shotnoise_JnIG_trajectory_Rcpp_(Tmax,N_time,x0,alpha,mu,tilde.sigma2,lambda))
