#' @useDynLib shotnoise
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname OU_Gamma
#'@title OU_Gamma
#'@description Simulation of Nsim values from Y(Tmax), with Y(t) being non-Gaussian OU-Gamma process (13) of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021),
#' using the algorithm in Y. Qu, A. Dassios, H. Zhao, J. Oper. Res. Soc. (2019)
#' dY(t) = -alpha Y(t)dt + 1*dZ(t), with Z(t) being a gamma process with
#' Lévy measure nu(dx)=tilde.alpha x^(-1)+exp(-beta x)dx,
#' with (tilde.alpha,beta)=
#' (mu/2,1/2) if the OU-Gamma process has been obtained as limit of shot noise with chi-square distributed jumps,
#' or (mu^2/tilde.sigma^2, mu/tilde.sigma^2) if the OU-Gamma process has been obtained as limit of shot noise with Gamma distributed jumps;
#'  (see Section 3.2 and Table 1 of TL2021)
#'@param Tmax time at which we want to sample from;
#'@param Ntime total number of time discretisation points of the OU-Gamma;
#'@param Nsim number of simulated values;
#'@param x0 starting condition of the processes;
#'@param alpha positive exponential decay rate;
#'@param alpha positive exponential decay rate of the shot noise process
#'@param mu positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (23), see also Table 1
#'@param scenario 1 if beta=0.5 (OU-Gamma obtained as limiting process of a shot noise with Chi-square distributed jumps) or 2 if beta=mu/tilde.sigma^2 (OU-Gamma obtained as limiting process of a shot noise with Gamma distributed jumps)
#'@return Nsim simulated values of the OU-Gamma process Y at time Tmax
#'@export
OU_Gamma <- function(Tmax,Ntime,Nsim,x0,alpha,mu,tilde.sigma2,scenario)   return(OU_Gamma_Rcpp_(Tmax,Ntime,Nsim,x0,alpha,mu,tilde.sigma2,scenario))


#'@rdname OU_Gamma_trajectory
#'@title OU_Gamma_trajectory
#'@description Simulation of a trajectory of the OU-Gamma process Y(t) for t in [0,Tmax].
#' Y(t) is the non-Gaussian OU-Gamma process (13) of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021),
#' using the algorithm in Y. Qu, A. Dassios, H. Zhao, J. Oper. Res. Soc. (2019)
#' dY(t) = -alpha Y(t)dt + 1*dZ(t), with Z(t) being a gamma process with
#' Lévy measure nu(dx)=tilde.alpha x^(-1)+exp(-beta x)dx,
#' with (tilde.alpha,beta)=(mu^2/tilde.sigma^2,beta),
#' with beta=0.5  if the OU-Gamma process has been obtained as limit of shot noise with chi-square distributed jumps,
#' or beta= mu/tilde.sigma^2 if the OU-Gamma process has been obtained as limit of shot noise with Gamma distributed jumps;
#'  (see Section 3.2 and Table 1 of TL2021)
#'@param Tmax time at which we want to sample from;
#'@param Ntime total number of time discretisation points of the OU-Gamma;
#'@param x0 non-negative starting condition of the processes;
#'@param alpha positive exponential decay rate;
#'@param alpha positive exponential decay rate of the shot noise process
#'@param mu positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (23), see also Table 1
#'@param scenario 1 if beta=0.5 (OU-Gamma obtained as limiting process of a shot noise with Chi-square distributed jumps) or 2 if beta=mu/tilde.sigma^2 (OU-Gamma obtained as limiting process of a shot noise with Gamma distributed jumps)
#' The length of each time-interval is given by Tmax/Ntime
#' @return Trajectory of the OU-Gamma process in [0,Tmax] with Ntime points
#' @export
OU_Gamma_trajectory <- function(Tmax,Ntime,x0,alpha,mu,tilde.sigma2,scenario)   return(OU_Gamma_trajectory_Rcpp_(Tmax,Ntime,x0,alpha,mu,tilde.sigma2,scenario))

#'@rdname OU_IG
#'@title OU_IG
#'@description Simulation of Nsim values from Y(Tmax), with Y(t) being non-Gaussian OU-Inverse Gaussian (OU-IG) process (13) of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021),
#' using the algorithm in Y. Qu, A. Dassios, H. Zhao, J. Appl. Probab. (2021) [QDZ2021].
#' dY(t) = -alpha Y(t)dt + 1*dZ(t), with Z(t) being an Inverse-Gaussian process with
#' Lévy measure nu(dx)=exp(-c^2 x/2)/sqrt(2*pi*x^3)dx, with c=sqrt(mu/tilde.sigma^2), rho=s=(from [QDZ2021])= sqrt(mu^3/tilde.sigma^2)
#'  (see Section 3.2 and Table 1 of TL2021)
#'@param Tmax time at which we want to sample from;
#'@param Ntime number of equidistant points to simulate in [0,Tmax]
#'@param Nsim number of simulated values from Y(t)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the mean and shape parameters of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the shape parameter  of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#'@export
OU_IG <- function(Tmax,Ntime,Nsim,x0,alpha,mu,tilde.sigma2)   return(OU_IG_Rcpp_(Tmax,Ntime,Nsim,x0,alpha,mu,tilde.sigma2))


#'@rdname OU_IG_trajectory
#'@title OU_IG_trajectory
#'@description Simulation of a trajectory of the OU-IG process Y(t) for t in [0,Tmax],
#'with Y(t) being non-Gaussian OU-Inverse Gaussian (OU-IG) process (13) of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021),
#' using the algorithm in Y. Qu, A. Dassios, H. Zhao, J. Appl. Probab. (2021) [QDZ2021].
#' dY(t) = -alpha Y(t)dt + 1*dZ(t), with Z(t) being an Inverse-Gaussian process with
#' Lévy measure nu(dx)=exp(-c^2 x/2)/sqrt(2*pi*x^3)dx, with c=sqrt(mu/tilde.sigma^2), rho=s=(from [QDZ2021])= sqrt(mu^3/tilde.sigma^2)
#'  (see Section 3.2 and Table 1 of TL2021)
#'@param Tmax time at which we want to sample from;
#'@param Ntime number of equidistant points to simulate in [0,Tmax]
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu  positive quantity entering into the mean and shape parameters of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit first infinitesimal moment, eq. (23)
#'@param tilde.sigma2  positive quantity entering into the shape parameter  of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
#' This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1
#' #' The length of each time-interval is given by Tmax/Ntime
#'@export
OU_IG_trajectory <- function(Tmax,Ntime,x0,alpha,mu,tilde.sigma2)   return(OU_IG_trajectory_Rcpp_(Tmax,Ntime,x0,alpha,mu,tilde.sigma2))


#'@rdname OU_Poisson
#'@title OU_Poisson
#'@description Simulation of Nsim values from Y(Tmax), with Y(t) being non-Gaussian OU-Poisson (see Tamborrino and Lansky, Physica D, 418 (2021) 132845, TL2021),
#' i.e. an OU process driven by a Poisson process with jump 1 and intensity mu>0
#'@param Tmax time at which we want to sample from;
#'@param Nsim number of simulated values from Y(Tmax)
#'@param x0 non-negative starting condition
#'@param alpha positive exponential decay rate
#'@param mu positive intensity of the underlying Poisson process
#'@export
OU_Poisson <- function(Tmax,Nsim,x0,alpha,mu)   return(OU_Poisson_Rcpp_(Tmax,Nsim,x0,alpha,mu))
